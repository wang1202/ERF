#include "ERF_SDWallLossHistory.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <random>
#include <sstream>
#include <stdexcept>
#include <system_error>
#include <utility>

namespace sdwl {
namespace {

namespace fs = std::filesystem;
using amrex::Real;
using Mass = std::array<Real,3>;
constexpr int max_bins = 64;
constexpr char sidecar[] = "SDWallLossHistory";
constexpr char version[] = "ERF SD wall loss checkpoint v1";
constexpr std::uint64_t basis = UINT64_C(14695981039346656037);

void require (bool ok, const std::string& why)
{
    if (!ok) { throw std::runtime_error(why); }
}

std::ostringstream output ()
{
    std::ostringstream os;
    os.exceptions(std::ios::badbit | std::ios::failbit);
    os.imbue(std::locale::classic());
    os << std::setprecision(std::numeric_limits<Real>::max_digits10);
    return os;
}

bool single_line (const std::string& s)
{
    return !s.empty() && s.find_first_of("\r\n\0", 0, 3) == std::string::npos;
}

void check_identity (const std::string& name, const std::string& signature)
{
    require(single_line(name) && name != "." && name != ".." &&
            name.find_first_of("/\\:") == std::string::npos, "expected a particle container basename");
    require(single_line(signature), "expected a nonempty single-line configuration signature");
}

std::string new_lineage ()
{
    // Host entropy on the IOProcessor only. Never touch AMReX's particle RNG.
    std::random_device entropy;
    std::uniform_int_distribution<unsigned> byte(0, 255);
    auto os = output();
    os << std::hex << std::setfill('0');
    for (int i = 0; i < 16; ++i) { os << std::setw(2) << byte(entropy); }
    return os.str();
}

bool nonnegative (Real x) { return std::isfinite(x) && x >= 0; }

bool valid_interval (Real time, Real dt)
{
    return std::isfinite(time) && std::isfinite(dt) && dt > 0 &&
           std::isfinite(time + dt) && time + dt > time;
}

// As in ERF's SBM level clock, bound roundoff by a fraction of a physical step.
bool same_time (Real a, Real b, Real width)
{
    // ERF advances in double. Rounding start and dt to particle-call Real
    // before addition can put the result one float away from the rounded end.
    if constexpr (sizeof(Real) < sizeof(double)) {
        if (width > 0 && std::nextafter(a, b) == b) { return true; }
    }
    const Real scale = std::max({Real(1), std::abs(a), std::abs(b)});
    return std::abs(a-b) <= std::min(width/32, 16*std::numeric_limits<Real>::epsilon()*scale);
}

// Do not follow symlinks (including dangling ones), open a FIFO, or mistake a
// permission error for a missing file. Filesystem exceptions carry the path.
bool regular_file_present (const std::string& path)
{
    const auto status = fs::symlink_status(path);
    if (status.type() == fs::file_type::not_found) { return false; }
    require(fs::is_regular_file(status), "not a regular file: " + path);
    return true;
}

// ERF::WriteCheckpointFile writes doubles at precision(17). istep[0] counts
// completed steps; t_new[0] is their END time, not the sampled start time.
// A standalone particle checkpoint has no such Header.
bool native_clock (const std::string& dir, int& step, double& time)
{
    const auto path = dir + "/Header";
    if (!regular_file_present(path)) { return false; }
    std::ifstream in(path);
    std::string line;
    auto next = [&] {
        require(bool(std::getline(in, line)), "missing level-zero clock in " + path);
        if (!line.empty() && line.back() == '\r') { line.pop_back(); }
    };
    next();
    require(line == "Checkpoint file for ERF", "invalid " + path);
    for (int i = 0; i < 5; ++i) { next(); } // finest, conserved, x/y/z counts
    auto first = [&](auto& value) {
        next();
        std::istringstream field(line);
        field.imbue(std::locale::classic());
        require(bool(field >> value) && (field.eof() || field.peek() == ' ' || field.peek() == '\t'),
                "invalid level-zero clock in " + path);
    };
    first(step);
    next(); // native dt need not be the last particle dt
    first(time);
    require(step >= 0 && std::abs(time) <= std::numeric_limits<Real>::max(), "invalid level-zero step/time in " + path);
    return true;
}

// Hash the entire physical segment (including its header), not just its tail.
// FNV-1a detects accidental edits; it is not an authentication mechanism.
struct Fingerprint {
    std::uint64_t bytes = 0, hash = basis;
    void add (const char* p, std::size_t n) {
        bytes += n;
        for (std::size_t i = 0; i < n; ++i) {
            hash = (hash ^ static_cast<unsigned char>(p[i])) * UINT64_C(1099511628211);
        }
    }
};

std::string read_sidecar (const std::string& path, const std::string& name, const std::string& signature)
{
    require(regular_file_present(path), "missing required " + path + " (no legacy history reset)");
    // At most 64 bin values, plus fixed metadata and two quoted identifiers.
    // A mistakenly supplied plotfile/log must not be slurped into memory.
    const auto n = fs::file_size(path);
    require(n > 0 && n <= 16384 + 2*(name.size() + signature.size()), "empty or oversized sidecar: " + path);
    std::ifstream in(path, std::ios::binary);
    require(in.is_open(), "cannot open required " + path + " (no legacy history reset)");
    std::string contents(static_cast<std::size_t>(n), '\0');
    require(in.read(contents.data(), static_cast<std::streamsize>(n)) && in.peek() == std::istream::traits_type::eof()
            && !in.bad(), "incomplete or changing sidecar: " + path);
    return contents;
}

// Binary mode keeps the checkpointed bytes identical across platforms. 'x'
// is essential: a fresh writer must never truncate an existing path.
void write_new (const std::string& path, const std::string& text)
{
    auto* f = std::fopen(path.c_str(), "wbx");
    if (!f) {
        const int error = errno;
        throw std::runtime_error("cannot exclusively create " + path + ": " + std::strerror(error));
    }
    const auto n = std::fwrite(text.data(), 1, text.size(), f);
    const int status = std::fclose(f); // includes flushing buffered writes
    require(n == text.size() && status == 0, "write or close failed: " + path);
}

// Serialize cooperating writers across check + append/publish. A process
// killed mid-write leaves an explicit lock, never an invitation to overwrite.
template <typename F>
void with_writer_lock (const std::string& path, F action)
{
    const auto lock = path + ".lock";
    require(fs::create_directory(lock), "writer lock exists: " + lock + " (check for an active writer before removal)");
    try {
        action(lock);
    } catch (...) {
        std::error_code ignored;
        fs::remove(lock, ignored); // leave a nonempty lock (unpublished sidecar) for inspection
        throw;
    }
    require(fs::remove(lock), "writer lock disappeared: " + lock);
}

bool matches (std::istream& file, const std::string& expected, bool at_end, const std::string& path)
{
    const auto n = static_cast<std::streamsize>(expected.size());
    file.seekg(at_end ? -n : 0, at_end ? std::ios::end : std::ios::beg);
    std::string actual(expected.size(), '\0');
    require(bool(file.read(actual.data(), n)), "cannot read history: " + path);
    return actual == expected;
}

Fingerprint fingerprint (std::istream& file, const std::string& path)
{
    file.seekg(0);
    Fingerprint fp;
    std::array<char, 32768> block;
    do {
        file.read(block.data(), block.size());
        fp.add(block.data(), static_cast<std::size_t>(file.gcount()));
    } while (file);
    require(file.eof() && !file.bad(), "cannot fingerprint history: " + path);
    file.clear(); // allow a subsequent seek/write on the same stream
    return fp;
}

} // namespace

struct History::State {
    std::string name, signature, lineage;
    int cadence = -1, bins = 0, step = 0;
    bool clock_known = false;
    Real time = 0, dt = 0; // step zero: initial time if supplied by native Header
    Mass kg{};
    int row_step = 0;
    std::string row; // exact numeric text, without newline, even across segments
    int origin_step = 0; // accepted step at which THIS physical segment starts
    Fingerprint segment;
    fs::file_time_type verified_at{}; // append fast path only; never checkpointed

    std::string log_path () const { return name + "_wall_loss.txt"; }

    void check_restart_clock (const std::string& dir) {
        int native_step;
        double native_end;
        if (native_clock(dir, native_step, native_end)) {
            const Real native_time = static_cast<Real>(native_end);
            require(native_step == step, "accepted step differs from " + dir + "/Header");
            if (clock_known) {
                require(same_time(time+dt, native_time, dt), "accepted time differs from " + dir + "/Header");
            } else {
                time = native_time; // an initially unbound step-zero clock
                clock_known = true;
            }
        }
    }

    std::string encode () const {
        auto os = output();
        os << version << '\n' << std::quoted(name) << '\n' << std::quoted(signature) << '\n'
           << lineage << '\n' << cadence << ' ' << bins << '\n'
           << step << ' ' << clock_known << ' ' << time << ' ' << dt << '\n'
           << kg[0] << ' ' << kg[1] << ' ' << kg[2] << '\n'
           << row_step << ' ' << std::quoted(row) << '\n'
           << origin_step << ' ' << segment.bytes << ' ' << segment.hash << '\n';
        auto body = os.str();
        Fingerprint fp;
        fp.add(body.data(), body.size());
        return body + "checksum " + std::to_string(fp.hash) + '\n';
    }

    void decode (const std::string& text) {
        std::istringstream in(text);
        in.imbue(std::locale::classic());
        std::string title;
        require(std::getline(in, title) && title == version, "unsupported wall-loss sidecar version");
        require(bool(in >> std::quoted(name) >> std::quoted(signature) >> lineage >> cadence >> bins
                        >> step >> clock_known >> time >> dt >> kg[0] >> kg[1] >> kg[2]
                        >> row_step >> std::quoted(row) >> origin_step >> segment.bytes >> segment.hash),
                "malformed wall-loss sidecar");
        check_identity(name, signature);
        // Also checks mass/row corruption between samples or in a new directory
        // where there is no previous log available for comparison.
        require(encode() == text, "sidecar checksum, canonical encoding or length mismatch");
        require(lineage.size() == 32 && lineage.find_first_not_of("0123456789abcdef") == std::string::npos,
                "invalid wall-loss lineage");
        require(step >= 0 && origin_step >= 0 && origin_step <= step &&
                std::all_of(kg.begin(), kg.end(), nonnegative), "invalid accepted step or cumulative mass");
        if (step == 0) {
            require(cadence >= -1 && bins == 0 && kg == Mass{} && dt == 0 &&
                    std::isfinite(time) && (clock_known || time == 0), "invalid step-zero checkpoint");
        } else {
            require(cadence >= 0 && bins >= 1 && bins <= max_bins && clock_known && valid_interval(time, dt),
                    "invalid accepted interval, cadence or bin count");
        }
        require(row_step == (cadence > 0 ? step - step%cadence : 0), "last sampled step disagrees with cadence");
        require(segment.bytes ? cadence > 0 : (origin_step == 0 && row_step == 0 && segment.hash == basis),
                "invalid segment metadata");
        if (row_step == 0) {
            require(row.empty(), "unexpected row before first sample");
            return;
        }
        std::istringstream fields(row);
        fields.imbue(std::locale::classic());
        Real t, width;
        require(fields >> t >> width && valid_interval(t, width), "bad last-row time/dt");
        auto canonical = output();
        canonical << t << ' ' << width;
        if (row_step == step) {
            require(t == time && width == dt, "last row does not describe the accepted step");
        } else if (row_step == step-1) {
            require(same_time(t+width, time, std::min(width, dt)), "last row is not adjacent to accepted step");
        } else {
            require(t+width < time, "last row overlaps accepted step");
        }
        for (int i = 0; i < 6+bins; ++i) {
            Real x;
            require(fields >> x && nonnegative(x), "invalid last-row count, mass or rate");
            canonical << ' ' << x;
            if (i >= 3 && i < 6) {
                require(row_step == step ? x == kg[i-3] : x <= kg[i-3], "last-row cumulative mass mismatch");
            }
        }
        require(canonical.str() == row, "noncanonical or trailing last-row columns");
    }

    std::string prefix () const {
        auto os = output();
        os << "# ERF SD wall loss log v1\n# name " << std::quoted(name)
           << "\n# signature " << std::quoted(signature) << "\n# lineage " << lineage
           << "\n# cadence " << cadence << '\n';
        return os.str();
    }

    std::string header (bool continuation) const {
        auto os = output();
        os << prefix();
        if (continuation) {
            // Deterministic: a second restart from this very checkpoint may
            // find just this header, even though its parent segment was absent.
            os << "# continuation accepted_step " << step << " clock_known " << clock_known
               << " end_time " << time+dt << " cumulative_kg " << kg[0] << ' ' << kg[1] << ' ' << kg[2]
               << " bins " << bins << " parent_segment " << segment.bytes << ' ' << segment.hash
               << " last_row_step " << row_step << " last_row " << std::quoted(row) << '\n';
        } else {
            os << "# fresh\n";
        }
        os << "# step counts and bin rates are sampled, not interval means; mass includes all steps.\n"
           << "# groups: four side walls, top, bottom; measured whole-chamber losses are attributed by wall area.\n"
           << "# time_start_s dt_s side_number top_number bottom_number side_cumulative_kg"
              " top_cumulative_kg bottom_cumulative_kg bin_number_rate_per_s[0..num_bins-1]\n";
        return os.str();
    }

    // For ordinary appends always check size, identity and sampled tail. Only
    // rescan the middle on mtime change. Restart and checkpoint always rescan,
    // even on filesystems with coarse timestamps or after a timestamp-preserving copy.
    void verify (std::istream& file, const std::string& path, bool full) const {
        const auto size = fs::file_size(path);
        const std::string id = prefix();
        require(size == segment.bytes && size >= id.size() && matches(file, id, false, path),
                "stale, ahead or foreign history: " + path);
        if (row_step > origin_step) {
            const auto tail = '\n' + row + '\n';
            require(size >= tail.size() && matches(file, tail, true, path), "last sampled row differs: " + path);
        }
        if (full || row_step <= origin_step || fs::last_write_time(path) != verified_at) {
            const auto fp = fingerprint(file, path);
            require(fp.bytes == segment.bytes && fp.hash == segment.hash, "history content differs: " + path);
        }
    }

    void account_write (const std::string& path, const std::string& text) {
        segment.add(text.data(), text.size());
        require(fs::file_size(path) == segment.bytes, "history size changed during write: " + path);
        verified_at = fs::last_write_time(path); // after a checked close/flush
    }

    void append_row (const std::string& next_row) {
        const auto path = log_path();
        with_writer_lock(path, [&](const auto&) {
            if (!segment.bytes) {
                const auto text = header(false) + next_row + '\n';
                write_new(path, text);
                account_write(path, text);
                return;
            }
            require(regular_file_present(path), "missing history: " + path);
            // in|out cannot create or truncate. Verify, then seek/write on the
            // SAME open stream under the writer lock (no verify/reopen gap).
            std::fstream file(path, std::ios::in | std::ios::out | std::ios::binary);
            require(file.is_open(), "cannot open history for append: " + path);
            verify(file, path, false);
            file.seekp(0, std::ios::end);
            const auto text = next_row + '\n';
            file << text;
            file.close();
            require(bool(file), "append or close failed: " + path);
            account_write(path, text);
        });
    }

    void attach_restart_log () {
        const auto path = log_path();
        const auto continuation = header(true); // uses the PARENT fingerprint
        with_writer_lock(path, [&](const auto&) {
            if (regular_file_present(path)) {
                std::ifstream file(path, std::ios::binary);
                require(file.is_open(), "cannot open history: " + path);
                if (fs::file_size(path) != continuation.size() || !matches(file, continuation, false, path)) {
                    verify(file, path, true);
                    verified_at = fs::last_write_time(path);
                    return;
                }
                // Exact, still-empty continuation from this very checkpoint.
            } else {
                write_new(path, continuation);
            }
            origin_step = step;
            segment = {};
            account_write(path, continuation);
        });
    }
};

History::History () = default;
History::~History () = default;

History::State& History::bind (const std::string& name, const std::string& signature) const
{
    if (!m_state) {
        check_identity(name, signature);
        m_state = std::make_unique<State>();
        m_state->name = name;
        m_state->signature = signature;
        m_state->lineage = new_lineage();
    }
    require(m_state->name == name && m_state->signature == signature, "container or configuration changed");
    return *m_state;
}

void History::record (const std::string& name, const std::string& signature, int cadence,
                      int step, Real time, Real dt, const amrex::Vector<Real>& global)
{
    if (!amrex::ParallelDescriptor::IOProcessor()) { return; }
    try {
        auto& s = bind(name, signature);
        require(cadence >= 0 && (s.cadence == -1 || cadence == s.cadence), "invalid or changed cadence");
        require(step > 0 && step-1 == s.step && valid_interval(time, dt), "nonconsecutive step or invalid interval");
        if (s.clock_known) {
            require(same_time(time, s.time+s.dt, std::min(dt, s.dt)), "step start differs from the accepted clock");
        }
        require(global.size() >= 7 && global.size() <= 6+max_bins, "expected 6 + (1..64) global values");
        const int bins = static_cast<int>(global.size())-6;
        require((s.bins == 0 || s.bins == bins) && std::all_of(global.begin(), global.end(), nonnegative),
                "changed bin count or nonfinite/negative removal");
        if (cadence > 0 && s.step == 0 && !s.segment.bytes) {
            require(!regular_file_present(s.log_path()), "fresh run refuses existing " + s.log_path());
        }
        s.bins = bins;
        s.cadence = cadence;
        for (int i = 0; i < 3; ++i) {
            s.kg[i] += global[3+i];
            require(std::isfinite(s.kg[i]), "cumulative mass overflow");
        }
        if (cadence > 0 && step%cadence == 0) {
            auto os = output();
            os << time << ' ' << dt;
            for (int i = 0; i < 3; ++i) { os << ' ' << global[i]; }
            for (Real kg : s.kg) { os << ' ' << kg; }
            for (int i = 0; i < bins; ++i) {
                const Real rate = global[6+i]/dt;
                require(std::isfinite(rate), "bin rate overflow");
                os << ' ' << rate;
            }
            auto row = os.str();
            s.append_row(row);
            s.row = std::move(row);
            s.row_step = step;
        }
        s.step = step;
        s.time = time;
        s.dt = dt;
        s.clock_known = true;
    } catch (const std::exception& e) {
        amrex::Abort("sdwl::History::record(" + name + ", step " + std::to_string(step) + "): " + e.what());
    }
}

void History::checkpoint (const std::string& dir, const std::string& name, const std::string& signature) const
{
    if (!amrex::ParallelDescriptor::IOProcessor()) { return; }
    try {
        auto& s = bind(name, signature);
        if (s.cadence > 0) {
            const auto log = s.log_path();
            with_writer_lock(log, [&](const auto&) {
                require(regular_file_present(log) == (s.segment.bytes != 0), "missing or unexpected history: " + log);
                if (s.segment.bytes) {
                    std::ifstream file(log, std::ios::binary);
                    require(file.is_open(), "cannot open history: " + log);
                    s.verify(file, log, true);
                }
            });
        }
        const auto path = dir + '/' + name + '/' + sidecar;
        const auto text = s.encode();
        with_writer_lock(path, [&](const std::string& lock) {
            if (regular_file_present(path)) {
                require(read_sidecar(path, name, signature) == text, "refusing to replace different checkpoint sidecar: " + path);
            } else {
                // Publish only the complete, closed file on the same filesystem.
                const auto tmp = lock + "/pending";
                write_new(tmp, text);
                fs::rename(tmp, path);
            }
        });
    } catch (const std::exception& e) {
        amrex::Abort("sdwl::History::checkpoint(" + dir + ", " + name + "): " + e.what());
    }
}

void History::restart (const std::string& dir, const std::string& name, const std::string& signature, int cadence)
{
    if (!amrex::ParallelDescriptor::IOProcessor()) { return; }
    try {
        check_identity(name, signature);
        auto s = std::make_unique<State>();
        s->decode(read_sidecar(dir + '/' + name + '/' + sidecar, name, signature));
        require(s->name == name && s->signature == signature, "checkpoint name or configuration signature differs");
        require(cadence >= 0 && (s->cadence == -1 || s->cadence == cadence), "checkpoint cadence differs");
        s->cadence = cadence;
        s->check_restart_clock(dir);
        if (cadence > 0) { s->attach_restart_log(); }
        m_state = std::move(s);
    } catch (const std::exception& e) {
        amrex::Abort("sdwl::History::restart(" + dir + ", " + name + "): " + e.what());
    }
}

} // namespace sdwl
