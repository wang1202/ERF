// CPU-only, single reader rank. This oracle uses AMReX's native plot reader,
// not the wall-loss kernel, accumulators, sidecar or diagnostic parser.
#include <AMReX.H>
#include <AMReX_BoxIterator.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
namespace fs = std::filesystem;
constexpr int nx = 16, ny = 8, nz = 8, cells = nx*ny*nz, bins = 4;
constexpr double dx = 2./nx, dy = 2./ny, dz = 1./nz, volume = dx*dy*dz;
constexpr double dt = 1./8192, n0 = 32768.*volume, mdrop = 1.5e-12 + 0.5e-12;
constexpr bool single = sizeof(amrex::Real) == 4 || sizeof(amrex::ParticleReal) == 4;
constexpr double log_tol = single ? 6.e-5 : 2.e-8;
constexpr double mesh_tol = single ? 6.e-5 : 2.e-6; // also allow float-encoded plots
const std::string container = "super_droplets_moisture";
using Groups = std::array<double,3>;     // side, top, bottom
using Row = std::array<double,2+6+bins>; // step start, dt, counts, cumulative kg, bin rates

void require (bool condition, const std::string& message)
{
    if (!condition) { throw std::runtime_error(message); }
}

void near (double got, double want, double tol, const std::string& what)
{
    // No dimensional absolute floor: zero-loss faces/bins must be zero and a
    // missing ~1e-9 kg signal is just as much an error as a missing count.
    if (!std::isfinite(got) || !(std::abs(got-want) <= tol*std::abs(want))) {
        std::ostringstream msg;
        msg << std::setprecision(16) << what << ": got " << got << ", expected " << want;
        throw std::runtime_error(msg.str());
    }
}

std::vector<Row> read_history (const fs::path& dir, int first, int last, int cadence)
{
    const auto path = dir/(container + "_wall_loss.txt");
    if (cadence == 0) {
        require(!fs::exists(path), "unexpected " + path.string());
        return {};
    }
    std::ifstream in(path);
    require(in.good(), "cannot read " + path.string());
    std::vector<Row> rows;
    std::string header, line;
    while (std::getline(in, line)) {
        std::istringstream s(line);
        s >> std::ws;
        if (s.eof()) { continue; }
        if (s.peek() == '#') {
            require(rows.empty(), "unexpected new header after data in " + path.string());
            header += line + '\n';
            continue;
        }
        Row row{};
        for (auto& x : row) {
            require(bool(s >> x) && std::isfinite(x) && x >= 0, "bad/short history row: " + line);
        }
        s >> std::ws;
        require(s.eof(), "extra history columns: " + line);
        rows.push_back(row);
    }
    require(in.eof(), "history read error");
    for (const char* word : {"time", "dt", "side", "top", "bottom", "kg", "rate", "sample"}) {
        require(header.find(word) != std::string::npos, std::string("header omits ") + word);
    }
    require((header.find("# continuation") != std::string::npos) == (first > 0),
            "wrong fresh/continuation history header");
    require(int(rows.size()) == last/cadence - first/cadence, "wrong number of history rows");
    return rows;
}

enum class Model { Off, Constant, Measured, Periodic, MeasuredPeriodic, Diffusion, Stochastic,
                   FaceControl, Settling, ReducedRoof, Threshold };

Model model_for (const std::string& name)
{
    if (name == "off") { return Model::Off; }
    if (name == "constant") { return Model::Constant; }
    if (name == "measured") { return Model::Measured; }
    if (name == "constant_periodic") { return Model::Periodic; }
    if (name == "measured_periodic") { return Model::MeasuredPeriodic; }
    if (name == "face_control") { return Model::FaceControl; }
    if (name == "settling") { return Model::Settling; }
    if (name == "reduced_roof") { return Model::ReducedRoof; }
    if (name == "threshold") { return Model::Threshold; }
    if (name == "diffusion") { return Model::Diffusion; }
    if (name == "stochastic") { return Model::Stochastic; }
    throw std::runtime_error("unknown model: " + name);
}

bool mixed_faces (Model model) {
    return model == Model::FaceControl || model == Model::Settling || model == Model::ReducedRoof;
}
double drop_mass (Model model) { return mixed_faces(model) ? 2.e-7 + 0.5e-12 : mdrop; }

Groups rates (Model model, int i, int j, int k)
{
    if (model == Model::Off) { return {}; }
    if (model == Model::Measured || model == Model::Stochastic) {
        // All cells, including the center. Chamber wall areas are 8:4:4.
        return {512., 256., 256.};
    }
    if (model == Model::MeasuredPeriodic) { return {0., 512., 512.}; }
    double side = 128., top = 192., bottom = 64.;
    // NaCl is soluble: water and solute contribute to the wet radius.
    const double pi = std::acos(-1.);
    const double radius = std::cbrt(drop_mass(model)/(4.*pi/3.*1000.));
    if (mixed_faces(model)) {
        // Atlas-Ulbrich uses diameter in mm. At this radius vt ~ 3.05 m/s:
        // the 2 m/s roof clamps, the 8 m/s roof still deposits, the floor gains vt.
        const double vt = model == Model::FaceControl ? 0. : 3.778*std::pow(2000.*radius, 0.67);
        side = 8.; top = std::max((model == Model::ReducedRoof ? 8. : 2.)-vt, 0.); bottom = 4.+vt;
    } else if (model == Model::Diffusion) {
        side = top = bottom = 1.380649e-23*293./(6.*pi*1.8e-8*radius*1.e-11);
    }
    if (model == Model::Periodic) { side = 0.; }
    // Each stationary cell retains its own population. Add both wall rates
    // at edges (three at a corner); interior cells have exactly zero loss.
    const bool xwall = i == nx-1 || (i == 0 && !mixed_faces(model)); // mixed xlo is Open
    return {side*(xwall/dx + (j == 0 || j == ny-1)/dy),
            top*(k == nz-1)/dz, bottom*(k == 0)/dz};
}

struct Snapshot {
    std::array<double,cells> n{}, mass{}; // actual represented count and kg, NOT SD count
    double number () const { return std::accumulate(n.begin(), n.end(), 0.); }
    double kg () const { return std::accumulate(mass.begin(), mass.end(), 0.); }
};

Snapshot read_plot (const fs::path& dir, int step, bool check_sd = false)
{
    std::ostringstream label;
    label << "plt" << std::setfill('0') << std::setw(5) << step;
    amrex::PlotFileData pf((dir/label.str()).string());
    require(pf.finestLevel() == 0 && pf.levelStep(0) == step, "wrong native plot level/step");
    near(pf.time(), step*dt, mesh_tol, "plot time");
    const auto domain = pf.probDomain(0);
    require(domain.smallEnd() == amrex::IntVect(0) && domain.bigEnd() == amrex::IntVect(nx-1,ny-1,nz-1),
            "wrong fixture domain");
    const auto spacing = pf.cellSize(0);
    require(spacing[0] == dx && spacing[1] == dy && spacing[2] == dz, "wrong fixture spacing");
    const auto& names = pf.varNames();
    auto component = [&](const std::string& name) {
        const auto it = std::find(names.begin(), names.end(), container + '_' + name);
        require(it != names.end(), "native plot lacks " + name);
        return static_cast<int>(it-names.begin());
    };
    const int num = component("number_density"), mass = component("mass_density");
    const int sd = check_sd ? component("sd_number_density") : -1;
    const auto fields = pf.get(0);
    Snapshot snap;
    int count = 0;
    for (amrex::MFIter mfi(fields); mfi.isValid(); ++mfi) {
        for (amrex::BoxIterator b(mfi.validbox()); b.ok(); ++b) {
            const auto& iv = b();
            const int c = iv[0] + nx*(iv[1] + ny*iv[2]);
            snap.n[c] = volume*fields[mfi](iv, num);
            snap.mass[c] = volume*fields[mfi](iv, mass);
            if (check_sd) { // distinguish terminal removal from an active zero-weight SD
                near(volume*fields[mfi](iv, sd), snap.n[c] > 0 ? 1. : 0., mesh_tol,
                     label.str() + ": active SD count, cell " + std::to_string(c));
            }
            ++count;
        }
    }
    require(count == cells, "checker did not read every valid cell (launch checker with one rank)");
    return snap;
}

// Verify the native mesh state *and* compute an independent reference for
// each of the six diagnostic quantities without evaluating production helpers.
void check_step (Model model, int step, const Snapshot& now, const Snapshot* before, const Row* row)
{
    Groups sample{}, kg{};
    const double mass = drop_mass(model);
    const std::string context = "step " + std::to_string(step) + ": ";
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                const int c = i + nx*(j + ny*k);
                const auto rate = rates(model, i, j, k);
                const double total = rate[0] + rate[1] + rate[2];
                auto population = [&](int s) {
                    const double left = n0*std::exp(-total*s*dt);
                    return model == Model::Threshold && left < 103. ? 0. : left;
                };
                const double surviving = population(step);
                const double n = now.n[c];
                if (model != Model::Stochastic || step == 0) {
                    near(n, surviving, mesh_tol, context + "represented N, cell " + std::to_string(c));
                } else {
                    require(before && n >= 0 && n <= before->n[c], context + "unbounded stochastic thinning");
                }
                near(now.mass[c], n*mass, mesh_tol, context + "total particle mass, cell " + std::to_string(c));
                if (!row || total == 0) { continue; }
                // A threshold crossing deposits the entire previous population,
                // including the nonzero residual, once and in the same face shares.
                const double lost = model == Model::Stochastic ? before->n[c]-n : population(step-1)-surviving;
                const double lost_to_date = n0 - (model == Model::Stochastic ? n : surviving);
                for (int g = 0; g < 3; ++g) {
                    sample[g] += lost*rate[g]/total;
                    kg[g] += lost_to_date*mass*rate[g]/total;
                }
            }
        }
    }
    if (model == Model::Stochastic && step > 0) {
        const auto bounds = std::minmax_element(now.n.begin(), now.n.end());
        require(*bounds.second - *bounds.first > 1., context + "no stochastic spread across identical cells");
        const double survival = std::exp(-1024.*step*dt);
        // Smoke bound on the ensemble, not a prescribed seedwise trajectory.
        require(std::abs(now.number()-cells*n0*survival) < 10.*std::sqrt(cells*n0*(1.-survival)),
                context + "stochastic mean inconsistent with uniform thinning");
    }
    if (!row) { return; }
    require(before != nullptr, "sampled step has no preceding particle state");
    const auto& r = *row;
    near(r[0], (step-1)*dt, log_tol, context + "row time (start of sampled step)");
    near(r[1], dt, log_tol, context + "row dt (not output interval)");
    const double tol = model == Model::Stochastic ? mesh_tol : log_tol;
    for (int g = 0; g < 3; ++g) {
        near(r[2+g], sample[g], tol, context + "sampled count, group " + std::to_string(g));
        near(r[5+g], kg[g], tol, context + "cumulative kg, group " + std::to_string(g));
    }
    const double removed = r[2]+r[3]+r[4];
    near(removed, before->number()-now.number(), mesh_tol, context + "mesh number loss in sampled step");
    // For weak/zero loss don't subtract almost identical float-encoded mesh
    // masses. Each deposited group is independently checked at log precision above.
    const double deposited = r[5]+r[6]+r[7], initial = cells*n0*mass;
    if (model == Model::Off || mixed_faces(model)) {
        near(now.kg()+deposited, initial, mesh_tol, context + "all-step suspended + deposited mass");
    } else {
        near(deposited, initial-now.kg(), mesh_tol, context + "all-step mesh mass loss");
    }
    // Wet radii: 7.8 um in bin 1 (not the diameter's bin 2); 363 um in bin 3.
    const int wet_bin = mixed_faces(model) ? 3 : 1;
    for (int b = 0; b < bins; ++b) {
        near(r[8+b], b == wet_bin ? removed/dt : 0., log_tol, context + "#/s in bin " + std::to_string(b));
    }
}

void check (const fs::path& dir, Model model, int first, int last, int cadence, const fs::path& reference)
{
    auto rows = read_history(dir, first, last, cadence);
    std::vector<Row> straight;
    if (!reference.empty()) { straight = read_history(reference, 0, last, cadence); }
    Snapshot prev;
    int sampled = 0;
    for (int step = first; step <= last; ++step) {
        const auto now = read_plot(dir, step, model == Model::Threshold);
        const bool emit = cadence && step > first && step%cadence == 0;
        check_step(model, step, now, step == first ? nullptr : &prev, emit ? &rows.at(sampled) : nullptr);
        if (!reference.empty()) {
            const auto other = read_plot(reference, step);
            for (int c = 0; c < cells; ++c) {
                near(now.n[c], other.n[c], mesh_tol, "straight-run N at step " + std::to_string(step));
                near(now.mass[c], other.mass[c], mesh_tol, "straight-run mass at step " + std::to_string(step));
            }
            if (emit) {
                const auto& want = straight.at(step/cadence-1);
                for (std::size_t col = 0; col < want.size(); ++col) {
                    near(rows.at(sampled)[col], want[col], log_tol, "straight-run history column " + std::to_string(col));
                }
            }
        }
        if (emit) { ++sampled; }
        prev = now;
    }
}
} // namespace

int main (int argc, char** argv)
{
    if (argc != 6 && argc != 7) {
        std::cerr << "usage: checker directory model first_step last_step cadence [straight_dir]\n";
        return 2;
    }
    amrex::Initialize(argc, argv, false);
    int rc = 0;
    try {
        check(argv[1], model_for(argv[2]), std::stoi(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]),
              argc == 7 ? argv[6] : "");
    } catch (const std::exception& e) {
        std::cerr << "SDWallLossCheck: " << argv[1] << ": " << e.what() << '\n';
        rc = 1;
    }
    amrex::Finalize();
    return rc;
}
