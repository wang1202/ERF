// Checker for the interval-mean Cloud Chamber wall-flux output
// (erf.wallflux_diagnostics, wallflux_lev0.txt).  See
// Tests/RunCloudChamberWallFlux.cmake for the runs it inspects.
//
//   dry    <run_dir> <one_box_dir> <disabled_dir> <plotfile> <rows> <row_dt>
//   satadj <run_dir> <plotfile> <rows> <row_dt> <dry_faces>
//
// <dry_faces> lists the dry walls of the SatAdj deck (e.g. "xhi,yhi").
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Reduce.H>
#include "ERF_Constants.H"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

namespace {

using amrex::MultiFab;
using amrex::PlotFileData;
using amrex::Real;

constexpr int NumFaces = 6;
constexpr int NumColumns = 19;
const char* const face_names[NumFaces] = {"xlo", "xhi", "ylo", "yhi", "zlo", "zhi"};

int column (int face, int channel) { return 1 + 3 * face + channel; }

int fail (const std::string& message)
{
    std::cerr << "CloudChamberWallFluxCheck: " << message << "\n";
    return 1;
}

// The comment header of wang1202/ERF_01272026 ERF_WriteWallFluxes.cpp, with
// the comment marker first: an independent statement of the output contract.
std::string expected_header ()
{
    std::ostringstream out;
    out << std::left << std::setw(15) << "# time[s]" << std::right;
    for (const char* face : face_names) {
        out << std::setw(18) << (std::string(face) + "_tau[Pa]")
            << std::setw(18) << (std::string(face) + "_H[W/m2]")
            << std::setw(18) << (std::string(face) + "_E[kg/m2/s]");
    }
    return out.str();
}

struct WallFluxFile {
    std::vector<std::array<double, NumColumns>> rows;
};

bool read_wallflux (const std::string& path, WallFluxFile& file, std::string& error)
{
    std::ifstream in(path);
    if (!in) { error = "cannot open " + path; return false; }
    std::string line;
    if (!std::getline(in, line)) { error = path + " is empty"; return false; }
    if (line != expected_header()) {
        error = path + ": header differs from the legacy layout:\n  got      '" + line +
                "'\n  expected '" + expected_header() + "'";
        return false;
    }
    const std::regex number(R"(-?[0-9]\.[0-9]{8}e[+-][0-9]{2,3})");
    int line_number = 1;
    while (std::getline(in, line)) {
        ++line_number;
        if (line.find('#') != std::string::npos) {
            error = path + ":" + std::to_string(line_number) + ": unexpected comment";
            return false;
        }
        std::istringstream tokens(line);
        std::array<double, NumColumns> row{};
        std::string token;
        int n = 0;
        while (tokens >> token) {
            if (n >= NumColumns) {
                error = path + ":" + std::to_string(line_number) + ": more than 19 columns";
                return false;
            }
            if (!std::regex_match(token, number)) {
                error = path + ":" + std::to_string(line_number) + ": '" + token +
                        "' is not scientific notation with eight decimals";
                return false;
            }
            row[n++] = std::stod(token);
        }
        if (n != NumColumns) {
            error = path + ":" + std::to_string(line_number) + ": " + std::to_string(n) +
                    " columns instead of 19";
            return false;
        }
        file.rows.push_back(row);
    }
    return true;
}

bool check_rows (const WallFluxFile& file, int rows, double row_dt, std::string& error)
{
    if (static_cast<int>(file.rows.size()) != rows) {
        error = "expected " + std::to_string(rows) + " rows, found " +
                std::to_string(file.rows.size());
        return false;
    }
    for (int r = 0; r < rows; ++r) {
        const double expected = (r + 1) * row_dt;
        if (std::abs(file.rows[r][0] - expected) > 1.0e-9 * std::max(1.0, expected)) {
            error = "row " + std::to_string(r) + " time " + std::to_string(file.rows[r][0]) +
                    " is not the window end " + std::to_string(expected);
            return false;
        }
        for (int f = 0; f < NumFaces; ++f) {
            const double tau = file.rows[r][column(f, 0)];
            if (!(tau >= 0.0)) {
                error = std::string(face_names[f]) + "_tau is negative or not finite";
                return false;
            }
        }
    }
    return true;
}

struct BudgetRow {
    double start_time = 0.0;
    double end_time = 0.0;
    std::string scalar;
    std::array<double, NumFaces> faces{};
};

bool read_budget (const std::string& path, std::vector<BudgetRow>& rows, std::string& error)
{
    std::ifstream in(path);
    if (!in) { error = "cannot open " + path; return false; }
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') { continue; }
        std::istringstream tokens(line);
        BudgetRow row;
        int start_step = 0;
        int end_step = 0;
        std::string units;
        tokens >> start_step >> row.start_time >> end_step >> row.end_time >> row.scalar >> units;
        for (auto& face : row.faces) { tokens >> face; }
        if (!tokens) { error = "malformed budget row: " + line; return false; }
        rows.push_back(row);
    }
    return true;
}

const BudgetRow* find_budget (const std::vector<BudgetRow>& rows,
                              const std::string& scalar, double end_time)
{
    for (const auto& row : rows) {
        if (row.scalar == scalar &&
            std::abs(row.end_time - end_time) <= 1.0e-9 * std::max(1.0, end_time)) {
            return &row;
        }
    }
    return nullptr;
}

std::array<double, NumFaces> face_areas (const PlotFileData& plotfile)
{
    const auto lo = plotfile.probLo();
    const auto hi = plotfile.probHi();
    const double lx = hi[0] - lo[0];
    const double ly = hi[1] - lo[1];
    const double lz = hi[2] - lo[2];
    return {ly * lz, ly * lz, lx * lz, lx * lz, lx * ly, lx * ly};
}

// Min and max of a plotfile variable over the cells of @p region.
std::pair<Real, Real> min_max (const MultiFab& field, const amrex::Box& region)
{
    amrex::ReduceOps<amrex::ReduceOpMin, amrex::ReduceOpMax> ops;
    amrex::ReduceData<Real, Real> data(ops);
    for (amrex::MFIter mfi(field); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox() & region;
        if (bx.isEmpty()) { continue; }
        const auto arr = field.const_array(mfi);
        ops.eval(bx, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            -> amrex::GpuTuple<Real, Real> { return {arr(i,j,k), arr(i,j,k)}; });
    }
    const auto result = data.value(ops);
    Real lo = amrex::get<0>(result);
    Real hi = amrex::get<1>(result);
    amrex::ParallelDescriptor::ReduceRealMin(lo);
    amrex::ParallelDescriptor::ReduceRealMax(hi);
    return {lo, hi};
}

bool close (double value, double reference, double rel, double abs_floor)
{
    return std::abs(value - reference) <= rel * std::abs(reference) + abs_floor;
}

int check_dry (const std::string& run_dir, const std::string& one_box_dir,
               const std::string& disabled_dir, const std::string& plotfile,
               int rows, double row_dt)
{
    std::string error;
    WallFluxFile multi;
    WallFluxFile one_box;
    if (!read_wallflux(run_dir + "/wallflux_lev0.txt", multi, error) ||
        !check_rows(multi, rows, row_dt, error) ||
        !read_wallflux(one_box_dir + "/wallflux_lev0.txt", one_box, error) ||
        !check_rows(one_box, rows, row_dt, error)) {
        return fail(error);
    }

    // Dry walls exchange no vapor, and the dry deck transports none.
    for (const auto& row : multi.rows) {
        for (int f = 0; f < NumFaces; ++f) {
            if (row[column(f, 2)] != 0.0) {
                return fail(std::string(face_names[f]) + "_E is not exactly zero in a dry run");
            }
        }
    }
    // Sheared walls: the uniform initial x-flow drags on the bulk zlo wall and
    // on the resolved no-slip y and z walls.
    for (int f : {2, 3, 4, 5}) {
        for (const auto& row : multi.rows) {
            if (!(row[column(f, 0)] > 0.0)) {
                return fail(std::string(face_names[f]) + "_tau is not positive under shear");
            }
        }
    }

    // Decomposition independence: one box on one rank against many boxes on
    // many ranks, at the printed precision.
    for (int channel = 0; channel < 3; ++channel) {
        double scale = 0.0;
        for (const auto& row : multi.rows) {
            for (int f = 0; f < NumFaces; ++f) {
                scale = std::max(scale, std::abs(row[column(f, channel)]));
            }
        }
        for (int r = 0; r < rows; ++r) {
            for (int f = 0; f < NumFaces; ++f) {
                const double a = multi.rows[r][column(f, channel)];
                const double b = one_box.rows[r][column(f, channel)];
                if (!close(a, b, 1.0e-7, 1.0e-9 * scale)) {
                    std::ostringstream message;
                    message << "decomposition mismatch row " << r << " column "
                            << column(f, channel) << ": " << a << " vs " << b;
                    return fail(message.str());
                }
            }
        }
    }

    // Enabling the diagnostic leaves the simulated fields unchanged.
    PlotFileData enabled(run_dir + "/" + plotfile);
    PlotFileData disabled(disabled_dir + "/" + plotfile);
    if (enabled.varNames() != disabled.varNames()) {
        return fail("enabled and disabled plotfiles hold different variables");
    }
    double max_relative = 0.0;
    for (const auto& name : enabled.varNames()) {
        MultiFab difference = enabled.get(0, name);
        const MultiFab other = disabled.get(0, name);
        const Real scale = std::max(Real(1.0), other.norm0(0, 0, false));
        MultiFab::Subtract(difference, other, 0, 0, 1, 0);
        const double relative = difference.norm0(0, 0, false) / scale;
        max_relative = std::max(max_relative, relative);
        if (!(relative <= 1.0e-10)) {
            return fail("diagnostic changed field " + name + ": max relative difference " +
                        std::to_string(relative));
        }
    }
    std::cout << "dry: rows=" << rows << " enabled_disabled_max_relative_difference="
              << max_relative << "\n";
    return 0;
}

int check_satadj (const std::string& run_dir, const std::string& plotfile_name,
                  int rows, double row_dt, const std::string& dry_faces)
{
    std::string error;
    WallFluxFile file;
    if (!read_wallflux(run_dir + "/wallflux_lev0.txt", file, error) ||
        !check_rows(file, rows, row_dt, error)) {
        return fail(error);
    }
    std::vector<BudgetRow> budget;
    if (!read_budget(run_dir + "/cloud_chamber_budget.dat", budget, error)) {
        return fail(error);
    }

    PlotFileData plotfile(run_dir + "/" + plotfile_name);
    const auto area = face_areas(plotfile);
    const MultiFab p_hse = plotfile.get(0, "pres_hse");
    const amrex::Box domain = plotfile.probDomain(0);
    amrex::Box bottom = domain; bottom.setBig(2, domain.smallEnd(2));
    amrex::Box top = domain; top.setSmall(2, domain.bigEnd(2));
    const auto p_bottom = min_max(p_hse, bottom);
    const auto p_top = min_max(p_hse, top);
    const auto exner = [] (Real p) { return std::pow(static_cast<double>(p) / p_0, R_d / Cp_d); };
    if (p_bottom.first != p_bottom.second || p_top.first != p_top.second) {
        return fail("hydrostatic pressure is not uniform on the bottom or top cell layer");
    }
    const double exner_face[NumFaces] = {0.0, 0.0, 0.0, 0.0,
                                         exner(p_bottom.first), exner(p_top.first)};

    double previous_time = 0.0;
    int wet_faces_seen = 0;
    for (int r = 0; r < rows; ++r) {
        const auto& row = file.rows[r];
        const BudgetRow* vapor = find_budget(budget, "water_vapor", row[0]);
        const BudgetRow* theta = find_budget(budget, "rhoTheta", row[0]);
        if (!vapor || !theta) {
            return fail("no budget window ends at " + std::to_string(row[0]));
        }
        const double elapsed = vapor->end_time - vapor->start_time;
        if (std::abs((row[0] - previous_time) - elapsed) > 1.0e-9 * std::max(1.0, elapsed)) {
            return fail("wall-flux and budget windows differ");
        }
        previous_time = row[0];

        double vapor_scale = 0.0;
        double theta_scale = 0.0;
        for (int f = 0; f < NumFaces; ++f) {
            vapor_scale = std::max(vapor_scale, std::abs(vapor->faces[f]));
            theta_scale = std::max(theta_scale, std::abs(theta->faces[f]));
        }
        for (int f = 0; f < NumFaces; ++f) {
            const std::string face = face_names[f];
            const bool dry = dry_faces.find(face) != std::string::npos;
            // E: the native vapor face integral over face area and window time,
            // in coordinate signs.
            const double E = row[column(f, 2)];
            const double E_native = vapor->faces[f] / (area[f] * elapsed);
            if (dry) {
                if (E != 0.0 || vapor->faces[f] != 0.0) {
                    return fail("dry wall " + face + " exchanges vapor");
                }
            } else {
                if (E == 0.0) { return fail("wet wall " + face + " exchanges no vapor"); }
                ++wet_faces_seen;
            }
            if (!close(E, E_native, 2.0e-8, 1.0e-12 * vapor_scale / (area[f] * elapsed))) {
                std::ostringstream message;
                message << face << "_E=" << std::setprecision(17) << E
                        << " differs from the budget's " << E_native;
                return fail(message.str());
            }

            // H: Cp_d * Exner * F(rho theta), converted per patch.  The
            // horizontal walls see one hydrostatic pressure, so there the
            // window mean must equal Cp_d * Exner times the budget's mean flux.
            // (On a sidewall Exner varies with height and the flux may change
            // sign along it, so no such identity or bound holds there.)
            if (f >= 4) {
                const double H = row[column(f, 1)];
                const double F_theta = theta->faces[f] / (area[f] * elapsed);
                const double H_expected = Cp_d * exner_face[f] * F_theta;
                if (!close(H, H_expected, 2.0e-8, 1.0e-12 * Cp_d * theta_scale / (area[f] * elapsed))) {
                    std::ostringstream message;
                    message << face << "_H=" << std::setprecision(17) << H
                            << " differs from Cp_d*Exner*F=" << H_expected;
                    return fail(message.str());
                }
                if (H == 0.0) { return fail(face + "_H is zero on a heated wall"); }
            }
        }
        std::cout << "satadj: t=" << row[0] << " elapsed=" << elapsed;
        for (int f = 0; f < NumFaces; ++f) {
            std::cout << " " << face_names[f] << "_E=" << row[column(f, 2)];
        }
        std::cout << "\n";
    }
    if (wet_faces_seen == 0) { return fail("no wet wall exchanged vapor"); }
    return 0;
}

} // namespace

int main (int argc, char** argv)
{
    amrex::Initialize(argc, argv, false);
    int status = 1;
    const std::string mode = argc > 1 ? argv[1] : std::string();
    if (mode == "dry" && argc == 8) {
        status = check_dry(argv[2], argv[3], argv[4], argv[5],
                           std::stoi(argv[6]), std::stod(argv[7]));
    } else if (mode == "satadj" && argc == 7) {
        status = check_satadj(argv[2], argv[3], std::stoi(argv[4]),
                              std::stod(argv[5]), argv[6]);
    } else {
        status = fail("usage: dry <run> <one_box> <disabled> <plotfile> <rows> <row_dt> | "
                      "satadj <run> <plotfile> <rows> <row_dt> <dry_faces>");
    }
    amrex::Finalize();
    return status;
}
