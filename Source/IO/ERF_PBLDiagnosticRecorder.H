#ifndef ERF_PBL_DIAGNOSTIC_RECORDER_H_
#define ERF_PBL_DIAGNOSTIC_RECORDER_H_

#include "ERF_InstrumentSimulators.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include <fstream>
#include <memory>
#include <string>

struct MoistureComponentIndices;
class ShocDriver;

class PBLDiagnosticRecorder
{
public:
    static constexpr amrex::Real missing_value = amrex::Real(-999.0);

    static std::string make_filename (const std::string& prefix, int site,
                                      bool profile);
    void init (const amrex::Geometry& geom, int max_level, bool restarting,
               const ShocDriver& driver);
    int level () const noexcept { return m_level; }
    bool summary_step (int step) const noexcept;
    bool profile_step (int step) const noexcept;

    void write (amrex::Real time, int step, const amrex::MultiFab& cons,
                const amrex::MultiFab& xvel, const amrex::MultiFab& yvel,
                const amrex::MultiFab& zvel, const amrex::MultiFab& z_phys_cc,
                const amrex::Geometry& geom, const ShocDriver& driver,
                const MoistureComponentIndices& moisture_indices,
                const amrex::MultiFab* qheating_rates);

private:
    int m_level = 0;
    int m_summary_interval = 1;
    int m_profile_interval = 60;
    bool m_write_profiles = true;
    bool m_restarting = false;
    std::string m_output_prefix = "pbl_diag_";
    amrex::Vector<int> m_iloc, m_jloc;
    InstrumentSimUtil::MoveSchedule m_move;
    amrex::Vector<std::unique_ptr<std::fstream>> m_summary_files, m_profile_files;
    amrex::Vector<amrex::Real> m_last_summary_time, m_last_profile_time;

    static amrex::Real read_last_time (const std::string& filename);
    static bool file_is_empty (const std::string& filename);
    void open_files (bool restarting);
    void write_summary_header (std::fstream& os, int site) const;
    void write_profile_header (std::fstream& os, int site) const;
};

#endif
