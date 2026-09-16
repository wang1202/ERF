#include "ERF_PBLDiagnosticRecorder.H"

#include "ERF_Constants.H"
#include "ERF_IndexDefines.H"
#include "ERF_ShocDriver.H"
#include "ERF_SurfaceFluxDiagnostics.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace {
constexpr amrex::Real missing = PBLDiagnosticRecorder::missing_value;

bool valid (amrex::Real x) { return std::isfinite(x) && x != missing; }
amrex::Real safe (amrex::Real x) { return valid(x) ? x : missing; }

amrex::Vector<amrex::Real> missing_column (int klo, int khi)
{ return amrex::Vector<amrex::Real>(khi-klo+1, missing); }

amrex::Vector<amrex::Real> column_or_missing (const amrex::MultiFab* mf, int comp,
                                               int i, int j, int klo, int khi)
{
    if (mf == nullptr || !mf->isDefined() || comp < 0 || comp >= mf->nComp()) {
        return missing_column(klo, khi);
    }
    return InstrumentSimUtil::extract_column(*mf, comp, i, j, klo, khi);
}

int nearest_level (const amrex::Vector<amrex::Real>& z, amrex::Real h)
{
    int best = -1;
    amrex::Real dmin = std::numeric_limits<amrex::Real>::max();
    for (int k = 0; k < static_cast<int>(z.size()); ++k) {
        if (valid(z[k]) && std::abs(z[k]-h) < dmin) { best = k; dmin = std::abs(z[k]-h); }
    }
    return best;
}

amrex::Real layer_dz (const amrex::Vector<amrex::Real>& z, int k)
{
    if (z.size() < 2) return missing;
    return k+1 < static_cast<int>(z.size()) ? z[k+1]-z[k] : z[k]-z[k-1];
}

void write_value (std::ostream& os, amrex::Real x)
{ os << std::setw(16) << std::setprecision(8) << std::scientific << safe(x); }
}

std::string PBLDiagnosticRecorder::make_filename (const std::string& prefix, int site,
                                                   bool profile)
{
    std::ostringstream os;
    os << prefix << std::setw(2) << std::setfill('0') << site
       << (profile ? "_profile.txt" : "_summary.txt");
    return os.str();
}

bool PBLDiagnosticRecorder::file_is_empty (const std::string& filename)
{
    std::ifstream is(filename, std::ios::binary);
    return !is.good() || is.peek() == std::ifstream::traits_type::eof();
}

amrex::Real PBLDiagnosticRecorder::read_last_time (const std::string& filename)
{
    std::ifstream is(filename);
    amrex::Real last = -std::numeric_limits<amrex::Real>::max();
    std::string line;
    while (std::getline(is, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream row(line);
        amrex::Real t;
        if (row >> t) last = t;
    }
    return last;
}

void PBLDiagnosticRecorder::open_files (bool restarting)
{
    const int nsites = static_cast<int>(m_iloc.size());
    m_summary_files.resize(nsites);
    m_profile_files.resize(nsites);
    m_last_summary_time.assign(nsites, -std::numeric_limits<amrex::Real>::max());
    m_last_profile_time.assign(nsites, -std::numeric_limits<amrex::Real>::max());
    if (!amrex::ParallelDescriptor::IOProcessor()) return;

    for (int s = 0; s < nsites; ++s) {
        const std::string sf = make_filename(m_output_prefix, s, false);
        const std::string pf = make_filename(m_output_prefix, s, true);
        const bool se = file_is_empty(sf), pe = file_is_empty(pf);
        const auto mode = std::ios::out | (restarting ? std::ios::app : std::ios::trunc);
        m_summary_files[s] = std::make_unique<std::fstream>(sf, mode);
        m_profile_files[s] = std::make_unique<std::fstream>(pf, mode);
        if (!m_summary_files[s]->good() || !m_profile_files[s]->good()) amrex::FileOpenFailed(sf.c_str());
        if (restarting && !se) m_last_summary_time[s] = read_last_time(sf);
        if (restarting && !pe) m_last_profile_time[s] = read_last_time(pf);
        if (!restarting || se) write_summary_header(*m_summary_files[s], s);
        if (!restarting || pe) write_profile_header(*m_profile_files[s], s);
    }
}

void PBLDiagnosticRecorder::init (const amrex::Geometry& geom, int max_level,
                                  bool restarting, const ShocDriver& driver)
{
    amrex::ParmParse pp("erf.pbl_recorder");
    pp.query("level", m_level);
    pp.getarr("i_loc", m_iloc);
    pp.getarr("j_loc", m_jloc);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!m_iloc.empty(), "erf.pbl_recorder requires at least one site");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_iloc.size() == m_jloc.size(), "erf.pbl_recorder i_loc and j_loc must match");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_level >= 0 && m_level <= max_level, "erf.pbl_recorder.level is outside the AMR hierarchy");
    for (int s = 0; s < static_cast<int>(m_iloc.size()); ++s) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            geom.Domain().contains(amrex::IntVect(m_iloc[s], m_jloc[s], geom.Domain().smallEnd(2))),
            "erf.pbl_recorder site is outside the selected level domain");
    }
    pp.query("summary_output_interval", m_summary_interval);
    pp.query("write_profiles", m_write_profiles);
    pp.query("profile_output_interval", m_profile_interval);
    pp.query("output_file", m_output_prefix);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_summary_interval > 0, "erf.pbl_recorder.summary_output_interval must be positive");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_profile_interval > 0, "erf.pbl_recorder.profile_output_interval must be positive");
    m_move.init(pp, "erf.pbl_recorder");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(driver.entrainment_diagnostics_enabled(),
        "erf.pbl_recorder requires erf.shoc.diagnose_entrainment = true");
    m_restarting = restarting;
    open_files(restarting);
    amrex::ParallelDescriptor::Barrier("PBLDiagnosticRecorder::init");
}

bool PBLDiagnosticRecorder::summary_step (int step) const noexcept
{ return step >= 0 && step % m_summary_interval == 0; }
bool PBLDiagnosticRecorder::profile_step (int step) const noexcept
{ return m_write_profiles && step >= 0 && step % m_profile_interval == 0; }

void PBLDiagnosticRecorder::write_summary_header (std::fstream& os, int site) const
{
    os << "# ERF native-SHOC moving-column PBL diagnostic recorder\n"
       << "# site=" << site << " initial_cell_i=" << m_iloc[site]
       << " initial_cell_j=" << m_jloc[site] << " moving=" << (m_move.moving()?1:0)
       << " restart_append=" << (m_restarting?1:0) << "\n"
       << "# missing_value=-999; heights=z_phys_cc; provenance=ERF build metadata\n"
       << "# jump_layers=erf.shoc.entrainment_jump_ncell; min_delta_theta_v=erf.shoc.entrainment_min_delta_theta_v\n"
       << "# columns: time[s] step x[m] y[m] i_loc j_loc pblh_agl[m] pbl_top_msl[m] we_kinematic[m/s] we_flux_jump[m/s] pblh_tendency[m/s] pblh_hadv[m/s] w_at_pblh[m/s] delta_theta_v[K] wthv_at_pblh[K m/s] entrainment_flux_ratio[-] wthv_sfc[K m/s] u_star[m/s] Olen[m] convective_velocity_scale[m/s] sensible_heat_flux[W/m2] latent_heat_flux[W/m2] tke_mean_pbl[m2/s2] shear_prod_int[m3/s3] buoy_prod_int[m3/s3] diss_tke_int[m3/s3] tke_net_int[m3/s3] brunt_max_inversion[s-2] cloud_fraction_max_pbl[-] lwp[kg/m2] iwp[kg/m2] qsrc_lw_min_near_top[K/s] qsrc_sw_mean_pbl[K/s]\n";
}

void PBLDiagnosticRecorder::write_profile_header (std::fstream& os, int site) const
{
    os << "# ERF native-SHOC profile; site=" << site << " initial_cell_i=" << m_iloc[site]
       << " initial_cell_j=" << m_jloc[site] << " moving=" << (m_move.moving()?1:0)
       << " missing_value=-999; z_msl=z_phys_cc; optional moisture fields use the active map\n"
       << "# provenance=ERF build metadata\n"
       << "# columns: time[s] step x[m] y[m] i_loc j_loc pblh_agl[m] we_kinematic[m/s] we_flux_jump[m/s] k z_msl[m] z_agl[m] theta[K] theta_v[K] qv[kg/kg] qc[kg/kg] qi[kg/kg] u[m/s] v[m/s] w[m/s] KE[m2/s2] wthv_sec[K m/s] wthl_sec[K m/s] wqw_sec[kg/kg m/s] w_sec[m/s] w3[m3/s3] brunt[s-2] shear_prod[m2/s3] buoy_prod[m2/s3] diss_tke[m2/s3] Lturb[m] Kmv[m2/s] Khv[m2/s] shoc_cldfrac[-] shoc_ql[kg/kg] shoc_cond[kg/kg] qsrc_sw[K/s] qsrc_lw[K/s]\n";
}

void PBLDiagnosticRecorder::write (amrex::Real time, int step, const amrex::MultiFab& cons,
                                   const amrex::MultiFab& xvel, const amrex::MultiFab& yvel,
                                   const amrex::MultiFab& zvel, const amrex::MultiFab& z_phys_cc,
                                   const amrex::Geometry& geom, const ShocDriver& driver,
                                   const MoistureComponentIndices& mi,
                                   const amrex::MultiFab* qheating_rates)
{
    const bool summary = summary_step(step), profile = profile_step(step);
    if (!summary && !profile) return;
    const int klo = geom.Domain().smallEnd(2), khi = geom.Domain().bigEnd(2), nz = khi-klo+1;
    for (int site = 0; site < static_cast<int>(m_iloc.size()); ++site) {
        const auto loc = InstrumentSimUtil::site_location_at_time(m_iloc[site], m_jloc[site], time, geom, m_move);
        const int i=loc.i, j=loc.j;
        auto z = InstrumentSimUtil::extract_column(z_phys_cc,0,i,j,klo,khi);
        auto rho = InstrumentSimUtil::extract_column(cons,Rho_comp,i,j,klo,khi);
        auto rt = InstrumentSimUtil::extract_column(cons,RhoTheta_comp,i,j,klo,khi);
        auto rho_ke = column_or_missing(&cons,RhoKE_comp,i,j,klo,khi);
        auto qv_mass = column_or_missing(&cons,mi.qv,i,j,klo,khi);
        auto qc_mass = column_or_missing(&cons,mi.qc,i,j,klo,khi);
        auto qi_mass = column_or_missing(&cons,mi.qi,i,j,klo,khi);
        auto theta = missing_column(klo,khi), theta_v = missing_column(klo,khi);
        auto qv = missing_column(klo,khi), qc = missing_column(klo,khi), qi = missing_column(klo,khi);
        auto ke = missing_column(klo,khi), u = missing_column(klo,khi), v = missing_column(klo,khi), w = missing_column(klo,khi);
        auto u0=InstrumentSimUtil::extract_column(xvel,0,i,j,klo,khi), u1=InstrumentSimUtil::extract_column(xvel,0,i+1,j,klo,khi);
        auto v0=InstrumentSimUtil::extract_column(yvel,0,i,j,klo,khi), v1=InstrumentSimUtil::extract_column(yvel,0,i,j+1,klo,khi);
        auto wf=InstrumentSimUtil::extract_face_column(zvel,0,i,j,klo,khi+1);
        auto pblh=column_or_missing(&driver.pblh_diagnostics(),0,i,j,klo,khi);
        auto wek=column_or_missing(&driver.we_kinematic_diagnostics(),0,i,j,klo,khi), wef=column_or_missing(&driver.we_flux_jump_diagnostics(),0,i,j,klo,khi);
        auto pblht=column_or_missing(&driver.pblh_tendency_diagnostics(),0,i,j,klo,khi), pblhh=column_or_missing(&driver.pblh_hadv_diagnostics(),0,i,j,klo,khi);
        auto wp=column_or_missing(&driver.w_at_pblh_diagnostics(),0,i,j,klo,khi), dtv=column_or_missing(&driver.delta_theta_v_diagnostics(),0,i,j,klo,khi), wtvp=column_or_missing(&driver.wthv_at_pblh_diagnostics(),0,i,j,klo,khi);
        auto wtv=column_or_missing(&driver.wthv_sec_diagnostics(),0,i,j,klo,khi), wtl=column_or_missing(&driver.wthl_sec_diagnostics(),0,i,j,klo,khi), wqw=column_or_missing(&driver.wqw_sec_diagnostics(),0,i,j,klo,khi), ws=column_or_missing(&driver.w_sec_diagnostics(),0,i,j,klo,khi), w3=column_or_missing(&driver.w3_diagnostics(),0,i,j,klo,khi);
        auto br=column_or_missing(&driver.brunt_diagnostics(),0,i,j,klo,khi), shear=column_or_missing(&driver.shear_prod_diagnostics(),0,i,j,klo,khi), buoy=column_or_missing(&driver.buoy_prod_diagnostics(),0,i,j,klo,khi), diss=column_or_missing(&driver.diss_tke_diagnostics(),0,i,j,klo,khi);
        auto eddy=column_or_missing(&driver.native_diagnostics(),EddyDiff::Turb_lengthscale,i,j,klo,khi), kmv=column_or_missing(&driver.native_diagnostics(),EddyDiff::Mom_v,i,j,klo,khi), khv=column_or_missing(&driver.native_diagnostics(),EddyDiff::Theta_v,i,j,klo,khi);
        auto ust=column_or_missing(&driver.shoc_ustar_diagnostics(),0,i,j,klo,khi), olen=column_or_missing(&driver.shoc_olen_diagnostics(),0,i,j,klo,khi), cf=column_or_missing(&driver.shoc_cldfrac_diagnostics(),0,i,j,klo,khi), qlc=column_or_missing(&driver.shoc_ql_diagnostics(),0,i,j,klo,khi), cond=column_or_missing(&driver.shoc_cond_diagnostics(),0,i,j,klo,khi);
        auto sens=column_or_missing(&driver.consumed_sens_flux_diagnostics(),0,i,j,klo,khi), latent=column_or_missing(&driver.consumed_laten_flux_diagnostics(),0,i,j,klo,khi), qsw=column_or_missing(qheating_rates,0,i,j,klo,khi), qlw=column_or_missing(qheating_rates,1,i,j,klo,khi);
        for (int n=0; n<nz; ++n) {
            if (valid(rho[n]) && rho[n] > 0.0) {
                theta[n]=rt[n]/rho[n];
                if (mi.qv >= 0 && valid(qv_mass[n])) qv[n]=qv_mass[n]/rho[n];
                if (mi.qc >= 0 && valid(qc_mass[n])) qc[n]=qc_mass[n]/rho[n];
                if (mi.qi >= 0 && valid(qi_mass[n])) qi[n]=qi_mass[n]/rho[n];
                if (valid(theta[n]) && valid(qv[n])) theta_v[n]=theta[n]*(1.0+0.61*qv[n]-(valid(qc[n])?qc[n]:0.0)-(valid(qi[n])?qi[n]:0.0));
                if (cons.nComp() > RhoKE_comp && valid(rho_ke[n])) ke[n]=rho_ke[n]/rho[n];
            }
            if (valid(u0[n]) && valid(u1[n])) u[n]=0.5*(u0[n]+u1[n]);
            if (valid(v0[n]) && valid(v1[n])) v[n]=0.5*(v0[n]+v1[n]);
            if (n+1 < static_cast<int>(wf.size()) && valid(wf[n]) && valid(wf[n+1])) w[n]=0.5*(wf[n]+wf[n+1]);
        }
        const int k0=0;
        const amrex::Real dz0=layer_dz(z,k0), z_sfc=valid(z[0])&&valid(dz0)?z[0]-0.5*dz0:missing, ph=safe(pblh[0]);
        const amrex::Real top=valid(z_sfc)&&valid(ph)?z_sfc+ph:missing;
        const int kt=nearest_level(z,top);
        amrex::Real tkei=0, shear_i=0, buoy_i=0, diss_i=0, depth=0, lwp=mi.qc>=0?0:missing, iwp=mi.qi>=0?0:missing, swsum=0, brmax=missing, cfmax=missing, lwmin=missing;
        int swcount=0;
        for (int n=0; n<nz; ++n) {
            const amrex::Real dz=layer_dz(z,n);
            if (valid(dz)&&dz>0&&valid(z[n])&&valid(top)&&z[n]<=top) {
                if (valid(ke[n])) { tkei+=ke[n]*dz; depth+=dz; }
                if (valid(shear[n])) shear_i+=shear[n]*dz;
                if (valid(buoy[n])) buoy_i+=buoy[n]*dz;
                if (valid(diss[n])) diss_i+=diss[n]*dz;
                if (valid(qc[n])&&valid(rho[n])) lwp+=rho[n]*qc[n]*dz;
                if (valid(qi[n])&&valid(rho[n])) iwp+=rho[n]*qi[n]*dz;
                if (valid(cf[n])) cfmax=valid(cfmax)?std::max(cfmax,cf[n]):cf[n];
                if (valid(qsw[n])) { swsum+=qsw[n]; ++swcount; }
            }
            const int window=std::max(1, kt>=0?std::abs(kt-k0):1);
            if (kt>=0&&std::abs(n-kt)<=window) {
                if(valid(br[n])) brmax=valid(brmax)?std::max(brmax,br[n]):br[n];
                if(valid(qlw[n])) lwmin=valid(lwmin)?std::min(lwmin,qlw[n]):qlw[n];
            }
        }
        const amrex::Real wtv0=valid(wtv[0])?wtv[0]:missing;
        const amrex::Real ratio=valid(wtv0)&&wtv0>1.e-12&&valid(wtvp[0])?-wtvp[0]/wtv0:missing;
        const amrex::Real wcube=valid(theta_v[0])&&valid(wtv0)&&wtv0>0&&valid(ph)&&ph>0?(CONST_GRAV/theta_v[0])*wtv0*ph:missing;
        const amrex::Real wstar=valid(wcube)&&wcube>0?std::cbrt(wcube):missing;
        if (summary && amrex::ParallelDescriptor::IOProcessor() && time>m_last_summary_time[site]) {
            auto& os=*m_summary_files[site];
            write_value(os,time); os<<' '<<step<<' '; write_value(os,loc.x); os<<' '; write_value(os,loc.y)<<' '<<i<<' '<<j<<' '; write_value(os,ph); os<<' '; write_value(os,top)<<' ';
            write_value(os,wek[0]);os<<' ';write_value(os,wef[0]);os<<' ';write_value(os,pblht[0]);os<<' ';write_value(os,pblhh[0]);os<<' ';write_value(os,wp[0]);os<<' ';write_value(os,dtv[0]);os<<' ';write_value(os,wtvp[0]);os<<' ';write_value(os,ratio);os<<' ';write_value(os,wtv0);os<<' ';write_value(os,ust[0]);os<<' ';write_value(os,olen[0]);os<<' ';write_value(os,wstar);os<<' ';
            write_value(os,valid(sens[0])?surface_flux_diagnostics::sensible_heat_flux_wm2_from_rhotheta_flux(sens[0]):missing);os<<' ';write_value(os,valid(latent[0])?surface_flux_diagnostics::latent_heat_flux_wm2_from_rhoqv_flux(latent[0]):missing);os<<' ';write_value(os,depth>0?tkei/depth:missing);os<<' ';write_value(os,shear_i);os<<' ';write_value(os,buoy_i);os<<' ';write_value(os,diss_i);os<<' ';write_value(os,shear_i+buoy_i-diss_i);os<<' ';write_value(os,brmax);os<<' ';write_value(os,cfmax);os<<' ';write_value(os,lwp);os<<' ';write_value(os,iwp);os<<' ';write_value(os,lwmin);os<<' ';write_value(os,swcount>0?swsum/swcount:missing)<<'\n';
            os.flush(); m_last_summary_time[site]=time;
        }
        if (profile && amrex::ParallelDescriptor::IOProcessor() && time>m_last_profile_time[site]) {
            auto& os=*m_profile_files[site];
            for (int n=0;n<nz;++n) {
                write_value(os,time);os<<' '<<step<<' ';write_value(os,loc.x);os<<' ';write_value(os,loc.y)<<' '<<i<<' '<<j<<' ';write_value(os,pblh[n]);os<<' ';write_value(os,wek[n]);os<<' ';write_value(os,wef[n]);os<<' '<<(klo+n)<<' ';write_value(os,z[n]);os<<' ';write_value(os,valid(z[n])&&valid(z_sfc)?z[n]-z_sfc:missing);os<<' ';write_value(os,theta[n]);os<<' ';write_value(os,theta_v[n]);os<<' ';write_value(os,qv[n]);os<<' ';write_value(os,qc[n]);os<<' ';write_value(os,qi[n]);os<<' ';write_value(os,u[n]);os<<' ';write_value(os,v[n]);os<<' ';write_value(os,w[n]);os<<' ';write_value(os,ke[n]);os<<' ';write_value(os,wtv[n]);os<<' ';write_value(os,wtl[n]);os<<' ';write_value(os,wqw[n]);os<<' ';write_value(os,ws[n]);os<<' ';write_value(os,w3[n]);os<<' ';write_value(os,br[n]);os<<' ';write_value(os,shear[n]);os<<' ';write_value(os,buoy[n]);os<<' ';write_value(os,diss[n]);os<<' ';write_value(os,eddy[n]);os<<' ';write_value(os,kmv[n]);os<<' ';write_value(os,khv[n]);os<<' ';write_value(os,cf[n]);os<<' ';write_value(os,qlc[n]);os<<' ';write_value(os,cond[n]);os<<' ';write_value(os,qsw[n]);os<<' ';write_value(os,qlw[n])<<'\n';
            }
            os.flush(); m_last_profile_time[site]=time;
        }
    }
}
