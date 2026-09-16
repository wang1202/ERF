.nr PS 10
.nr VS 12
.nr LL 6.5i
.nr PO 1.0i
.nr HM 0.65i
.nr FM 0.65i
.ds CH Aerosol-Cloud Interactions in ERF
.ds CF - % -
.
.TL
Aerosol-Cloud Interactions in ERF
.AU
A practical experiment and development plan using bulk microphysics and SHOC
.AI
Technical note based on the ERF source tree reviewed 31 August 2026
.AB
This note describes how to investigate aerosol-cloud interactions with the
microphysics schemes currently available in ERF. Morrison and WDM6 provide the
most direct bulk-microphysics pathways: Morrison accepts a prescribed cloud
droplet number, while WDM6 carries cloud condensation nuclei (CCN), cloud
droplet number, and rain number. Kessler, SAM, and WSM6 are useful for
controlled process-sensitivity studies but do not currently represent aerosol
activation explicitly. SHOC adds the turbulence and cloud-macrophysics
dimension, although its native state-update coupling presently lacks the
number closure required by Morrison and WDM6. A staged program is recommended
to separate aerosol effects, uncertain warm-rain process rates, turbulent
mixing, and radiation coupling.
.AE
.
.NH 1
Main conclusions
.PP
The scientifically cleanest near-term program is to use Morrison-NoIce and
WDM6 for warm-cloud aerosol perturbations, with identical dynamics, forcing,
resolution, and surface conditions. Morrison supplies a controlled
prescribed-droplet-number experiment. WDM6 is the more literal bulk aerosol
experiment because it predicts CCN, cloud-droplet number, and rain number.
.PP
Kessler, SAM, and WSM6 can still answer important questions about cloud and
precipitation susceptibility. Their autoconversion thresholds, collection
efficiencies, and accretion rates should be interpreted as microphysical
process uncertainty or aerosol proxies, not as explicit aerosol activation.
.PP
SHOC is not an aerosol model. It changes the thermodynamic variance,
subgrid-cloud diagnosis, mixing, and entrainment through which aerosol
perturbations operate. The ideal target is therefore a factorial experiment
combining aerosol level with SHOC closure settings. In the current native
SHOC implementation, however, number-aware microphysics requires an additional
closure before that full experiment is possible.
.
.NH 1
Physical framing
.PP
The chain to resolve is:
.DS C
.ft B
Aerosol/CCN -> activation -> droplet number and size
     -> autoconversion/accretion -> precipitation
     -> cloud water, cloud fraction, lifetime, and radiation
.ft R
.DE
.PP
SHOC acts mainly on activation conditions, cloud fraction, condensate, and
entrainment through its representation of subgrid variability and turbulent
mixing. Bulk-microphysics tuning acts mainly on the conversion from condensate
to precipitation. Separating these links is essential: changing an accretion
coefficient can reproduce a precipitation response, but it is not equivalent
to changing aerosol concentration.
.
.NH 1
Capabilities of the available schemes
.na
.TS
center box tab(|);
lw(0.9i) | lw(2.25i) | lw(1.85i).
T{
\fBScheme\fP
T}|T{
\fBCurrent aerosol-related representation\fP
T}|T{
\fBRecommended role\fP
T}
_
Kessler|T{
Warm-rain, single moment; predicts water vapor, cloud water, and rainwater only.
T}|T{
Autoconversion and accretion process sensitivity; aerosol proxy only.
T}
SAM|T{
Single-moment mixed-phase scheme with fixed thresholds, intercepts, and collection efficiencies.
T}|T{
Process uncertainty for liquid, ice, snow, and graupel conversion.
T}
WSM6|T{
Single-moment mixed phase; cloud-droplet number is a fixed internal assumption, xncr = 3e8 m^-3 (300 cm^-3).
T}|T{
Expose xncr as a runtime input for prescribed-droplet-number proxy experiments.
T}
Morrison|T{
Two-moment framework, but ERF currently uses constant cloud-droplet number selected by erf.morrison_ndcnst.
T}|T{
Best no-code prescribed-Nc experiment; use Morrison-NoIce first for warm clouds.
T}
WDM6|T{
Predicts CCN aerosol (nn), cloud-droplet number (nc), and rain number (nr); ice remains single moment.
T}|T{
Best existing bulk scheme for explicit aerosol-cloud interaction experiments.
T}
SuperDroplets|T{
Represents aerosol size/composition, activation, condensation, stochastic collision-coalescence, and sedimentation with particles.
T}|T{
Higher-fidelity benchmark for a limited subset of bulk experiments.
T}
.TE
.ad
.
.NH 2
Kessler
.PP
Kessler has no aerosol or number-concentration variable. Its warm-rain
conversion currently uses a cloud-water autoconversion threshold of
1e-3 kg kg^-1, an autoconversion coefficient of 1e-3 s^-1, and an accretion
term proportional to 2.2 qc qr^0.875. These values are compiled constants.
.PP
Useful additions would be scheme-local runtime inputs for the autoconversion
threshold, autoconversion multiplier, and accretion multiplier. Defaults should
exactly preserve the present solution. Such experiments quantify structural
sensitivity of precipitation formation; they should not be labeled direct
aerosol experiments.
.
.NH 2
SAM
.PP
SAM is a single-moment mixed-phase scheme. It shares fixed liquid and ice
autoconversion thresholds and uses fixed rain, snow, and graupel distribution
intercepts and collection efficiencies. Recommended tuners are liquid and ice
autoconversion thresholds, conversion-rate multipliers, and collection-
efficiency multipliers. These are particularly useful for assessing whether an
apparent aerosol signal is robust to uncertainty in ice and precipitation
processes.
.
.NH 2
WSM6
.PP
WSM6 assumes xncr = 3e8 m^-3, equivalent to 300 cm^-3, inside the scheme. A
runtime wsm6.xncr option would enable clean, intermediate, and polluted proxy
experiments. Coefficients derived from xncr must be recomputed consistently in
both the native C++ and Fortran paths. Reasonable first levels are 50, 300, and
1000 cm^-3, followed by a denser sweep if the response is nonlinear.
.PP
Higher prescribed droplet number generally decreases droplet size and delays
warm-rain formation in this class of parameterization. It remains a prescribed
single-moment assumption: aerosol transport, activation, depletion, and
regeneration are absent.
.
.NH 2
Morrison
.PP
ERF exposes the constant droplet concentration through
erf.morrison_ndcnst, with a default of 250 cm^-3. This is immediately suitable
for a prescribed-Nc ensemble. Begin with Morrison-NoIce to isolate the warm-
cloud pathway, using 50, 250, and 1000 cm^-3 as clean, intermediate, and
polluted conditions. Repeat selected cases with full Morrison to evaluate
mixed-phase dependence.
.PP
This experiment alters droplet size and the number-dependent warm-rain
pathway, but it still does not calculate activation from a prognostic aerosol
population. The interpretation is therefore a cloud-droplet-number
susceptibility, not a complete aerosol life-cycle simulation.
.
.NH 2
WDM6
.PP
WDM6 is the strongest current bulk option for aerosol-cloud interaction. It
uses a runtime background CCN concentration, wdm6.ccn0, and carries nn, nc, and
nr along with the hydrometeor mass fields. A useful first sweep is 1e8, 5e8,
and 1e9 m^-3, equivalent to 100, 500, and 1000 cm^-3.
.PP
There is an important implementation caveat. At every macro timestep, nn is
clamped into the interval 1e8 to 2e10 m^-3. Values below 1e8 m^-3 are therefore
ineffective, and activation depletion below the floor is replenished. The
current ccn0 should be described as a maintained background reservoir rather
than a fully conserved aerosol population. The land/water selection also
affects the autoconversion regime, so erf.is_land and wdm6.hail_opt must remain
fixed across the aerosol ensemble.
.
.NH 2
SuperDroplets as a benchmark
.PP
The SuperDroplets option can specify dry aerosol number, size distribution,
and composition and then represent activation and collision-coalescence
explicitly. It is more expensive and statistically noisy, but a small matched
set is valuable for determining where Morrison and WDM6 bulk assumptions
control the result. Match aerosol number and dry-size distribution as closely
as possible, and use multiple particle realizations.
.
.NH 1
Using SHOC to study turbulence-aerosol coupling
.PP
The most ACI-relevant shared SHOC controls are listed below. Change them one at
a time by approximately plus or minus 25 percent around the default before
attempting a multi-parameter ensemble.
.na
.TS
center box tab(|);
lw(1.7i) | lw(0.6i) | lw(2.7i).
T{
\fBInput\fP
T}|T{
\fBDefault\fP
T}|T{
\fBACI interpretation\fP
T}
_
erf.shoc.qw2tune|1.0|T{
Controls total-water variance used by the subgrid closure; influences cloud occurrence and condensate.
T}
erf.shoc.qwthl2tune|1.0|T{
Controls water-temperature covariance and therefore the joint thermodynamic PDF.
T}
erf.shoc.thl2tune|1.0|T{
Controls liquid-water-potential-temperature variance and saturation variability.
T}
erf.shoc.w2tune|1.0|T{
Controls vertical-velocity variance, relevant to updraft variability and activation conditions.
T}
erf.shoc.length_fac|0.5|T{
Changes the turbulent length scale and consequently mixing and entrainment.
T}
erf.shoc.coeff_kh|0.1|T{
Changes scalar and heat diffusivity; useful as a secondary mixing sensitivity.
T}
.TE
.ad
.PP
Enable erf.shoc.extra_shoc_diags and save shoc_cldfrac, shoc_ql, shoc_ql2,
shoc_cond, wqls_sec, w_sec, thl_sec, qw_sec, and qwthl_sec where available.
These diagnostics distinguish a turbulence-mediated condensate response from a
change produced directly by precipitation conversion.
.
.NH 2
Current coupling restriction
.PP
Native SHOC in moist simulations must use state_update because moist
host_diffusion does not own the required cloud macrophysics. The native
state-update driver currently rejects layouts containing cloud- or ice-number
concentrations: changing cloud mass without a corresponding number tendency
would leave the size distribution inconsistent. Thus native SHOC currently
supports Kessler, SAM, and WSM6 more naturally than Morrison or WDM6.
.PP
For Morrison with prescribed Nc, a first closure could consistently reimpose
the prescribed cloudy-cell Nc and clear it when cloud water vanishes. For WDM6,
a defensible closure must transfer number between nn and nc during activation
and return or remove number consistently during evaporation. This should be
implemented and tested before interpreting a WDM6-by-SHOC factorial ensemble.
.
.NH 1
Recommended experiment program
.NH 2
Cases and controls
.PP
Use RICO as the primary precipitating warm-cumulus case. It exercises the
aerosol-to-warm-rain pathway while remaining interpretable. Use DYCOMS as a
secondary stratocumulus case for cloud-water and lifetime responses, and BOMEX
as a weakly precipitating control. ERF already contains canonical cases for
these configurations.
.PP
For every comparison, hold the grid, timestep, forcing, sounding, radiation,
surface fluxes, land/water mask, microphysics substepping, and output cadence
fixed. Use at least three initial-perturbation realizations when feasible; a
single deterministic LES difference can be dominated by chaotic sampling.
.
.NH 2
Stage 1: no-code aerosol and droplet-number experiments
.na
.TS
center box tab(|);
lw(1.45i) | lw(2.0i) | lw(1.55i).
T{
\fBConfiguration\fP
T}|T{
\fBPerturbation levels\fP
T}|T{
\fBPrimary interpretation\fP
T}
_
Morrison-NoIce|T{
morrison_ndcnst = 50, 250, 1000 cm^-3
T}|T{
Prescribed cloud-droplet-number susceptibility.
T}
WDM6|T{
ccn0 = 1e8, 5e8, 1e9 m^-3
T}|T{
Bulk CCN activation and warm-rain response, subject to the nn floor.
T}
SuperDroplets subset|T{
Matched clean and polluted aerosol distributions
T}|T{
Benchmark for activation and collision-coalescence assumptions.
T}
.TE
.ad
.
.NH 2
Stage 2: process attribution
.PP
Add runtime multipliers for autoconversion and accretion, initially preserving
all present defaults. Perturb one rate at a time, then use a compact factorial
design around the aerosol levels. The goal is to determine whether the aerosol
response exceeds uncertainty in the bulk conversion formulae and whether the
sign of the response is robust.
.PP
Add budget diagnostics for activation (nn to nc), autoconversion (qc to qr),
accretion of cloud water by rain, rain evaporation, sedimentation, and any
ice-phase conversions used in the selected case. Record mass tendencies in
kg kg^-1 s^-1 and number tendencies in kg^-1 s^-1, with column-integrated
counterparts when possible.
.
.NH 2
Stage 3: SHOC interaction experiment
.PP
Before number closure is available, use native SHOC with WSM6 after xncr is
made configurable, treating xncr as an aerosol proxy. Once number closure is
implemented, repeat the Morrison and WDM6 clean/polluted endpoints with the
baseline and selected perturbed SHOC settings.
.PP
For any outcome Y, quantify the turbulence modulation with the
difference-in-differences interaction:
.EQ
I = (Y sub "polluted,SHOC" - Y sub "clean,SHOC")
    - (Y sub "polluted,noSHOC" - Y sub "clean,noSHOC") .
.EN
A nonzero I shows that the aerosol response depends on the turbulence and
subgrid-cloud treatment rather than simply adding to it.
.
.NH 2
Stage 4: radiation coupling
.PP
The present RRTMGP interface diagnoses cloud fraction as either zero or one in
each grid cell and assigns a fixed liquid effective radius of 10 micrometers
(and ice radius of 25 micrometers). It does not use predicted Nc or native
SHOC cloud fraction. Therefore existing experiments can capture cloud amount,
water-path, lifetime, and precipitation responses, but not a complete Twomey
effect.
.PP
Before claiming aerosol-cloud-radiation interaction, pass a microphysics-
derived liquid effective radius and the SHOC subgrid cloud fraction/in-cloud
water to radiation. The mass-number estimate
.EQ
r sub v = left [ {3 q sub c} over {4 pi rho sub w N sub c} right ] sup {1/3}
.EN
provides a volume-mean radius when qc and Nc use compatible per-dry-air-mass
units. Effective radius additionally requires a spectral-shape factor
consistent with the microphysics size distribution.
.
.NH 1
Diagnostics and evaluation metrics
.PP
Save three-dimensional qv, qc, qi, qrain, qsnow, and qgraup as applicable;
number concentrations nc, nr, nn, ni, ns, and ng where available; accumulated
surface precipitation; and SHOC thermodynamic and turbulence diagnostics.
Construct the following analysis products:
.IP \(bu 3
Liquid and ice water paths, cloud fraction, cloud-base and cloud-top height,
and cloud lifetime.
.IP \(bu 3
Surface precipitation rate, accumulated precipitation, onset time, and
precipitating area fraction.
.IP \(bu 3
Cloud-droplet and rain-number distributions, droplet-size estimates, and
joint qc-Nc histograms.
.IP \(bu 3
Domain-integrated microphysical process budgets and their vertical profiles.
.IP \(bu 3
Aerosol susceptibility, such as d ln(LWP)/d ln(N) and
d ln(precipitation)/d ln(N), with ensemble uncertainty.
.IP \(bu 3
Energy and water closure, including precipitation and sedimentation fluxes.
.PP
Compare time windows only after a common spin-up and report ensemble means,
spread, and temporal autocorrelation-aware uncertainty. Avoid interpreting a
single domain-mean difference without inspecting changes in the cloud
population distribution.
.
.NH 1
Implementation priorities
.IP 1. 4
Expose scheme-local Kessler and SAM conversion/collection parameters and WSM6
xncr as runtime inputs while preserving current defaults exactly.
.IP 2. 4
Add named microphysical process-rate diagnostics and mass/number budget tests.
.IP 3. 4
Implement and validate number closure for native SHOC with Morrison and WDM6.
.IP 4. 4
Revisit WDM6's per-timestep CCN floor so conserved and maintained-reservoir
experiments can be selected explicitly.
.IP 5. 4
Couple microphysics-derived effective radius and SHOC cloud fraction to
RRTMGP before evaluating the first aerosol indirect radiative effect.
.
.NH 1
Example ERF inputs
.NH 2
Morrison prescribed-droplet-number ensemble
.DS L
.ft CW
.ps 9
.nf
erf.moisture_model = Morrison_NoIce
erf.morrison_ndcnst = 50.0   # repeat: 250.0, 1000.0 cm^-3
erf.plot_vars_1 = qv qc qrain nc nr rain_accum
.fi
.ps 10
.ft R
.DE
.NH 2
WDM6 background-CCN ensemble
.DS L
.ft CW
.ps 9
.nf
erf.moisture_model = WDM6
wdm6.ccn0 = 1.0e8            # repeat: 5.0e8, 1.0e9 m^-3
wdm6.hail_opt = 0
erf.is_land = 0               # maritime idealized case; keep fixed
erf.plot_vars_1 = qv qc qrain nc nn nr rain_accum
.fi
.ps 10
.ft R
.DE
.NH 2
Native SHOC with a single-moment proxy scheme
.DS L
.ft CW
.ps 9
.nf
zlo.type = surface_layer
erf.pbl_type = NATIVE_SHOC
erf.moisture_model = WSM6
erf.shoc.transport_mode = state_update
erf.shoc.extra_shoc_diags = true
.fi
.ps 10
.ft R
.DE
.PP
This SHOC configuration is not an aerosol ensemble until WSM6 xncr is exposed
and varied. Request the SHOC diagnostic plot variables explicitly in addition
to enabling their calculation.
.
.NH 1
Code and documentation references
.ps 9
.vs 11p
.na
.PP
Paths are relative to the ERF repository root.
.PP
\fBConstants and single-moment schemes:\fP Source/ERF_Constants.H;
Source/Microphysics/Kessler/ERF_KesslerUtils.H; and
Source/Microphysics/WSM6/ERF_WSM6.H.
.PP
\fBMorrison and WDM6:\fP
Source/Microphysics/Morrison/ERF_AdvanceMorrison.cpp;
Source/Microphysics/WDM6/ERF_InitWDM6.cpp; and
Source/Microphysics/WDM6/ERF_AdvanceWDM6.cpp.
.PP
\fBSHOC and radiation:\fP Source/PBL/Shoc/ERF_ShocDriver.H and
Source/PhysicsInterfaces/Radiation/ERF_Radiation.cpp.
.PP
\fBDocumentation:\fP Docs/sphinx_doc/Inputs.rst;
Docs/sphinx_doc/theory/Microphysics.rst; and
Docs/sphinx_doc/plotfiles/Plotfile3DReference.rst.
.ad
.ps 10
.vs 12p
