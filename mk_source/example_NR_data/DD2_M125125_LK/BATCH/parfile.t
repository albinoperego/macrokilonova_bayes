ActiveThorns = "
    ADMBase
    ADMCoupling
    ADMDerivatives
    ADMMass
    ADMMacros
    AEILocalInterp
    AHFinderDirect
    BNSAnalysis
    BNSTrackerGen
    Boundary
    Carpet
    CarpetInterp
    CarpetInterp2
    CarpetIOASCII
    CarpetIOBasic
    CarpetIOHDF5
    CarpetIOScalar
    CarpetMask
    CarpetLib
    CarpetRegrid2
    CarpetReduce
    CarpetSlab
    CartesianCoordinates
    CartGrid3D
    Constants
    Composition
    CoordBase
    CoordGauge
    CTGBase
    Dissipation
    Formaline
    EOS_Barotropic
    EOS_Thermal
    EOS_Thermal_Extable
    EOS_Thermal_Table3d
    Fortran
    GlobalDerivative
    HydroBase
    HRSCCore
    ID_Switch_EOS
    InitBase
    IOUtil
    LocalInterp
    LoreneID
    LoopControl
    CTGEvolution
    CTGGauge
    CTGMatter
    CTGConstraints
    CTGRadiativeBC
    MoL
    Multipole
    NaNChecker
    NewRad
    Outflow
    PizzaBase
    PizzaIDBase
    PizzaNumUtils
    QuasiLocalMeasures
    ReflectionSymmetry
    Refluxing
    Slab
    SpaceMask
    SphericalSurface
    SphericalHarmonicDecomp
    StaticConformal
    SummationByParts
    SymBase
    SystemTopology
    TerminationTrigger
    THC_Core
    THC_LeakageBase
    THC_Refluxing
    THC_Tracer
    Time
    TimerReport
    TmunuBase
    Volomnia
    WeakRates
    WeylScal4
    WorldTube
    WatchDog
"

Cactus::terminate			= "never"

Formaline::output_source_subdirectory		= "../../cactus-source"

TerminationTrigger::max_walltime		= @WALLTIME_HOURS@
TerminationTrigger::on_remaining_walltime	= 20
TerminationTrigger::create_termination_file	= "yes"
TerminationTrigger::termination_from_file	= "yes"
TerminationTrigger::termination_file		= "../TERMINATE"

TimerReport::out_every		= 4096
TimerReport::output_all_timers_readable		= "yes"

# =============================================================================
# Grid
# =============================================================================
Grid::avoid_origin			= "no"
Grid::domain			= "full"
Grid::type			= "coordbase"

ReflectionSymmetry::reflection_x		= "no"
ReflectionSymmetry::reflection_y		= "no"
ReflectionSymmetry::reflection_z		= "yes"
ReflectionSymmetry::avoid_origin_x		= "yes"
ReflectionSymmetry::avoid_origin_y		= "yes"
ReflectionSymmetry::avoid_origin_z		= "yes"

Volomnia::symm_weight		= 2

CoordBase::xmin			= -1024
CoordBase::xmax			=  1024
CoordBase::ymin			= -1024
CoordBase::ymax			=  1024
CoordBase::zmin			=  0
CoordBase::zmax			=  1024

CoordBase::spacing			= "numcells"
CoordBase::ncells_x			= 256
CoordBase::ncells_y			= 256
CoordBase::ncells_z			= 128

CoordBase::boundary_size_x_lower		= 3
CoordBase::boundary_size_x_upper		= 3
CoordBase::boundary_shiftout_x_lower		= 0
CoordBase::boundary_shiftout_x_upper		= 0
CoordBase::boundary_staggered_x_lower		= "yes"
CoordBase::boundary_staggered_x_upper		= "yes"

CoordBase::boundary_size_y_lower		= 3
CoordBase::boundary_size_y_upper		= 3
CoordBase::boundary_shiftout_y_lower		= 0
CoordBase::boundary_shiftout_y_upper		= 0
CoordBase::boundary_staggered_y_lower		= "yes"
CoordBase::boundary_staggered_y_upper		= "yes"

CoordBase::boundary_size_z_lower		= 3
CoordBase::boundary_size_z_upper		= 3
CoordBase::boundary_shiftout_z_lower		= 0
CoordBase::boundary_shiftout_z_upper		= 0
CoordBase::boundary_staggered_z_lower		= "yes"
CoordBase::boundary_staggered_z_upper		= "yes"

Driver::ghost_size			= 3
Driver::ghost_size_x		= 3
Driver::ghost_size_y		= 3
Driver::ghost_size_z		= 3

Carpet::refinement_centering		= "cell"
Carpet::domain_from_coordbase		= "yes"

InitBase::initial_data_setup_method		= "init_some_levels"

Carpet::max_refinement_levels		= 8
Carpet::prolongation_order_space		= 3
Carpet::prolongation_order_time		= 2
Carpet::use_buffer_zones		= "yes"
Carpet::enable_all_storage		= "no"
Carpet::init_fill_timelevels		= "yes"
Carpet::time_refinement_factors		= "[1,1,2,4,8,16,32,64,128,256,512]"

Carpet::grid_coordinates_filename		= "grid.carpet"

CarpetLib::poison_new_memory		= "yes"

CarpetRegrid2::num_centres		= 3
CarpetRegrid2::regrid_every		= 128
CarpetRegrid2::snap_to_coarse		= "yes"
CarpetRegrid2::freeze_unaligned_levels		= "yes"
CarpetRegrid2::freeze_unaligned_parent_levels	= "yes"

# ----------- Region 1 --------------------
CarpetRegrid2::active_1		= "yes"
CarpetRegrid2::num_levels_1		= 7
CarpetRegrid2::position_x_1		= 13.54
CarpetRegrid2::position_y_1		= 0.0
CarpetRegrid2::radius_1[1]		= 300.0
CarpetRegrid2::radius_1[2]		= 160.0
CarpetRegrid2::radius_1[3]		= 80.0
CarpetRegrid2::radius_1[4]		= 40.0
CarpetRegrid2::radius_1[5]		= 24.0
CarpetRegrid2::radius_1[6]		= 12.0

# ----------- Region 2 --------------------
CarpetRegrid2::active_2		= "yes"
CarpetRegrid2::num_levels_2		= 7
CarpetRegrid2::position_x_2		= -13.54
CarpetRegrid2::position_y_2		= -0.0
CarpetRegrid2::radius_2[1]		= 300.0
CarpetRegrid2::radius_2[2]		= 160.0
CarpetRegrid2::radius_2[3]		= 80.0
CarpetRegrid2::radius_2[4]		= 40.0
CarpetRegrid2::radius_2[5]		= 24.0
CarpetRegrid2::radius_2[6]		= 12.0

# ----------- Region 3 --------------------
CarpetRegrid2::active_3		= "yes"
CarpetRegrid2::num_levels_3		= 3
CarpetRegrid2::position_x_3		= 0

CarpetRegrid2::radius_3[1]		= 320.0
CarpetRegrid2::radius_3[2]		= 160.0
CarpetRegrid2::radius_3[3]		= 80.0
CarpetRegrid2::radius_3[4]		= 40.0
CarpetRegrid2::radius_x_3[5]		= 24.0
CarpetRegrid2::radius_y_3[5]		= 24.0
CarpetRegrid2::radius_z_3[5]		= 16.0
CarpetRegrid2::radius_3[6]		= 12.0
CarpetRegrid2::radius_3[7]		= 5.0
# -----------------------------------------

BNSTrackerGen::sym_pi		= "no"
BNSTrackerGen::analysis_reflevel		= 6
BNSTrackerGen::analyze_every		= 128
BNSTrackerGen::merge_separation		= 1.5
BNSTrackerGen::collapse_separation		= -1
BNSTrackerGen::add_levels_post_merge		= 4
BNSTrackerGen::add_levels_post_collapse		= 1

CarpetRegrid2::symmetry_rotating180		= "no"
CarpetRegrid2::verbose		= "no"

NaNChecker::check_every		= 128
NaNChecker::check_vars		= "
    HydroBase::rho
    HydroBase::vel
    HydroBase::w_lorentz
    THC_Core::densxn
    THC_Core::densxp
    THC_Core::scon
    THC_Core::tau
    THC_Core::volform
"
NaNChecker::action_if_found		= "abort"

# =============================================================================
# Time integration
# =============================================================================
Carpet::num_integrator_substeps		= 3

MoL::ode_method			= "RK3"
MoL::MoL_Intermediate_Steps		= 3
MoL::MoL_Num_Scratch_Levels		= 0
MoL::verbose			= "register"

HydroBase::timelevels		= 3

Time::timestep_method		= "courant_static"
Time::dtfac			= 0.075

# =============================================================================
# Initial data
# =============================================================================
ADMBase::initial_data		= "LoreneBNS"
ADMBase::initial_lapse		= "LoreneBNS"
ADMBase::initial_shift		= "zero"
ADMBase::initial_dtlapse		= "zero"
ADMBase::initial_dtshift		= "zero"
HydroBase::initial_hydro		= "LoreneBNS"
HydroBase::initial_entropy		= "THCode"
HydroBase::initial_temperature		= "LoreneBNS"
HydroBase::initial_Y_e		= "LoreneBNS"
HydroBase::initial_Abar		= "zero"

# File describing a one-parametric EOS in Pizza format. Used only for initial data.
PizzaIDBase::eos_file		= "@HOME@/Data/EOS/DD2/0.01MeV/eos_dd2_adb.pizza"
LoreneID::lorene_bns_file		= "@HOME@/Data/LoreneData/DD2_T05_125125_40km/resu.d"
LoreneID::rho_cut			= 1.0e-5

# Geometric unit system for initial data, specified by length unit.
# use CACTUS units
PizzaBase::length_unit		= 1476.7161818921163

# Switch EOS
ID_Switch_EOS::sync_eps_temp		= "yes"
ID_Switch_EOS::temp_from_eps		= "no"
ID_Switch_EOS::limit_efrac		= "yes"

# =============================================================================
# Templated hydrodynamics code
# =============================================================================
HydroBase::evolution_method		= "THCode"
HydroBase::temperature_evolution_method		= "THCode"
HydroBase::entropy_evolution_method		= "THCode"
HydroBase::Y_e_evolution_method		= "THCode"
HydroBase::Abar_evolution_method		= "Composition"

THC_Core::eos_type			= "nuclear"
THC_Core::physics			= "GRHD"

TmunuBase::prolongation_type		= "none"
TmunuBase::stress_energy_storage		= "yes"
TmunuBase::stress_energy_at_RHS		= "yes"
TmunuBase::support_old_CalcTmunu_mechanism	= "no"

THC_Core::bc_type			= "none"

HRSCCore::scheme			= "FV"
HRSCCore::pplim			= "yes"
HRSCCore::refluxing			= "yes"
HRSCCore::reconstruction		= "MP5"
HRSCCore::riemann_solver		= "HLLE"
HRSCCore::flux_split		= "LLF"
HRSCCore::system_split		= "components"
HRSCCore::speed_eps			= 0.05

THC_Core::atmo_rho			= 1e-14
THC_Core::atmo_tolerance		= 0.0
THC_Core::atmo_temperature		= 0.02

THC_Core::c2a_BH_alp		= 0.15
THC_Core::c2a_rho_strict		= 2e-5
THC_Core::c2a_set_to_nan_on_failure		= "no"
THC_Core::c2a_fix_conservatives		= "yes"
THC_Core::c2a_kill_on_failure		= "no"

EOS_Thermal::evol_eos_name		= "Extable"
EOS_Thermal_Extable::rho_max		= 1e10
EOS_Thermal_Extable::temp_max		= 1000
EOS_Thermal_Extable::extend_ye		= "yes"

EOS_Thermal_Table3d::eos_db_loc		= "@HOME@/Data/EOS"
EOS_Thermal_Table3d::eos_folder		= "DD2"
EOS_Thermal_Table3d::eos_filename		= "DD2_DD2_hydro_30-Mar-2015.h5"

THC_LeakageBase::neu_abs_type		= "None"
THC_LeakageBase::num_sph_grids		= 2
THC_LeakageBase::center_grid1[0]		= 13.54
THC_LeakageBase::center_grid2[0]		= -13.54
THC_LeakageBase::store_free_rates		= "yes"
THC_LeakageBase::store_neu_luminosity		= "yes"

THC_Refluxing::nvars		= 6

WeakRates::table_filename		= "@HOME@/Data/EOS/DD2/DD2_DD2_weak_30-Mar-2015.h5"
WeakRates::use_rho_max_ext		= "yes"
WeakRates::use_rho_min_ext		= "yes"
WeakRates::use_temp_max_ext		= "yes"
WeakRates::use_temp_min_ext		= "yes"
WeakRates::use_ye_max_ext		= "yes"
WeakRates::use_ye_min_ext		= "yes"

# =============================================================================
# Spacetime evolution
# =============================================================================
ADMBase::evolution_method		= "CTGamma"
ADMBase::lapse_evolution_method		= "1+log"
ADMBase::shift_evolution_method		= "gamma-driver"
ADMBase::dtlapse_evolution_method		= "1+log"
ADMBase::dtshift_evolution_method		= "gamma-driver"

CTGBase::timelevels			= 3
CTGBase::conformal_factor_type		= "w"
CTGBase::use_matter			= "yes"

CTGEvolution::bc			= radiative
CTGGauge::eta			= 1.0
CTGGauge::damping_factor_method		= prescribed
CTGGauge::damping_factor_type		= "Schnetter-Simple"
CTGGauge::eta_damping_radius		= 250.0
CTGEvolution::force_lndetg_zero		= yes
CTGBase::evolution_system_type		= Z4c
CTGEvolution::kappa1		= 0.0
CTGEvolution::kappa2		= 0.0
CTGEvolution::MaxNumEvolvedVars		= 18
CTGEvolution::MaxNumConstrainedVars		= 13

CTGConstraints::constraints_persist		= yes

SummationByParts::order		= 4
SummationByParts::sbp_upwind_deriv		= no
SummationByParts::onesided_outer_boundaries	= yes
SummationByParts::onesided_interpatch_boundaries	= no
SummationByParts::sbp_1st_deriv		= yes
SummationByParts::sbp_2nd_deriv		= no

SummationByParts::use_dissipation		= no
GlobalDerivative::use_dissipation		= yes
SummationByParts::scale_with_h		= yes
SummationByParts::dissipation_type		= "Kreiss-Oliger"
SummationByParts::epsdis		= 0.1
GlobalDerivative::epsdis_for_level	[0]	= 0.1

SummationByParts::vars		= "
ADMBase::lapse
ADMBase::shift
CTGBase::conformal_factor
CTGBase::conformal_metric
CTGBase::curvature_scalar_a
CTGBase::curvature_scalar_b
CTGBase::curvature_tensor
CTGBase::Gamma"

THC_Core::fd_order			= 4

# =============================================================================
# Analysis
# =============================================================================
ADMMass::ADMMass_volume_radius[0]		= 450.0

#----- Wave Extraction using "WeylScal4" and "Multipole" --------
WeylScal4::offset			= 1.0e-8
WeylScal4::fd_order			= "4th" # old style
# WeylScal4::fdOrder		= 4  # New feature

WeylScal4::verbose			= 0

Multipole::nradii			= 4

Multipole::radius[0]		= 200
Multipole::radius[1]		= 300
Multipole::radius[2]		= 400
Multipole::radius[3]		= 500
Multipole::ntheta			= 120
Multipole::nphi			= 240

Multipole::variables		= "
    WeylScal4::Psi4r{sw=-2 cmplx='WeylScal4::Psi4i' name='Psi4'}"
Multipole::l_max			= 8

#------ Spheres for characteristic extraction --------------------
SphericalSlice::nslices		= 2
SphericalSlice::precalc_sYlms		= no
SphericalSlice::use_carpet_interp1		= yes

SphericalSlice::set_spherical	[0]	= yes
SphericalSlice::radius	[0]	= 200.0
SphericalSlice::type	[0]	= 1patch
SphericalSlice::use_Llama	[0]	= no
SphericalSlice::nghostzones	[0]	= 0
SphericalSlice::ntheta	[0]	= 120
SphericalSlice::nphi	[0]	= 240

SphericalSlice::set_spherical	[1]	= yes
SphericalSlice::radius	[1]	= 500.0
SphericalSlice::type	[1]	= 1patch
SphericalSlice::use_Llama	[1]	= no
SphericalSlice::nghostzones	[1]	= 0
SphericalSlice::ntheta	[1]	= 120
SphericalSlice::nphi	[1]	= 240

#----- CCE data -------------------------------

WorldTube::boundary_behavior		= CCE
WorldTube::lmax			= 8
WorldTube::ntubes			= 2
WorldTube::which_slice_to_take	[0]	= 0
WorldTube::which_slice_to_take	[1]	= 1
ADMDerivatives::store_time_derivatives		= yes
ADMDerivatives::store_radial_derivatives	= yes
ADMDerivatives::timelevels		= 3
ADMBase::lapse_timelevels		= 3
ADMBase::shift_timelevels		= 3
ADMBase::metric_timelevels		= 3

#------SphericalSurface-----------------------
SphericalSurface::nsurfaces		= 6
SphericalSurface::maxntheta		= 140
SphericalSurface::maxnphi		= 240

SphericalSurface::ntheta	[0]	= 55
SphericalSurface::nphi	[0]	= 96
SphericalSurface::nghoststheta	[0]	= 2
SphericalSurface::nghostsphi	[0]	= 2

SphericalSurface::ntheta	[1]	= 55
SphericalSurface::nphi	[1]	= 96
SphericalSurface::nghoststheta	[1]	= 2
SphericalSurface::nghostsphi	[1]	= 2

SphericalSurface::ntheta	[2]	= 55
SphericalSurface::nphi	[2]	= 96
SphericalSurface::nghoststheta	[2]	= 2
SphericalSurface::nghostsphi	[2]	= 2
SphericalSurface::set_spherical	[2]	= "yes"
SphericalSurface::radius	[2]	= 200

SphericalSurface::ntheta	[3]	= 55
SphericalSurface::nphi	[3]	= 96
SphericalSurface::nghoststheta	[3]	= 2
SphericalSurface::nghostsphi	[3]	= 2
SphericalSurface::set_spherical	[3]	= "yes"
SphericalSurface::radius	[3]	= 300

SphericalSurface::ntheta	[4]	= 55
SphericalSurface::nphi	[4]	= 96
SphericalSurface::nghoststheta	[4]	= 2
SphericalSurface::nghostsphi	[4]	= 2
SphericalSurface::set_spherical	[4]	= "yes"
SphericalSurface::radius	[4]	= 400

SphericalSurface::ntheta	[5]	= 55
SphericalSurface::nphi	[5]	= 96
SphericalSurface::nghoststheta	[5]	= 2
SphericalSurface::nghostsphi	[5]	= 2
SphericalSurface::set_spherical	[5]	= "yes"
SphericalSurface::radius	[5]	= 500

#------AHFinderDirect-----------------------
AHFinderDirect::N_horizons		= 2
AHFinderDirect::find_every		= 128
AHFinderDirect::output_h_every		= 0
AHFinderDirect::max_Newton_iterations__initial	= 50
AHFinderDirect::max_Newton_iterations__subsequent	= 50
AHFinderDirect::max_allowable_Theta_growth_iterations	= 10
AHFinderDirect::max_allowable_Theta_nonshrink_iterations	= 10
AHFinderDirect::geometry_interpolator_name	= "Lagrange polynomial interpolation"
AHFinderDirect::geometry_interpolator_pars	= "order=4"
AHFinderDirect::surface_interpolator_name	= "Lagrange polynomial interpolation"
AHFinderDirect::surface_interpolator_pars	= "order=4"
AHFinderDirect::verbose_level		= "physics highlights"
AHFinderDirect::move_origins		= "yes"

AHFinderDirect::origin_x[1]		= 0.0
AHFinderDirect::initial_guess__coord_sphere__x_center[1]	= 0.0
AHFinderDirect::initial_guess__coord_sphere__y_center[1]	= 0.0
AHFinderDirect::initial_guess__coord_sphere__z_center[1]	= 0.0
AHFinderDirect::initial_guess__coord_sphere__radius[1]	= 1.0
AHFinderDirect::which_surface_to_store_info[1]	= 0
AHFinderDirect::set_mask_for_individual_horizon[1]	= "no"
AHFinderDirect::reset_horizon_after_not_finding[1]	= "no"
AHFinderDirect::find_after_individual_time[1]	= 1000.0
AHFinderDirect::max_allowable_horizon_radius[1]	= 15.0

AHFinderDirect::origin_x[2]		= 0.0
AHFinderDirect::initial_guess__coord_sphere__x_center[2]	= 0.0
AHFinderDirect::initial_guess__coord_sphere__y_center[2]	= 0.0
AHFinderDirect::initial_guess__coord_sphere__z_center[2]	= 0.0
AHFinderDirect::initial_guess__coord_sphere__radius[2]	= 3.0
AHFinderDirect::which_surface_to_store_info[2]	= 1
AHFinderDirect::set_mask_for_individual_horizon[2]	= "no"
AHFinderDirect::reset_horizon_after_not_finding[2]	= "no"
AHFinderDirect::find_after_individual_time[2]	= 1000.0
AHFinderDirect::max_allowable_horizon_radius[2]	= 15.0

# new parameters suggested by Erik, also for stability at recovery
AHFinderDirect::reshape_while_moving		= "yes"
AHFinderDirect::predict_origin_movement		= "yes"

#--------QuasiLocalMeasures-----------------
#QuasiLocalMeasures::verbose		= yes
QuasiLocalMeasures::interpolator		= "Lagrange polynomial interpolation"
QuasiLocalMeasures::interpolator_options	= "order=4"
QuasiLocalMeasures::spatial_order		= 4

QuasiLocalMeasures::num_surfaces		= 6
QuasiLocalMeasures::surface_index	[0]	= 0
QuasiLocalMeasures::surface_index	[1]	= 1
QuasiLocalMeasures::surface_index	[2]	= 2
QuasiLocalMeasures::surface_index	[3]	= 3
QuasiLocalMeasures::surface_index	[4]	= 4
QuasiLocalMeasures::surface_index	[5]	= 5


#--------Outflow----------------------------
Outflow::num_detectors		= 4
Outflow::surface_index[0]		= 2
Outflow::surface_index[1]		= 3
Outflow::surface_index[2]		= 4
Outflow::surface_index[3]		= 5
Outflow::threshold_on_var		= "eninf"
Outflow::num_thresholds		= 1
Outflow::extra_variables		= "
    ADMBase::lapse
    HydroBase::rho
    HydroBase::vel
    HydroBase::Y_e
    HydroBase::entropy
    HydroBase::temperature
"
Outflow::output_2d_data		= "yes"
Outflow::out_format			= ".9g"

#--------CarpetMask--------------------------
CarpetMask::excluded_surface	[0]	= 0
CarpetMask::excluded_surface_factor	[0]	= 1
CarpetMask::excluded_surface	[1]	= 1
CarpetMask::excluded_surface_factor	[1]	= 1

# =============================================================================
# Checkpoint
# =============================================================================
CarpetIOHDF5::checkpoint		= "yes"
CarpetIOHDF5::use_reflevels_from_checkpoint	= "yes"

IOUtil::abort_on_io_errors		= "yes"
IOUtil::checkpoint_on_terminate		= "yes"
IOUtil::checkpoint_every		= 8192
IOUtil::checkpoint_keep		= 2
IOUtil::recover			= "autoprobe"
IOUtil::checkpoint_dir		= "./checkpoint"
IOUtil::recover_dir			= "./parent/checkpoint"

# =============================================================================
# Output
# =============================================================================
IOUtil::out_dir			= "./data/"
IOUtil::strict_io_parameter_check		= "yes"
IOUtil::parfile_write		= "copy"

CarpetIOBasic::outinfo_vars		= "
    ADMBase::lapse
    HydroBase::rho
"

CarpetIOScalar::outscalar_vars		= "
    ADMBase::lapse
    ADMBase::shift
    ADMBase::curv
    ADMBase::metric
    BNSAnalysis::dens_unbnd
    BNSAnalysis::dens_unbnd_bernoulli
    BNSAnalysis::dens_unbnd_garching
    HydroBase::rho
    HydroBase::vel
    HydroBase::w_lorentz
    HydroBase::entropy
    HydroBase::eps
    HydroBase::press
    HydroBase::temperature
    HydroBase::Y_e
    CTGConstraints::H
    THC_Core::dens
    THC_Core::densgain
    THC_Core::volform
    THC_LeakageBase::thc_leakage_eff_rates
    THC_LeakageBase::thc_leakage_luminosity
    THC_LeakageBase::thc_leakage_free_rates
    THC_LeakageBase::thc_leakage_opacity
    THC_LeakageBase::thc_leakage_optd
"

CarpetIOASCII::out0D_vars		= "
    ADMMass::ADMMass_Masses
    BNSTrackerGen::bns_positions
    Carpet::physical_time_per_hour
    QuasiLocalMeasures::qlm_scalars
    Volomnia::grid_volume
    Volomnia::cell_volume
"

CarpetIOASCII::out1d_vars		= "
    ADMBase::lapse
    ADMBase::shift
    ADMBase::curv
    ADMBase::metric
    BNSAnalysis::dens_unbnd
    BNSAnalysis::dens_unbnd_bernoulli
    BNSAnalysis::dens_unbnd_garching
    HydroBase::rho
    HydroBase::vel
    HydroBase::w_lorentz
    HydroBase::eps
    HydroBase::press
    HydroBase::entropy
    HydroBase::temperature
    HydroBase::Y_e
    CTGConstraints::H
    THC_Core::bitmask
    THC_Core::dens
    THC_Core::densgain
    THC_Core::csound
    THC_Core::scon
    THC_Core::tau
    THC_Core::volform
    THC_LeakageBase::thc_leakage_eff_rates
    THC_LeakageBase::thc_leakage_luminosity
    THC_LeakageBase::thc_leakage_free_rates
    THC_LeakageBase::thc_leakage_opacity
    THC_LeakageBase::thc_leakage_optd
"

CarpetIOASCII::out2d_vars		= "
    WorldTube::extracted_vars
"

CarpetIOHDF5::out2d_vars		= "
    ADMBase::lapse
    ADMBase::shift
    BNSAnalysis::dens_unbnd
    BNSAnalysis::s_phi
    HydroBase::rho
    HydroBase::vel
    HydroBase::w_lorentz
    HydroBase::eps
    HydroBase::press
    HydroBase::entropy
    HydroBase::temperature
    HydroBase::Y_e
    CTGConstraints::H
    THC_Core::bitmask
    THC_Core::volform
    THC_LeakageBase::thc_leakage_eff_rates
    THC_LeakageBase::thc_leakage_optd
"

CarpetIOHDF5::out_vars		= "
    ADMBase::lapse
    ADMBase::shift
    ADMBase::metric
    HydroBase::rho
    HydroBase::vel
    HydroBase::w_lorentz
    HydroBase::temperature
    HydroBase::Y_e
    THC_Core::volform
"

CarpetIOASCII::one_file_per_group		= "yes"
IO::out_single_precision		= "yes"

CarpetIOASCII::out1D_x		= "yes"
CarpetIOASCII::out1D_y		= "yes"
CarpetIOASCII::out1D_z		= "yes"
CarpetIOASCII::out1D_d		= "yes"

CarpetIOBasic::outinfo_every		= 256
CarpetIOScalar::outscalar_criterion		= "divisor"
CarpetIOASCII::out0d_criterion		= "divisor"
CarpetIOASCII::out0d_every		= 256
CarpetIOASCII::out1d_every		= 1024
CarpetIOASCII::out2d_every		= 256
CarpetIOScalar::outscalar_every		= 256
CarpetIOHDF5::out2d_every		= 2048
CarpetIOHDF5::out_every		= 16384

Outflow::compute_every		= 512
Multipole::out_every		= 256

# vim: set ft=sh tabstop=20 noexpandtab :
# /* -*- Mode: C; indent-tabs-mode: t; tab-width: 20 -*- */
