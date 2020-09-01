! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2015,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!

!
! Interfaces to RRTMGP.
!  Two interfaces are implemented:
!  -- rrtmgp_lw_opt uses optical properties
!  (optical depth, single-scattering albedo, and phase function) directly
!  -- rrtmgp_lw_phys converts from a physical description to optical properties
!  before calling the other interface
!
! The interfaces are overloaded via rrtmgp_lw()
!
! For both interfaces information is supplied as derived types following the convention
!   ty_TYPE from module mo_TYPE.
!   -- gas_desc encapsulates gas concentrations provided by users
!      concentrations should be set before calling rrtmgp_lw()
!   -- gas_optics_specification describes the spectral discretization;
!      it's read from a file using a load() procedure.
!
! If optical properties are provided these are packaged together using
!   -- optical_props; variables are expected for clouds and optional for aerosols
!
! If physical properties are provided users must provide two (or three, if using aerosols) derived types
!   -- rng: encapsulates the random number state or seed (depending on the implementation)
!   -- cloud_desc: encapsulates the physical description of clouds
!   -- aerosol_desc: encapsulates the physical description of aerosols
!
! Each of the types has procedures associated with them, and defined in the abstract interface,
!    that RRTMGP will call
!
! Output is via type ty_fluxes defined in this module
!   (this means the argument list to rrtmgp_lw() need not change).
!   Fluxes by g-point are available within the code. Most applications will expect some
!   reduced version of this detail e.g. broadband flux profiles.
! This implementation provides broadband fluxes but users are free to add others in
!   rrtmgp_lw_opt()
!
!
module mo_rrtmgp_lw
  use mo_rrtmgp_kind,   only: wp
  use mo_gas_optics_specification, &
                        only: ty_gas_optics_specification
  use mo_gas_concentrations, & 
                        only: ty_gas_concs
  use mo_optical_props, only: ty_optical_props, ty_optical_props_2str, ty_optical_props_nstr
  use mo_fluxes,        only: ty_fluxes
  use mo_lw_solver,     only: lw_solver_init, lw_solver
  ! Following modules are needed only by rrtmgp_lw_phys
  use mo_rng
  use mo_cloud_optics
  use mo_aerosol_optics
  implicit none
  private

  interface rrtmgp_lw
    module procedure rrtmgp_lw_opt, rrtmgp_lw_phys
  end interface rrtmgp_lw

  public :: rrtmgp_lw_init, rrtmgp_lw

  !
  ! Configuration information
  !
  integer :: nstreams = 0 ! Number of streams: 0 means no scattering
                          ! Require even, positive value when setting

contains
  ! --------------------------------------------------
  ! Initialization function
  ! --------------------------------------------------
  function rrtmgp_lw_init(nlwstreams, nangles) result(error_msg)
    integer,           optional, intent( in) :: nlwstreams ! Scattering/no scattering
    integer,           optional, intent( in) :: nangles    ! number of quadrature angles for 
                                                           ! no-scattering calculation
    character(len=128)                       :: error_msg

    error_msg = ""
    if(present(nlwstreams)) then
      if(nlwstreams >= 0) then  ! Check for an even number?
        nstreams = nlwstreams
      else
        error_msg = "rrtmgp_lw_init: nlwstreams provided is less than 0"
        return 
      end if
    end if
    if(present(nangles)) then
      error_msg = lw_solver_init(n_angles = nangles)
      if(error_msg /= "") return 
    end if
  end function rrtmgp_lw_init

  ! --------------------------------------------------
  ! Public interfaces to rrtmgp_lw()
  ! --------------------------------------------------
  !
  ! Physical description of cloud, surface, potentially aerosols
  !
  function rrtmgp_lw_phys(k_dist, gas_concs, p_lay, t_lay, p_lev, t_sfc, emis_sfc,   &
                          clouds, rngs, allsky_fluxes, clrsky_fluxes,                &
                          aerosols, col_dry, t_lev, inc_flux) result(error_msg)
    type(ty_gas_optics_specification), &
                              intent(in ) :: k_dist    !< derived type with spectral information
    type(ty_gas_concs),        intent(in ) :: gas_concs !< derived type encapsulating gas concentrations
    real(wp), dimension(:,:), intent(in ) :: p_lay     !< pressure at layer centers     [Pa] (ncol, nlay)
    real(wp), dimension(:,:), intent(in ) :: t_lay     !< temperature at layer centers  [K]  (ncol, nlay)
    real(wp), dimension(:,:), intent(in ) :: p_lev     !< pressure at levels/interfaces [Pa] (ncol, nlay+1)
    real(wp), dimension(:),   intent(in ) :: t_sfc     !< surface temperature   [K] (ncol)
    real(wp), dimension(:,:), intent(in ) :: emis_sfc  !< emissivity at surface []  (nbands, ncol)
    class(ty_cloud_desc),     intent(in ) :: clouds    !< cloud physical description
    class(ty_rng), &
        dimension(:),       intent(inout) :: rngs      !< random number generator states
    class(ty_fluxes),       intent(inout) :: allsky_fluxes, clrsky_fluxes
    class(ty_aerosol_desc),  &
              optional,       intent(in ) :: aerosols  !< aerosol physical description
    real(wp), dimension(:,:), target, &
              optional,       intent(in ) :: col_dry   !< Molecular number density 
    real(wp), dimension(:,:), target, &
              optional,       intent(in ) :: t_lev     !< temperature at levels [K] (ncol, nlay+1)
    real(wp), dimension(:,:),         &
              optional,       intent(in ) :: inc_flux   !< incident flux at domain top [W/m2] (ncol, ngpts)
    character(len=128)                    :: error_msg
    !---------------------------------
    ! Local variables

    integer :: ncol, nlay, ngpt, nband
    ! Variables to represent clouds, aerosol optical props
    class(ty_optical_props), allocatable :: cloud_props, aer_props
    !---------------------------------

    if(.not. k_dist%is_initialized()) then
      error_msg = "rrtmgp_lw: k-distribution isn't initialized"
      return
    end if

    ! Problem sizes
    ncol = size(p_lay,dim=1)
    nlay = size(p_lay,dim=2)
    ngpt = k_dist%get_ngpt()
    nband = k_dist%get_nband()

    !
    ! Initialize optical properties objects: is LW scattering to be included?
    !
    select case(nstreams)
      case(0)       ! No scattering
        allocate(ty_optical_props::cloud_props)
      case(2)       ! two-stream calculation
        allocate(ty_optical_props_2str::cloud_props)
      case default  ! n-stream calculation
        allocate(ty_optical_props_nstr::cloud_props)
    end select

    select type (cloud_props)
      class is (ty_optical_props)      ! No scattering
        error_msg = cloud_props%init_1scalar(ncol, nlay, ngpt)
      class is (ty_optical_props_2str) ! two-stream calculation
        error_msg = cloud_props%init_2stream(ncol, nlay, ngpt)
      class is (ty_optical_props_nstr) ! n-stream calculation
        error_msg = cloud_props%init_nstream(nstreams/2, ncol, nlay, ngpt)
    end select
    if (error_msg /= '') return

    !
    ! No need to allocate memory for aerosols if they aren't provided.
    !
    if(present(aerosols)) then
      allocate(aer_props, source=cloud_props) 
      select type (aer_props)
        class is (ty_optical_props)      ! No scattering
          error_msg = aer_props%init_1scalar(ncol, nlay, nband)
        class is (ty_optical_props_2str) ! two-stream calculation
          error_msg = aer_props%init_2stream(ncol, nlay, nband)
        class is (ty_optical_props_nstr) ! n-stream calculation
          error_msg = aer_props%init_nstream(nstreams/2, ncol, nlay, nband)
      end select
    end if
    if (error_msg /= '') return

    !
    ! Map physical properties to optical properties
    !
    error_msg = clouds%sample_and_optics(rngs, k_dist, cloud_props)
    if(error_msg /= '') return
    if(present(aerosols)) call aerosols%optics(k_dist, aer_props)

    if(present(aerosols)) then
      error_msg = rrtmgp_lw_opt(k_dist, gas_concs, p_lay, t_lay, p_lev, t_sfc, emis_sfc, &
                                cloud_props, allsky_fluxes, clrsky_fluxes,               &
                                aer_props, col_dry, t_lev, inc_flux)
    else
      error_msg = rrtmgp_lw_opt(k_dist, gas_concs, p_lay, t_lay, p_lev, t_sfc, emis_sfc, &
                                cloud_props, allsky_fluxes, clrsky_fluxes,               &
                                col_dry = col_dry,       & 
                                t_lev=t_lev,             &
                                inc_flux=inc_flux)
    end if
  end function rrtmgp_lw_phys
  ! --------------------------------------------------
  !
  ! Optical description of atmosphere and surface
  ! All state variables except gas amounts have been converted to optical properties
  !   This function does error checking on array sizes before calling the core function
  !
  function rrtmgp_lw_opt(k_dist, gas_concs, p_lay, t_lay, p_lev, &
                         t_sfc, emis_sfc, cloud_props,           &
                         allsky_fluxes, clrsky_fluxes,           &
                         aer_props, col_dry, t_lev, inc_flux) result(error_msg)
    type(ty_gas_optics_specification), intent(in   ) :: k_dist    !< derived type with spectral information
    type(ty_gas_concs),                 intent(in   ) :: gas_concs !< derived type encapsulating gas concentrations
    real(wp), dimension(:,:),  intent(in   ) :: p_lay     !< pressure at layer centers     [Pa] (ncol,nlay)
    real(wp), dimension(:,:),  intent(in   ) :: t_lay     !< temperature at layer centers  [K]  (ncol,nlay)
    real(wp), dimension(:,:),  intent(in   ) :: p_lev     !< pressure at levels/interfaces [Pa] (ncol,nlay+1)
    real(wp), dimension(:),    intent(in   ) :: t_sfc     !< surface temperature           [K]  (ncol)
    real(wp), dimension(:,:),  intent(in   ) :: emis_sfc  !< emissivity at surface         []   (nband, ncol)
    class(ty_optical_props),   intent(in   ) :: cloud_props !< cloud optical properties (ncol,nlay,ngpt)

    class(ty_fluxes),          intent(inout) :: allsky_fluxes, clrsky_fluxes

    ! Optional inputs
    class(ty_optical_props),  &
              optional,       intent(in ) :: aer_props !< aerosol optical properties (ncol,nlay,nband?)
    real(wp), dimension(:,:), target, &
              optional,       intent(in ) :: col_dry   !< Molecular number density (ncol, nlay)
    real(wp), dimension(:,:), target, &
              optional,       intent(in ) :: t_lev     !< temperature at levels [K] (ncol, nlay+1)
    real(wp), dimension(:,:), target, &
              optional,       intent(in ) :: inc_flux   !< incident flux at domain top [W/m2] (ncol, ngpts)

    character(len=128)                    :: error_msg

    ! --------------------------------
    integer :: ncol, nlay, ngpt, nband
    ! --------------------------------
    error_msg = ""
    !
    ! Has the k-distribution been initialized?
    !
    if(.not. k_dist%is_initialized()) then
      error_msg = "rrtmgp_lw: k-distribution isn't initialized"
      return
    end if

    ! Problem sizes
    ncol = size(p_lay,dim=1)
    nlay = size(p_lay,dim=2)
    ngpt = k_dist%get_ngpt()
    nband = k_dist%get_nband()

    !
    ! Check array sizes: input arrays and also those within derived types
    !   Wait until the end to fail to catch as many errors as possible.
    !
    if(.not. all([size(t_lay, dim=1) == ncol,   &
                  size(t_lay, dim=2) == nlay,   &
                  size(p_lev, dim=1) == ncol,   &
                  size(p_lev, dim=2) == nlay+1, &
                  size(t_sfc       ) == ncol,   &
                  size(emis_sfc, dim=1) == nband, &
                  size(emis_sfc, dim=2) == ncol])) then
      error_msg = "rrtmgp_lw: input arrays are inconsistently sized"
      return
    end if
    if(.not. all([size(cloud_props%tau, dim=1) == ncol,   &
                  size(cloud_props%tau, dim=2) == nlay,   &
                  size(cloud_props%tau, dim=3) == ngpt])) then
      error_msg = "rrtmgp_lw: cloud optical properties sizes don't match input arrays"
      return
    end if

    ! Optional arguments
    if(present(aer_props)) then
      if(.not. all([size(aer_props%tau, dim=1) == ncol,   &
                    size(aer_props%tau, dim=2) == nlay,   &
                    size(aer_props%tau, dim=3) == nband])) then
        error_msg = "rrtmgp_lw: aerosol optical properties sizes don't match input arrays"
        return
      end if
    end if

    if(present(col_dry)) then
      if(.not. all([size(col_dry, dim=1) == ncol,   &
                    size(col_dry, dim=2) == nlay])) then
        error_msg = "rrtmgp_lw: array col_dry is inconsistently sized"
        return
      end if
    end if

    if(present(t_lev)) then
      if(.not. all([size(t_lev, dim=1) == ncol,   &
                    size(t_lev, dim=2) == nlay+1])) then
        error_msg = "rrtmgp_lw: array t_lev is inconsistently sized"
        return
      end if
    end if

    if(present(inc_flux)) then
      if(.not. all([size(inc_flux, dim=1) == ncol,   &
                    size(inc_flux, dim=2) == ngpt])) then
        error_msg = "rrtmgp_lw: array inc_flux is inconsistently sized"
        return
      end if
    end if
    if(len_trim(error_msg) > 0) return

    error_msg = &
      rrtmgp_lw_core(ncol, nlay, nband, ngpt,                                              &
                     k_dist, gas_concs, p_lay, t_lay, p_lev, t_sfc, emis_sfc, cloud_props, &
                     allsky_fluxes, clrsky_fluxes,                                         &
                     aer_props, col_dry, t_lev, inc_flux)

  end function rrtmgp_lw_opt
  ! --------------------------------------------------
  ! Private interface to rrtmpg_lw()
  !   Inputs are free of errors and sizes are known.
  ! --------------------------------------------------
  function rrtmgp_lw_core(ncol, nlay, nband, ngpt,                                              &
                          k_dist, gas_concs, p_lay, t_lay, p_lev, t_sfc, emis_sfc, cloud_props, &
                          allsky_fluxes, clrsky_fluxes,                                         &
                          aer_props, col_dry, t_lev, inc_flux) result(error_msg)
    !
    ! All state variables except gas amounts have been converted to optical properties
    !
    integer,                           intent(in   ) :: ncol, nlay, nband, ngpt
    type(ty_gas_optics_specification), intent(in   ) :: k_dist  !< derived type with spectral information
    type(ty_gas_concs),                 intent(in   ) :: gas_concs !< derived type encapsulating gas concentrations
    real(wp), dimension(ncol,nlay  ),  intent(in   ) :: p_lay     !< pressure at layer centers     [Pa]
    real(wp), dimension(ncol,nlay  ),  intent(in   ) :: t_lay     !< temperature at layer centers  [K]
    real(wp), dimension(ncol,nlay+1),  intent(in   ) :: p_lev     !< pressure at levels/interfaces [Pa]
    real(wp), dimension(ncol       ),  intent(in   ) :: t_sfc     !< surface temperature   [K]
    real(wp), dimension(nband,ncol ),  intent(in   ) :: emis_sfc  !< emissivity at surface []
    class(ty_optical_props),           intent(in   ) :: cloud_props !< cloud optical properties (ncol,nlay,ngpt)

    class(ty_fluxes),                  intent(inout) :: allsky_fluxes, clrsky_fluxes

    class(ty_optical_props),  &
              optional,       intent(in ) :: aer_props !< aerosol optical properties (ncol,nlay,nband)
    real(wp), dimension(ncol,nlay),   target, &
              optional,       intent(in ) :: col_dry   !< Molecular number density (ncol, nlay)
    real(wp), dimension(ncol,nlay+1), target, &
              optional,       intent(in ) :: t_lev     !< temperature at levels [K]
    real(wp), dimension(ncol,ngpt  ), target, &
              optional,       intent(in ) :: inc_flux   !< incident flux at domain top [W/m2] (ncol, ngpts)

    character(len=128)                    :: error_msg
    ! --------------------------------
    !
    ! Local variables
    !
    integer :: icol, ilay
    class(ty_optical_props), allocatable :: all_sky   !

    REAL(wp), DIMENSION(ncol,nlay+1,ngpt) :: gpt_flux_up, gpt_flux_dn
                                                        ! Fluxes by gpoint  [W/m2]
    REAL(wp), DIMENSION(ncol,nlay,  ngpt) :: lay_source
                                                        ! Planck source at layer average temperature
                                                        ! [W/m2] (ncol, nlay, ngpt)
    REAL(wp), DIMENSION(ncol,nlay+1,ngpt) :: lev_source_inc, lev_source_dec
                                      ! Planck source at layer edge for radiation, [W/m2] (ncol, nlay+1, ngpt)
                                      !   in increasing/decreasing ilay direction
                                      ! Includes spectral weighting that accounts for state-dependent frequency to g-space mapping
    REAL(wp), DIMENSION(ncol,       ngpt) :: sfc_source
                                      ! Planck source for radiation from surface, [W/m2] (ncol, ngpt)
                                      ! Includes spectral weighting 
    REAL(wp), DIMENSION(ncol,       ngpt) :: emis_sfc_gpt

    LOGICAL :: top_at_1
    integer :: top_lev
    ! ------------------------------------------------------------------------------------
    ! Executable statements start here
    ! True if layer 1 is the top of model (pressure increases with layer index)
    top_at_1 = p_lay(1, 1) < p_lay(1, nlay)
    top_lev = MERGE(1, nlay+1, top_at_1)
    
    select case(nstreams)
      case(0)
        allocate(ty_optical_props::all_sky)
      case(2)
        allocate(ty_optical_props_2str::all_sky)
      case default
        allocate(ty_optical_props_nstr::all_sky)
    end select

    select type (all_sky)
      class is (ty_optical_props) ! No scattering
        error_msg = all_sky%init_1scalar(ncol, nlay, ngpt)
      class is (ty_optical_props_2str)
        error_msg =  all_sky%init_2stream(ncol, nlay, ngpt)
      class is (ty_optical_props_nstr)
        error_msg =  all_sky%init_nstream(nstreams/2, ncol, nlay, ngpt)
    end select
    if (error_msg /= '') return
    ! ----------------------------------------------------
    ! Optical properties of atmospheric components
    !
    ! Gas optical depth -- pressure need to be expressed as hPa 

    error_msg = k_dist%gas_optics(p_lay/100._wp, p_lev/100._wp, t_lay, t_sfc, gas_concs,  &
                                  all_sky,                                                &
                                  lay_source, lev_source_inc, lev_source_dec, sfc_source, &
                                  col_dry, t_lev)
    if (error_msg /= '') return
    ! ----------------------------------------------------
    ! Boundary conditions
    ! Upper boundary condition
    if(present(inc_flux)) then 
      gpt_flux_dn(:,top_lev,:) = inc_flux(:,:)  
    else 
      gpt_flux_dn(:,top_lev,:) = 0._wp 
    end if 

    ! Lower boundary condition -- expand surface emissivity by band to gpoints
    do icol = 1, ncol
      emis_sfc_gpt(icol, 1:ngpt) = k_dist%expand(emis_sfc(:,icol))
    end do

    ! ----------------------------------------------------
    ! Clear sky is gases + aerosols (if they're supplied)
    if(present(aer_props)) & 
       error_msg = all_sky%increment_by(aer_props, k_dist%get_band_gpoint_limits()) 
    if (error_msg /= '') return

    if(clrsky_fluxes%are_desired()) then
      error_msg = lw_solver(ncol, nlay, ngpt, top_at_1,                          &
                            all_sky, lay_source, lev_source_inc, lev_source_dec, &
                            emis_sfc_gpt, sfc_source, gpt_flux_up, gpt_flux_dn)
      if(len_trim(error_msg) > 0) return
      ! Reduce spectral fluxes to desired output quantities
      error_msg = clrsky_fluxes%reduce(gpt_flux_up, gpt_flux_dn, k_dist, top_at_1)
      if(len_trim(error_msg) > 0) return
    end if

    ! ----------------------------------------------------
    ! All sky is clear sky plus clouds
    ! all_sky = all_sky + cloud_props
    if(allsky_fluxes%are_desired()) then
      error_msg = all_sky%increment_by(cloud_props)
      if (error_msg /= '') return
      error_msg = lw_solver(ncol, nlay, ngpt, top_at_1,                          &
                            all_sky, lay_source, lev_source_inc, lev_source_dec, &
                            emis_sfc_gpt, sfc_source, gpt_flux_up, gpt_flux_dn)
      if (error_msg /= '') return

      ! Reduce spectral fluxes to desired output quantities
      error_msg = allsky_fluxes%reduce(gpt_flux_up, gpt_flux_dn, k_dist, top_at_1)
      if(len_trim(error_msg) > 0 ) return
    end if
  end function rrtmgp_lw_core

end module mo_rrtmgp_lw
