! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2015-2016,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!

!
! Interfaces to RRTMGP.
!  Two interfaces are implemented:
!  -- rrtmgp_sw_opt uses optical properties
!  (optical depth, single-scattering albedo, and phase function) directly
!  -- rrtmgp_sw_phys converts from a physical description to optical properties
!  before calling the other interface
!
! The interfaces are overloaded via rrtmgp_sw()
!
! For both interfaces information is supplied as derived types following the convention
!   ty_TYPE from module mo_TYPE.
!   -- gas_desc encapsulates gas concentrations provided by users
!      concentrations should be set before calling rrtmgp_sw()
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
!   (this means the argument list to rrtmgp_sw() need not change).
!   Fluxes by g-point are available within the code. Most applications will expect some
!   reduced version of this detail e.g. broadband flux profiles.
! This implementation provides broadband fluxes but users are free to add others in
!   rrtmgp_sw_opt()
!
!
module mo_rrtmgp_sw
  use mo_rrtmgp_kind,   only: wp
  use mo_gas_optics_specification, &
                        only: ty_gas_optics_specification
  use mo_gas_concentrations, & 
                        only: ty_gas_concs
  use mo_optical_props, only: ty_optical_props, ty_optical_props_2str, ty_optical_props_nstr
  use mo_fluxes,        only: ty_fluxes
  use mo_sw_solver,     only: sw_solver
  ! Following modules are needed only by rrtmgp_sw_phys
  use mo_rng
  use mo_cloud_optics
  use mo_aerosol_optics
  implicit none
  private

  interface rrtmgp_sw
    module procedure rrtmgp_sw_opt, rrtmgp_sw_phys
  end interface rrtmgp_sw

  public :: rrtmgp_sw_init, rrtmgp_sw

  !
  ! Configuration information
  !
  integer :: nstreams = 2 ! Number of streams: 0 means no scattering
                          ! Require even, positive value when setting

contains
  ! --------------------------------------------------
  ! Initialization function
  ! --------------------------------------------------
  function rrtmgp_sw_init(nswstreams) result(error_msg)
    integer,           optional, intent( in) :: nswstreams
    character(len=128)                       :: error_msg

    error_msg = ""
    if(present(nswstreams)) then
      if(nswstreams >= 0) then  ! Check for an even number?
        nstreams = nswstreams
      else
        error_msg = "rrtmgp_sw_init: nswstreams provided is less than 0"
      end if
    end if
  end function rrtmgp_sw_init

  ! --------------------------------------------------
  ! Public interfaces to rrtmgp_sw()
  ! --------------------------------------------------
  !
  ! Physical description of cloud, surface, potentially aerosols
  !
  function rrtmgp_sw_phys(k_dist, gas_concs, p_lay, t_lay, p_lev, &
                          mu0, sfc_alb_dir, sfc_alb_dif,          &
                          clouds, rngs, allsky_fluxes, clrsky_fluxes, &
                          aerosols, col_dry, inc_flux) result(error_msg)
    type(ty_gas_optics_specification), &
                              intent(in ) :: k_dist  !< derived type with spectral information
    type(ty_gas_concs),        intent(in ) :: gas_concs !< derived type encapsulating gas concentrations
    real(wp), dimension(:,:), intent(in ) :: p_lay     !< pressure at layer centers     [Pa] (ncol, nlay)
    real(wp), dimension(:,:), intent(in ) :: t_lay     !< temperature at layer centers  [K]  (ncol, nlay)
    real(wp), dimension(:,:), intent(in ) :: p_lev     !< pressure at levels/interfaces [Pa] (ncol, nlay+1)
    real(wp), dimension(:),   intent(in ) :: mu0       !< cosine of solar zenith angle [] (ncol,nband) 
    real(wp), dimension(:,:), intent(in ) :: sfc_alb_dir, &  !< surface albedo for direct and diffuse radiation 
                                             sfc_alb_dif     !< [] (nband, ncol) 
    class(ty_cloud_desc),     intent(in ) :: clouds    !< cloud physical description
    class(ty_rng), &
        dimension(:),       intent(inout) :: rngs      !< random number generator states
    class(ty_fluxes),       intent(inout) :: allsky_fluxes, clrsky_fluxes
    class(ty_aerosol_desc),  &
              optional,       intent(in ) :: aerosols !< aerosol physical description
    real(wp), dimension(:,:), target, &
              optional,       intent(in ) :: col_dry, & !< Molecular number density (ncol, nlay)
                                             inc_flux   !< incident flux at domain top [W/m2] (ncol, ngpts)
    character(len=128)                    :: error_msg
    !---------------------------------
    ! Local variables

    integer :: ncol, nlay, ngpt, nband
    ! Variables to represent clouds, aerosol optical props
    class(ty_optical_props), allocatable :: cloud_props, aer_props
    !---------------------------------

    error_msg = ""
    if(.not. k_dist%is_initialized()) then
      error_msg = "rrtmgp_sw: k-distribution isn't initialized"
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
        error_msg =  cloud_props%init_2stream(ncol, nlay, ngpt)
      class is (ty_optical_props_nstr) ! n-stream calculation
        error_msg =  cloud_props%init_nstream(nstreams/2, ncol, nlay, ngpt)
    end select
    if(error_msg /= '') return

    !
    ! No need to allocate memory for aerosols if they aren't provided.
    !
    if(present(aerosols)) then
      allocate(aer_props, source=cloud_props) 
      select type (aer_props)
        class is (ty_optical_props)      ! No scattering
          error_msg =  aer_props%init_1scalar(ncol, nlay, nband)
        class is (ty_optical_props_2str) ! two-stream calculation
          error_msg =  aer_props%init_2stream(ncol, nlay, nband)
        class is (ty_optical_props_nstr) ! n-stream calculation
          error_msg =  aer_props%init_nstream(nstreams/2, ncol, nlay, nband)
      end select
    end if
    if(error_msg /= '') return

    !
    ! Map physical properties to optical properties
    !
    error_msg = clouds%sample_and_optics(rngs, k_dist, cloud_props)
    if(len_trim(error_msg) > 0) return
    if(present(aerosols)) call aerosols%optics(k_dist, aer_props)

    if(present(aerosols)) then
      error_msg = rrtmgp_sw_opt(k_dist, gas_concs, p_lay, t_lay, p_lev, & 
                                mu0, sfc_alb_dir, sfc_alb_dif, cloud_props, & 
                                allsky_fluxes, clrsky_fluxes,           &
                                aer_props, col_dry, inc_flux)
    else
      error_msg = rrtmgp_sw_opt(k_dist, gas_concs, p_lay, t_lay, p_lev, & 
                                mu0, sfc_alb_dir, sfc_alb_dif, cloud_props, & 
                                allsky_fluxes, clrsky_fluxes,           &
                                col_dry = col_dry, inc_flux=inc_flux)
    end if
  end function rrtmgp_sw_phys
  ! --------------------------------------------------
  !
  ! Optical description of atmosphere and surface
  ! All state variables except gas amounts have been converted to optical properties
  !   This function does error checking on array sizes before calling the core function
  !
  function rrtmgp_sw_opt(k_dist, gas_concs, p_lay, t_lay, p_lev, &
                         mu0, sfc_alb_dir, sfc_alb_dif, cloud_props, &
                         allsky_fluxes, clrsky_fluxes,           &
                         aer_props, col_dry, inc_flux) result(error_msg)
    type(ty_gas_optics_specification), intent(in   ) :: k_dist       !< derived type with spectral information
    type(ty_gas_concs),                intent(in   ) :: gas_concs    !< derived type encapsulating gas concentrations
    real(wp), dimension(:,:),          intent(in   ) :: p_lay, t_lay !< pressure [Pa], temperature [K] at layer centers (ncol,nlay)
    real(wp), dimension(:,:),          intent(in   ) :: p_lev        !< pressure at levels/interfaces [Pa] (ncol,nlay+1)
    real(wp), dimension(:  ),          intent(in   ) :: mu0          !< cosine of solar zenith angle 
    real(wp), dimension(:,:),          intent(in   ) :: sfc_alb_dir, sfc_alb_dif  
                                                        !  surface albedo for direct and diffuse radiation (band, col) 
    class(ty_optical_props),           intent(in   ) :: cloud_props !< cloud optical properties (ncol,nlay,ngpt)

    class(ty_fluxes),                  intent(inout) :: allsky_fluxes, clrsky_fluxes

    ! Optional inputs
    class(ty_optical_props),  &
              optional,       intent(in ) :: aer_props !< aerosol optical properties (ncol,nlay,nband?)
    real(wp), dimension(:,:), target, &
              optional,       intent(in ) :: col_dry, & !< Molecular number density (ncol, nlay)
                                             inc_flux   !< incident flux at domain top [W/m2] (ncol, ngpts)

    character(len=128)                    :: error_msg

    integer :: ncol, nlay, ngpt, nband
    ! --------------------------------
    error_msg = ""
    !
    ! Has the k-distribution been initialized?
    !
    if(.not. k_dist%is_initialized()) then
      error_msg = "rrtmgp_sw: k-distribution isn't initialized"
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
                  size(mu0,   dim=1) == ncol,   &
                  size(sfc_alb_dir, dim=1) == nband, &
                  size(sfc_alb_dir, dim=2) == ncol,  & 
                  size(sfc_alb_dif, dim=1) == nband, &
                  size(sfc_alb_dif, dim=2) == ncol])) then
      error_msg = "rrtmgp_sw: input arrays are inconsistently sized"
      return
    end if
    if(.not. all([size(cloud_props%tau, dim=1) == ncol,   &
                  size(cloud_props%tau, dim=2) == nlay,   &
                  size(cloud_props%tau, dim=3) == ngpt])) then
      error_msg = "rrtmgp_sw: cloud optical properties sizes don't match input arrays"
      return
    end if

    ! Optional arguments
    if(present(aer_props)) then
      if(.not. all([size(aer_props%tau, dim=1) == ncol,   &
                    size(aer_props%tau, dim=2) == nlay,   &
                    size(aer_props%tau, dim=3) == nband])) then
        error_msg = "rrtmgp_sw: aerosol optical properties sizes don't match input arrays"
        return
      end if
    end if

    if(present(col_dry)) then
      if(.not. all([size(col_dry, dim=1) == ncol,   &
                    size(col_dry, dim=2) == nlay])) then
        error_msg = "rrtmgp_sw: array col_dry is inconsistently sized"
        return
      end if
    end if

    if(present(inc_flux)) then
      if(.not. all([size(inc_flux, dim=1) == ncol,   &
                    size(inc_flux, dim=2) == ngpt])) then
        error_msg = "rrtmgp_sw: array inc_flux is inconsistently sized"
        return
      end if
    end if
    error_msg = &
      rrtmgp_sw_core(ncol, nlay, nband, ngpt,                    &
                     k_dist, gas_concs, p_lay, t_lay, p_lev,     & 
                     mu0, sfc_alb_dir, sfc_alb_dif, cloud_props, &
                     allsky_fluxes, clrsky_fluxes,               &
                     aer_props, col_dry, inc_flux)

  end function rrtmgp_sw_opt
  ! --------------------------------------------------
  ! Private interface to rrtmpg_sw()
  !   Inputs are free of errors and sizes are known.
  ! --------------------------------------------------

  function rrtmgp_sw_core(ncol, nlay, nband, ngpt,                    &
                          k_dist, gas_concs, p_lay, t_lay, p_lev,     & 
                          mu0, sfc_alb_dir, sfc_alb_dif, cloud_props, &
                          allsky_fluxes, clrsky_fluxes,               &
                          aer_props, col_dry, inc_flux) result(error_msg)
    !
    ! All state variables except gas amounts have been converted to optical properties
    !
    integer,                           intent(in   ) :: ncol, nlay, nband, ngpt
    type(ty_gas_optics_specification), intent(in   ) :: k_dist  !< derived type with spectral information
    type(ty_gas_concs),                 intent(in   ) :: gas_concs !< derived type encapsulating gas concentrations
    real(wp), dimension(ncol,nlay  ),  intent(in   ) :: p_lay, t_lay !< pressure [Pa], temperature [K] at layer centers
    real(wp), dimension(ncol,nlay+1),  intent(in   ) :: p_lev        !< pressure at levels/interfaces [Pa]
    real(wp), dimension(ncol),         intent(in   ) :: mu0          !< cosine of solar zenith angle 
    real(wp), dimension(nband, ncol),  intent(in   ) :: sfc_alb_dir, sfc_alb_dif !< surface albedo for direct and diffuse radiation 
    class(ty_optical_props),           intent(in   ) :: cloud_props  !< cloud optical properties (ncol,nlay,ngpt)

    class(ty_fluxes),                  intent(inout) :: allsky_fluxes, clrsky_fluxes

    class(ty_optical_props),  &
              optional,       intent(in ) :: aer_props !< aerosol optical properties (ncol,nlay,nband)
    real(wp), dimension(ncol,nlay),   target, &
              optional,       intent(in ) :: col_dry   !< Molecular number density (ncol, nlay)
    real(wp), dimension(ncol,ngpt  ), target, &
              optional,       intent(in ) :: inc_flux   !< incident flux at domain top [W/m2] (ncol, ngpts)

    character(len=128)                    :: error_msg
    ! --------------------------------
    !
    ! Local variables
    !
    integer :: icol, ilay
    class(ty_optical_props), allocatable :: all_sky   !

    REAL(wp), DIMENSION(ncol,nlay+1,ngpt) :: gpt_flux_up, gpt_flux_dn, gpt_flux_dir
                                                        ! Fluxes by gpoint  [W/m2]
    real(wp), dimension(ncol,       ngpt) :: inc_flux_def 
                                                        ! Incident fluxes by g-point at TOA [W/m2] 
                                                        ! May be overridden by user-provided values 
    REAL(wp), DIMENSION(ncol,       ngpt) :: sfc_alb_dir_gpt, sfc_alb_dif_gpt
    LOGICAL :: top_at_1
    integer :: top_lev
!++dbg
    character(len=*), parameter :: sub='rrtmgp_sw_core'
!--dbg
    ! ------------------------------------------------------------------------------------
    ! Executable statements start here
    ! True if layer 1 is the top of model (pressure increases with layer index)
    error_msg = ""

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
        error_msg = all_sky%init_2stream(ncol, nlay, ngpt)
      class is (ty_optical_props_nstr)
        error_msg = all_sky%init_nstream(nstreams/2, ncol, nlay, ngpt)
    end select
!++dbg
    if (error_msg /= '') then
       error_msg = sub//': init all_sky: '//trim(error_msg)
       return
    end if
!--dbg
    
    ! ----------------------------------------------------
    ! Optical properties of atmospheric components
    !
    ! Gas optical depth -- pressure need to be expressed as hPa 
    error_msg = k_dist%gas_optics(p_lay/100._wp, p_lev/100._wp, t_lay, gas_concs,  &
                                  all_sky, inc_flux_def,                           &
                                  col_dry)
!++dbg
!    if (error_msg /= '') return
    if (error_msg /= '') then
       error_msg = sub//': k_dist%gas_optics: '//trim(error_msg)
       return
    end if

!--dbg
    ! ----------------------------------------------------
    ! Boundary conditions
    ! Upper boundary condition
    top_at_1 = p_lay(1, 1) < p_lay(1, nlay)
    top_lev = MERGE(1, nlay+1, top_at_1)
    
    if(present(inc_flux)) then 
      gpt_flux_dir(:,top_lev,:) = inc_flux(:,:)  
    else 
      gpt_flux_dir(:,top_lev,:) = inc_flux_def(:,:) 
    end if 

!++dbg, apply solar zenith angle
    gpt_flux_dir(:,top_lev,:) = gpt_flux_dir(:,top_lev,:) * spread(mu0(:), dim = 2, ncopies = ngpt)
!--dbg

    ! Should allow for diffuse upper BC for generality
    ! 
    gpt_flux_dn(:,top_lev,:) = 0._wp 

    ! Lower boundary condition -- expand surface albedos by band to gpoints 
    !   and switch dimension ordering
    do icol = 1, ncol
      sfc_alb_dir_gpt(icol, 1:ngpt) = k_dist%expand(sfc_alb_dir(:,icol))
    end do
    do icol = 1, ncol
      sfc_alb_dif_gpt(icol, 1:ngpt) = k_dist%expand(sfc_alb_dif(:,icol))
    end do

    ! ----------------------------------------------------
    ! Clear sky is gases + aerosols (if they're supplied)
    if(present(aer_props)) & 
      error_msg = all_sky%increment_by(aer_props, k_dist%get_band_gpoint_limits()) 
!++dbg
!    if(error_msg /= "") return 
    if (error_msg /= '') then
       error_msg = sub//': all_sky%increment_by(aer_props): '//trim(error_msg)
       return
    end if
!--dbg
    
    if(clrsky_fluxes%are_desired()) then
      error_msg = sw_solver(ncol, nlay, ngpt, top_at_1,              &
                            all_sky, mu0, sfc_alb_dir_gpt, sfc_alb_dif_gpt,  &
                            gpt_flux_up, gpt_flux_dn, gpt_flux_dir)

!++dbg
!      if(len_trim(error_msg) > 0) return
    if (error_msg /= '') then
       error_msg = sub//': clear sky sw_solver: '//trim(error_msg)
       return
    end if
!--dbg
      ! Reduce spectral fluxes to desired output quantities
      error_msg = clrsky_fluxes%reduce(gpt_flux_up, gpt_flux_dn, k_dist, top_at_1, gpt_flux_dir)
!++dbg
!      if(len_trim(error_msg) > 0) return
    if (error_msg /= '') then
       error_msg = sub//': clear sky reduce: '//trim(error_msg)
       return
    end if
!--dbg
    end if
    ! ----------------------------------------------------
    ! All sky is clear sky plus clouds
    ! all_sky = all_sky + cloud_props
    if(allsky_fluxes%are_desired()) then
      error_msg = all_sky%increment_by(cloud_props)
      if(error_msg /= "") return 
      error_msg = sw_solver(ncol, nlay, ngpt, top_at_1,              &
                            all_sky, mu0, sfc_alb_dir_gpt, sfc_alb_dif_gpt,  &
                            gpt_flux_up, gpt_flux_dn, gpt_flux_dir)
      if(len_trim(error_msg) > 0) return
      ! Reduce spectral fluxes to desired output quantities
      error_msg = allsky_fluxes%reduce(gpt_flux_up, gpt_flux_dn, k_dist, top_at_1, gpt_flux_dir)
      if(len_trim(error_msg) > 0 ) return
    end if
  end function rrtmgp_sw_core

end module mo_rrtmgp_sw
