! Module: mo_gas_optics_specification

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
! Description:  Specifies all properties related to the k-distributions.  This includes each band's spectral properties and
! absorbing gases.

module mo_gas_optics_specification
  use mo_rrtmgp_kind,        only: wp
  use mo_rrtmgp_constants,   only: avogad, pi, m_dry, m_h2o, grav
  use mo_gas_optics_kernels, only: interpolation, gas_optical_depths_major, & 
                                   gas_optical_depths_continuum, gas_optical_depths_minor, & 
                                   gas_optical_depths_rayleigh, source
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_optical_props,      only: ty_optical_props, ty_optical_props_2str, ty_optical_props_nstr
  use mo_util_reorder
  implicit none
  private

  ! -----------------------------------------------------------------------------------
  type, public :: ty_gas_optics_specification
    private
    character(32), & 
              dimension(:),   allocatable :: gas_names  ! gas names

    integer,  dimension(:,:), allocatable :: flavor        ! major species pair; (2,nflav)
    integer,  dimension(:,:), allocatable :: gpoint_flavor ! flavor = gpoint_flavor(lower or upper atmosphere, g-point)

    integer,  dimension(:,:), allocatable :: band2gpt       ! (begin g-point, end g-point) = band2gpt(2,band)
    integer,  dimension(:),   allocatable :: gpt2band       ! band = gpt2band(g-point)

    real(wp), dimension(:,:), allocatable :: band_lims_wavenum  ! (upper and lower wavenumber by band) = band_lims_wavenum(2,band)
    ! -----------------------------------------------------------------------------------
    ! Temperature and pressure interpolation grids 
    real(wp), dimension(:),  allocatable :: press_ref,  press_ref_log, temp_ref
    ! volume mixing ratios for reference atmosphere; vmr_ref(lower or upper atmosphere, gas, temp)
    real(wp), dimension(:,:,:), allocatable :: vmr_ref

    real(wp) :: press_ref_min, press_ref_max, &  ! min and max pressure of interpolation grids
                temp_ref_min,  temp_ref_max      ! min and max temperature 
    real(wp) :: press_ref_log_delta, & ! difference in ln pressure between consecutive reference levels
                temp_ref_delta,      & ! Temperature difference between consecutive reference levels
                press_ref_trop_log     ! log of reference pressure separating the lower and upper atmosphere
    real(wp) :: stpfac                 ! standard pressure:temperature ratio
    ! -----------------------------------------------------------------------------------
    ! Absorption coefficients
      ! ----- major gas absorption coefficients ; kmajor(g-point,eta,pressure,temperature)
    real(wp), dimension(:,:,:,:), allocatable :: kmajor
    ! ----- water vapor continuum
      ! stored absorption coefficients due to water vapor self continuum; selfrefin(eta,temperature,g-point)
    real(wp), dimension(:,:,:),   allocatable :: selfrefin, forrefin
    ! ----- minor species
      ! stored absorption coefficients due to minor absorbing gases in lower/upper part of atmosphere;
      ! kminor_lower(minor_gas,g-point,eta,temperature)
    real(wp), dimension(:,:,:,:), allocatable :: kminor_lower, kminor_upper
    ! -----------------------------------------------------------------------------------
    ! ----- Rayleigh scattering
      ! stored scattering coefficients due to molecules in atmosphere;
      ! krayl(g-point,eta,temperature,upper/lower atmosphere)
    real(wp), dimension(:,:,:,:), allocatable :: krayl
    ! -----------------------------------------------------------------------------------
    ! Planck function spectral mapping
    !   Allocated only when gas optics object is internal-source "flavor" 
    !
    real(wp), dimension(:,:,:,:), allocatable :: planck_frac   ! input Planck fractions
                                                               ! planck_frac(eta,temperature,pressure,g-point)
    real(wp), dimension(:,:),     allocatable :: totplnk       ! integrated Planck function by band; (reference temperatures,band)
    real(wp)                                  :: totplnk_delta ! temperature steps in totplnk
    ! -----------------------------------------------------------------------------------
    ! Solar source function spectral mapping 
    !   Allocated only when gas optics object is external-source "flavor" 
    !
    real(wp), dimension(:), allocatable :: solar_src ! solar source
    ! -----------------------------------------------------------------------------------
    ! Ancillary
    ! computed by fill_kminor_activity(); g-points where minor gases are active; (g-point, minor gas) = kminor_activity(2,:)
    integer, dimension(:,:), allocatable :: kminor_activity
    ! -----------------------------------------------------------------------------------

  contains
    ! Type-bound procedures
    ! Public procedures
    ! public interface
    generic, public :: init       => init_int,       init_ext
    generic, public :: gas_optics => gas_optics_int, gas_optics_ext
    procedure, public :: is_initialized
    procedure, public :: is_internal_source_present
    procedure, public :: is_external_source_present
    procedure, public :: weight_bandvals_by_gpoint
    procedure, public :: expand
    procedure, public :: get_ngas
    procedure, public :: get_gases
    procedure, public :: get_nband
    procedure, public :: get_ngpt
    procedure, public :: get_press_ref_min
    procedure, public :: get_press_ref_max
    procedure, public :: get_temp_ref_min
    procedure, public :: get_temp_ref_max
    procedure, public :: get_band_gpoint_limits
    procedure, public :: convert_band2gpt
    procedure, public :: convert_gpt2band
    procedure, public :: get_band_lims_wavenumber
    procedure, public :: get_band_lims_wavelength
    ! Internal procedures
    procedure, public :: init_int
    procedure, public :: init_ext
    procedure, public :: gas_optics_int
    procedure, public :: gas_optics_ext
    procedure, private :: fill_kminor_activity
    procedure, private :: get_nflav
    procedure, private :: get_nlay_ref
    procedure, private :: get_neta
    procedure, private :: gas2id
    procedure, private :: compute_gas_tau_core
  end type
  ! -----------------------------------------------------------------------------------
  public :: get_col_dry ! Utility function, not type-bound

  interface check_range
    module procedure check_range_1D, check_range_2D, check_range_3D 
  end interface check_range

  interface check_extent
    module procedure check_extent_1D, check_extent_2D, check_extent_3D 
  end interface check_extent
contains
  ! --------------------------------------------------------------------------------------
  !
  ! Public procedures
  !
  ! --------------------------------------------------------------------------------------
  !
  ! Two functions to define array sizes needed by gas_optics() 
  !
  pure function get_ngas(this)
    ! return the number of gases registered in the spectral configuration
    class(ty_gas_optics_specification), intent(in) :: this
    integer                                        :: get_ngas

    get_ngas = size(this%gas_names)
  end function get_ngas
  !--------------------------------------------------------------------------------------------------------------------
  pure function get_nflav(this)
    ! return the number of major species pairs
    class(ty_gas_optics_specification), intent(in) :: this
    integer                                        :: get_nflav

    get_nflav = size(this%flavor,dim=2)
  end function get_nflav
  !--------------------------------------------------------------------------------------------------------------------

  ! Compute gas optical depth and, optionally, Planck source functions,
  !  given temperature, pressure, and composition
  function gas_optics_int(this,                                   &
                      play, plev, tlay, tsfc, gas_desc,           & ! mandatory inputs
                      optical_props,                              & ! mandatory outputs
                      lay_src, lev_src_inc, lev_src_dec, sfc_src, & ! internal-source specific outputs
                      col_dry, tlev)                              & ! optional inputs
                      result(error_msg)
    ! inputs
    class(ty_gas_optics_specification), intent(in) :: this
    real(wp), dimension(:,:), intent(in   ) :: play, &   ! layer pressures [hPa, mb]; (ncol,nlay)
                                               plev, &   ! level pressures [hPa, mb]; (ncol,nlay+1)
                                               tlay      ! layer temperatures [K]; (ncol,nlay)
    real(wp), dimension(:),   intent(in   ) :: tsfc      ! surface skin temperatures [K]; (ncol)
    type(ty_gas_concs),        intent(in   ) :: gas_desc  ! Gas volume mixing ratios
    ! output
    class(ty_optical_props),  intent(inout) :: optical_props
    character(len=128)                      :: error_msg
    ! source functions (LW only)
    ! These include spectral weighting that accounts for state-dependent frequency to k-distribution mapping
    ! [W/m2] 
    real(wp), dimension(:,:,:), intent(  out) :: lay_src, &  ! source for average layer temperature; (ncol,nlay,ngpt)
                                                 lev_src_inc, lev_src_dec  
                                                             ! level source radiances in increasing/decreasing 
                                                             ! ilay direction (ncol,nlay+1,ngpt)
    real(wp), dimension(:,:),   intent(  out) :: sfc_src     ! Surface Planck source; (ncol,ngpts)
    ! Optional inputs
    real(wp), dimension(:,:),   intent(in   ), &
                           optional, target :: col_dry, &  ! Column dry amount; dim(ncol,nlay)
                                               tlev        ! level temperatures [K]l (ncol,nlay+1)
    ! ----------------------------------------------------------
    ! Local variables
    ! Interpolation coefficients to save for use in source function 
    integer,  dimension(size(play,dim=1), size(play,dim=2)) :: jtemp, jpress
    logical,  dimension(size(play,dim=1), size(play,dim=2)) :: tropo
    real(wp), dimension(2,2,2,this%get_nflav(),size(play,dim=1), size(play,dim=2)) :: fmajor
    integer,  dimension(2,    this%get_nflav(),size(play,dim=1), size(play,dim=2)) :: jeta
    ! ----------------------------------------------------------
    ! dimensions - determined from problem size
    integer :: ncol, nlay ! number of columns, layers
    ! dimensions - provided by k-distribution
    integer :: ngpt, nband, ngas, nflav ! Number of g-points, bands, gas, gas "flavors" (major species combinations) 
    ! index
    integer :: icol, ilay, igpt, igas

    ! Variables for temperature at layer edges
    ! [K] (ncol, nlay+1)
    real(wp), dimension(size(play,dim=1),size(play,dim=2)+1), target  :: tlev_arr
    real(wp), dimension(:,:),                                 pointer :: tlev_wk => NULL()

    ! ----------------------------------------------------------
    ! Code starts
    !
    error_msg = ""
    error_msg = compute_gas_taus(this,                       &
                                 play, plev, tlay, gas_desc, & 
                                 optical_props,              & 
                                 jtemp, jpress, jeta, tropo, fmajor, & 
                                 col_dry)    
    if(error_msg  /= '') return

    ! init from array dimensions
    ncol = size(play,dim=1)
    nlay = size(play,dim=2)
    ngpt = this%get_ngpt()
    nband = this%get_nband()
    ngas = this%get_ngas()
    nflav = this%get_nflav()

    !
    ! Planck source function 
    !   Check input data sizes and values
    !
    error_msg = check_extent(tsfc, ncol, 'tsfc') 
    if(error_msg  /= '') return
    error_msg = check_range(tsfc, this%temp_ref_min,  this%temp_ref_max,  'tsfc')
    if(error_msg  /= '') return
    if(present(tlev)) then 
      error_msg = check_extent(tlev, ncol, nlay+1, 'tlev') 
      if(error_msg  /= '') return
      error_msg = check_range(tlev, this%temp_ref_min, this%temp_ref_max, 'tlev')
      if(error_msg  /= '') return
    end if
    
    !
    !   output extents 
    !
    error_msg = check_extent(sfc_src,     ncol,         ngpt, 'sfc_src') 
    if(error_msg  /= '') return
    error_msg = check_extent(lay_src,     ncol, nlay,   ngpt, 'lay_src') 
    if(error_msg  /= '') return
    error_msg = check_extent(lev_src_inc, ncol, nlay+1, ngpt, 'lev_src_inc') 
    if(error_msg  /= '') return
    error_msg = check_extent(lev_src_dec, ncol, nlay+1, ngpt, 'lev_src_dec') 
    if(error_msg  /= '') return

    !
    ! Source function needs temperature at interfaces/levels and at layer centers
    !
    if (present(tlev)) then
      !   Users might have provided these
      tlev_wk => tlev
    else
       tlev_wk => tlev_arr
       !
       ! Interpolate temperature to levels if not provided
       !   Interpolation and extrapolation at boundaries is weighted by pressure
       !
       do icol = 1, ncol
         tlev_arr(icol,1) = tlay(icol,1) &
                           + (plev(icol,1)-play(icol,1))*(tlay(icol,2)-tlay(icol,1))  &
              &                                           / (play(icol,2)-play(icol,1))
       end do
       do ilay = 2, nlay
         do icol = 1, ncol
           tlev_arr(icol,ilay) = (play(icol,ilay-1)*tlay(icol,ilay-1)*(plev(icol,ilay  )-play(icol,ilay)) &
                                +  play(icol,ilay  )*tlay(icol,ilay  )*(play(icol,ilay-1)-plev(icol,ilay))) /  &
                                  (plev(icol,ilay)*(play(icol,ilay-1) - play(icol,ilay)))
         end do
       end do
       do icol = 1, ncol
         tlev_arr(icol,nlay+1) = tlay(icol,nlay)                                                             &
                                + (plev(icol,nlay+1)-play(icol,nlay))*(tlay(icol,nlay)-tlay(icol,nlay-1))  &
                                                                      / (play(icol,nlay)-play(icol,nlay-1))
       end do
     end if

    call source(ncol, nlay, ngpt, nband, ngas, nflav, &
                tlay, tlev_wk, tsfc, merge(1,nlay,play(1,1) > play(1,nlay)), & 
                fmajor, jeta, tropo, jtemp, jpress,                    & 
                this%band2gpt, this%planck_frac, this%temp_ref_min,    & 
                this%totplnk_delta, this%totplnk, this%gpoint_flavor,  &  
                sfc_src, lay_src, lev_src_inc, lev_src_dec)
   
   
  end function gas_optics_int
  !------------------------------------------------------------------------------------------

  ! Compute gas optical depth
  !  given temperature, pressure, and composition
  function gas_optics_ext(this,        &
    play, plev, tlay, gas_desc,        & ! mandatory inputs
    optical_props, toa_src,            & ! mandatory outputs
    col_dry) result(error_msg)           ! optional input 
    
    class(ty_gas_optics_specification), intent(in) :: this
    real(wp), dimension(:,:), intent(in   ) :: play, &   ! layer pressures [hPa, mb]; (ncol,nlay)
                                               plev, &   ! level pressures [hPa, mb]; (ncol,nlay+1)
                                               tlay      ! layer temperatures [K]; (ncol,nlay)
    type(ty_gas_concs),        intent(in   ) :: gas_desc  ! Gas volume mixing ratios
    ! output
    class(ty_optical_props),  intent(inout) :: optical_props
    real(wp), dimension(:,:), intent(  out) :: toa_src     ! Top-of-atmosphere flux
    character(len=128)                      :: error_msg

    ! Optional inputs
    real(wp), dimension(:,:), intent(in   ), &
                           optional, target :: col_dry ! Column dry amount; dim(ncol,nlay)
    ! ----------------------------------------------------------
    ! Local variables
    ! Interpolation coefficients to save for use in source function 
    integer,  dimension(size(play,dim=1), size(play,dim=2)) :: jtemp, jpress
    logical,  dimension(size(play,dim=1), size(play,dim=2)) :: tropo
    real(wp), dimension(2,2,2,this%get_nflav(),size(play,dim=1), size(play,dim=2)) :: fmajor
    integer,  dimension(2,    this%get_nflav(),size(play,dim=1), size(play,dim=2)) :: jeta
    ! ----------------------------------------------------------
    integer :: ncol, ngpt 

    ! ----------------------------------------------------------
    ! Code starts
    !
    error_msg = ""
    error_msg = compute_gas_taus(this,                       &
                                 play, plev, tlay, gas_desc, & 
                                 optical_props,              & 
                                 jtemp, jpress, jeta, tropo, fmajor, & 
                                 col_dry)    
    if(error_msg  /= '') return

    ncol = size(play,dim=1)
    ngpt = this%get_ngpt()
    error_msg = check_extent(toa_src,     ncol,         ngpt, 'toa_src') 
    if(error_msg  /= '') return
    toa_src(:,:) = spread(this%solar_src(:), dim=1, ncopies=ncol) 
   
  end function gas_optics_ext
  !------------------------------------------------------------------------------------------
  !
  ! Returns optical properties and interpolation coefficients 
  !
  function compute_gas_taus(this,                       &
                            play, plev, tlay, gas_desc, & 
                            optical_props,              & 
                            jtemp, jpress, jeta, tropo, fmajor, & 
                            col_dry) result(error_msg)    
    class(ty_gas_optics_specification), & 
                                      intent(in   ) :: this
    real(wp), dimension(:,:),         intent(in   ) :: play, &   ! layer pressures [hPa, mb]; (ncol,nlay)
                                                       plev, &   ! level pressures [hPa, mb]; (ncol,nlay+1)
                                                       tlay      ! layer temperatures [K]; (ncol,nlay)
    type(ty_gas_concs),                intent(in   ) :: gas_desc  ! Gas volume mixing ratios

    class(ty_optical_props),          intent(inout) :: optical_props
    ! Interpolation coefficients  for use in source function 
    integer,  dimension(:,:),         intent(  out) :: jtemp, jpress
    integer,  dimension(:,:,:,:),     intent(  out) :: jeta
    logical,  dimension(:,:),         intent(  out) :: tropo
    real(wp), dimension(:,:,:,:,:,:), intent(  out) :: fmajor
    character(len=128)                            :: error_msg

    ! Optional inputs
    real(wp), dimension(:,:), intent(in   ), &
                           optional, target :: col_dry ! Column dry amount; dim(ncol,nlay)
    ! ----------------------------------------------------------
    ! Local variables
    ! gas amounts
    real(wp), dimension(size(play,dim=1), size(play,dim=2), this%get_ngas()) :: vmr     ! volume mixing ratios; (nlay,ncol,ngas)
    real(wp), dimension(size(play,dim=1), size(play,dim=2))                  :: one_vmr ! a single volume mixing ratio, (ncol, nlay)

    real(wp), dimension(size(optical_props%tau,dim=3), &
                        size(optical_props%tau,dim=2), &
                        size(optical_props%tau,dim=1)) :: tau  ! optical depth; (ngpt, nlay, ncol)
    real(wp), dimension(size(optical_props%tau,dim=3), &
                        size(optical_props%tau,dim=2), &
                        size(optical_props%tau,dim=1)) :: tau_rayleigh ! optical depth; (ngpt, nlay, ncol)

    ! ----------------------------------------------------------
    ! dimensions - determined from problem size
    integer :: ncol, nlay ! number of columns, layers
    ! dimensions - provided by k-distribution
    integer :: ngpt, nband, ngas, nflav ! Number of g-points, bands, gas, gas "flavors" (major species combinations) 
    ! index
    integer :: igas, idx_h2o
!++dbg
    integer :: icol
    character(len=*), parameter :: sub='compute_gas_taus'
!--dbg
    ! Number of molecules per cm^2
    real(wp), dimension(size(play,dim=1), size(play,dim=2)), target  :: col_dry_arr
    real(wp), dimension(:,:),                                pointer :: col_dry_wk => NULL()
    ! ----------------------------------------------------------
    ! Code starts
    !
    error_msg = ''
    ! Check for initialization
    if (.not. this%is_initialized()) then
      error_msg = 'ERROR: spectral configuration not loaded'
      return
    end if

    ! init from array dimensions
    ncol = size(play,dim=1)
    nlay = size(play,dim=2)
    ngpt = this%get_ngpt()
    nband = this%get_nband()
    ngas = this%get_ngas()
    nflav = this%get_nflav()
    idx_h2o = this%gas2id('h2o')

    !
    ! Check input data sizes and values
    !
    error_msg = check_extent(play, ncol, nlay,   'play') 
    if(error_msg  /= '') return
    error_msg = check_extent(plev, ncol, nlay+1, 'plev') 
    if(error_msg  /= '') return
    error_msg = check_extent(tlay, ncol, nlay,   'tlay') 
    if(error_msg  /= '') return

!++dbg, disable this check... plev is sufficient
!    error_msg = check_range(play, this%press_ref_min, this%press_ref_max, 'play')
!    if(error_msg  /= '') return
!--dbg

    error_msg = check_range(plev, this%press_ref_min, this%press_ref_max, 'plev')
!++dbg
!    if(error_msg  /= '') return
    if(error_msg  /= '') then
!       do icol = 1, ncol
!          if (plev(icol,1) > this%press_ref_max) then
!             print*,sub//': INFO: plev > press_ref_max=', plev(icol,1), this%press_ref_max
!          end if
!       end do
       print*,sub//': INFO: plev > press_ref_max'
       error_msg = ''
    end if
!--dbg

    error_msg = check_range(tlay, this%temp_ref_min,  this%temp_ref_max,  'tlay')
    if(error_msg  /= '') return
    if(present(col_dry)) then 
      error_msg = check_extent(col_dry, ncol, nlay, 'col_dry')
      if(error_msg  /= '') return
      error_msg = check_range(col_dry, 0._wp, huge(col_dry), 'col_dry')
      if(error_msg  /= '') return
    end if
    
    ! Expand volume mixing ratio to 3D fields
    do igas = 1, ngas
      error_msg = gas_desc%get_vmr(this%gas_names(igas),one_vmr)
      if (error_msg /= '') return
      vmr(:,:,igas) = one_vmr
    end do

    ! Compute column amounts (number of molecule per cm^2) if user hasn't provided them
    if (present(col_dry)) then
      col_dry_wk => col_dry
    else
      col_dry_arr = get_col_dry(vmr(:,:,idx_h2o), plev, tlay) ! column dry amounts computation
      col_dry_wk => col_dry_arr
    end if

    error_msg = this%compute_gas_tau_core(play, tlay, vmr, col_dry_wk, &
                                     ncol, nlay, ngpt, nband, ngas, nflav, &
                                     tau, tau_rayleigh, & 
                                     fmajor, jeta, tropo, jtemp, jpress)
    if (error_msg /= '') return

    call combine_and_reorder(tau, tau_rayleigh, allocated(this%krayl), optical_props) 
   
  end function compute_gas_taus
  !------------------------------------------------------------------------------------------
  function compute_gas_tau_core(this,            &
    play, tlay, vmr, col_dry,               & !  inputs
    ncol, nlay, ngpt, nband, ngas, nflav,   &
    tau, tau_rayleigh,                      & ! mandatory outputs
    fmajor_out, jeta_out, tropo_out, jtemp_out, jpress_out) result(error_msg)
    
    class(ty_gas_optics_specification), intent(in) :: this
    ! dimensions
    integer, intent(in) :: ncol  ! Number of columns
    integer, intent(in) :: nlay  ! Number of layers
    integer, intent(in) :: ngpt  ! Number of gpts
    integer, intent(in) :: nband ! Number of bands
    integer, intent(in) :: ngas  ! Number of gases
    integer, intent(in) :: nflav ! Number of gas flavors

    real(wp), dimension(ncol,nlay  ), intent(in) :: play   ! Layer pressures [hPa, mb]
    real(wp), dimension(ncol,nlay  ), intent(in) :: tlay   ! Layer temperatures [K]
    real(wp), dimension(ncol,nlay  ,ngas), &
                                      intent(in) :: vmr ! volume mixing ratios
    real(wp), dimension(ncol,nlay  ), intent(in) :: col_dry ! Column amount of dry air

    ! output
    real(wp), dimension(ngpt,nlay,ncol), intent(out) :: tau          ! optical depth, will be transposed
    real(wp), dimension(ngpt,nlay,ncol), intent(out) :: tau_rayleigh ! optical depth, will be transposed

    real(wp), dimension(2,2,2,nflav,ncol,nlay), optional, intent(out) :: fmajor_out 
    integer,  dimension(2,    nflav,ncol,nlay), optional, intent(out) :: jeta_out  
    logical,  dimension(            ncol,nlay), optional, intent(out) :: tropo_out
    integer,  dimension(            ncol,nlay), optional, intent(out) :: jtemp_out
    integer,  dimension(            ncol,nlay), optional, intent(out) :: jpress_out

    ! result
    character(len=128) :: error_msg
    ! ----------------------------------------------------------
    ! Local variables
    ! index
    integer :: igas
    ! Planck fractions
    ! gas amounts
    real(wp), dimension(ncol,nlay,ngas) :: col_gas ! column amounts for each gas

    ! temperature variables
    integer,  dimension(ncol,nlay) :: jtemp ! interpolation index for temperature
    ! pressure variables
    integer,  dimension(ncol,nlay) :: jpress ! interpolation index for pressure
    logical,  dimension(ncol,nlay) :: tropo ! true lower atmosphere; false upper atmosphere

    integer, dimension(2,     nflav,ncol,nlay) :: jeta ! interpolation index for binary species parameter (eta)
                                                     ! index(1) : reference temperature level
                                                     ! index(2) : flavor
                                                     ! index(3) : layer

    real(wp), dimension(2,    nflav,ncol,nlay) :: col_mix ! combination of major species's column amounts
                                                         ! index(1) : reference temperature level
                                                         ! index(2) : flavor
                                                         ! index(3) : layer

    real(wp), dimension(2,2,2,nflav,ncol,nlay) :: fmajor ! interpolation fractions for major species
                                                            ! index(1) : reference eta level (temperature dependent)
                                                            ! index(2) : reference pressure level
                                                            ! index(3) : reference temperature level
                                                            ! index(4) : flavor
                                                            ! index(5) : layer

    real(wp), dimension(2,2,  nflav,ncol,nlay) :: fminor ! interpolation fractions for minor species and continuum
                                                          ! index(1) : reference eta level (temperature dependent)
                                                          ! index(2) : reference temperature level
                                                          ! index(3) : flavor
                                                          ! index(4) : layer

    integer :: idx_h2o, idx_o2, idx_n2 ! index of some gases

    ! ----------------------------------------------------------
    ! Code starts
    !
    error_msg = ''

    ! special gases
    idx_h2o = this%gas2id('h2o')
    idx_o2 = this%gas2id('o2')
    idx_n2 = this%gas2id('n2')

    ! vmr and column gas amounts
    do igas = 1, ngas
      col_gas(:,:,igas) = vmr(:,:,igas) * col_dry(:,:)
    end do

    tau(:,:,:) = 0._wp
    ! ---- calculate gas optical depths ----
    call interpolation( &
      ncol,nlay,nflav,this%get_neta(), & ! dimensions
      this%flavor,this%press_ref_log,this%temp_ref,this%press_ref_log_delta,this%temp_ref_min, & ! inputs from object
      this%temp_ref_delta, this%press_ref_trop_log,this%vmr_ref,this%get_nlay_ref(), &
      play,tlay,col_gas, & ! local input
      jtemp,fmajor,fminor,col_mix,tropo,jeta,jpress) ! output 
      
    call gas_optical_depths_major( &
      ncol,nlay,ngpt,nflav, & ! dimensions
      this%gpoint_flavor,this%kmajor, & ! inputs from object
      col_mix,fmajor,& 
      jeta,tropo,jtemp,jpress, & ! local input
      tau)
      
    call gas_optical_depths_continuum( &
      ncol,nlay,ngpt,ngas,nflav, & ! dimensions
      this%flavor,this%gpoint_flavor,this%selfrefin,this%forrefin,this%stpfac, & ! inputs from object
      idx_h2o ,play,tlay,vmr,col_gas,& 
      fminor,jeta,tropo,jtemp, & ! local input
      tau)
      
    call gas_optical_depths_minor( &
      ncol,nlay,ngpt,ngas,nflav, & ! dimensions
      this%gpoint_flavor,this%band2gpt,this%kminor_lower,this%kminor_upper,this%kminor_activity, & ! inputs from object
      idx_h2o,idx_o2,idx_n2,play,tlay,col_dry,col_gas,& 
      fminor,jeta,tropo,jtemp, & ! local input
      tau)
      
    if (allocated(this%krayl)) then
      call gas_optical_depths_rayleigh( &
        ncol,nlay,ngpt,ngas,nflav, & ! dimensions
        this%gpoint_flavor,this%krayl, & ! inputs from object
        idx_h2o,idx_o2,idx_n2,play,tlay,col_dry,col_gas,& 
        fminor,jeta,tropo,jtemp, & ! local input
        tau_rayleigh)
    end if

    ! This is an internal function -- we can assume that all or none of these are present 
    if(present(fmajor_out)) then 
      fmajor_out = fmajor
      jeta_out   = jeta 
      tropo_out  = tropo 
      jtemp_out  = jtemp 
      jpress_out = jpress
    end if 
     
  end function compute_gas_tau_core
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Initialization 
  !
  !--------------------------------------------------------------------------------------------------------------------
  ! Initialize object based on data read from netCDF file however the user desires. 
  !  Rayleigh scattering tables may or may not be present; this is indicated with allocation status 
  ! This interface is for the internal-sources object -- includes Plank functions and fractions
  ! 
  function init_int(this, gas_names, key_species,        & 
                    band2gpt, band_lims_wavenum,            & 
                    press_ref, press_ref_trop, temp_ref, & 
                    temp_ref_p, temp_ref_t, vmr_ref,     & 
                    kmajor, selfrefin, forrefin, kminor_lower, kminor_upper, & 
                    totplnk, planck_frac, rayl_lower, rayl_upper) result(err_message) 
    class(ty_gas_optics_specification), intent(inout) :: this
    character(len=*), dimension(:), intent(in) :: gas_names
    integer,  dimension(:,:,:),   intent(in) :: key_species
    integer,  dimension(:,:),     intent(in) :: band2gpt 
    real(wp), dimension(:,:),     intent(in) :: band_lims_wavenum
    real(wp), dimension(:),       intent(in) :: press_ref, temp_ref
    real(wp),                     intent(in) :: press_ref_trop, temp_ref_p, temp_ref_t
    real(wp), dimension(:,:,:),   intent(in) :: vmr_ref
    real(wp), dimension(:,:,:,:), intent(in) :: kmajor
    real(wp), dimension(:,:,:),   intent(in) :: selfrefin, forrefin
    real(wp), dimension(:,:,:,:), intent(in) :: kminor_lower, kminor_upper
    real(wp), dimension(:,:),     intent(in) :: totplnk
    real(wp), dimension(:,:,:,:), intent(in) :: planck_frac
    real(wp), dimension(:,:,:),   intent(in), & 
                                 allocatable :: rayl_lower, rayl_upper
    character(len = 128) err_message
    ! ---- 
    err_message = init_abs_coeffs(this, & 
                                  gas_names, key_species,    & 
                                  band2gpt, band_lims_wavenum, &
                                  press_ref, temp_ref,       & 
                                  press_ref_trop, temp_ref_p, temp_ref_t, &
                                  vmr_ref,                   & 
                                  kmajor, selfrefin, forrefin, kminor_lower, kminor_upper, & 
                                  rayl_lower, rayl_upper) 
    ! Planck function tables 
    ! 
    this%totplnk = totplnk
    this%planck_frac = planck_frac
    ! Temperature steps for Planck function interpolation 
    !   Assumes that temperature minimum is the same on both scales and that Planck scale
    !   is equally spaced 
    this%totplnk_delta =  (this%temp_ref_max-this%temp_ref_min) / (size(this%totplnk,dim=1)-1)
  end function init_int

  !--------------------------------------------------------------------------------------------------------------------
  ! Initialize object based on data read from netCDF file however the user desires. 
  !  Rayleigh scattering tables may or may not be present; this is indicated with allocation status 
  ! This interface is for the external-sources object -- includes TOA source function table 
  ! 
  function init_ext(this, gas_names, key_species,        & 
                    band2gpt, band_lims_wavenum,           & 
                    press_ref, press_ref_trop, temp_ref, & 
                    temp_ref_p, temp_ref_t, vmr_ref,     & 
                    kmajor, selfrefin, forrefin, kminor_lower, kminor_upper, & 
                    solar_src, rayl_lower, rayl_upper)  result(err_message)
    class(ty_gas_optics_specification), intent(inout) :: this
    character(len=*), & 
              dimension(:),       intent(in) :: gas_names
    integer,  dimension(:,:,:),   intent(in) :: key_species
    integer,  dimension(:,:),     intent(in) :: band2gpt 
    real(wp), dimension(:,:),     intent(in) :: band_lims_wavenum
    real(wp), dimension(:),       intent(in) :: press_ref, temp_ref
    real(wp),                     intent(in) :: press_ref_trop, temp_ref_p, temp_ref_t
    real(wp), dimension(:,:,:),   intent(in) :: vmr_ref
    real(wp), dimension(:,:,:,:), intent(in) :: kmajor
    real(wp), dimension(:,:,:),   intent(in) :: selfrefin, forrefin
    real(wp), dimension(:,:,:,:), intent(in) :: kminor_lower, kminor_upper
    real(wp), dimension(:),       intent(in), allocatable :: solar_src 
                                                            ! allocatable status to change when solar source is present in file 
    real(wp), dimension(:,:,:), intent(in), allocatable :: rayl_lower, rayl_upper
    character(len = 128) err_message
    ! ---- 
    err_message = init_abs_coeffs(this, & 
                                  gas_names, key_species,    & 
                                  band2gpt, band_lims_wavenum, &
                                  press_ref, temp_ref,       & 
                                  press_ref_trop, temp_ref_p, temp_ref_t, &
                                  vmr_ref,                   & 
                                  kmajor, selfrefin, forrefin, kminor_lower, kminor_upper, & 
                                  rayl_lower, rayl_upper) 
    !
    ! Something something solar source table init here 
    ! 
    this%solar_src = solar_src
    
  end function init_ext
  !--------------------------------------------------------------------------------------------------------------------
  ! Initialize absorption coefficient arrays,  
  !   including Rayleigh scattering tables if provided (allocated) 
  ! 
  function init_abs_coeffs(this, & 
                           gas_names, key_species,    & 
                           band2gpt, band_lims_wavenum, &
                           press_ref, temp_ref,       & 
                           press_ref_trop, temp_ref_p, temp_ref_t, &
                           vmr_ref,                   & 
                           kmajor, selfrefin, forrefin, kminor_lower, kminor_upper, & 
                           rayl_lower, rayl_upper) result(err_message) 
    class(ty_gas_optics_specification), intent(inout) :: this
    character(len=*), & 
              dimension(:),       intent(in) :: gas_names
    integer,  dimension(:,:,:),   intent(in) :: key_species
    integer,  dimension(:,:),     intent(in) :: band2gpt 
    real(wp), dimension(:,:),     intent(in) :: band_lims_wavenum
    real(wp), dimension(:),       intent(in) :: press_ref, temp_ref
    real(wp),                     intent(in) :: press_ref_trop, temp_ref_p, temp_ref_t
    real(wp), dimension(:,:,:),   intent(in) :: vmr_ref
    real(wp), dimension(:,:,:,:), intent(in) :: kmajor
    real(wp), dimension(:,:,:),   intent(in) :: selfrefin, forrefin
    real(wp), dimension(:,:,:,:), intent(in) :: kminor_lower, kminor_upper
    
    real(wp), dimension(:,:,:),   intent(in), & 
                                 allocatable :: rayl_lower, rayl_upper
    character(len=128)                       :: err_message 
    ! --------------------------------------
    err_message = "" 
    ! Assignment 
    !   includes allocation 
    this%gas_names = gas_names
    this%press_ref = press_ref
    this%temp_ref  = temp_ref 
    this%vmr_ref   = vmr_ref
    this%kmajor       = kmajor
    this%selfrefin    = selfrefin
    this%forrefin     = forrefin 
    this%kminor_lower = kminor_lower
    this%kminor_upper = kminor_upper
    this%band2gpt        = band2gpt       
    this%band_lims_wavenum = band_lims_wavenum

    if(allocated(rayl_lower) .neqv. allocated(rayl_upper)) then 
      err_message = "rayl_lower and rayl_upper must have the same allocation status"
      return 
    end if 
    if (allocated(rayl_lower)) then
      allocate(this%krayl(size(rayl_lower,dim=1),size(rayl_lower,dim=2),size(rayl_lower,dim=3),2))
      this%krayl(:,:,:,1) = rayl_lower
      this%krayl(:,:,:,2) = rayl_upper
    end if


    ! ---- post processing ----

    this%press_ref(:) = this%press_ref(:) * 0.01_wp ! convert reference pressure from hPa to Pa
    ! creates log reference pressure
    allocate(this%press_ref_log(size(this%press_ref)))
    this%press_ref_log(:) = log(this%press_ref(:))

    ! log scale of reference pressure
    this%press_ref_trop_log = log(press_ref_trop)
    ! factor needed for continuum optical depth
    this%stpfac = temp_ref_t/temp_ref_p*100._wp

    ! create flavor list
    call create_flavor(key_species, this%flavor)
    ! create gpt2band
    call create_gpt2band(this%band2gpt, this%gpt2band)
    ! create gpoint_flavor list
    call create_gpoint_flavor(key_species, this%gpt2band, this%flavor, this%gpoint_flavor)

    ! minimum, maximum reference temperature, pressure -- assumes low-to-high ordering 
    !   for T, high-to-low ordering for p 
    this%temp_ref_min  = this%temp_ref (1)
    this%temp_ref_max  = this%temp_ref (size(this%temp_ref))
    this%press_ref_min = this%press_ref(size(this%press_ref))
    this%press_ref_max = this%press_ref(1)

    ! creates press_ref_log, temp_ref_delta
    this%press_ref_log_delta = (log(this%press_ref_min)-log(this%press_ref_max))/(size(this%press_ref)-1)
    this%temp_ref_delta      = (this%temp_ref_max-this%temp_ref_min)/(size(this%temp_ref)-1)

    ! fills kminor_activity; a list where minor species are active
    call this%fill_kminor_activity()
  end function init_abs_coeffs
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Inquiry functions 
  !
  !--------------------------------------------------------------------------------------------------------------------
  ! return true if initialized, false otherwise
  pure function is_initialized(this)
    class(ty_gas_optics_specification), intent(in) :: this
    logical                                        :: is_initialized
    is_initialized = allocated(this%gas_names)
  end function is_initialized
  !--------------------------------------------------------------------------------------------------------------------

  ! return true if initialized for internal sources, false otherwise
  pure function is_internal_source_present(this)
    class(ty_gas_optics_specification), intent(in) :: this
    logical                                        :: is_internal_source_present
    is_internal_source_present = allocated(this%totplnk).and.allocated(this%planck_frac)
  end function is_internal_source_present
  !--------------------------------------------------------------------------------------------------------------------

  ! return true if initialized for external sources, false otherwise
  pure function is_external_source_present(this)
    class(ty_gas_optics_specification), intent(in) :: this
    logical                                        :: is_external_source_present
    is_external_source_present = allocated(this%solar_src)
  end function is_external_source_present

  !--------------------------------------------------------------------------------------------------------------------
  ! return the gas names
  pure function get_gases(this)
    class(ty_gas_optics_specification), intent(in) :: this
    character(32), dimension(this%get_ngas())     :: get_gases

    get_gases = this%gas_names
  end function get_gases
  !--------------------------------------------------------------------------------------------------------------------
  ! return the number of bands
  pure function get_nband(this)
    class(ty_gas_optics_specification), intent(in) :: this
    integer                                        :: get_nband

    get_nband = size(this%band2gpt,dim=2)
  end function get_nband

  !--------------------------------------------------------------------------------------------------------------------
  ! return the number of g-points
  pure function get_ngpt(this)
    class(ty_gas_optics_specification), intent(in) :: this
    integer                                        :: get_ngpt

    get_ngpt = size(this%gpoint_flavor,dim=2)
  end function get_ngpt

  !--------------------------------------------------------------------------------------------------------------------
  ! return the minimum pressure on the interpolation grids
  pure function get_press_ref_min(this)
    class(ty_gas_optics_specification), intent(in) :: this
    real(wp)                                       :: get_press_ref_min

    get_press_ref_min = this%press_ref_min
  end function get_press_ref_min

  !--------------------------------------------------------------------------------------------------------------------
  ! return the maximum pressure on the interpolation grids
  pure function get_press_ref_max(this)
    class(ty_gas_optics_specification), intent(in) :: this
    real(wp)                                       :: get_press_ref_max

    get_press_ref_max = this%press_ref_max
  end function get_press_ref_max

  !--------------------------------------------------------------------------------------------------------------------
  ! return the minimum temparature on the interpolation grids
  pure function get_temp_ref_min(this)
    class(ty_gas_optics_specification), intent(in) :: this
    real(wp)                                       :: get_temp_ref_min

    get_temp_ref_min = this%temp_ref_min
  end function get_temp_ref_min

  !--------------------------------------------------------------------------------------------------------------------
  ! return the maximum temparature on the interpolation grids
  pure function get_temp_ref_max(this)
    class(ty_gas_optics_specification), intent(in) :: this
    real(wp)                                       :: get_temp_ref_max

    get_temp_ref_max = this%temp_ref_max
  end function get_temp_ref_max

  !--------------------------------------------------------------------------------------------------------------------
  ! return the first and last g-point of all the bands at once
  ! dimension (2, nbands) 
  ! 
  pure function get_band_gpoint_limits(this)
    class(ty_gas_optics_specification), intent(in) :: this
    integer, dimension(size(this%band2gpt,dim=1), size(this%band2gpt,dim=2)) & 
                                                   :: get_band_gpoint_limits

    get_band_gpoint_limits = this%band2gpt
  end function get_band_gpoint_limits

  !--------------------------------------------------------------------------------------------------------------------
  ! return the first and last g-point of a band
  pure function convert_band2gpt(this, band)
    class(ty_gas_optics_specification), intent(in) :: this
    integer,                            intent(in) :: band
    integer, dimension(2)                         :: convert_band2gpt

    convert_band2gpt(:) = this%band2gpt(:,band)
  end function convert_band2gpt

  !--------------------------------------------------------------------------------------------------------------------
  ! return the lower and upper wavenumber of a band
  ! (upper and lower wavenumber by band) = band_lims_wavenum(2,band)
  pure function get_band_lims_wavenumber(this)
    class(ty_gas_optics_specification), intent(in) :: this
    real(wp), dimension(size(this%band_lims_wavenum,1), size(this%band_lims_wavenum,2)) & 
                                                   :: get_band_lims_wavenumber

    get_band_lims_wavenumber(:,:) = this%band_lims_wavenum(:,:)
  end function get_band_lims_wavenumber

  !--------------------------------------------------------------------------------------------------------------------
  ! return the lower and upper wavelength of a band
  pure function get_band_lims_wavelength(this)
    class(ty_gas_optics_specification), intent(in) :: this
    real(wp), dimension(size(this%band_lims_wavenum,1), size(this%band_lims_wavenum,2)) & 
                                                   :: get_band_lims_wavelength

    get_band_lims_wavelength(:,:) = 1._wp/this%band_lims_wavenum(:,:)
  end function get_band_lims_wavelength

  !--------------------------------------------------------------------------------------------------------------------
  ! return the band of a g-point
  pure function convert_gpt2band(this, gpt)
    class(ty_gas_optics_specification), intent(in) :: this
    integer,                            intent(in) :: gpt
    integer                                        :: convert_gpt2band

    convert_gpt2band = this%gpt2band(gpt)
  end function convert_gpt2band

  !--------------------------------------------------------------------------------------------------------------------
  ! The result is the band_value multiplied by the corresponding g-point weight.
  ! We use knowledge of the starting and ending g-point for each band to know the values of i and j in
  ! weight_bandvals_by_gpoint(j) = gpt_weights(j) * band_vals(i)
  pure function weight_bandvals_by_gpoint(this, gpt_weights)
    class(ty_gas_optics_specification),     intent(in) :: this
    real(wp), dimension(:),                 intent(in) :: gpt_weights ! dim(# of g-points)
    real(wp), dimension(size(this%gpoint_flavor,dim=2)) :: weight_bandvals_by_gpoint

    integer :: iband, igpt

    do iband=1,this%get_nband()
      do igpt=this%band2gpt(1,iband), this%band2gpt(2,iband)
        weight_bandvals_by_gpoint = gpt_weights(igpt) * iband
      end do
    end do
  end function weight_bandvals_by_gpoint

  !--------------------------------------------------------------------------------------------------------------------
  ! expands an array of dimension arr_in(nband) to dimension arr_out(ngpt)
  pure function expand(this, arr_in) result(arr_out)
    class(ty_gas_optics_specification), intent(in) :: this
    real(wp), dimension(:), intent(in) :: arr_in ! (nband)
    real(wp), dimension(size(this%gpoint_flavor,dim=2)) :: arr_out
    integer :: iband
    do iband=1,this%get_nband()
      arr_out(this%band2gpt(1,iband):this%band2gpt(2,iband)) = arr_in(iband)
    end do
  end function expand


  !--------------------------------------------------------------------------------------------------------------------
  ! --- gas optical depth calculations
  !--------------------------------------------------------------------------------------------------------------------
  ! Utility function, provided for user convenience
  ! computes column amounts of dry air using hydrostatic equation
  function get_col_dry(vmr_h2o, plev, tlay, latitude) result(col_dry)
    ! input
    real(wp), dimension(:,:), intent(in) :: vmr_h2o  ! volume mixing ratio of all gases excluding water; (ncol,nlay)
    real(wp), dimension(:,:), intent(in) :: plev     ! Layer boundary pressures [hPa, mb] (ncol,nlay+1)
    real(wp), dimension(:,:), intent(in) :: tlay     ! Layer temperatures [K] (ncol,nlay)
    real(wp), dimension(:),   optional, &
                              intent(in) :: latitude ! Latitude [degrees] (ncol)
    ! output
    real(wp), dimension(size(tlay,dim=1),size(tlay,dim=2)) :: col_dry ! Column dry amount (ncol,nlay)
    ! ------------------------------------------------
    ! first and second term of Helmert formula
    real(wp), parameter :: helmert1 = 9.80665_wp
    real(wp), parameter :: helmert2 = 0.02586_wp
    ! local variables
    real(wp), dimension(size(tlay,dim=1)                 ) :: g0 ! (ncol)
    real(wp), dimension(size(tlay,dim=1),size(tlay,dim=2)) :: delta_plev ! (ncol,nlay)
    real(wp), dimension(size(tlay,dim=1),size(tlay,dim=2)) :: m_air ! average mass of air; (ncol,nlay)
    integer :: nlev, nlay
    ! ------------------------------------------------
    nlay = size(tlay, dim=2)
    nlev = size(plev, dim=2)

    if(present(latitude)) then
      g0(:) = helmert1 - helmert2 * cos(2.0_wp * pi * latitude(:) / 180.0_wp) ! acceleration due to gravity [m/s^2]
    else
      g0(:) = grav
    end if
    delta_plev(:,:) = plev(:,1:nlev-1) - plev(:,2:nlev)

    ! Get average mass of air
    m_air(:,:) = (m_dry+m_h2o*vmr_h2o(:,:))/(1.+vmr_h2o(:,:))

    ! Hydrostatic equation
    col_dry(:,:) = 1000._wp*delta_plev(:,:)*avogad/(1000._wp*m_air(:,:)*100._wp*spread(g0(:),dim=2,ncopies=nlay))
    col_dry(:,:) = col_dry(:,:)/(1._wp+vmr_h2o(:,:))
  end function get_col_dry
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Internal procedures 
  ! 
  !--------------------------------------------------------------------------------------------------------------------
  pure function rewrite_key_species_pair(key_species_pair)
    ! (x,0) becomes (x,x)
    ! (0,0) becomes (2,2) -- because absorption coefficients for these g-points will be 0. 
    integer, dimension(2) :: rewrite_key_species_pair
    integer, dimension(2), intent(in) :: key_species_pair
    rewrite_key_species_pair = key_species_pair
    if (key_species_pair(2).eq.0) then
      rewrite_key_species_pair(2) = key_species_pair(1)
    end if
    if (all(key_species_pair(:).eq.(/0,0/))) then
      rewrite_key_species_pair(:) = (/2,2/)
    end if
  end function

  ! ---------------------------------------------------------------------------------------
  ! true is key_species_pair exists in key_species_list
  pure function key_species_pair_exists(key_species_list, key_species_pair)
    logical :: key_species_pair_exists
    integer, dimension(:,:), intent(in) :: key_species_list
    integer, dimension(2), intent(in) :: key_species_pair
    integer :: i
    do i=1,size(key_species_list,dim=2)
      if (all(key_species_list(:,i).eq.key_species_pair(:))) then
        key_species_pair_exists = .true.
        return
      end if
    end do
    key_species_pair_exists = .false.
  end function key_species_pair_exists

  ! ---------------------------------------------------------------------------------------
  ! create flavor list -- 
  !   an unordered array of extent (2,:) containing all possible pairs of key species 
  !   used in either upper or lower atmos
  !
  subroutine create_flavor(key_species, flavor)
    integer, dimension(:,:,:), intent(in) :: key_species
    integer, dimension(:,:), allocatable, intent(out) :: flavor
    integer, dimension(2,size(key_species,3)*2) :: key_species_list

    integer :: ibnd, iatm, i, iflavor
    ! prepare list of key_species
    i = 1
    do ibnd=1,size(key_species,3)
      do iatm=1,size(key_species,1)
        key_species_list(:,i) = key_species(:,iatm,ibnd)
        i = i + 1
      end do
    end do
    ! rewrite single key_species pairs
    do i=1,size(key_species_list,2)
        key_species_list(:,i) = rewrite_key_species_pair(key_species_list(:,i))
    end do
    ! count unique key species pairs
    iflavor = 0
    do i=1,size(key_species_list,2)
      if (.not.key_species_pair_exists(key_species_list(:,1:i-1),key_species_list(:,i))) then
        iflavor = iflavor + 1
      end if
    end do
    ! fill flavors
    allocate(flavor(2,iflavor))
    iflavor = 0
    do i=1,size(key_species_list,2)
      if (.not.key_species_pair_exists(key_species_list(:,1:i-1),key_species_list(:,i))) then
        iflavor = iflavor + 1
        flavor(:,iflavor) = key_species_list(:,i)
      end if
    end do
  end subroutine create_flavor
! ---------------------------------------------------------------------------------------

  ! returns flavor index; -1 if not found
  pure function key_species_pair2flavor(flavor, key_species_pair)
    integer :: key_species_pair2flavor
    integer, dimension(:,:), intent(in) :: flavor
    integer, dimension(2), intent(in) :: key_species_pair
    integer :: iflav
    do iflav=1,size(flavor,2)
      if (all(key_species_pair(:).eq.flavor(:,iflav))) then
        key_species_pair2flavor = iflav
        return
      end if
    end do
    key_species_pair2flavor = -1
  end function key_species_pair2flavor

  ! ---------------------------------------------------------------------------------------
  ! create gpoint_flavor list
  !   a map pointing from each g-point to the corresponding entry in the "flavor list" 
  !
  subroutine create_gpoint_flavor(key_species, gpt2band, flavor, gpoint_flavor)
    integer, dimension(:,:,:), intent(in) :: key_species
    integer, dimension(:), intent(in) :: gpt2band
    integer, dimension(:,:), intent(in) :: flavor
    integer, dimension(:,:), intent(out), allocatable :: gpoint_flavor
    integer :: ngpt, igpt, iatm
    ngpt = size(gpt2band)
    allocate(gpoint_flavor(2,ngpt))
    do igpt=1,ngpt
      do iatm=1,2
        gpoint_flavor(iatm,igpt) = key_species_pair2flavor( &
          flavor, &
          rewrite_key_species_pair(key_species(:,iatm,gpt2band(igpt))) &
        )
      end do
    end do
  end subroutine create_gpoint_flavor

  !--------------------------------------------------------------------------------------------------------------------
  ! create gpt2band -- a map from g-points to bands 
  !   includes allocating the result 
  ! 
  subroutine create_gpt2band(band2gpt, gpt2band)
    integer, dimension(:,:), intent(in) :: band2gpt
    integer, dimension(:), allocatable, intent(out) :: gpt2band
    integer :: iband, ngpt
    ngpt = maxval(band2gpt)
    allocate(gpt2band(ngpt))
    do iband=1,size(band2gpt,dim=2)
      gpt2band(band2gpt(1,iband):band2gpt(2,iband)) = iband
    end do
  end subroutine create_gpt2band

 !--------------------------------------------------------------------------------------------------------------------
 !
 ! Utility function to combine optical properties for absorption and Rayleigh scattering 
 !   (and reorder them for convenience, while we're at it) 
 !
 subroutine combine_and_reorder(tau, tau_rayleigh, has_rayleigh, optical_props) 
    real(wp), dimension(:,:,:), intent(in) :: tau
    real(wp), dimension(:,:,:), intent(in) :: tau_rayleigh
    logical,                    intent(in) :: has_rayleigh
    class(ty_optical_props),    intent(inout) :: optical_props
    
    integer :: icol, ilay, igpt, ncol, nlay, ngpt
    
    ncol = size(tau, 3) 
    nlay = size(tau, 2) 
    ngpt = size(tau, 1)
     
    if (.not. has_rayleigh) then
      ! index reorder (ngpt, nlay, ncol) -> (ncol,nlay,gpt)
      optical_props%tau = reorder123x321(tau)
    else
      ! combine optical depth and rayleigh scattering
      select type(optical_props)
        type is (ty_optical_props)
          ! User is asking for absorption optical depth 
          optical_props%tau = reorder123x321(tau)
        type is (ty_optical_props_2str) 
          do icol = 1, ncol
            do ilay = 1, nlay
              do igpt = 1, ngpt 
                optical_props%tau(icol,ilay,igpt) = tau(igpt,ilay,icol) + tau_rayleigh(igpt,ilay,icol)
                optical_props%ssa(icol,ilay,igpt) = tau_rayleigh(igpt,ilay,icol) / optical_props%tau(icol,ilay,igpt)
              end do 
            end do 
          end do 
          optical_props%g = 0._wp 
        type is (ty_optical_props_nstr) ! We ought to be able to combine this with above
          do icol = 1, ncol
            do ilay = 1, nlay
              do igpt = 1, ngpt 
                optical_props%tau(icol,ilay,igpt) = tau(igpt,ilay,icol) + tau_rayleigh(igpt,ilay,icol)
                optical_props%ssa(icol,ilay,igpt) = tau_rayleigh(igpt,ilay,icol) / optical_props%tau(icol,ilay,igpt)
              end do 
            end do 
          end do 
          optical_props%p = 0._wp
          optical_props%p(2,:,:,:) = 0.1_wp 
      end select
      
    end if
  end subroutine combine_and_reorder

  !--------------------------------------------------------------------------------------------------------------------
  ! return the number of reference pressure layers
  pure function get_nlay_ref(this)
    class(ty_gas_optics_specification), intent(in) :: this
    integer                                        :: get_nlay_ref

    get_nlay_ref = size(this%kmajor,dim=3)
  end function get_nlay_ref

  !--------------------------------------------------------------------------------------------------------------------
  ! return eta dimension
  pure function get_neta(this)
    class(ty_gas_optics_specification), intent(in) :: this
    integer                                        :: get_neta

    get_neta = size(this%selfrefin,dim=2)
  end function

  !--------------------------------------------------------------------------------------------------------------------
  ! translate a gas string to an index
  ! slow, simple, but sufficient implementation
  pure function gas2id(this, gas)
    use mo_util_string, only : lower_case
    class(ty_gas_optics_specification), intent(in) :: this
    character(len=*),                   intent(in) :: gas
    integer                                        :: gas2id

    integer :: i

    do i = 1, size(this%gas_names)
      if (lower_case(trim(this%gas_names(i))) == lower_case(trim(gas))) then
        gas2id = i
        return
      end if
    end do
    gas2id = -1
  end function gas2id

  !--------------------------------------------------------------------------------------------------------------------
  ! fills kminor_activity; a list where minor species are active
  !   looks to see for which entires the absorption coefficients are zero at all p, eta
  !   compiles a list of extent (2, N) where each entry is (gpt, minor_gas_index) 
  ! 
  subroutine fill_kminor_activity(this)
    class(ty_gas_optics_specification) :: this

    integer :: igpt, imnr, ilist, ngpt, nmnr, nlist

    ngpt = this%get_ngpt()
    nmnr = size(this%kminor_lower, dim=1)
    nlist = 0
    ! determine list length
    do igpt = 1, ngpt
      do imnr = 1, nmnr
        if (sum(this%kminor_lower(imnr,igpt,:,:))+sum(this%kminor_upper(imnr,igpt,:,:)) .gt. 1.E-40_wp) then
          nlist = nlist + 1
        end if
      end do
    end do
    ! fill list with (g-point, minor-gas) tuples
    ilist = 0
    if (allocated(this%kminor_activity)) deallocate(this%kminor_activity)
    allocate(this%kminor_activity(2,nlist))
    do igpt = 1, ngpt
      do imnr = 1, nmnr
        if (sum(this%kminor_lower(imnr,igpt,:,:))+sum(this%kminor_upper(imnr,igpt,:,:)) .gt. 1.E-40_wp) then
          ilist = ilist + 1
          this%kminor_activity(:,ilist) = (/ igpt,imnr /)
        end if
      end do
    end do
  end subroutine fill_kminor_activity
  !--------------------------------------------------------------------------------------------------------------------
  ! Generic procedures for checking sizes, limits
  !--------------------------------------------------------------------------------------------------------------------
  !
  ! Extents 
  !
  function check_extent_1d(array, nx, label) 
    real(wp), dimension(:),     intent(in) :: array
    integer,                    intent(in) :: nx
    character(len=*),           intent(in) :: label 
    character(len=128)                     :: check_extent_1d
  
    check_extent_1d = "" 
    if(size(array,1) /= nx) & 
      check_extent_1d = trim(label) // ' has incorrect size.' 
  end function check_extent_1d
  
  function check_extent_2d(array, nx, ny, label) 
    real(wp), dimension(:,:),   intent(in) :: array
    integer,                    intent(in) :: nx, ny 
    character(len=*),           intent(in) :: label 
    character(len=128)                     :: check_extent_2d
  
    check_extent_2d = "" 
    if(size(array,1) /= nx .or. size(array,2) /= ny) & 
      check_extent_2d = trim(label) // ' has incorrect size.' 
  end function check_extent_2d
  
  function check_extent_3d(array, nx, ny, nz, label) 
    real(wp), dimension(:,:,:), intent(in) :: array
    integer,                    intent(in) :: nx, ny, nz
    character(len=*),           intent(in) :: label 
    character(len=128)                     :: check_extent_3d
  
    check_extent_3d = "" 
    if(size(array,1) /= nx .or. size(array,2) /= ny .or. size(array,3) /= nz) & 
      check_extent_3d = trim(label) // ' has incorrect size.' 
  end function check_extent_3d
  
  !
  ! Values 
  !
  function check_range_1D(val, minV, maxV, label)
    real(wp), dimension(:),     intent(in) :: val
    real(wp),                   intent(in) :: minV, maxV
    character(len=*),           intent(in) :: label 
    character(len=128)                     :: check_range_1D
    
    check_range_1D = ""
    if(any(val < minV) .or. any(val > maxV)) & 
      check_range_1D = trim(label) // ' values out of range.' 
  end function check_range_1D

  function check_range_2D(val, minV, maxV, label)
    real(wp), dimension(:,:),   intent(in) :: val
    real(wp),                   intent(in) :: minV, maxV
    character(len=*),           intent(in) :: label 
    character(len=128)                     :: check_range_2D
    
    check_range_2D = ""
    if(any(val < minV) .or. any(val > maxV)) & 
      check_range_2D = trim(label) // ' values out of range.' 
  end function check_range_2D
  
  function check_range_3D(val, minV, maxV, label)
    real(wp), dimension(:,:,:), intent(in) :: val
    real(wp),                   intent(in) :: minV, maxV
    character(len=*),           intent(in) :: label 
    character(len=128)                     :: check_range_3D
    
    check_range_3D = ""
    if(any(val < minV) .or. any(val > maxV)) & 
      check_range_3D = trim(label) // ' values out of range.' 
  end function check_range_3D
  !------------------------------------------------------------------------------------------


end module mo_gas_optics_specification
