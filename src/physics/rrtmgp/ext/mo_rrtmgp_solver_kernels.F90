! This code is part of
! RRTM for GCM Applications - Parallel (RRTMGP)
!
! Eli Mlawer and Robert Pincus
! Andre Wehe and Jennifer Delamere
! email:  rrtmgp@aer.com
!
! Copyright 2016,  Atmospheric and Environmental Research and
! Regents of the University of Colorado.  All right reserved.
!
! Use and duplication is permitted under the terms of the
!    BSD 3-clause license, see http://opensource.org/licenses/BSD-3-Clause
!
! Description:  Numeric calculations for radiative transfer solvers. 

module mo_rrtmgp_solver_kernels
  use mo_rrtmgp_kind,          only: wp
  use mo_rrtmgp_constants, &
                        only: cp_dry, grav ! Only needed for heating rate calculation 
  implicit none
  private 
  public :: lw_solver_noscat,  & 
            sw_solver_noscat, sw_solver_2stream, & 
            compute_heating_rate
            
  ! These routines don't really need to be visible but making them so is useful for testing.          
  public :: two_stream, &  
            adding_sw 
          
contains 
  ! Compute heating rate from fluxes 
  ! heating rate H [K/sec] = 1/(rho cp) d f_net/d z
  ! Here use hydrostatic equation for density and heat capacity of dry air
  function compute_heating_rate(flux_up, flux_dn, plev, heating_rate) result(error_msg) 
    real(wp), dimension(:,:), intent(in ) :: flux_up, flux_dn, & !< fluxes at interfaces [W/m2]
                                             plev                !< pressure at interfaces [Pa] 
    real(wp), dimension(:,:), intent(out) :: heating_rate        !< heating rate within layer [K/sec] 
    character(len=128)                    :: error_msg
    ! ---------
    integer :: ncol, nlay, ilay 
    ! ---------
    error_msg = "" 
    ncol = size(flux_up, 1) 
    nlay = size(flux_up, 2) - 1 
    
    if(any([size(flux_dn, 1), size(flux_dn, 2)] /= [ncol, nlay+1])) & 
      error_msg = "heating_rate: flux_dn array inconsistently sized." 
    if(any([size(plev,    1), size(plev,     2)] /= [ncol, nlay+1])) & 
      error_msg = "heating_rate: plev array inconsistently sized." 
    if(any([size(heating_rate, 1), size(heating_rate, 2)] /= [ncol, nlay])) & 
      error_msg = "heating_rate: heating_rate array inconsistently sized." 
    if(error_msg /= "") return 
    
    do ilay = 1, nlay
      heating_rate(1:ncol,ilay) = (flux_up(1:ncol,ilay+1) - flux_up(1:ncol,ilay) - &
                                   flux_dn(1:ncol,ilay+1) + flux_dn(1:ncol,ilay)) * &
                                  grav / (cp_dry * (plev(1:ncol,ilay+1) - plev(1:ncol,ilay)))
    end do

  end function compute_heating_rate 
  ! ---------------------------------------------------------------
  !
  ! Description: Direct beam calculation
  !
  subroutine lw_solver_noscat(ncol, nlay, ngpt, &
                              top_is_1, tau, mu, lay_source, lev_source_inc, lev_source_dec, sfc_emis, sfc_src, radn_up, radn_dn)
    integer,                    intent( in) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical,                    intent( in) :: top_is_1
    real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: tau          ! Absorption optical thickness []
    real(wp), dimension(ncol,       ngpt), intent( in) :: mu           ! cosine of diffusivity angle, dimension []
    real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: lay_source   ! Planck source at layer average temperature [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), target, &
                                           intent( in) :: lev_source_inc
                                        ! Planck source at layer edge for radiation in increasing ilay direction [W/m2]
                                        ! Includes spectral weighting that accounts for state-dependent frequency to g-space mapping
    real(wp), dimension(ncol,nlay+1,ngpt), target, &
                                           intent( in) :: lev_source_dec
                                               ! Planck source at layer edge for radiation in decreasing ilay direction [W/m2]
    real(wp), dimension(ncol,       ngpt), intent( in) :: sfc_emis         ! Surface emissivity      []
    real(wp), dimension(ncol,       ngpt), intent( in) :: sfc_src          ! Surface source function [W/m2]
    real(wp), dimension(ncol,nlay+1,ngpt), intent(out) :: radn_up, radn_dn ! Radiances [W/m2-str]
                                                                           ! Top level must contain incident flux boundary condition

    ! Local variables
    integer :: sfcLev
    real(wp), dimension(ncol, nlay) :: tau_loc, & ! path length (tau/mu)
                                       trans      ! transmittance, exp(-tau/mu) []

    ! Level Planck sources for upward and downward radiation
    real(wp), dimension(:,:,:), pointer :: lev_source_up, lev_source_dn

    integer :: icol, ilev, igpt
    ! ------------------------------------
    ! Which way is up?
    if(top_is_1) then
      sfcLev = nlay+1
      lev_source_up => lev_source_dec
      lev_source_dn => lev_source_inc
    else
      sfcLev = 1
      lev_source_up => lev_source_inc
      lev_source_dn => lev_source_dec
    end if

    do igpt = 1, ngpt
      ! Floor on optical depth to avoid ill-conditioning in layer emission function
      do ilev = 1, nlay
        tau_loc(:,ilev) = max(tau(:,ilev,igpt)/mu(:,igpt), 2._wp * TINY(1._wp))
      end do
      ! Transmittance: exp(-tau/mu)
      trans(:,:)  = exp(-tau_loc(:,:))

      ! Indexing into arrays for upward and downward propagation depends on the vertical
      !   orientation of the arrays (whether the domain top is at the first or last index)
      ! We write the loops out explicitly so compilers will have no trouble optimizing them.

      ! Downward propagation
      if(top_is_1) then
        ! For the flux at this level, what was the previous level, and which layer has the
        !   radiation just passed through?
        ! layer index = level index - 1
        ! previous level is up (-1)
        do ilev = 2, nlay+1
          radn_dn(:,ilev,igpt) = trans(:,ilev-1) * radn_dn(:,ilev-1,igpt) +  &
                                (1._wp - trans(:,ilev-1)) *                  &
                                lay_emission(lay_source(:,ilev-1,igpt),  &
                                             lev_source_dn(:,ilev,igpt), &
                                             tau_loc(:,ilev-1), trans(:,ilev-1))
        end do
      else
        ! layer index = level index
        ! previous level is up (+1)
        do ilev = nlay, 1, -1
          radn_dn(:,ilev,igpt) = trans(:,ilev  ) * radn_dn(:,ilev+1,igpt) +  &
                                (1._wp - trans(:,ilev  )) *                  &
                                lay_emission(lay_source(:,ilev  ,igpt),  &
                                             lev_source_dn(:,ilev,igpt), &
                                             tau_loc(:,ilev  ), trans(:,ilev  ))
        end do
      end if

      ! Surface reflection and emission
      radn_up(:,sfcLev,igpt) = radn_dn(:, sfcLev,igpt) * (1._wp - sfc_emis(:,igpt)) + &
                               sfc_src(:        ,igpt)         *  sfc_emis(:,igpt)

      ! Upward propagation
      if(top_is_1) then
        ! layer index = level index
        ! previous level is down (+1)
        do ilev = nlay, 1, -1
          radn_up(:,ilev,igpt) = trans(:,ilev  ) * radn_up(:,ilev+1,igpt) +  &
                                 (1._wp - trans(:,ilev  )) *                 &
                                 lay_emission(lay_source(:,ilev  ,igpt),  &
                                              lev_source_up(:,ilev,igpt), &
                                              tau_loc(:,ilev  ), trans(:,ilev  ))
        end do
      else
        ! layer index = level index - 1
        ! previous level is down (-1)
        do ilev = 2, nlay+1
          radn_up(:,ilev,igpt) = trans(:,ilev-1) * radn_up(:,ilev-1,igpt) +  &
                                 (1._wp - trans(:,ilev-1)) *                 &
                                 lay_emission(lay_source(:,ilev-1,igpt),  &
                                              lev_source_up(:,ilev,igpt), &
                                              tau_loc(:,ilev-1), trans(:,ilev-1))
        end do
      end if
    end do  ! g point loop

  end subroutine lw_solver_noscat
! ---------------------------------------------------------------
  elemental function lay_emission(lay_src, lev_src, tau, trans)
    real(wp), intent(in) :: lay_src, lev_src, tau, trans
    real(wp)             :: lay_emission

    ! Effective Planck function
    ! Provide the amount of emission from a layer given layer average source,
    !   layer edge source, layer transmittance and optical depth
    ! Here we're using "linear in tau" approximation: source function is assumed
    !   to vary linearly with optical depth
    ! See Clough et al., 1992, doi: 10.1029/92JD01419, Eq 13

    lay_emission = merge(lev_src + 2._wp * (lay_src - lev_src) * &
                         (1._wp/tau - trans/(1._wp - trans)),    &
                          0._wp, trans < 1._wp - spacing(1._wp))

  end function lay_emission
! ---------------------------------------------------------------
!   Shortwave kernels 
! ---------------------------------------------------------------
  pure subroutine sw_solver_noscat(ncol, nlay, ngpt, &
                              top_is_1, tau, mu0, flux_dir)  
    integer,                    intent( in) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical,                    intent( in) :: top_is_1
    real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: tau          ! Absorption optical thickness []
    real(wp), dimension(ncol            ), intent( in) :: mu0          ! cosine of solar zenith angle 
    real(wp), dimension(ncol,nlay+1,ngpt), intent(out) :: flux_dir     ! Direct-beam flux, spectral [W/m2]
                                                                       ! Top level must contain incident flux boundary condition
    integer :: icol, ilev, igpt
    ! ------------------------------------
    ! Indexing into arrays for upward and downward propagation depends on the vertical
    !   orientation of the arrays (whether the domain top is at the first or last index)
    ! We write the loops out explicitly so compilers will have no trouble optimizing them.

    ! Downward propagation
    if(top_is_1) then
      ! For the flux at this level, what was the previous level, and which layer has the
      !   radiation just passed through?
      ! layer index = level index - 1
      ! previous level is up (-1)
      do igpt = 1, ngpt
        do ilev = 2, nlay+1
          flux_dir(:,ilev,igpt) = flux_dir(:,ilev-1,igpt) * exp(-tau(:,ilev,igpt)/mu0(:)) 
        end do
      end do 
    else
      ! layer index = level index
      ! previous level is up (+1)
      do igpt = 1, ngpt
        do ilev = nlay, 1, -1
          flux_dir(:,ilev,igpt) = flux_dir(:,ilev+1,igpt) * exp(-tau(:,ilev,igpt)/mu0(:)) 
        end do
      end do 
    end if
  end subroutine sw_solver_noscat
  ! ---------------------------------------------------------------
   subroutine sw_solver_2stream (ncol, nlay, ngpt, top_is_1, &
                                 tau, ssa, g, mu0,           &
                                 sfc_alb_dir, sfc_alb_dif,   & 
!++dbg
!                                 flux_up, flux_dn, flux_dir)   
                                 flux_up, flux_dn, flux_dir, ierr)
!--dbg
    integer,                    intent( in) :: ncol, nlay, ngpt ! Number of columns, layers, g-points
    logical,                    intent( in) :: top_is_1
    real(wp), dimension(ncol,nlay,  ngpt), intent( in) :: tau, &  ! Optical thickness, 
                                                          ssa, &  ! single-scattering albedo, 
                                                          g       ! asymmetry parameter []
    real(wp), dimension(ncol            ), intent( in) :: mu0     ! cosine of solar zenith angle 
    real(wp), dimension(ncol,       ngpt), intent( in) :: sfc_alb_dir, sfc_alb_dif 
                                                                  ! Spectral albedo of surface to direct and diffuse radiation 
    real(wp), dimension(ncol,nlay+1,ngpt), & 
                                           intent(out) :: flux_up, flux_dn, &  ! Fluxes [W/m2]
                                                          flux_dir             ! Downward direct 
                                                                               ! Top level (= merge(1, nlay+1, top_is_1)
                                                                               ! must contain incident flux boundary condition
!++dbg
    logical, intent(out) :: ierr
!--dbg
    integer :: igpt
    real(wp), dimension(ncol,nlay) :: Rdif, Tdif, Rdir, Tdir, Tnoscat
    ! ------------------------------------
    !
    ! It's not clear how to be most efficient here. The two-stream calculation is atomic 
    !   but the transport (adding) is not. Combining the loops reduces the memory footprint 
    !   for the transmission and reflection arrays by a factor of ngpt, normally a lot, 
    !   but then the size for the two-stream problem is greatly reduced. 
    !
    do igpt = 1, ngpt
      !
      ! Compute cell properties 
      ! 
      call two_stream(ncol, nlay, mu0,                                &
                      tau (:,:,igpt), ssa (:,:,igpt), g   (:,:,igpt), & 
!++dbg
!                      Rdif, Tdif, Rdir, Tdir, Tnoscat) 
                      Rdif, Tdif, Rdir, Tdir, Tnoscat, ierr) 
      if (ierr) then
         print*, 'ERROR return from two_stream at igpt=', igpt
         return
      end if
!--dbg
      call adding_sw(ncol, nlay, top_is_1,                     & 
                     Rdif, Tdif, Rdir, Tdir, Tnoscat,          & 
                     sfc_alb_dif(:,igpt), sfc_alb_dir(:,igpt), & 
                     flux_up(:,:,igpt), flux_dn(:,:,igpt), flux_dir(:,:,igpt))
    end do 
    
  end subroutine sw_solver_2stream   

! ---------------------------------------------------------------
! Two-stream solutions to direct and diffuse reflectance and transmittance for a layer
!    with optical depth tau, single scattering albedo w0, and asymmetery parameter g.  
!
! Equations are developed in Meador and Weaver, 1980, 
!    doi:10.1175/1520-0469(1980)037<0630:TSATRT>2.0.CO;2
!
!++dbg
!  pure subroutine two_stream(ncol, nlay, mu0, tau, w0, g, & 
        subroutine two_stream(ncol, nlay, mu0, tau, w0, g, & 
!                             Rdif, Tdif, Rdir, Tdir, Tnoscat) 
                             Rdif, Tdif, Rdir, Tdir, Tnoscat, ierr) 
!--dbg
    integer,                        intent(in)  :: ncol, nlay
    real(wp), dimension(ncol),      intent(in)  :: mu0
    real(wp), dimension(ncol,nlay), intent(in)  :: tau, w0, g 
    real(wp), dimension(ncol,nlay), intent(out) :: Rdif, Tdif, Rdir, Tdir, Tnoscat
!++dbg
    logical, intent(out) :: ierr
!--dbg    
    
    ! -----------------------
    integer  :: i, j 
    
    ! Variables used in Meador and Weaver
    real(wp) :: gamma1(ncol), gamma2(ncol), gamma3(ncol), gamma4(ncol)
    real(wp) :: alpha1(ncol), alpha2(ncol), k(ncol)
    
    ! Ancillary variables
    real(wp) :: RT_term(ncol)
    real(wp) :: exp_minusktau(ncol), exp_minus2ktau(ncol)
    real(WP) :: k_mu, k_gamma3, k_gamma4
    
    ! ---------------------------------
!++dbg
    ierr = .false.
!--dbg    
    do j = 1, nlay
      do i = 1, ncol  
        ! Zdunkowski Practical Improved Flux Method "PIFM" 
        !  (Zdunkowski et al., 1980;  Contributions to Atmospheric Physics 53, 147-66)
        !
        gamma1(i)= (8._wp - w0(i,j) * (5._wp + 3._wp * g(i,j))) * .25_wp
        gamma2(i)=  3._wp *(w0(i,j) * (1._wp -         g(i,j))) * .25_wp
        gamma3(i)= (2._wp - 3._wp * mu0(i) *           g(i,j) ) * .25_wp
        gamma4(i)=  1._wp - gamma3(i)

        alpha1(i) = gamma1(i) * gamma4(i) + gamma2(i) * gamma3(i)           ! Eq. 16
        alpha2(i) = gamma1(i) * gamma3(i) + gamma2(i) * gamma4(i)           ! Eq. 17
!++dbg
!        if ( gamma2(i) > gamma1(i)) then
!           print*, 'ERROR in calc of k at i,j=',i,j
!           print*, 'gamma1, gamma2=', gamma1(i), gamma2(i)
!           print*, 'tau, w0, g=', tau(i,j), w0(i,j), g(i,j)
!           ierr = .true.
!           return
!        end if
!--dbg
      end do
    
      ! Written to encourage vectorization of exponential, square root 
      ! Eq 18;  k = SQRT(gamma1**2 - gamma2**2), limited below to avoid div by 0. 
      !   k = 0 for isotropic, conservative scattering; this lower limit on k 
      !   gives relative error with respect to conservative solution  
      !   of < 0.1% in Rdif down to tau = 10^-9
!++dbg
!      k(1:ncol) = max(sqrt((gamma1(1:ncol) - gamma2(1:ncol)) * & 
!                           (gamma1(1:ncol) + gamma2(1:ncol))), & 
!                      1.e-6_wp)
      k(1:ncol) = sqrt(max((gamma1(1:ncol) - gamma2(1:ncol))* & 
                           (gamma1(1:ncol) + gamma2(1:ncol)), & 
                           1.e-12_wp))
!--dbg
      exp_minusktau(1:ncol) = exp(-tau(1:ncol,j)*k(1:ncol)) 
        
      !
      ! Diffuse reflection and transmission
      ! 
      do i = 1, ncol
        exp_minus2ktau(i) = exp_minusktau(i) * exp_minusktau(i) 

        ! Refactored to avoid rounding errors when k, gamma1 are of very different magnitudes 
        RT_term(i) = 1._wp / (k     (i) * (1._wp + exp_minus2ktau(i))  + & 
                              gamma1(i) * (1._wp - exp_minus2ktau(i)) ) 

        ! Equation 25
        Rdif(i,j) = RT_term(i) * gamma2(i) * (1._wp - exp_minus2ktau(i)) 

        ! Equation 26
        Tdif(i,j) = RT_term(i) * 2._wp * k(i) * exp_minusktau(i)
      end do 

      !
      ! Transmittance of direct, unscattered beam. Also used below 
      !
      Tnoscat(1:ncol,j) = exp(-tau(1:ncol,j)/mu0(1:ncol))

      !
      ! Direct reflect and transmission
      !
      do i = 1, ncol 
        k_mu     = k(i) * mu0(i)
        k_gamma3 = k(i) * gamma3(i)
        k_gamma4 = k(i) * gamma4(i)
      
        !
        ! Equation 14, multiplying top and bottom by exp(-k*tau) 
        !   and rearranging to avoid div by 0. 
        !
        RT_term(i) =  w0(i,j) * RT_term(i)/merge(1._wp - k_mu*k_mu, & 
                                                 epsilon(1._wp),    &
                                                 abs(1._wp - k_mu*k_mu) >= epsilon(1._wp))

        Rdir(i,j) = RT_term(i)  *                                        &
            ((1._wp - k_mu) * (alpha2(i) + k_gamma3)                     - &
             (1._wp + k_mu) * (alpha2(i) - k_gamma3) * exp_minus2ktau(i) - &
             2.0_wp * (k_gamma3 - alpha2(i) * k_mu)  * exp_minusktau (i) * Tnoscat(i,j)) 

        !
        ! Equation 15, multiplying top and bottom by exp(-k*tau), 
        !   multiplying through by exp(-tau/mu0) to 
        !   prefer underflow to overflow
        !
        Tdir(i,j) = Tnoscat(i,j) - & 
                  RT_term(i) * ((1._wp + k_mu) * (alpha1(i) + k_gamma4)                     * Tnoscat(i,j) - &
                                (1._wp - k_mu) * (alpha1(i) - k_gamma4) * exp_minus2ktau(i) * Tnoscat(i,j) - &
                                2.0_wp * (k_gamma4 + alpha1(i) * k_mu)  * exp_minusktau (i))

        Tdir(i,j) = Tdir(i,j) - Tnoscat(i,j)
      end do 
    end do 
  end subroutine two_stream
! ---------------------------------------------------------------
! Transport of solar radiation through a vertically layered atmosphere.
!   Equations are after Shonk and Hogan 2008, doi:10.1175/2007JCLI1940.1 (SH08)
!
 ! Indexing into arrays for upward and downward propagation depends on the vertical
 !   orientation of the arrays (whether the domain top is at the first or last index)
 ! We write the loops out explicitly so compilers will have no trouble optimizing them.

  pure subroutine adding_sw(ncol, nlay, top_is_1,            & 
                            rdif, tdif, Rdir, Tdir, Tnoscat, & 
                            sfc_alb_dif, sfc_alb_dir,        & 
                            flux_up, flux_dn_dif, flux_dn_dir) 
    integer,                          intent(in   ) :: ncol, nlay
    logical,                          intent(in   ) :: top_is_1
    real(wp), dimension(ncol,nlay  ), intent(in   ) :: rdif, tdif, Rdir, Tdir, Tnoscat
    real(wp), dimension(ncol       ), intent(in   ) :: sfc_alb_dif, sfc_alb_dir
    real(wp), dimension(ncol,nlay+1), intent(  out) :: flux_up
    ! intent(inout) because top layer includes incident flux
    real(wp), dimension(ncol,nlay+1), intent(inout) :: flux_dn_dif, flux_dn_dir
    ! ------------------
    integer :: ilev 
    real(wp), dimension(ncol,nlay+1) :: albedo, &  ! reflectivity to diffuse radiation below this level
                                                   ! alpha in SH08 
                                        dir_src    ! source of diffuse radiation from direct radiation 
                                                   ! G in SH08
    real(wp), dimension(ncol,nlay  ) :: denom      ! beta in SH08
    ! ------------------
    !
    ! Direct beam flux
    !
    if(top_is_1) then
      ! For the flux at this level, what was the previous level, and which layer has the
      !   radiation just passed through?
      ! layer index = level index - 1
      ! previous level is up (-1)
      do ilev = 2, nlay+1
        flux_dn_dir(:,ilev) = flux_dn_dir(:,ilev-1) * Tnoscat(:,ilev)
      end do
    else
      ! layer index = level index
      ! previous level is up (+1)
      do ilev = nlay, 1, -1
        flux_dn_dir(:,ilev) = flux_dn_dir(:,ilev+1) * Tnoscat(:,ilev)
      end do
    end if
    
    if(top_is_1) then
      ! layer index = level index - 1
      ! previous level is up (-1)
      ilev = nlay + 1
      ! Albedo of lowest level is the surface albedo... 
      albedo(:,ilev)  = sfc_alb_dif(:)
      ! ... and source of diffuse radiation is direct beam 
      dir_src(:,ilev) = sfc_alb_dir(:) * flux_dn_dir(:,ilev) ! Should this be weighted by cos(mu)? 
      
      !
      ! From bottom to top of atmosphere -- 
      !   compute albedo and contribution of direct beam to diffuse radiation
      !
      do ilev = nlay, 1, -1
        denom(:, ilev) = 1._wp/(1._wp - rdif(:,ilev)*albedo(:,ilev+1))                 ! Eq 10
        albedo(:,ilev) = rdif(:,ilev) + & 
                         tdif(:,ilev)*tdif(:,ilev) * albedo(:,ilev+1) * denom(:,ilev) ! Equation 9
        !
        ! Equation 11 -- source is "scattering of the direct solar beam into the diffuse components"
        !   i.e. reflected flux from this layer (Rdir*flux_dn_dir) and 
        !   transmitted through the layer (Tdir*flux_dn_dir) and reflected from layers below (albedo)
        !
        dir_src(:,ilev) = Rdir(:,ilev) * flux_dn_dir(:,ilev) + &       
                          tdif(:,ilev) * denom(:,ilev) *       & 
                            (dir_src(:,ilev+1) + albedo(:,ilev+1)*Tdir(:,ilev)*flux_dn_dir(:,ilev))
      end do
      
      ! Eq 12, at the top of the domain upwelling diffuse is due to ...
      ilev = 1 
      flux_up(:,ilev) = flux_dn_dif(:,ilev) * albedo(:,ilev) + & ! ... reflection of incident diffuse and
                        dir_src(:,ilev)                          ! scattering by the direct beam below 
                        
      !
      ! From the top of the atmosphere downward -- compute fluxes
      !
      do ilev = 2, nlay+1
        flux_up(:,ilev) = flux_dn_dif(:,ilev) * albedo(:,ilev) + & ! Equation 12
                          dir_src(:,ilev) 
        ! ... The equation doesn't have the denominator in it but it seems like it should 
        flux_dn_dif(:,ilev) = (tdif(:,ilev-1)*flux_dn_dif(:,ilev-1) + &  ! Equation 13
                               rdif(:,ilev-1)*dir_src(:,ilev) +       & 
                               Tdir(:,ilev-1)*flux_dn_dir(:,ilev-1)) * denom(:,ilev-1)
      end do 
    else
      ! layer index = level index
      ! previous level is up (+1)
      ilev = 1
      ! Albedo of lowest level is the surface albedo... 
      albedo(:,ilev)  = sfc_alb_dif(:)
      ! ... and source of diffuse radiation is direct beam 
      dir_src(:,ilev) = sfc_alb_dir(:) * flux_dn_dir(:,ilev) ! Should this be weighted by cos(mu)? 
      
      !
      ! From bottom to top of atmosphere -- 
      !   compute albedo and contribution of direct beam to diffuse radiation
      !
      do ilev = 1, nlay
        denom(:, ilev  ) = 1._wp/(1._wp - rdif(:,ilev)*albedo(:,ilev))                ! Eq 10
        albedo(:,ilev+1) = rdif(:,ilev) + & 
                           tdif(:,ilev)*tdif(:,ilev) * albedo(:,ilev) * denom(:,ilev) ! Equation 9
        !
        ! Equation 11 -- source is "scattering of the direct solar beam into the diffuse components"
        !   i.e. reflected flux from this layer (Rdir*flux_dn_dir) and 
        !   transmitted through the layer (Tdir*flux_dn_dir) and reflected from layers below (albedo)
        !
        dir_src(:,ilev+1) = Rdir(:,ilev) * flux_dn_dir(:,ilev+1) + &       
                            tdif(:,ilev) * denom(:,ilev) *       & 
                              (dir_src(:,ilev) + albedo(:,ilev)*Tdir(:,ilev)*flux_dn_dir(:,ilev+1))
      end do
      
      ! Eq 12, at the top of the domain upwelling diffuse is due to ...
      ilev = nlay+1 
      flux_up(:,ilev) = flux_dn_dif(:,ilev) * albedo(:,ilev) + & ! ... reflection of incident diffuse and
                        dir_src(:,ilev)                          ! scattering by the direct beam below 
                        
      !
      ! From the top of the atmosphere downward -- compute fluxes
      !
      do ilev = nlay, 1, -1
        ! ... The equation doesn't have the denominator in it but it seems like it should 
        flux_dn_dif(:,ilev) = (tdif(:,ilev)*flux_dn_dif(:,ilev+1) + &  ! Equation 13
                               rdif(:,ilev)*dir_src(:,ilev) +       & 
                               Tdir(:,ilev)*flux_dn_dir(:,ilev+1)) * denom(:,ilev)
        flux_up(:,ilev) = flux_dn_dif(:,ilev) * albedo(:,ilev) + & ! Equation 12
                          dir_src(:,ilev) 
      end do 
    end if 
   
  end subroutine adding_sw
end module mo_rrtmgp_solver_kernels
