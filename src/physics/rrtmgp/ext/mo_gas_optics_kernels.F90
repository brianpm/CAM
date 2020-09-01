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
! Description: Numeric calculations for gas optics. Absorption and Rayleigh optical depths, 
!   source functions. 

module mo_gas_optics_kernels
  use mo_rrtmgp_kind,      only: wp
  implicit none 
  
  interface interpolate2D
    module procedure interpolate2D_1, interpolate2D_all
  end interface interpolate2D
contains
  ! --------------------------------------------------------------------------------------
  subroutine gas_optical_depths_major(ncol,nlay,ngpt,nflav, &
    gpoint_flavor,kmajor,col_mix,fmajor,&
    jeta,tropo,jtemp,jpress, &
    tau)
    ! input dimensions
    integer, intent(in) :: ncol, nlay, ngpt, nflav ! dimensions

    ! inputs from object
    integer,  dimension(2,ngpt),  intent(in) :: gpoint_flavor
    real(wp), dimension(:,:,:,:), intent(in) :: kmajor
    
    ! inputs from profile or parent function
    real(wp), dimension(2,    nflav,ncol,nlay), intent(in) :: col_mix
    real(wp), dimension(2,2,2,nflav,ncol,nlay), intent(in) :: fmajor
    integer,  dimension(2,    nflav,ncol,nlay), intent(in) :: jeta
    logical,  dimension(ncol,nlay), intent(in) :: tropo
    integer,  dimension(ncol,nlay), intent(in) :: jtemp, jpress

    ! outputs
    real(wp), dimension(ngpt,nlay,ncol), intent(inout) :: tau
    ! -----------------
    ! local variables
    real(wp) :: tau_major ! major species optical depth
    ! local index
    integer :: icol, ilay, iflav, igpt, itropo
    ! -----------------

    do ilay = 1, nlay
      do icol = 1, ncol 

        ! itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
        itropo = merge(1,2,tropo(icol,ilay))

        ! optical depth calculation for major species
        do igpt = 1, ngpt
          iflav = gpoint_flavor(itropo, igpt)
          tau_major = &
            ! interpolation in temperature, pressure, and eta
            interpolate3D(col_mix(:,iflav,icol,ilay), &
                          fmajor(:,:,:,iflav,icol,ilay), kmajor(igpt,:,:,:), &
                          jeta(:,iflav,icol,ilay), jtemp(icol,ilay),jpress(icol,ilay)+itropo)
          tau(igpt,ilay,icol) = tau(igpt,ilay,icol) + tau_major
        end do ! igpt
      end do 
    end do ! ilay
  end subroutine gas_optical_depths_major

  ! ----------------------------------------------------------
  ! compute water vapor continuum optical depths
  subroutine gas_optical_depths_continuum( &
    ncol,nlay,ngpt,ngas,nflav, &
    flavor,gpoint_flavor,selfrefin,forrefin,stpfac, &
    idx_h2o, play,tlay,vmr,col_gas,fminor,jeta,tropo,jtemp, &
    tau)
    ! input dimensions
    integer, intent(in) :: ncol,nlay,ngpt,ngas,nflav

    ! inputs from object
    integer,  dimension(:,:),   intent(in) :: flavor
    integer,  dimension(:,:),   intent(in) :: gpoint_flavor
    real(wp), dimension(:,:,:), intent(in) :: selfrefin, forrefin
    real(wp),                   intent(in) :: stpfac

    ! inputs from profile or parent function
    integer,                             intent(in) :: idx_h2o
    real(wp), dimension(ncol,nlay),      intent(in) :: play, tlay
    real(wp), dimension(ncol,nlay,ngas), intent(in) :: vmr
    real(wp), dimension(ncol,nlay,ngas), intent(in) :: col_gas
    real(wp), dimension(2,2,nflav,ncol,nlay), intent(in) :: fminor
    integer,  dimension(2,  nflav,ncol,nlay), intent(in) :: jeta
    logical,  dimension(ncol,nlay),      intent(in) :: tropo
    integer,  dimension(ncol,nlay),      intent(in) :: jtemp

    ! outputs
    real(wp), dimension(ngpt,nlay,ncol), intent(inout) :: tau
    ! -----------------
    ! local variables
    ! factor needed for foreign continuum optical depth calculation
    real(wp) :: forfac(ncol),  selffac(ncol) ! scaling variables for foreign and self-continuum 
    real(wp) :: tau_for(ngpt), tau_self(ngpt)

    ! local index
    integer :: icol, ilay, iflav, igpt, itropo
    logical :: water_is_key(ngpt, 2) ! Continuum only contributes if water is a key species
                                     ! Second dimension is lower/upper  
    ! -----------------
    do igpt = 1, ngpt
      water_is_key(igpt, 1) = any(flavor(:, gpoint_flavor(1, igpt)) == idx_h2o)
      water_is_key(igpt, 2) = any(flavor(:, gpoint_flavor(2, igpt)) == idx_h2o) 
    end do     
    
    do ilay = 1, nlay
      do icol = 1, ncol
        ! foreign continuum variables
        forfac(icol) = stpfac * play(icol,ilay)/tlay(icol,ilay) *  & 
                       col_gas(icol,ilay,idx_h2o) / (1._wp + vmr(icol,ilay,idx_h2o))
        ! factor that multiplies self continuum absorption coefficient
        selffac(icol)  = vmr(icol,ilay,idx_h2o) * forfac(icol)
      end do
      do icol = 1, ncol
        itropo = merge(1,2,tropo(icol,ilay)) ! itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
        do igpt = 1, ngpt
          iflav = gpoint_flavor(itropo, igpt)
          ! interpolation in temperature and eta
          tau_for(igpt)  = forfac (icol) * & 
                           interpolate2D(fminor(:,:,iflav,icol,ilay), forrefin(igpt,:,:), & 
                                         jeta(:,iflav,icol,ilay), jtemp(icol,ilay))
          tau_self(igpt) = selffac(icol) * & 
                           interpolate2D(fminor(:,:,iflav,icol,ilay), selfrefin(igpt,:,:), & 
                                         jeta(:,iflav,icol,ilay), jtemp(icol,ilay))
        end do 
        tau(1:ngpt,ilay,icol) = tau(1:ngpt,ilay,icol) + merge(tau_for(1:ngpt) + tau_self(1:ngpt), & 
                                                              0._wp,                              & 
                                                              water_is_key(1:ngpt, itropo))
      end do 
    end do ! ilay
  end subroutine gas_optical_depths_continuum

  ! ----------------------------------------------------------
  ! compute minor species optical depths
  subroutine gas_optical_depths_minor( &
    ncol,nlay,ngpt,ngas,nflav, & ! input dimensions
    gpoint_flavor,band2gpt,kminor_lower,kminor_upper,kminor_activity, &
    idx_h2o,idx_o2,idx_n2, play,tlay,col_dry,col_gas,fminor,jeta,tropo,jtemp, &
    tau )
    ! input dimensions
    integer, intent(in) :: ncol,nlay,ngpt,ngas,nflav

    ! inputs from object
    integer, dimension(:,:), intent(in) :: gpoint_flavor
    integer, dimension(:,:), intent(in) :: band2gpt
    real(wp), dimension(:,:,:,:), intent(in) :: kminor_lower, kminor_upper
    integer,  dimension(:,:),     intent(in) :: kminor_activity

    ! inputs from profile or parent function
    integer, intent(in) :: idx_h2o,idx_o2,idx_n2
    real(wp), dimension(ncol,nlay),      intent(in) :: play, tlay
    real(wp), dimension(ncol,nlay),      intent(in) :: col_dry
    real(wp), dimension(ncol,nlay,ngas), intent(in) :: col_gas
    real(wp), dimension(2,2,nflav,ncol,nlay), intent(in) :: fminor
    integer,  dimension(2,  nflav,ncol,nlay), intent(in) :: jeta
    logical,  dimension(ncol,nlay),      intent(in) :: tropo
    integer,  dimension(ncol,nlay),      intent(in) :: jtemp

    ! outputs
    real(wp), dimension(ngpt,nlay,ncol), intent(inout) :: tau
    ! -----------------
    ! local variables
    real(wp) :: kminor, tau_minor ! minor species absorption coefficient, optical depth 

    ! local index
    integer ::  icol, ilay, iflav, igpt, itropo, imnr, ilist, nlist
    real(wp) :: scaling 
    ! -----------------
    do ilay = 1, nlay
      do icol = 1, ncol 
        itropo = merge(1,2,tropo(icol,ilay)) ! itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
        ! loop through the list of species/band combinations for which contributions of minor species are included
        nlist = size(kminor_activity, dim=2)
        do ilist = 1, nlist
          igpt = kminor_activity(1,ilist)
          imnr = kminor_activity(2,ilist)
          iflav = gpoint_flavor(itropo, igpt)
          
          ! interpolation in temperature and eta
          if (tropo(icol,ilay)) then ! lower atmosphere
            kminor = &
              interpolate2D(fminor(:,:,iflav,icol,ilay), kminor_lower(imnr,igpt,:,:), jeta(:,iflav,icol,ilay), jtemp(icol,ilay))
          else ! upper atmosphere
            kminor = &
              interpolate2D(fminor(:,:,iflav,icol,ilay), kminor_upper(imnr,igpt,:,:), jeta(:,iflav,icol,ilay), jtemp(icol,ilay))
          end if

          scaling = col_gas(icol,ilay,imnr) ! standard treatment of minor gases
          ! different treatment for collision-induced absorption
          ! oxygen
          if (imnr == idx_o2) & 
            scaling = scaling * play(icol,ilay) / tlay(icol,ilay)
          ! nitrogen
          if (imnr == idx_n2) then
            scaling = scaling * play(icol,ilay) / tlay(icol,ilay)
            if (igpt < band2gpt(2,1)) &  ! if band == 1 then
             ! absorption only due to N2-N2
              scaling = scaling * (col_gas(icol,ilay,imnr) / (col_dry(icol,ilay) + col_gas(icol,ilay,idx_h2o))) 
          end if 
          
          tau_minor = kminor * scaling
          tau(igpt,ilay,icol) = tau(igpt,ilay,icol) + tau_minor
        end do ! ilist
      end do 
    end do ! ilay
  end subroutine gas_optical_depths_minor

  ! ----------------------------------------------------------
  ! compute minor species optical depths
  subroutine gas_optical_depths_rayleigh( &
    ncol,nlay,ngpt,ngas,nflav, & ! input dimensions
    gpoint_flavor,krayl, &
    idx_h2o,idx_o2,idx_n2, play,tlay,col_dry,col_gas,fminor,jeta,tropo,jtemp, &
    tau_rayleigh)
    
    ! input dimensions
    integer, intent(in) :: ncol,nlay,ngpt,ngas,nflav

    ! inputs from object
    integer, dimension(:,:), intent(in) :: gpoint_flavor
    real(wp), dimension(:,:,:,:), intent(in) :: krayl

    ! inputs from profile or parent function
    integer, intent(in) :: idx_h2o,idx_o2,idx_n2
    real(wp), dimension(ncol,nlay),      intent(in) :: play, tlay
    real(wp), dimension(ncol,nlay),      intent(in) :: col_dry
    real(wp), dimension(ncol,nlay,ngas), intent(in) :: col_gas
    real(wp), dimension(2,2,nflav,ncol,nlay), intent(in) :: fminor
    integer,  dimension(2,  nflav,ncol,nlay), intent(in) :: jeta
    logical,  dimension(ncol,nlay),      intent(in) :: tropo
    integer,  dimension(ncol,nlay),      intent(in) :: jtemp

    ! outputs
    real(wp), dimension(ngpt,nlay,ncol), intent(inout) :: tau_rayleigh
    ! -----------------
    ! local variables
    real(wp) :: k, tau_rayl ! rayleigh scattering coefficient, optical depth 

    ! local index
    integer :: icol, ilay, iflav, igpt, itropo, imnr, ilist, nlist
    ! -----------------
    do ilay = 1, nlay
      do icol = 1, ncol
        itropo = merge(1,2,tropo(icol,ilay)) ! itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
        do igpt = 1, ngpt
          iflav = gpoint_flavor(itropo, igpt)
          k = interpolate2D(fminor(:,:,iflav,icol,ilay), krayl(igpt,:,:,itropo), jeta(:,iflav,icol,ilay), jtemp(icol,ilay))
          tau_rayleigh(igpt,ilay,icol) =  k * (col_gas(icol,ilay,idx_h2o)+col_dry(icol,ilay))
        end do ! igpt
      end do 
    end do ! ilay
  end subroutine gas_optical_depths_rayleigh
  
  ! ----------------------------------------------------------
  subroutine source(ncol, nlay, ngpt, nbnd, ngas, nflav,  & 
                    tlay, tlev, tsfc, sfc_lay,            & 
                    fmajor, jeta, tropo, jtemp, jpress,   & 
                    band2gpt, pfracin, temp_ref_min, totplnk_delta, totplnk, gpoint_flavor, & 
                    sfc_src, lay_src, lev_src_inc, lev_src_dec)
    integer,            intent(in) :: ncol, nlay, ngpt, nbnd, ngas, nflav, sfc_lay   
    real(wp), dimension(ncol,nlay  ) :: tlay 
    real(wp), dimension(ncol,nlay+1) :: tlev 
    real(wp), dimension(ncol       ) :: tsfc 
    ! Interpolation variables 
    real(wp), dimension(2,2,2,nflav,ncol,nlay), intent(in) :: fmajor
    integer,  dimension(2,    nflav,ncol,nlay), intent(in) :: jeta
    logical,  dimension(            ncol,nlay), intent(in) :: tropo
    integer,  dimension(            ncol,nlay), intent(in) :: jtemp, jpress
    ! Table-specific
    integer, dimension(2, nbnd),  intent(in) :: band2gpt! Starting and ending g-points of bands
    real(wp),                     intent(in) :: temp_ref_min, totplnk_delta
    real(wp), dimension(:,:,:,:), intent(in) :: pfracin       ! change to provide size 
    real(wp), dimension(:,:),     intent(in) :: totplnk       ! change to provide size 
    integer,  dimension(:,:),     intent(in) :: gpoint_flavor ! change to provide size 

    real(wp), dimension(ncol,       ngpt), intent(out) :: sfc_src 
    real(wp), dimension(ncol,nlay,  ngpt), intent(out) :: lay_src 
    real(wp), dimension(ncol,nlay+1,ngpt), intent(out) :: lev_src_inc, lev_src_dec
    ! -----------------
    ! local
    integer  :: ilay, icol, igpt, ibnd, itropo, iflav
    real(wp) :: pfrac(ngpt, ncol, nlay)
    real(wp) :: planck_function(nbnd)
    ! -----------------

    do ilay = 1, nlay
      do icol = 1, ncol
        ! itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
        itropo = merge(1,2,tropo(icol,ilay))
        ! Planck fraction calculation
        do igpt = 1, ngpt
          iflav = gpoint_flavor(itropo, igpt)
          pfrac(igpt,icol,ilay) = &
            ! interpolation in temperature, pressure, and eta
            interpolate3D((/1._wp,1._wp/), fmajor(:,:,:,iflav,icol,ilay), pfracin(igpt,:,:,:), &
                          jeta(:,iflav,icol,ilay), jtemp(icol,ilay),jpress(icol,ilay)+itropo)
        end do ! igpt
      end do ! col
    end do 

    ! compute surface source radiances for each g-point
    do icol = 1, ncol
      planck_function(:) = interpolate1D(tsfc(icol), temp_ref_min, totplnk_delta, totplnk)
      do ibnd = 1, nbnd
        sfc_src(icol,band2gpt(1,ibnd):band2gpt(2,ibnd)) = &  
               pfrac(band2gpt(1,ibnd):band2gpt(2,ibnd),icol,sfc_lay) * planck_function(ibnd) 
      end do 
    end do ! icol

    ! compute layer source radiances for each g-point
    do icol = 1, ncol
      do ilay = 1, nlay
        planck_function(:) = interpolate1D(tlay(icol,ilay), temp_ref_min, totplnk_delta, totplnk)
        do ibnd = 1, nbnd
          lay_src(icol,ilay,band2gpt(1,ibnd):band2gpt(2,ibnd)) = &
                      pfrac(band2gpt(1,ibnd):band2gpt(2,ibnd),icol,ilay) * planck_function(ibnd) 
        end do 
      end do ! ilay
    end do ! icol

    ! compute level Planck source function for each g-point in increasing ilay direction
    do icol = 1, ncol
      do ilay = 2, nlay
        planck_function(:) = interpolate1D(tlev(icol,ilay), temp_ref_min, totplnk_delta, totplnk)
        do ibnd = 1, nbnd
          lev_src_inc(icol,ilay,band2gpt(1,ibnd):band2gpt(2,ibnd)) = &
                          pfrac(band2gpt(1,ibnd):band2gpt(2,ibnd),icol,ilay-1) * planck_function(ibnd) 
          lev_src_dec(icol,ilay,band2gpt(1,ibnd):band2gpt(2,ibnd)) = &
                          pfrac(band2gpt(1,ibnd):band2gpt(2,ibnd),icol,ilay  ) * planck_function(ibnd) 
        end do 
      end do ! ilay
    end do ! icol

    ! Edge cases 
    lev_src_inc(:,1,     :) = 0 ! this value is padding
    lev_src_dec(:,nlay+1,:) = 0 ! this value is padding
    do icol = 1, ncol
        ilay = 1 
        planck_function(:) = interpolate1D(tlev(icol,ilay), temp_ref_min, totplnk_delta, totplnk)
        do ibnd = 1, nbnd
          lev_src_dec(icol,ilay,band2gpt(1,ibnd):band2gpt(2,ibnd)) = &
                          pfrac(band2gpt(1,ibnd):band2gpt(2,ibnd),icol,ilay) * planck_function(ibnd) 
        end do 
        ilay = nlay+1
        planck_function(:) = interpolate1D(tlev(icol,ilay), temp_ref_min, totplnk_delta, totplnk)
        do ibnd = 1, nbnd
          lev_src_inc(icol,ilay,band2gpt(1,ibnd):band2gpt(2,ibnd)) = &
                          pfrac(band2gpt(1,ibnd):band2gpt(2,ibnd),icol,ilay-1) * planck_function(ibnd) 
        end do         
    end do ! icol

  end subroutine source
  ! ----------------------------------------------------------
  !
  ! One dimensional interpolation -- return all values along second table dimension
  ! 
  pure function interpolate1D(val, offset, delta, table) result(res)
    ! input
    real(wp), intent(in) :: val,    & ! axis value at which to evaluate table 
                            offset, & ! minimum of table axis
                            delta     ! step size of table axis 
    real(wp), dimension(:,:), & 
              intent(in) :: table ! dimensions (axis, values)
    ! output
    real(wp), dimension(size(table,dim=2)) :: res 
    
    ! local
    real(wp) :: val0 ! fraction index adjusted by offset and delta
    integer :: index ! index term
    real(wp) :: frac ! fractional term
    ! -------------------------------------
    val0 = (val - offset) / delta
    frac = val0 - int(val0) ! get fractional part
    index = min(size(table,dim=1)-1, max(1, int(val0)+1)) ! limit the index range
    res(:) = table(index,:) + frac * (table(index+1,:) - table(index,:))
  end function interpolate1D
 ! ----------------------------------------------------------
 !
 ! interpolation in temperature and eta
 !   First function returns all values along first axis of k (absorption coefficent) table 
 !
  pure function interpolate2D_all(fminor, k, jeta, jtemp) result(res)
    real(wp), dimension(2,2),    intent(in) :: fminor ! interpolation fractions for minor species and continuum
                                       ! index(1) : reference eta level (temperature dependent)
                                       ! index(2) : reference temperature level
    real(wp), dimension(:, :,:), intent(in) :: k ! (gpoint, eta, temp)
    integer,                     intent(in) :: jtemp ! interpolation index for temperature
    integer, dimension(2),       intent(in) :: jeta ! interpolation index for binary species parameter (eta)
    real(wp), dimension(size(k,1))          :: res ! the result

    res(:) =  &
      fminor(1,1) * k(:, jeta(1)  , jtemp  ) + &
      fminor(2,1) * k(:, jeta(1)+1, jtemp  ) + &
      fminor(1,2) * k(:, jeta(2)  , jtemp+1) + &
      fminor(2,2) * k(:, jeta(2)+1, jtemp+1)
  end function interpolate2D_all
  ! ------------
 !   This function returns a single value from a subset (in gpoint) of the k table 
 !
  pure function interpolate2D_1(fminor, k, jeta, jtemp) result(res)
    real(wp), dimension(2,2), intent(in) :: fminor ! interpolation fractions for minor species and continuum
                                       ! index(1) : reference eta level (temperature dependent)
                                       ! index(2) : reference temperature level
    real(wp), dimension(:,:), intent(in) :: k ! (eta, temp)
    integer,                  intent(in) :: jtemp ! interpolation index for temperature
    integer, dimension(2),    intent(in) :: jeta ! interpolation index for binary species parameter (eta)
    real(wp)                             :: res ! the result

    res =  &
      fminor(1,1) * k(jeta(1)  , jtemp  ) + &
      fminor(2,1) * k(jeta(1)+1, jtemp  ) + &
      fminor(1,2) * k(jeta(2)  , jtemp+1) + &
      fminor(2,2) * k(jeta(2)+1, jtemp+1)
  end function interpolate2D_1

  ! ----------------------------------------------------------
  ! interpolation in temperature, pressure, and eta
  pure function interpolate3D(scaling, fmajor, k, jeta, jtemp, jpress) result(res)
    real(wp), dimension(2),     intent(in) :: scaling
    real(wp), dimension(2,2,2), intent(in) :: fmajor ! interpolation fractions for major species
                                                     ! index(1) : reference eta level (temperature dependent)
                                                     ! index(2) : reference pressure level
                                                     ! index(3) : reference temperature level
    real(wp), dimension(:,:,:),  intent(in) :: k ! (eta,temp,press)
    integer,                     intent(in) :: jpress ! interpolation index for pressure
    integer, dimension(2),       intent(in) :: jeta ! interpolation index for binary species parameter (eta)
    integer,                     intent(in) :: jtemp ! interpolation index for temperature
    real(wp)                                :: res ! the result
    ! each code block is for a different reference temperature
    res =  &
      scaling(1) * &
      ( &
        fmajor(1,1,1) * k(jeta(1)  , jpress-1, jtemp  ) + &
        fmajor(2,1,1) * k(jeta(1)+1, jpress-1, jtemp  ) + &
        fmajor(1,2,1) * k(jeta(1)  , jpress  , jtemp  ) + &
        fmajor(2,2,1) * k(jeta(1)+1, jpress  , jtemp  ) &
      ) + &
      scaling(2) * &
      ( &
        fmajor(1,1,2) * k(jeta(2)  , jpress-1, jtemp+1) + &
        fmajor(2,1,2) * k(jeta(2)+1, jpress-1, jtemp+1) + &
        fmajor(1,2,2) * k(jeta(2)  , jpress  , jtemp+1) + &
        fmajor(2,2,2) * k(jeta(2)+1, jpress  , jtemp+1) &
      )
  end function interpolate3D

  ! ----------------------------------------------------------
  ! Compute interpolation coefficients
  ! for calculations of major optical depths, minor optical depths,
  ! continuum optical depths, and Planck fractions
  subroutine interpolation(ncol,nlay,nflav,neta, &
    flavor,press_ref_log,temp_ref,press_ref_log_delta,temp_ref_min,temp_ref_delta,press_ref_trop_log,vmr_ref,nlay_ref, &
    play,tlay,col_gas, &
    jtemp,fmajor,fminor,col_mix,tropo,jeta,jpress)
    ! input dimensions
    integer, intent(in) :: ncol,nlay,nflav,neta

    ! inputs from object
    integer,  dimension(:,:),   intent(in) :: flavor
    real(wp), dimension(:),     intent(in) :: press_ref_log
    real(wp), dimension(:),     intent(in) :: temp_ref
    real(wp), intent(in)                 :: press_ref_log_delta, & 
                                            temp_ref_min, temp_ref_delta, & 
                                            press_ref_trop_log
    real(wp), dimension(:,:,:), intent(in) :: vmr_ref
    integer,                    intent(in) :: nlay_ref

    ! inputs from profile or parent function
    real(wp), dimension(:,:),   intent(in) :: play, tlay
    real(wp), dimension(:,:,:), intent(in) :: col_gas

    ! outputs
    integer,  dimension(ncol,nlay), intent(out) :: jtemp, jpress
    logical,  dimension(ncol,nlay), intent(out) :: tropo
    integer,  dimension(2,    nflav,ncol,nlay), intent(out) :: jeta
    real(wp), dimension(2,    nflav,ncol,nlay), intent(out) :: col_mix
    real(wp), dimension(2,2,2,nflav,ncol,nlay), intent(out) :: fmajor
    real(wp), dimension(2,2,  nflav,ncol,nlay), intent(out) :: fminor
    ! -----------------
    ! local
    real(wp), dimension(ncol,nlay) :: ftemp, fpress ! interpolation fraction for temperature, pressure 
    real(wp) :: locpress ! needed to find location in pressure grid
    real(wp) :: ratio_eta_half ! ratio of vmrs of major species that defines eta=0.5 
                               ! for given flavor and reference temperature level
    real(wp) :: eta, feta   ! binary_species_parameter, interpolation variable for eta
    real(wp) :: loceta ! needed to find location in eta grid
    ! -----------------
    ! local indexes
    integer :: icol, ilay, iflav, igases(2), itropo

    do ilay = 1, nlay
      do icol = 1, ncol
        ! index and factor for temperature interpolation
        jtemp(icol,ilay) = int((tlay(icol,ilay) - (temp_ref_min - temp_ref_delta)) / temp_ref_delta)
        jtemp(icol,ilay) = min(size(temp_ref) - 1, max(1, jtemp(icol,ilay))) ! limit the index range
        ftemp(icol,ilay) = (tlay(icol,ilay) - temp_ref(jtemp(icol,ilay))) / temp_ref_delta

        ! index and factor for pressure interpolation
        locpress = 1._wp + (log(play(icol,ilay)) - press_ref_log(1)) / press_ref_log_delta
        jpress(icol,ilay) = min(nlay_ref-2, max(1, int(locpress)))
        fpress(icol,ilay) = locpress - float(jpress(icol,ilay))

        ! determine if in lower or upper part of atmosphere
        tropo(icol,ilay) = log(play(icol,ilay)) > press_ref_trop_log

        ! itropo = 1 lower atmosphere; itropo = 2 upper atmosphere
        itropo = merge(1,2,tropo(icol,ilay))

        ! loop over implemented combinations of major species
        do iflav = 1, nflav
          igases(:) = flavor(:,iflav)

          ! compute interpolation fractions needed for lower reference temperature level
          ! compute binary species parameter (eta) for flavor and temperature, and associated interpolation index and factor
          ratio_eta_half = vmr_ref(itropo,igases(1),jtemp(icol,ilay)) / vmr_ref(itropo,igases(2),jtemp(icol,ilay))
          col_mix(1,iflav,icol,ilay) = col_gas(icol,ilay,igases(1)) + ratio_eta_half * col_gas(icol,ilay,igases(2))
          eta = merge(col_gas(icol,ilay,igases(1)) / col_mix(1,iflav,icol,ilay), & 
                      0.5_wp, col_mix(1,iflav,icol,ilay) > 2._wp * tiny(col_mix))
          loceta = eta * float(neta-1)
          jeta(1,iflav,icol,ilay) = min(int(loceta)+1, neta-1)
          feta = mod(loceta, 1.0_wp)

          ! compute interpolation fractions needed for minor species and continuum
          fminor(1,1,iflav,icol,ilay) = (1._wp-feta) * (1._wp-ftemp(icol,ilay))
          fminor(2,1,iflav,icol,ilay) =        feta  * (1._wp-ftemp(icol,ilay))
          ! compute interpolation fractions needed for major species
          fmajor(1,1,1,iflav,icol,ilay) = (1._wp-fpress(icol,ilay)) * fminor(1,1,iflav,icol,ilay)
          fmajor(2,1,1,iflav,icol,ilay) = (1._wp-fpress(icol,ilay)) * fminor(2,1,iflav,icol,ilay)
          fmajor(1,2,1,iflav,icol,ilay) =        fpress(icol,ilay)  * fminor(1,1,iflav,icol,ilay)
          fmajor(2,2,1,iflav,icol,ilay) =        fpress(icol,ilay)  * fminor(2,1,iflav,icol,ilay)

          ! compute interpolation fractions needed for lower reference temperature level
          ! compute binary species parameter (eta) for flavor and temperature, and associated interpolation index and factor
          ratio_eta_half = vmr_ref(itropo,igases(1),jtemp(icol,ilay)+1) / vmr_ref(itropo,igases(2),jtemp(icol,ilay)+1)
          col_mix(2,iflav,icol,ilay) = col_gas(icol,ilay,igases(1)) + ratio_eta_half * col_gas(icol,ilay,igases(2))
          eta = merge(col_gas(icol,ilay,igases(1)) / col_mix(2,iflav,icol,ilay), & 
                      0.5_wp, col_mix(2,iflav,icol,ilay) > 2._wp * tiny(col_mix))
          loceta = eta * float(neta-1)
          jeta(2,iflav,icol,ilay) = min(int(loceta)+1, neta-1)
          feta = mod(loceta, 1.0_wp)

          ! compute interpolation fractions needed for minor species and continuum
          fminor(1,2,iflav,icol,ilay) = (1._wp-feta) * ftemp(icol,ilay)
          fminor(2,2,iflav,icol,ilay) =        feta  * ftemp(icol,ilay)
          ! compute interpolation fractions needed for major species
          fmajor(1,1,2,iflav,icol,ilay) = (1._wp-fpress(icol,ilay)) * fminor(1,2,iflav,icol,ilay)
          fmajor(2,1,2,iflav,icol,ilay) = (1._wp-fpress(icol,ilay)) * fminor(2,2,iflav,icol,ilay)
          fmajor(1,2,2,iflav,icol,ilay) =        fpress(icol,ilay)  * fminor(1,2,iflav,icol,ilay)
          fmajor(2,2,2,iflav,icol,ilay) =        fpress(icol,ilay)  * fminor(2,2,iflav,icol,ilay)

        end do ! iflav
      end do ! icol,ilay
    end do 
  end subroutine interpolation

end module mo_gas_optics_kernels