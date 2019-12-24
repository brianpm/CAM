module check_energy

!---------------------------------------------------------------------------------
!
! Energy and mass conservation checking
!
!   1. vertically integrated total energy and water conservation for each
!      column within the physical parameterizations
!
!   2. global mean total energy conservation between the physics output state
!      and the input state on the next time step.
!
!   3. add a globally uniform heating term to account for any change of total energy in 2.
!
!---------------------------------------------------------------------------------

use shr_kind_mod,    only: r8 => shr_kind_r8
use spmd_utils,      only: masterproc
use physconst,       only: gravit, latvap, latice, cpair, cpairv, rearth, omega
use phys_control,    only: phys_getopts
use constituents,    only: pcnst, cnst_get_ind, cnst_name, cnst_get_type_byind

use ppgrid,          only: pcols, pver, begchunk, endchunk
use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_ptend_init

use physics_buffer,  only: physics_buffer_desc, pbuf_add_field, dtype_r8, dyn_time_lvls,   &
                           pbuf_get_index, pbuf_get_field, pbuf_get_chunk, pbuf_set_field, &
                           pbuf_old_tim_idx, pbuf_register_subcol

use gmean_mod,       only: gmean
use time_manager,    only: is_first_step, get_nstep
use subcol_utils,    only: is_subcol_on

use cam_history,     only: addfld, add_default, horiz_only, outfld, hist_fld_active
use cam_logfile,     only: iulog

use cam_abortutils,  only: endrun

implicit none
private
save


public :: &
   check_energy_readnl,       &! read namelist values
   check_energy_register,     &! register fields in physics buffer
   check_energy_init,         &! initialization of module
   check_energy_timestep_init,&! timestep initialization of energy integrals and cumulative 
                               ! boundary fluxes
   check_energy_chng,         &! check changes in integrals against cumulative boundary fluxes
   check_energy_gmean,        &! global means of physics input and output total energy
   check_energy_gmean_diags,  &! global means of energy, water, angular momentum diagnostics
   check_energy_fix,          &! add global mean energy difference as a heating
   check_energy_write,        &! outfld calls for diagnostics
   check_energy_after_phys,   &! save state
   check_energy_get_integrals,&! get energy integrals computed in check_energy_gmean
   check_tracers_init,        &! initialize tracer integrals and cumulative boundary fluxes
   check_tracers_chng,        &! check changes in integrals against cumulative boundary fluxes
   calc_te_and_aam_budgets     ! calculate and output total energy and axial angular momentum
                               ! budget diagnostics

public :: check_tracers_data
type check_tracers_data
   real(r8) :: tracer(pcols,pcnst)     ! initial vertically integrated total (kinetic + static) energy
   real(r8) :: tracer_tnd(pcols,pcnst) ! cumulative boundary flux of total energy
   integer  :: count(pcnst)            ! count of values with significant imbalances
end type check_tracers_data

! Private module data

logical :: print_energy_errors = .false.
logical :: check_energy_diags  = .false.
logical :: check_water_diags   = .false.
logical :: check_tracer_diags  = .false.
logical :: check_aam_diags     = .false.

real(r8) :: teout_glob           ! global mean energy of state output by physics
real(r8) :: teinp_glob           ! global mean energy of state input to physics
real(r8) :: tedif_glob           ! global mean energy difference
real(r8) :: psurf_glob           ! global mean surface pressure
real(r8) :: ptopb_glob           ! global mean top boundary pressure
real(r8) :: heat_glob            ! global mean heating rate

! Physics buffer indices
integer :: &
   teout_idx  = 0,  &
   dtcore_idx = 0,  &
   qini_idx   = 0

! constituent indices
integer :: &
   cldliq_idx = 0,   &
   cldice_idx = 0,   &
   rainqm_idx = 0,   &
   snowqm_idx = 0,   &
   grauqm_idx = 0,   &
   tt_lw_idx = 0

! global mean diagnostics
integer, parameter :: n_energy = 8
integer, parameter :: n_water  = 6
integer, parameter :: n_tracer = 2
integer, parameter :: n_aam    = 8
real(r8), allocatable :: col_energy(:,:,:), gm_energy(:)
real(r8), allocatable :: col_water(:,:,:), gm_water(:)
real(r8), allocatable :: col_tracer(:,:,:), gm_tracer(:)
real(r8), allocatable :: col_aam(:,:,:), gm_aam(:)


! constants for aam calc
real(r8) :: mr_cnst, mo_cnst

!=========================================================================================
contains
!=========================================================================================

subroutine check_energy_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use spmd_utils,      only: mpicom, mstrid=>masterprocid, mpi_logical
   use cam_abortutils,  only: endrun

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: sub = 'check_energy_readnl'

   namelist /check_energy_nl/ print_energy_errors, check_energy_diags, check_water_diags, &
                              check_tracer_diags, check_aam_diags
   !----------------------------------------------------------------------------

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'check_energy_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, check_energy_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(sub//': FATAL: reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

   call mpi_bcast(print_energy_errors, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: print_energy_errors")

   call mpi_bcast(check_energy_diags, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: check_energy_diags")

   call mpi_bcast(check_water_diags, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: check_water_diags")

   call mpi_bcast(check_tracer_diags, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: check_tracer_diags")

   call mpi_bcast(check_aam_diags, 1, mpi_logical, mstrid, mpicom, ierr)
   if (ierr /= 0) call endrun(sub//": FATAL: mpi_bcast: check_aam_diags")

   if (masterproc) then
      write(iulog,*) 'check_energy options:'
      write(iulog,*) '  print_energy_errors =', print_energy_errors
      write(iulog,*) '  check_energy_diags  =', check_energy_diags
      write(iulog,*) '  check_water_diags   =', check_water_diags
      write(iulog,*) '  check_tracer_diags  =', check_tracer_diags
      write(iulog,*) '  check_aam_diags     =', check_aam_diags
   end if

end subroutine check_energy_readnl

!===============================================================================

subroutine check_energy_register()

   call pbuf_add_field('TEOUT', 'global',dtype_r8 , (/pcols,dyn_time_lvls/),      teout_idx)
   call pbuf_add_field('DTCORE','global',dtype_r8,  (/pcols,pver,dyn_time_lvls/),dtcore_idx)

   if (is_subcol_on()) then
      call pbuf_register_subcol('TEOUT', 'phys_register', teout_idx)
      call pbuf_register_subcol('DTCORE', 'phys_register', dtcore_idx)
   end if

end subroutine check_energy_register

!===============================================================================

subroutine check_energy_init()

   ! local variables
   integer :: ierr
   logical :: history_budget, history_waccm
   integer :: history_budget_histfile_num   ! history file number for budget fields
   !----------------------------------------------------------------------------

   mr_cnst = rearth**3/gravit
   mo_cnst = omega*rearth**4/gravit

   call cnst_get_ind('CLDICE', cldice_idx, abort=.false.)
   call cnst_get_ind('CLDLIQ', cldliq_idx, abort=.false.)
   call cnst_get_ind('RAINQM', rainqm_idx, abort=.false.)
   call cnst_get_ind('SNOWQM', snowqm_idx, abort=.false.)
   call cnst_get_ind('GRAUQM', grauqm_idx, abort=.false.)
   call cnst_get_ind('TT_LW' , tt_lw_idx,  abort=.false.)

   ! Q at begining of physics package
   qini_idx = pbuf_get_index('QINI', errcode=ierr)

   ! register history variables
   call addfld('TEINP',  horiz_only,  'A', 'J/m2', 'Total energy of physics input')
   call addfld('TEOUT',  horiz_only,  'A', 'J/m2', 'Total energy of physics output')
   call addfld('TEFIX',  horiz_only,  'A', 'J/m2', 'Total energy after fixer')
   call addfld('EFIX',   horiz_only,  'A', 'W/m2', 'Effective sensible heat flux due to energy fixer')
   call addfld('DTCORE', (/ 'lev' /), 'A', 'K/s' , 'T tendency due to dynamical core')

   ! energy budget diagnostics
   if (check_energy_diags) then
      call addfld('SE_pBF', horiz_only, 'A', 'J/m2','Dry Static Energy before energy fixer')
      call addfld('SE_pBP', horiz_only, 'A', 'J/m2','Dry Static Energy before parameterizations')
      call addfld('SE_pAP', horiz_only, 'A', 'J/m2','Dry Static Energy after parameterizations')
      call addfld('SE_pAM', horiz_only, 'A', 'J/m2','Dry Static Energy after mass correction')

      call addfld('KE_pBF', horiz_only, 'A', 'J/m2','Kinetic Energy before energy fixer')
      call addfld('KE_pBP', horiz_only, 'A', 'J/m2','Kinetic Energy before parameterizations')
      call addfld('KE_pAP', horiz_only, 'A', 'J/m2','Kinetic Energy after parameterizations')
      call addfld('KE_pAM', horiz_only, 'A', 'J/m2','Kinetic Energy after mass correction')

      call add_default('SE_pBF', 1, ' ')
      call add_default('SE_pBP', 1, ' ')
      call add_default('SE_pAP', 1, ' ')
      call add_default('SE_pAM', 1, ' ')
                         
      call add_default('KE_pBF', 1, ' ')
      call add_default('KE_pBP', 1, ' ')
      call add_default('KE_pAP', 1, ' ')
      call add_default('KE_pAM', 1, ' ')

      allocate(col_energy(pcols,begchunk:endchunk,n_energy), gm_energy(n_energy))
   end if

   if (check_water_diags) then
      call addfld('WV_pBP', horiz_only, 'A', 'kg/m2','Total column water vapor before parameterizations')
      call addfld('WV_pAP', horiz_only, 'A', 'kg/m2','Total column water vapor after parameterizations')

      call addfld('WL_pBP', horiz_only, 'A', 'kg/m2','Total column cloud water before parameterizations')
      call addfld('WL_pAP', horiz_only, 'A', 'kg/m2','Total column cloud water after parameterizations')

      call addfld('WI_pBP', horiz_only, 'A', 'kg/m2','Total column cloud ice before parameterizations')
      call addfld('WI_pAP', horiz_only, 'A', 'kg/m2','Total column cloud ice after parameterizations')

      call add_default('WV_pBP', 1, ' ')
      call add_default('WV_pAP', 1, ' ')
      call add_default('WL_pBP', 1, ' ')
      call add_default('WL_pAP', 1, ' ')
      call add_default('WI_pBP', 1, ' ')
      call add_default('WI_pAP', 1, ' ')

      allocate(col_water(pcols,begchunk:endchunk,n_water), gm_water(n_water))
   end if

   if (check_tracer_diags) then
      call addfld('TT_pBP', horiz_only, 'A', 'kg/m2','Total column test tracer before parameterizations')
      call addfld('TT_pAP', horiz_only, 'A', 'kg/m2','Total column test tracer after parameterizations')

      call add_default('TT_pBP', 1, ' ')
      call add_default('TT_pAP', 1, ' ')

      allocate(col_tracer(pcols,begchunk:endchunk,n_tracer), gm_tracer(n_tracer))
   end if

   if (check_aam_diags) then
      ! Axial Angular Momentum diagnostics
      call addfld('MR_pBF', horiz_only, 'A', 'kg*m2/s*rad2',&
               'Total column wind axial angular momentum before energy fixer')
      call addfld('MR_pBP', horiz_only, 'A', 'kg*m2/s*rad2',&
               'Total column wind axial angular momentum before parameterizations')
      call addfld('MR_pAP', horiz_only, 'A', 'kg*m2/s*rad2',&
               'Total column wind axial angular momentum after parameterizations')
      call addfld('MR_pAM', horiz_only, 'A', 'kg*m2/s*rad2',&
               'Total column wind axial angular momentum after mass correction')

      call addfld('MO_pBF', horiz_only, 'A', 'kg*m2/s*rad2',&
               'Total column mass axial angular momentum before energy fixer')
      call addfld('MO_pBP', horiz_only, 'A', 'kg*m2/s*rad2',&
               'Total column mass axial angular momentum before parameterizations')
      call addfld('MO_pAP', horiz_only, 'A', 'kg*m2/s*rad2',&
               'Total column mass axial angular momentum after parameterizations')
      call addfld('MO_pAM', horiz_only, 'A', 'kg*m2/s*rad2',&
               'Total column mass axial angular momentum after mass correction')

      call add_default('MR_pBF', 1, ' ')
      call add_default('MR_pBP', 1, ' ')
      call add_default('MR_pAP', 1, ' ')
      call add_default('MR_pAM', 1, ' ')

      call add_default('MO_pBF', 1, ' ')
      call add_default('MO_pBP', 1, ' ')
      call add_default('MO_pAP', 1, ' ')
      call add_default('MO_pAM', 1, ' ')

      allocate(col_aam(pcols,begchunk:endchunk,n_aam), gm_aam(n_aam))
   end if

   call phys_getopts(history_budget_out = history_budget, &
                     history_budget_histfile_num_out = history_budget_histfile_num, &
                     history_waccm_out = history_waccm )

   if ( history_budget ) then
      call add_default ('DTCORE', history_budget_histfile_num, ' ')
   end if
   if ( history_waccm ) then
      call add_default ('DTCORE', 1, ' ')
   end if

end subroutine check_energy_init

!=========================================================================================

subroutine check_energy_timestep_init(state, pbuf, col_type)

   ! At begining of physics package run, compute values of energy and water integrals.

   ! Arguments
   type(physics_state),   intent(inout)    :: state
   type(physics_buffer_desc), pointer      :: pbuf(:)
   integer, optional                       :: col_type  ! Flag inidicating whether using grid or subcolumns

   ! Local variables
   integer :: lchnk              ! chunk index
   integer :: ncol               ! number of columns
   integer :: i, k               ! column, level indices

   real(r8) :: mass(state%ncol,pver) ! air mass (wet)
   real(r8) :: ke(state%ncol)    ! vertical integral of kinetic energy
   real(r8) :: se(state%ncol)    ! vertical integral of static energy
   real(r8) :: wv(state%ncol)    ! vertical integral of water (vapor)
   real(r8) :: wl(state%ncol)    ! vertical integral of water (liquid)
   real(r8) :: wi(state%ncol)    ! vertical integral of water (ice)
   real(r8) :: tmp(state%ncol)
   !----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   do k = 1, pver
      do i = 1, ncol
         mass(i,k) = state%pdel(i,k)/gravit
      end do
   end do

   call column_energy(state%psetcols, lchnk, ncol, state%u(:ncol,:), state%v(:ncol,:), &
      state%t(:ncol,:), mass, state%phis(:ncol), state%ps(:ncol), ke, se)

   call column_mass(ncol, state%q(:ncol,:,1), mass, wv)

   wl = 0._r8
   if (cldliq_idx > 0) then
      call column_mass(ncol, state%q(:ncol,:,cldliq_idx), mass, wl)
   end if
   if (rainqm_idx > 0) then
      call column_mass(ncol, state%q(:ncol,:,rainqm_idx), mass, tmp)
      wl = wl + tmp
   end if

   wi = 0._r8
   if (cldice_idx > 0) then
      call column_mass(ncol, state%q(:ncol,:,cldice_idx), mass, wi)
   end if
   if (snowqm_idx > 0) then
      call column_mass(ncol, state%q(:ncol,:,snowqm_idx), mass, tmp)
      wi = wi + tmp
   end if
   if (grauqm_idx > 0) then
      call column_mass(ncol, state%q(:ncol,:,grauqm_idx), mass, tmp)
      wi = wi + tmp
   end if

   ! Compute vertical integrals of frozen static energy and total water.
   do i = 1, ncol
      state%te_ini(i) = se(i) + ke(i) + (latvap+latice)*wv(i) + latice*wl(i)
      state%tw_ini(i) = wv(i) + wl(i) + wi(i)

      state%te_cur(i) = state%te_ini(i)
      state%tw_cur(i) = state%tw_ini(i)
   end do

   state%count = 0

   ! initialize physics buffer
   if (is_first_step()) then
      call pbuf_set_field(pbuf, teout_idx, state%te_ini, col_type=col_type)
   end if

end subroutine check_energy_timestep_init

!=========================================================================================

subroutine check_energy_chng(state, tend, name, ztodt,        &
                             flx_vap, flx_cnd, flx_ice, flx_sen)

   ! Check that the energy and water change matches the boundary fluxes
   ! ** Note that the precip and ice fluxes are in precip units (m/s). **

   ! Arguments
   type(physics_state), intent(inout) :: state
   type(physics_tend ), intent(inout) :: tend
   character*(*),       intent(in)    :: name    ! parameterization name for fluxes
   real(r8), intent(in) :: ztodt       ! 2 delta t (model time increment)
   real(r8), intent(in) :: flx_vap(:)  ! (pcols) - boundary flux of vapor (kg/m2/s)
   real(r8), intent(in) :: flx_cnd(:)  ! (pcols) -boundary flux of liquid+ice (m/s) (precip?)
   real(r8), intent(in) :: flx_ice(:)  ! (pcols) -boundary flux of ice        (m/s) (snow?)
   real(r8), intent(in) :: flx_sen(:)  ! (pcols) -boundary flux of sensible heat (w/m2)

   ! Local variables
   integer :: nstep                  ! current timestep number (for log output)
   integer :: lchnk                  ! chunk identifier
   integer :: ncol                   ! number of atmospheric columns
   integer :: i, k                   ! column, level indices

   real(r8) :: mass(state%ncol,pver) ! air mass (wet)
   real(r8) :: te_xpd(state%ncol)    ! expected value (f0 + dt*boundary_flux)
   real(r8) :: te_dif(state%ncol)    ! energy of input state - original energy
   real(r8) :: te_tnd(state%ncol)    ! tendency from last process
   real(r8) :: te_rer(state%ncol)    ! relative error in energy column

   real(r8) :: tw_xpd(state%ncol)    ! expected value (w0 + dt*boundary_flux)
   real(r8) :: tw_dif(state%ncol)    ! tw_inp - original water
   real(r8) :: tw_tnd(state%ncol)    ! tendency from last process
   real(r8) :: tw_rer(state%ncol)    ! relative error in water column

   real(r8) :: ke(state%ncol)        ! vertical integral of kinetic energy
   real(r8) :: se(state%ncol)        ! vertical integral of static energy
   real(r8) :: wv(state%ncol)        ! vertical integral of water (vapor)
   real(r8) :: wl(state%ncol)        ! vertical integral of water (liquid)
   real(r8) :: wi(state%ncol)        ! vertical integral of water (ice)
   real(r8) :: tmp(state%ncol)

   real(r8) :: te(state%ncol)        ! vertical integral of total energy
   real(r8) :: tw(state%ncol)        ! vertical integral of total water
   !----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   do k = 1, pver
      do i = 1, ncol
         mass(i,k) = state%pdel(i,k)/gravit
      end do
   end do

   call column_energy(state%psetcols, lchnk, ncol, state%u(:ncol,:), state%v(:ncol,:), &
      state%t(:ncol,:), mass, state%phis(:ncol), state%ps(:ncol), ke, se)

   call column_mass(ncol, state%q(:ncol,:,1), mass, wv)

   wl = 0._r8
   if (cldliq_idx > 0) then
      call column_mass(ncol, state%q(:ncol,:,cldliq_idx), mass, wl)
   end if
   if (rainqm_idx > 0) then
      call column_mass(ncol, state%q(:ncol,:,rainqm_idx), mass, tmp)
      wl = wl + tmp
   end if

   wi = 0._r8
   if (cldice_idx > 0) then
      call column_mass(ncol, state%q(:ncol,:,cldice_idx), mass, wi)
   end if
   if (snowqm_idx > 0) then
      call column_mass(ncol, state%q(:ncol,:,snowqm_idx), mass, tmp)
      wi = wi + tmp
   end if
   if (grauqm_idx > 0) then
      call column_mass(ncol, state%q(:ncol,:,grauqm_idx), mass, tmp)
      wi = wi + tmp
   end if

   ! Compute vertical integrals of frozen static energy and total water.
   do i = 1, ncol
      te(i) = se(i) + ke(i) + (latvap+latice)*wv(i) + latice*wl(i)
      tw(i) = wv(i) + wl(i) + wi(i)
   end do

   ! compute expected values and tendencies
   do i = 1, ncol
      ! change in static energy and total water
      te_dif(i) = te(i) - state%te_cur(i)
      tw_dif(i) = tw(i) - state%tw_cur(i)

      ! expected tendencies from boundary fluxes for last process
      te_tnd(i) = flx_vap(i)*(latvap+latice) - (flx_cnd(i) - flx_ice(i))*1000._r8*latice + flx_sen(i)
      tw_tnd(i) = flx_vap(i) - flx_cnd(i) *1000._r8

      ! cummulative tendencies from boundary fluxes
      tend%te_tnd(i) = tend%te_tnd(i) + te_tnd(i)
      tend%tw_tnd(i) = tend%tw_tnd(i) + tw_tnd(i)

      ! expected new values from previous state plus boundary fluxes
      te_xpd(i) = state%te_cur(i) + te_tnd(i)*ztodt
      tw_xpd(i) = state%tw_cur(i) + tw_tnd(i)*ztodt

      ! relative error, expected value - input state / previous state
      te_rer(i) = (te_xpd(i) - te(i)) / state%te_cur(i)
   end do

   ! relative error for total water (allow for dry atmosphere)
   tw_rer = 0._r8
   where (state%tw_cur(:ncol) > 0._r8)
      tw_rer(:ncol) = (tw_xpd(:ncol) - tw(:ncol)) / state%tw_cur(:ncol)
   end where

   ! error checking
   if (print_energy_errors) then
      nstep = get_nstep()
      if (any(abs(te_rer(1:ncol)) > 1.E-14_r8 .or. abs(tw_rer(1:ncol)) > 1.E-10_r8)) then
         do i = 1, ncol
            ! the relative error threshold for the water budget has been reduced to 1.e-10
            ! to avoid messages generated by QNEG3 calls
            if (abs(te_rer(i)) > 1.E-14_r8 ) then
               state%count = state%count + 1
               write(iulog,*) "significant energy conservation error after ", name,        &
                      " count", state%count, " nstep", nstep, "chunk", lchnk, "col", i
               write(iulog,*) te(i),te_xpd(i),te_dif(i),tend%te_tnd(i)*ztodt,  &
                      te_tnd(i)*ztodt,te_rer(i)
            endif
            if ( abs(tw_rer(i)) > 1.E-10_r8) then
               state%count = state%count + 1
               write(iulog,*) "significant water conservation error after ", name,        &
                      " count", state%count, " nstep", nstep, "chunk", lchnk, "col", i
               write(iulog,*) tw(i),tw_xpd(i),tw_dif(i),tend%tw_tnd(i)*ztodt,  &
                      tw_tnd(i)*ztodt,tw_rer(i)
            end if
         end do
      end if
   end if

   ! copy new value to state
   do i = 1, ncol
      state%te_cur(i) = te(i)
      state%tw_cur(i) = tw(i)
   end do

end subroutine check_energy_chng

!=========================================================================================

subroutine check_energy_gmean(state, pbuf2d, dtime)

   ! Compute global mean total energy of physics input and output states

   ! Arguments
   type(physics_state),       intent(in) :: state(begchunk:endchunk)
   type(physics_buffer_desc), pointer    :: pbuf2d(:,:)
   real(r8),                  intent(in) :: dtime

   ! Local variables
   integer :: ncol                      ! number of active columns
   integer :: lchnk                     ! chunk index
   integer :: nstep                     ! current timestep number

   real(r8) :: te(pcols,begchunk:endchunk,3)
                                         ! total energy of input/output states (copy)
   real(r8) :: te_glob(3)               ! global means of total energy
   real(r8), pointer :: teout(:)
   !----------------------------------------------------------------------------

   nstep = get_nstep()

   ! Copy total energy out of input and output states
   do lchnk = begchunk, endchunk
      ncol = state(lchnk)%ncol
      ! input energy
      te(:ncol,lchnk,1) = state(lchnk)%te_ini(:ncol)
      ! output energy
      call pbuf_get_field(pbuf_get_chunk(pbuf2d,lchnk), teout_idx, teout)

      te(:ncol,lchnk,2) = teout(1:ncol)
      ! surface pressure for heating rate
      te(:ncol,lchnk,3) = state(lchnk)%pint(:ncol,pver+1)
   end do

   ! Compute global means of input and output energies and of
   ! surface pressure for heating rate (assume uniform ptop)
   call gmean(te, te_glob, 3)

   if (begchunk .le. endchunk) then
      teinp_glob = te_glob(1)
      teout_glob = te_glob(2)
      psurf_glob = te_glob(3)
      ptopb_glob = state(begchunk)%pint(1,1)

      ! Global mean total energy difference
      tedif_glob =  teinp_glob - teout_glob
      heat_glob  = -tedif_glob/dtime * gravit / (psurf_glob - ptopb_glob)

      if (masterproc) then
         write(iulog,'(1x,a9,1x,i8,4(1x,e25.17))') "nstep, te", nstep, teinp_glob, teout_glob, heat_glob, psurf_glob
      end if
   else
      heat_glob = 0._r8
   end if

end subroutine check_energy_gmean

!===============================================================================

subroutine check_energy_gmean_diags()

   ! Compute global means for energy, water, and angular momentum diagnostics

   ! Local variables
   integer :: nstep                     ! current timestep number
   !----------------------------------------------------------------------------

   nstep = get_nstep()

   if (check_energy_diags) then
      call gmean(col_energy, gm_energy, n_energy)
      if (masterproc) then
         write(iulog,'(1x,a9,1x,i8,4(1x,e25.17))') "pBF nstep, SE, KE", nstep, gm_energy(1), gm_energy(2)
         write(iulog,'(1x,a9,1x,i8,4(1x,e25.17))') "pBP nstep, SE, KE", nstep, gm_energy(3), gm_energy(4)
         write(iulog,'(1x,a9,1x,i8,4(1x,e25.17))') "pAP nstep, SE, KE", nstep, gm_energy(5), gm_energy(6)
         write(iulog,'(1x,a9,1x,i8,4(1x,e25.17))') "pAM nstep, SE, KE", nstep, gm_energy(7), gm_energy(8)
      end if
   end if

end subroutine check_energy_gmean_diags

!===============================================================================

subroutine check_energy_fix(state, ptend, eshflx)

   ! Add heating rate required for global mean total energy conservation

   ! Arguments
   type(physics_state), intent(in)  :: state
   type(physics_ptend), intent(out) :: ptend
   real(r8),            intent(out) :: eshflx(pcols)  ! effective sensible heat flux

   ! Local variables
   integer :: i
   integer :: lchnk, ncol
   !----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol = state%ncol

   call physics_ptend_init(ptend, state%psetcols, 'chkenergyfix', ls=.true.)

#if ( defined OFFLINE_DYN )
   ! disable the energy fix for offline driver
   heat_glob = 0._r8
#endif
   ! add (-) global mean total energy difference as heating
   ptend%s(:ncol,:pver) = heat_glob

   ! compute effective sensible heat flux
   do i = 1, ncol
      eshflx(i) = heat_glob * (state%pint(i,pver+1) - state%pint(i,1)) / gravit
   end do

   call outfld('EFIX', eshflx, pcols, lchnk)

end subroutine check_energy_fix

!=========================================================================================

subroutine check_energy_write(state, pbuf, ztodt)

   ! Arguments
   type(physics_state),       intent(in) :: state
   type(physics_buffer_desc), pointer    :: pbuf(:)
   real(r8),                  intent(in) :: ztodt   ! physics package timestep

   ! local variables
   integer :: itim_old
   integer :: lchnk
   integer :: ncol
   integer :: nstep

   real(r8), pointer :: teout(:)
   real(r8), pointer :: dtcore(:,:)
   !----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol
   nstep = get_nstep()

   ! Associate pointers with physics buffer fields
   itim_old = pbuf_old_tim_idx()

   ! Associate pointers with physics buffer fields
   call pbuf_get_field(pbuf, teout_idx, teout, (/1,itim_old/), (/pcols,1/))
   call pbuf_get_field(pbuf, dtcore_idx, dtcore, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   call outfld('TEOUT', teout       , pcols, lchnk   )
   call outfld('TEINP', state%te_ini, pcols, lchnk   )
   call outfld('TEFIX', state%te_cur, pcols, lchnk   )

   ! T tendency due to dynamics
   if (nstep > dyn_time_lvls-1) then
      dtcore(:ncol,:pver) = (state%t(:ncol,:pver) - dtcore(:ncol,:pver))/ztodt
      call outfld( 'DTCORE', dtcore, pcols, lchnk )
   end if

end subroutine check_energy_write

!=========================================================================================

subroutine check_energy_after_phys(state, pbuf)

   ! Arguments
   type(physics_state),       intent(in) :: state
   type(physics_buffer_desc), pointer    :: pbuf(:)

   ! local variables
   integer :: itim_old
   integer :: lchnk
   integer :: ncol

   real(r8), pointer :: teout(:)
   real(r8), pointer :: dtcore(:,:)
   !----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   ! Associate pointers with physics buffer fields
   itim_old = pbuf_old_tim_idx()

   ! Associate pointers with physics buffer fields
   call pbuf_get_field(pbuf, teout_idx, teout, (/1,itim_old/), (/pcols,1/))
   call pbuf_get_field(pbuf, dtcore_idx, dtcore, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

   ! Save total enery after physics for energy conservation checks
   teout = state%te_cur

   ! Store T in buffer for use in computing dynamics T-tendency in next timestep
   dtcore(:ncol,:) = state%t(:ncol,:)

end subroutine check_energy_after_phys

!=========================================================================================

subroutine check_energy_get_integrals( tedif_glob_out, heat_glob_out )

   ! Return energy integrals

   real(r8), intent(out), optional :: tedif_glob_out
   real(r8), intent(out), optional :: heat_glob_out
   !----------------------------------------------------------------------------

   if (present(tedif_glob_out)) then
      tedif_glob_out = tedif_glob
   end if

   if (present(heat_glob_out)) then
      heat_glob_out = heat_glob
   end if

end subroutine check_energy_get_integrals

!=========================================================================================

subroutine check_tracers_init(state, tracerint)

   ! Compute initial values of tracer column burdens
   ! zero cumulative tendencies

   ! Arguments
   type(physics_state),      intent(in)  :: state
   type(check_tracers_data), intent(out) :: tracerint

   ! Local storage

   integer :: ncol
   integer :: i, k, m

   real(r8) :: tr(pcols)            ! vertical integral of tracer
   real(r8) :: trpdel(pcols, pver)  ! pdel for tracer
   !----------------------------------------------------------------------------

   ncol  = state%ncol

   do m = 1, pcnst
      ! dont process water substances.  They are checked in check_energy.
      if ( any(m == (/ 1, cldliq_idx, cldice_idx, &
                          rainqm_idx, snowqm_idx, grauqm_idx /)) ) cycle

      if (cnst_get_type_byind(m).eq.'dry') then
         trpdel(:ncol,:) = state%pdeldry(:ncol,:)
      else
         trpdel(:ncol,:) = state%pdel(:ncol,:)
      endif

      ! Compute vertical integrals of tracer and save in tracerint object
      tr = 0._r8
      do k = 1, pver
         do i = 1, ncol
            tr(i) = tr(i) + state%q(i,k,m)*trpdel(i,k)/gravit
         end do
      end do

      do i = 1, ncol
         tracerint%tracer(i,m) = tr(i)
      end do

      ! zero cummulative boundary fluxes
      tracerint%tracer_tnd(:ncol,m) = 0._r8

      tracerint%count(m) = 0

   end do

end subroutine check_tracers_init

!===============================================================================

subroutine check_tracers_chng(state, tracerint, name, ztodt, cflx)

   ! Check that the tracers and water change matches the boundary fluxes
   ! these checks are not save when there are tracers transformations, as
   ! they only check to see whether a mass change in the column is
   ! associated with a flux

   ! Arguments

   type(physics_state)    , intent(in   ) :: state
   type(check_tracers_data), intent(inout) :: tracerint! tracers integrals and boundary fluxes
   character*(*),intent(in) :: name               ! parameterization name for fluxes
   real(r8), intent(in   ) :: ztodt               ! 2 delta t (model time increment)
   real(r8), intent(in   ) :: cflx(pcols,pcnst)       ! boundary flux of tracers       (kg/m2/s)

   ! Local storage

   real(r8) :: tracer_inp(pcols,pcnst)                   ! total tracer of new (input) state
   real(r8) :: tracer_xpd(pcols,pcnst)                   ! expected value (w0 + dt*boundary_flux)
   real(r8) :: tracer_dif(pcols,pcnst)                   ! tracer_inp - original tracer
   real(r8) :: tracer_tnd(pcols,pcnst)                   ! tendency from last process
   real(r8) :: tracer_rer(pcols,pcnst)                   ! relative error in tracer column

   real(r8) :: tr(pcols)                           ! vertical integral of tracer
   real(r8) :: trpdel(pcols, pver)                       ! pdel for tracer

   integer lchnk                                  ! chunk identifier
   integer ncol                                   ! number of atmospheric columns
   integer  i,k                                   ! column, level indices
   integer :: m                            ! tracer index
   integer :: nstep               ! current timestep number
   character(len=8) :: tracname   ! tracername
   !----------------------------------------------------------------------------

   nstep = get_nstep()

   lchnk = state%lchnk
   ncol  = state%ncol

   do m = 1,pcnst

      if ( any(m == (/ 1, cldliq_idx, cldice_idx, &
                          rainqm_idx, snowqm_idx, grauqm_idx /)) ) cycle ! dont process water substances
                                                                 ! they are checked in check_energy

      tracname = cnst_name(m)
      if (cnst_get_type_byind(m).eq.'dry') then
         trpdel(:ncol,:) = state%pdeldry(:ncol,:)
      else
         trpdel(:ncol,:) = state%pdel(:ncol,:)
      endif

      ! Compute vertical integrals tracers
      tr = 0._r8
      do k = 1, pver
         do i = 1, ncol
            tr(i) = tr(i) + state%q(i,k,m)*trpdel(i,k)/gravit
         end do
      end do

      ! Compute vertical integrals of tracer
      do i = 1, ncol
         tracer_inp(i,m) = tr(i)
      end do

      ! compute expected values and tendencies
      do i = 1, ncol
         ! change in tracers
         tracer_dif(i,m) = tracer_inp(i,m) - tracerint%tracer(i,m)

         ! expected tendencies from boundary fluxes for last process
         tracer_tnd(i,m) = cflx(i,m)

         ! cummulative tendencies from boundary fluxes
         tracerint%tracer_tnd(i,m) = tracerint%tracer_tnd(i,m) + tracer_tnd(i,m)

         ! expected new values from original values plus boundary fluxes
         tracer_xpd(i,m) = tracerint%tracer(i,m) + tracerint%tracer_tnd(i,m)*ztodt

         ! relative error, expected value - input value / original
         tracer_rer(i,m) = (tracer_xpd(i,m) - tracer_inp(i,m)) / tracerint%tracer(i,m)
      end do

      ! final loop for error checking
      if ( maxval(tracer_rer) > 1.E-14_r8 ) then
         write(iulog,*) "CHECK_TRACERS TRACER large rel error"
         write(iulog,*) tracer_rer
      endif

      do i = 1, ncol
         ! error messages
         if (abs(tracer_rer(i,m)) > 1.E-14_r8 ) then
            tracerint%count = tracerint%count + 1
            write(iulog,*) "CHECK_TRACERS TRACER significant conservation error after ", name,        &
                  " count", tracerint%count, " nstep", nstep, "chunk", lchnk, "col",i
            write(iulog,*)' process name, tracname, index ',  name, tracname, m
            write(iulog,*)" input integral              ",tracer_inp(i,m)
            write(iulog,*)" expected integral           ", tracer_xpd(i,m)
            write(iulog,*)" input - inital integral     ",tracer_dif(i,m)
            write(iulog,*)" cumulative tend      ",tracerint%tracer_tnd(i,m)*ztodt
            write(iulog,*)" process tend         ",tracer_tnd(i,m)*ztodt
            write(iulog,*)" relative error       ",tracer_rer(i,m)
            call endrun()
         end if
      end do
   end do

end subroutine check_tracers_chng

!=========================================================================================

subroutine calc_te_and_aam_budgets(state, pbuf, suf)

   ! Compute column integrals for energy, water, and angular momentum budgets.  Output
   ! column integrals to history file and save values for later calculation of global
   ! mean values.

   ! Arguments
   type(physics_state),       intent(in) :: state
   type(physics_buffer_desc), pointer    :: pbuf(:)
   character*(*),             intent(in) :: suf     ! suffix for "outfld" names

   ! Local storage
   integer :: lchnk   ! chunk identifier
   integer :: ncol    ! number of atmospheric columns
   integer :: i, k    ! column, level indices
   integer :: e_idx   ! start index for storing energy diags
   integer :: w_idx   ! start index for storing water diags
   integer :: t_idx   ! start index for storing tracer diags
   integer :: m_idx   ! start index for storing aam diags

   real(r8) :: mass(state%ncol,pver)    ! air mass (wet)
   real(r8) :: drymass(state%ncol,pver) ! dry air mass
   real(r8) :: fdq(state%ncol)          ! moist mass conversion
   real(r8) :: pdeladj(state%ncol,pver) ! pdel adjusted for change in q by moist physics
   real(r8), pointer :: qini(:,:)       ! q at begining of physics package
   !-----------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   ! The following case selector sets the start index for storing column integrals in
   ! arrays to be passed to global mean calculations.
   ! Note that the 'pAM' diagnostics are handled during the call for 'pAP'

   select case (trim(suf))
   case('pBF')
      e_idx = 0
      w_idx = -1  ! no water masses before fixer
      t_idx = -1  ! no tracer masses before fixer
      m_idx = 0
   case('pBP')
      e_idx = 2   ! increment by 2 for SE and KE
      w_idx = 0
      t_idx = 0
      m_idx = 2   ! increment by 2 for MR and MO
   case('pAP')
      e_idx = 4   ! increment by 2 for SE and KE
      w_idx = 3   ! increment by 3 for WV, WL, WI
      t_idx = 1   ! increment by 1 for TT
      m_idx = 4   ! increment by 2 for MR and MO
   case default
      call endrun('calc_te_and_aam_budgets: bad suffix:'//trim(suf))
   end select

   if (check_energy_diags .or. check_water_diags .or. check_tracer_diags) then
      do k = 1, pver
         do i = 1, ncol
            mass(i,k) = state%pdel(i,k)/gravit
         end do
      end do
   end if

   if (check_energy_diags) then
      call column_energy(state%psetcols, lchnk, ncol, state%u(:ncol,:), state%v(:ncol,:), &
                         state%t(:ncol,:), mass, state%phis(:ncol), state%ps(:ncol),      &
                         col_energy(:ncol,lchnk,e_idx+1), col_energy(:ncol,lchnk,e_idx+2))

      call outfld('SE_'//trim(suf), col_energy(:ncol,lchnk,e_idx+1), ncol, lchnk)
      call outfld('KE_'//trim(suf), col_energy(:ncol,lchnk,e_idx+2), ncol, lchnk)
   end if

   if (check_water_diags .and. w_idx>-1) then

      call column_mass(ncol, state%q(:ncol,:,1), mass, col_water(:ncol,lchnk,w_idx+1))
      call outfld('WV_'//trim(suf), col_water(:ncol,lchnk,w_idx+1), ncol, lchnk)

      
      col_water(:ncol,lchnk,w_idx+2) = 0._r8
      if (cldliq_idx > 0) then
         call column_mass(ncol, state%q(:ncol,:,cldliq_idx), mass, col_water(:ncol,lchnk,w_idx+2))
      end if
      call outfld('WL_'//trim(suf), col_water(:ncol,lchnk,w_idx+2), ncol, lchnk)
      
      col_water(:ncol,lchnk,w_idx+3) = 0._r8
      if (cldice_idx > 0) then
         call column_mass(ncol, state%q(:ncol,:,cldice_idx), mass, col_water(:ncol,lchnk,w_idx+3))
      end if
      call outfld('WL_'//trim(suf), col_water(:ncol,lchnk,w_idx+3), ncol, lchnk)
   end if
      
   if (check_tracer_diags .and. t_idx>-1) then

      col_tracer(:ncol,lchnk,t_idx+1) = 0._r8

      ! Test tracers are dry, so need dry air mass
      do k = 1, pver
         do i = 1, ncol
            drymass(i,k) = state%pdeldry(i,k)/gravit
         end do
      end do

      if (tt_lw_idx > 0) then
         call column_mass(ncol, state%q(:ncol,:,tt_lw_idx), drymass, col_tracer(:ncol,lchnk,t_idx+1))
      end if
      call outfld('TT_'//trim(suf), col_tracer(:ncol,lchnk,t_idx+1), ncol, lchnk)
   end if

   ! For the call "after physics" also compute the SE and KE integrals after updating
   ! the moist air mass to account for water vapor changes due to physics package.
   if (trim(suf) == 'pAP') then

      e_idx = 6
      
      if (check_energy_diags) then

         call pbuf_get_field(pbuf, qini_idx, qini)

         ! update moist mass
         do k = 1, pver
            ! fdq is (moist air mass after physics)/(moist air mass before physics)
            fdq(:ncol) = 1._r8 + state%q(:ncol,k,1) - qini(:ncol,k)
            mass(:ncol,k) = state%pdel(:ncol,k)*fdq(:ncol)/gravit
         end do

         call column_energy(state%psetcols, lchnk, ncol, state%u(:ncol,:), state%v(:ncol,:), &
            state%t(:ncol,:), mass, state%phis(:ncol), state%ps(:ncol), &
            col_energy(:ncol,lchnk,e_idx+1), col_energy(:ncol,lchnk,e_idx+2))

         call outfld('SE_pAM', col_energy(:ncol,lchnk,e_idx+1), ncol, lchnk)
         call outfld('KE_pAM', col_energy(:ncol,lchnk,e_idx+2), ncol, lchnk)

      end if
   end if

   ! Axial angular momentum diagnostics
   !
   ! Code follows
   !
   ! Lauritzen et al., (2014): Held-Suarez simulations with the Community Atmosphere Model
   ! Spectral Element (CAM-SE) dynamical core: A global axial angularmomentum analysis using Eulerian
   ! and floating Lagrangian vertical coordinates. J. Adv. Model. Earth Syst. 6,129-140,
   ! doi:10.1002/2013MS000268
   !
   ! MR is equation (6) without \Delta A and sum over areas (areas are in units of radians**2)
   ! MO is equation (7) without \Delta A and sum over areas (areas are in units of radians**2)

   if (check_aam_diags) then

      call column_aam(ncol, state%lat(:ncol), state%u(:ncol,:), state%pdel(:ncol,:), &
         col_aam(:ncol,lchnk,m_idx+1), col_aam(:ncol,lchnk,m_idx+2))
      call outfld('MR_'//trim(suf), col_aam(:ncol,lchnk,m_idx+1), pcols, lchnk)
      call outfld('MO_'//trim(suf), col_aam(:ncol,lchnk,m_idx+2), pcols, lchnk)

   end if

   ! For the call "after physics" also compute the momentum integrals after updating
   ! the moist air mass to account for water vapor changes due to physics package.
   if (trim(suf) == 'pAP') then

      m_idx = 6
      
      if (check_aam_diags) then

         call pbuf_get_field(pbuf, qini_idx, qini)

         ! update pdel
         do k = 1, pver
            ! fdq is (moist air mass after physics)/(moist air mass before physics)
            fdq(:ncol) = 1._r8 + state%q(:ncol,k,1) - qini(:ncol,k)
            pdeladj(:ncol,k) = state%pdel(:ncol,k)*fdq(:ncol)
         end do

         call column_aam(ncol, state%lat(:ncol), state%u(:ncol,:), pdeladj(:ncol,:), &
            col_aam(:ncol,lchnk,m_idx+1), col_aam(:ncol,lchnk,m_idx+2))
         call outfld('MR_pAM', col_aam(:ncol,lchnk,m_idx+1), pcols, lchnk)
         call outfld('MO_pAM', col_aam(:ncol,lchnk,m_idx+2), pcols, lchnk)

      end if
   end if

end subroutine calc_te_and_aam_budgets

!=========================================================================================
! Private
!=========================================================================================

subroutine column_energy(psetcols, lchnk, ncol, u, v, t, mass, phis, ps, ke, se)

   ! arguments
   integer,   intent(in)  :: psetcols
   integer,   intent(in)  :: lchnk
   integer,   intent(in)  :: ncol
   real(r8),  intent(in)  :: u(ncol,pver)
   real(r8),  intent(in)  :: v(ncol,pver)
   real(r8),  intent(in)  :: t(ncol,pver)
   real(r8),  intent(in)  :: mass(ncol,pver) ! air mass (wet)
   real(r8),  intent(in)  :: phis(ncol)
   real(r8),  intent(in)  :: ps(ncol)
   real(r8),  intent(out) :: ke(ncol)        ! vertical integral of kinetic energy
   real(r8),  intent(out) :: se(ncol)        ! vertical integral of static energy

   ! local variables
   integer  :: i, k
   real(r8), allocatable :: cpairv_loc(:,:)
   !----------------------------------------------------------------------------

   ! cpairv_loc needs to be allocated to a size which matches state and ptend
   ! If psetcols == pcols, cpairv is the correct size and just copy into cpairv_loc
   ! If psetcols > pcols and all cpairv match cpair, then assign the constant cpair
   allocate(cpairv_loc(psetcols,pver))
   if (psetcols == pcols) then
      cpairv_loc(:,:) = cpairv(:,:,lchnk)
   else if (psetcols > pcols .and. all(cpairv(:,:,lchnk) == cpair)) then
      cpairv_loc(:,:) = cpair
   else
      call endrun('check_energy::column_integrals:'//&
                  ' cpairv is not allowed to vary when subcolumns are turned on')
   end if

   ke = 0._r8
   se = 0._r8
   do k = 1, pver
      do i = 1, ncol
         ke(i) = ke(i) + 0.5_r8*(u(i,k)**2 + v(i,k)**2)*mass(i,k)
         se(i) = se(i) + t(i,k)*cpairv_loc(i,k)*mass(i,k)
      end do
   end do
   do i = 1, ncol
      se(i) = se(i) + phis(i)*ps(i)/gravit
   end do

   deallocate(cpairv_loc)

end subroutine column_energy

!=========================================================================================

subroutine column_mass(ncol, q, mass, qmass)

   ! arguments
   integer,   intent(in)  :: ncol
   real(r8),  intent(in)  :: q(ncol,pver)    ! constituent mmr
   real(r8),  intent(in)  :: mass(ncol,pver) ! air mass (consistent w/ constituent mmr basis)
   real(r8),  intent(out) :: qmass(ncol)     ! constituent column mass

   ! local variables
   integer  :: i, k
   !----------------------------------------------------------------------------

   qmass = 0._r8
   do k = 1, pver
      do i = 1, ncol
         qmass(i) = qmass(i) + q(i,k)*mass(i,k)
      end do
   end do

end subroutine column_mass

!=========================================================================================

subroutine column_aam(ncol, lat, u, pdel, mr, mo)

   ! arguments
   integer,   intent(in)  :: ncol
   real(r8),  intent(in)  :: lat(ncol)      ! latitude (radians)
   real(r8),  intent(in)  :: u(ncol,pver)
   real(r8),  intent(in)  :: pdel(ncol,pver)
   real(r8),  intent(out) :: mr(ncol)
   real(r8),  intent(out) :: mo(ncol)

   ! local variables
   integer  :: i, k
   real(r8) :: cos_lat(ncol)
   real(r8) :: cos_lat2(ncol)
   !----------------------------------------------------------------------------

   cos_lat = cos(lat)
   cos_lat2 = cos_lat**2

   mr = 0._r8
   mo = 0._r8
   do k = 1, pver
      do i = 1, ncol
         mr(i) = mr(i) + u(i,k)*pdel(i,k)*cos_lat(i)
         mo(i) = mo(i) + pdel(i,k)*cos_lat2(i)
      end do
   end do
   mr = mr_cnst*mr
   mo = mo_cnst*mo

end subroutine column_aam

!=========================================================================================
!=========================================================================================

end module check_energy
