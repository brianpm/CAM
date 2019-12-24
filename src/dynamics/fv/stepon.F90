module stepon

!----------------------------------------------------------------------
! stepon provides a standard interface that is called from the higher level
! CAM component run methods while leaving non-standardized dycore interface
! methods to be called from this layer.  Plan is to move the
! stardardization into the dynamics interface layer.
!----------------------------------------------------------------------

use shr_kind_mod,       only: r8 => shr_kind_r8
use spmd_utils,         only: mpicom, iam, masterproc

use ppgrid,             only: begchunk, endchunk
use constituents,       only: pcnst
use physconst,          only: zvir, cappa, physconst_calc_kappav, rair, cpair

use camsrfexch,         only: cam_out_t    
use cam_control_mod,    only: initial_run, moist_physics
use inic_analytic,      only: analytic_ic_active

use dyn_comp,           only: dyn_import_t, dyn_export_t, initial_mr, dyn_run
use dynamics_vars,      only: t_fvdycore_state, t_fvdycore_grid
use dyn_internal_state, only: get_dyn_state, get_dyn_state_grid
use advect_tend,        only: compute_adv_tends_xyz

use dp_coupling,        only: d_p_coupling, p_d_coupling

use physics_types,      only: physics_state, physics_tend
use physics_buffer,     only: physics_buffer_desc

use fv_prints,          only: fv_out

use time_manager,       only: get_step_size, get_curr_date
use cam_logfile,        only: iulog
use cam_abortutils,     only: endrun
use perf_mod,           only: t_startf, t_stopf, t_barrierf

implicit none
private
save

public :: &
   stepon_init,  &! Initialization
   stepon_d_p,   &! dynamics to physics coupling
   stepon_p_d,   &! physics to dynamics coupling
   stepon_run,   &! run dycore
   stepon_final   ! Finalization

integer  :: pdt        ! Physics time step
real(r8) :: te0        ! Total energy before dynamics

! for fv_out
logical, parameter :: fv_monitor=.true.  ! Monitor Mean/Max/Min fields
                                         ! This is CPU-time comsuming;
                                         ! set it to false for production runs
real (r8) :: ptop

!=========================================================================================
contains
!=========================================================================================

subroutine stepon_init(dyn_in, dyn_out)

   type (dyn_import_t)   :: dyn_in             ! Dynamics import container
   type (dyn_export_t)   :: dyn_out            ! Dynamics export container

   ! local variables:
   type (t_fvdycore_grid), pointer :: grid

   integer :: im, km
   integer :: ifirstxy, ilastxy, jfirstxy, jlastxy
   integer :: i,k,j,m             ! longitude, level, latitude and tracer indices
   logical :: nlres = .false.  ! true => restart or branch run

   integer :: ks
   real (r8), pointer :: ak(:)
   real (r8), pointer :: bk(:)

   real(r8), allocatable :: delpdryxy(:,:,:)
   real(r8), allocatable :: cap3vi(:,:,:), cappa3v(:,:,:)
   !----------------------------------------------------------------------------

   if (.not. initial_run) nlres=.true.

   grid => get_dyn_state_grid()
   im      =  grid%im
   km      =  grid%km


   ifirstxy  =  grid%ifirstxy
   ilastxy  =  grid%ilastxy
   jfirstxy  =  grid%jfirstxy
   jlastxy  =  grid%jlastxy

   ks     =  grid%ks
   ptop   =  grid%ptop
   ak     => grid%ak
   bk     => grid%bk

   do j = jfirstxy, jlastxy
      do i=ifirstxy, ilastxy
         dyn_in%pe(i,1,j) = ptop
      enddo
   enddo

   if ( nlres) then ! restart or branch run
      !
      ! read_restart_dynamics delivers phis, ps, u3s, v3s, delp, pt
      ! in XY decomposition

      !
      ! Do not recalculate delta pressure (delp) if this is a restart run.
      ! Re. SJ Lin: The variable "delp" (pressure thikness for a Lagrangian
      ! layer) must be in the restart file. This is because delp will be
      ! modified "after" the physics update (to account for changes in water
      ! vapor), and it can not be reproduced by surface pressure and the
      ! ETA coordinate's a's and b's.

!$omp parallel do private(i,j,k)
      do j = jfirstxy, jlastxy
        do k=1, km
          do i=ifirstxy, ilastxy
            dyn_in%pe(i,k+1,j) = dyn_in%pe(i,k,j) + dyn_in%delp(i,j,k)
          enddo
        enddo
      enddo
   else
 
      ! Initial run --> generate pe and delp from the surface pressure
 
!$omp parallel do private(i,j,k)
         do j = jfirstxy, jlastxy
            do k=1,km+1
               do i=ifirstxy, ilastxy
                  dyn_in%pe(i,k,j) = ak(k) + bk(k) * dyn_in%ps(i,j)
               enddo
            enddo
         enddo

!$omp parallel do private(i,j,k)
         do k = 1, km
            do j = jfirstxy, jlastxy
               do i= ifirstxy, ilastxy
                  dyn_in%delp(i,j,k) = dyn_in%pe(i,k+1,j) - dyn_in%pe(i,k,j)
               enddo
            enddo
         enddo
   endif

   !----------------------------------------------------------
   ! Check total dry air mass; set to 982.22 mb if initial run
   ! Print out diagnostic message if restart run
   !----------------------------------------------------------

   if ( moist_physics .and. .not. analytic_ic_active()) then
      call dryairm( grid, .true., dyn_in%ps, dyn_in%tracer,  &
                    dyn_in%delp, dyn_in%pe, nlres )
   endif

   if (grid%iam < grid%npes_xy) then

      allocate( cappa3v(ifirstxy:ilastxy,jfirstxy:jlastxy,km) )
      allocate( cap3vi(ifirstxy:ilastxy,jfirstxy:jlastxy,km+1) )
      if (grid%high_alt) then
         call physconst_calc_kappav( ifirstxy,ilastxy,jfirstxy,jlastxy,1,km, grid%ntotq, dyn_in%tracer, cappa3v )

!$omp parallel do private(i,j,k)
         do k=2,km
            do j=jfirstxy,jlastxy
               do i=ifirstxy,ilastxy
                  cap3vi(i,j,k) = 0.5_r8*(cappa3v(i,j,k-1)+cappa3v(i,j,k))
               enddo
            enddo
         enddo
         cap3vi(:,:,1) = 1.5_r8 * cappa3v(:,:,1) - 0.5_r8 * cappa3v(:,:,2)
         cap3vi(:,:,km+1) = 1.5_r8 * cappa3v(:,:,km) - 0.5_r8 * cappa3v(:,:,km-1)
      else
         cappa3v = rair/cpair
         cap3vi = rair/cpair
      endif

      ! Initialize pk, edge pressure to the cappa power.  Do this with constituent dependent cappa

!$omp parallel do private(i,j,k)
      do k = 1, km+1
         do j = jfirstxy, jlastxy
            do i = ifirstxy, ilastxy
               dyn_in%pk(i,j,k) = dyn_in%pe(i,k,j)**cap3vi(i,j,k)
            enddo
         enddo
      enddo

      ! Generate pkz, the conversion factor betw pt and t3

      call pkez(1,      im,   km,       jfirstxy,  jlastxy,              &
                1,      km,   ifirstxy, ilastxy,    dyn_in%pe,    &
                dyn_in%pk, cappa3v,  ks, dyn_out%peln, dyn_out%pkz,  .false., grid%high_alt )

      deallocate( cappa3v, cap3vi )

   endif

   if (initial_run) then

      ! pt is output by the dycore as virtual temperature.  This is what d_p_coupling
      ! expects for input

!$omp parallel do private(i,j,k)
      do k = 1, km
         do j = jfirstxy, jlastxy
            do i = ifirstxy, ilastxy
               dyn_in%pt(i,j,k) = dyn_in%t3(i,j,k)*(1._r8 + zvir*dyn_in%tracer(i,j,k,1))
            enddo
         enddo
      enddo

      !----------------------------------------------------------------
      ! Convert mixing ratios initialized as dry to moist for dynamics
      !----------------------------------------------------------------

      ! on initial time step, dry mixing ratio advected constituents have been
      !    initialized to dry mixing ratios. dynpkg expects moist m.r. so convert here.

      ! first calculate delpdry. The set_pdel_state subroutine
      !   is called after the dynamics in d_p_coupling to set more variables.
      !   This is not in tracers.F90 because it is only used by LR dynamics.
      allocate (delpdryxy(ifirstxy:ilastxy,jfirstxy:jlastxy,1:km))
      do k = 1, km
         do j = jfirstxy, jlastxy
            do i = ifirstxy, ilastxy
               delpdryxy(i,j,k) = dyn_in%delp(i,j,k)*          &
                    (1._r8 - dyn_in%tracer(i,j,k,1))
            enddo
         enddo
      enddo
      do m = 1,pcnst
         if (initial_mr(m) == 'dry') then
            do k=1, km
               do j = jfirstxy, jlastxy
                  do i = ifirstxy, ilastxy
                     dyn_in%tracer(i,j,k,m) =               &
                             dyn_in%tracer(i,j,k,m)*        &
                             delpdryxy(i,j,k)/dyn_in%delp(i,j,k)
                  end do
               end do
            end do
         end if
      end do
      deallocate (delpdryxy)
      
   end if

end subroutine stepon_init

!=========================================================================================

subroutine stepon_d_p(phys_state, phys_tend, pbuf2d, &
                      dyn_in, dyn_out, dtime_phys)

   ! arguments
   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend),  intent(inout) :: phys_tend(begchunk:endchunk)
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   type(dyn_import_t)                 :: dyn_in       ! Dynamics import container
   type(dyn_export_t)                 :: dyn_out      ! Dynamics export container
   real(r8),            intent(out)   :: dtime_phys   ! physics package time-step

   ! local variables
   type(T_FVDYCORE_GRID),  pointer :: grid
   !----------------------------------------------------------------------------

   pdt = get_step_size()    ! physics step size
   dtime_phys = real(pdt, r8)

   grid => get_dyn_state_grid()

   ! Dump state variables to IC file
   call t_barrierf('sync_diag_dynvar_ic', mpicom)
   call t_startf ('diag_dynvar_ic')
   call diag_dynvar_ic (grid, dyn_out%phis, dyn_out%ps,             &
                        dyn_out%t3, dyn_out%u3s, dyn_out%v3s, dyn_out%tracer  )
   call t_stopf  ('diag_dynvar_ic')

   ! Couple from dynamics to physics
   call t_barrierf('sync_d_p_coupling', mpicom)
   call t_startf('d_p_coupling')
   call d_p_coupling(grid, phys_state, phys_tend,  pbuf2d, dyn_out)
   call t_stopf('d_p_coupling')

end subroutine stepon_d_p

!=========================================================================================

subroutine stepon_p_d(phys_state, phys_tend, pbuf2d, dyn_in, dyn_out)

   ! INPUT/OUTPUT PARAMETERS:
   type(physics_state), intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend),  intent(inout) :: phys_tend(begchunk:endchunk)
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   type(dyn_import_t),  intent(inout) :: dyn_in  ! Dynamics import container
   type(dyn_export_t),  intent(inout) :: dyn_out ! Dynamics export container

   ! local variables
   real(r8) :: dtime_phys ! Physics time step
   type(T_FVDYCORE_GRID),  pointer :: grid

   integer :: rc 
   !-----------------------------------------------------------------------

   dtime_phys = real(pdt, r8)
   grid => get_dyn_state_grid()

   call t_barrierf('sync_p_d_coupling', mpicom)
   call t_startf ('p_d_coupling')
   call p_d_coupling(grid, phys_state, phys_tend, pbuf2d, &
                     dyn_in, dtime_phys, zvir, cappa, ptop)
   call t_stopf  ('p_d_coupling')

end subroutine stepon_p_d

!=========================================================================================

subroutine stepon_run(phys_state, dyn_in, dyn_out, cam_out)

   ! arguments
   type(physics_state), intent(in)    :: phys_state(begchunk:endchunk)
   type (dyn_import_t), intent(inout) :: dyn_in  ! Dynamics import container
   type (dyn_export_t), intent(inout) :: dyn_out ! Dynamics export container
   type(cam_out_t),     intent(inout) :: cam_out(begchunk:endchunk)

   ! LOCAL VARIABLES:
   integer :: yr, mon, day      ! year, month, day components of a date
   integer :: ncsec             ! time of day relative to current date [seconds]
   integer :: ncdate            ! current date in integer format [yyyymmdd]

   type(t_fvdycore_state), pointer :: dyn_state
   type(t_fvdycore_grid),  pointer :: grid

   integer :: freq_diag
   integer :: rc
   !----------------------------------------------------------------------------

   ! Monitor max/min/mean of selected fields
   ! Note that fv_out uses both dynamics and physics instantiations.

   call get_curr_date(yr, mon, day, ncsec)
   ncdate = yr*10000 + mon*100 + day

   dyn_state => get_dyn_state()
   grid => dyn_state%grid

   freq_diag = dyn_state%check_dt

   if (fv_monitor .and. mod(ncsec, freq_diag) == 0) then

      call t_barrierf('sync_fv_out', mpicom)
      call t_startf('fv_out')
      call fv_out(grid, dyn_out%pk, dyn_out%pt,         &
                  ptop, dyn_out%ps, dyn_out%tracer,     &
                  dyn_out%delp, dyn_out%pe, cam_out,    &
                  phys_state, ncdate, ncsec, moist_physics)
      call t_stopf('fv_out')
   endif

   call t_startf('comp_adv_tends1')
   call compute_adv_tends_xyz(grid, dyn_in%tracer )
   call t_stopf('comp_adv_tends1')

   call t_barrierf('sync_dyn_run', mpicom)
   call t_startf('dyn_run')
   call dyn_run(ptop,      pdt,     te0,         &
                dyn_state, dyn_in,  dyn_out,  rc)
   if ( rc /= 0 ) then
     write(iulog,*) "STEPON_RUN: dyn_run returned bad error code", rc
     write(iulog,*) "Quitting."
     call endrun
   endif 
   call t_stopf('dyn_run')

   call t_startf('comp_adv_tends2')
   call compute_adv_tends_xyz(grid, dyn_out%tracer )
   call t_stopf('comp_adv_tends2')

end subroutine stepon_run

!=========================================================================================

subroutine stepon_final(dyn_in, dyn_out)

   ! arguments
   type (dyn_import_t), intent(out) :: dyn_in  ! Dynamics import container
   type (dyn_export_t), intent(out) :: dyn_out ! Dynamics export container
   !----------------------------------------------------------------------------

end subroutine stepon_final

!=========================================================================================

end module stepon
