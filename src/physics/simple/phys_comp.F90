module phys_comp

!-----------------------------------------------------------------------
!
! Interfaces to the physics package for CAM simple models.
!
!-----------------------------------------------------------------------

use shr_kind_mod,    only: r8 => shr_kind_r8
use shr_sys_mod,     only: shr_sys_flush
use spmd_utils,      only: masterproc, mpicom

! Note: ideal_phys is true for Held-Suarez (1994) physics
use cam_control_mod, only: moist_physics, adiabatic, ideal_phys, kessler_phys, tj2016_phys
use physconst,       only: physconst_init, mwh2o, cpwv
use constituents,    only: cnst_add, cnst_chk_dim

use ppgrid,          only: begchunk, endchunk, pcols, pver, pverp
use phys_grid,       only: get_ncols_p
use physics_types,   only: physics_state, physics_tend, physics_type_alloc, &
                           physics_state_set_grid

use physics_buffer,  only: physics_buffer_desc, dtype_r8, pbuf_add_field,   &
                           pbuf_initialize, pbuf_get_index, pbuf_get_chunk, &
                           pbuf_allocate, pbuf_deallocate,                  &
                           pbuf_init_time, pbuf_update_tim_idx

use camsrfexch,      only: cam_out_t, cam_in_t, cam_export

use physpkg,         only: physpkg_init, physpkg_run
use held_suarez_cam, only: held_suarez_init
use kessler_cam,     only: kessler_register, kessler_init
use tj2016_cam,      only: thatcher_jablonowski_register, thatcher_jablonowski_init
use chemistry,       only: chem_register, chem_init, chem_is_active
use tracers,         only: tracers_register, tracers_init

use wv_saturation,   only: wv_sat_init

use check_energy,    only: check_energy_register, check_energy_init, check_energy_gmean, &
                           check_energy_gmean_diags
use phys_gmean,      only: gmean_mass
use qneg_module,     only: qneg_init

use cam_grid_support,only: cam_grid_check, cam_grid_id, cam_grid_get_dim_names

use cam_diagnostics, only: diag_register, diag_allocate, diag_init, diag_deallocate
use phys_debug_util, only: phys_debug_init
use perf_mod,        only: t_barrierf, t_startf, t_stopf, t_adj_detailf
use cam_logfile,     only: iulog
use cam_abortutils,  only: endrun

implicit none
private
save

public :: &
   phys_register, &! Register constituents and physics buffers
   phys_init,     &! Initialize
   phys_run,      &! run physics package
   phys_final      ! Finalize

!=========================================================================================
contains
!=========================================================================================

subroutine phys_register

   ! Register constituents and physics buffer fields.

   ! Local variables
   integer  :: mm       ! constituent index
   !----------------------------------------------------------------------------

   ! Initialize dyn_time_lvls
   call pbuf_init_time()

   ! Register water vapor.
   ! ***** N.B. ***** This must be the first call to cnst_add so that
   !                  water vapor is constituent 1.
   if (moist_physics) then
      call cnst_add('Q', mwh2o, cpwv, 1.E-12_r8, mm, &
         longname='Specific humidity', readiv=.true.)
   else
      call cnst_add('Q', mwh2o, cpwv, 0.0_r8, mm, &
         longname='Specific humidity', readiv=.false.)
   end if

   if (kessler_phys) then
      call kessler_register()
   else if (tj2016_phys) then
      call thatcher_jablonowski_register()
   end if

   ! check energy package
   call check_energy_register()

   ! register chemical constituents including aerosols ...
   call chem_register()

   ! Register test tracers
   call tracers_register()

   ! All tracers registered, check that the dimensions are correct
   call cnst_chk_dim()

   ! ***NOTE*** No registering constituents after the call to cnst_chk_dim.

   ! Register diagnostics PBUF fields
   call diag_register()

end subroutine phys_register

!======================================================================================

subroutine phys_init(phys_state, phys_tend, pbuf2d, cam_out)

   ! Initialize physics parameterization and package infrastructure

   ! Input/output arguments
   type(physics_state), pointer       :: phys_state(:)
   type(physics_tend ), pointer       :: phys_tend(:)
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   type(cam_out_t),intent(inout)      :: cam_out(begchunk:endchunk)

   ! local variables
   integer :: lchnk
   !----------------------------------------------------------------------------

   call physics_type_alloc(phys_state, phys_tend, begchunk, endchunk, pcols)

   do lchnk = begchunk, endchunk
      call physics_state_set_grid(lchnk, phys_state(lchnk))
   end do

   call pbuf_initialize(pbuf2d)

   ! Initialize any variables in physconst which are not temporally and/or
   ! spatially constant
   call physconst_init()

   ! Initialize debugging a physics column
   call phys_debug_init()

   call diag_init(pbuf2d)

   call check_energy_init()

   if (moist_physics) then
      call wv_sat_init()
   end if

   call tracers_init()

   call physpkg_init()

   call phys_inidat(cam_out, pbuf2d)

   if (ideal_phys) then
      call held_suarez_init()
   else if (kessler_phys) then
      call kessler_init()
   else if (tj2016_phys) then
      call thatcher_jablonowski_init()
   end if

   if (chem_is_active()) then
      call chem_init(phys_state, pbuf2d)
   end if

   call qneg_init()

end subroutine phys_init

!=========================================================================================

subroutine phys_run(ztodt, phys_state, phys_tend, pbuf2d, cam_in, cam_out)

   ! Run physics package

   ! arguments
   real(r8),                  intent(in)    :: ztodt
   type(physics_state),       intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend ),       intent(inout) :: phys_tend(begchunk:endchunk)
   type(physics_buffer_desc), pointer       :: pbuf2d(:,:)
   type(cam_in_t),            intent(inout) :: cam_in(begchunk:endchunk)
   type(cam_out_t),           intent(inout) :: cam_out(begchunk:endchunk)

   ! Local variables
   integer                            :: c                    ! chunk index
   type(physics_buffer_desc), pointer :: phys_buffer_chunk(:)
   !----------------------------------------------------------------------------

   call t_startf('phys_run_st1')

   call pbuf_allocate(pbuf2d, 'physpkg')
   call diag_allocate()

   ! Compute total energy of input to physics and previous output from physics
   call t_startf('chk_en_gmean')
   call check_energy_gmean(phys_state, pbuf2d, ztodt)
   call t_stopf('chk_en_gmean')

   ! Advance time interpolated data, save off state, etc.
   call phys_timestep_init(phys_state, cam_in, cam_out, pbuf2d)

   call t_stopf('phys_run_st1')

#ifdef TRACER_CHECK
   call gmean_mass('before physpkg_run', phys_state)
#endif

   call t_barrierf('sync_physpkg_run', mpicom)
   call t_startf('physpkg_run')
   call t_adj_detailf(+1)

!$OMP PARALLEL DO PRIVATE (c, phys_buffer_chunk)
   do c = begchunk, endchunk

      phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)

      call physpkg_run(ztodt, cam_in(c), cam_out(c), phys_state(c), phys_tend(c), &
                       phys_buffer_chunk)
   end do

   call t_adj_detailf(-1)
   call t_stopf('physpkg_run')

   call t_startf('phys_run_st2')

#ifdef TRACER_CHECK
   call gmean_mass('after physpkg_run', phys_state)
#endif

   call check_energy_gmean_diags()

   call pbuf_deallocate(pbuf2d, 'physpkg')
   call pbuf_update_tim_idx()

   call diag_deallocate()
   call t_stopf('phys_run_st2')

end subroutine phys_run

!=========================================================================================

subroutine phys_final( phys_state, phys_tend, pbuf2d)

   ! Finalize physics package

   ! Input/output arguments
   type(physics_state), pointer :: phys_state(:)
   type(physics_tend ), pointer :: phys_tend(:)
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)
   !----------------------------------------------------------------------------

   if(associated(pbuf2d)) then
      call pbuf_deallocate(pbuf2d,'global')
      deallocate(pbuf2d)
   end if
   deallocate(phys_state)
   deallocate(phys_tend)

end subroutine phys_final

!=========================================================================================
! Private routines
!=========================================================================================

subroutine phys_inidat(cam_out, pbuf2d)

   ! Dynamics variables are handled in dyn_init - here we read variables
   ! needed for physics but not dynamics

   ! Dummy arguments
   type(cam_out_t), intent(inout)     :: cam_out(begchunk:endchunk)
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   ! Local variables
   character(len=8)                   :: dim1name, dim2name
   integer                            :: grid_id ! grid ID for data mapping
   character(len=*), parameter        :: subname='phys_inidat'
   !----------------------------------------------------------------------------

   grid_id = cam_grid_id('physgrid')
   if (.not. cam_grid_check(grid_id)) then
      call endrun(subname//': Internal error, no "physgrid" grid')
   end if
   call cam_grid_get_dim_names(grid_id, dim1name, dim2name)

end subroutine phys_inidat

!=========================================================================================

subroutine phys_timestep_init(phys_state, cam_in, cam_out, pbuf2d)

   ! Call initializations needed at beginning of each timestep.
   ! Generally this is used to update time interpolated fields from
   ! boundary datasets.

   type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
   type(cam_in_t),      intent(inout), dimension(begchunk:endchunk) :: cam_in
   type(cam_out_t),     intent(inout), dimension(begchunk:endchunk) :: cam_out

   type(physics_buffer_desc), pointer                               :: pbuf2d(:,:)
   !----------------------------------------------------------------------------

end subroutine phys_timestep_init

!=========================================================================================

end module phys_comp
