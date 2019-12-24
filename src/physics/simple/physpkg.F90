module physpkg

!-----------------------------------------------------------------------
!
! Physics package driver (run method) for CAM simple models.
!
! The current simplified physical parameterization packages are:
! . no physics forcing, adiabatic
! . dry Held-Suarez forcing, see Held and Suarez (BAMS, 1994)
! . Kessler warm-rain precipitation, see Klemp et al. (JAMES, 2015) and DCMIP-2016
! . "moist Held-Suarez" physics package by Thatcher and Jablonowski (GMD, 2016)
!
!-----------------------------------------------------------------------

use shr_kind_mod,    only: r8 => shr_kind_r8
use shr_sys_mod,     only: shr_sys_flush
use spmd_utils,      only: masterproc, mpicom

! Note: ideal_phys is true for Held-Suarez (1994) physics
use cam_control_mod, only: moist_physics, adiabatic, ideal_phys, kessler_phys, tj2016_phys
use dycore,          only: dycore_is
use phys_control,    only: phys_getopts

use camsrfexch,      only: cam_out_t, cam_in_t, cam_export

use physics_types,   only: physics_state, physics_tend, physics_ptend, &
                           physics_tend_init, physics_update,          &
                           physics_state_check

use ppgrid,          only: pcols

use physics_buffer,  only: physics_buffer_desc

use held_suarez_cam, only: held_suarez_tend
use kessler_cam,     only: kessler_tend
use tj2016_cam,      only: thatcher_jablonowski_sfc_pbl_hs_tend, thatcher_jablonowski_precip_tend
use chemistry,       only: chem_is_active, chem_timestep_tend

use check_energy,    only: check_energy_chng, check_energy_timestep_init, check_energy_fix, &
                           calc_te_and_aam_budgets, check_energy_after_phys,                &
                           check_energy_write,                                              &
                           check_tracers_data, check_tracers_init, check_tracers_chng
use cam_diagnostics, only: diag_surf, diag_before_phys, diag_after_phys

use perf_mod,        only: t_startf, t_stopf
use cam_logfile,     only: iulog
use cam_abortutils,  only: endrun

implicit none
private
save

public :: &
   physpkg_init, &
   physpkg_run

logical :: state_debug_checks

!=========================================================================================
contains
!=========================================================================================

subroutine physpkg_init()

   ! Get physics options
   call phys_getopts(state_debug_checks_out = state_debug_checks)

end subroutine physpkg_init

!=========================================================================================

subroutine physpkg_run(ztodt, cam_in, cam_out, state, tend, pbuf)

   ! Arguments
   real(r8),                  intent(in)    :: ztodt   ! physics package timestep
   type(cam_in_t),            intent(inout) :: cam_in
   type(cam_out_t),           intent(inout) :: cam_out
   type(physics_state),       intent(inout) :: state
   type(physics_tend ),       intent(inout) :: tend
   type(physics_buffer_desc), pointer       :: pbuf(:)

   ! Local variables
   type(physics_ptend) :: ptend           ! parameterization tendencies
   real(r8)            :: flx_heat(pcols) ! flux corresponding to global energy fixer
   real(r8), parameter :: zero(pcols) = 0._r8
   !----------------------------------------------------------------------------

   ! Verify state coming from the dynamics
   if (state_debug_checks) then
      call physics_state_check(state, name="before physics")
   end if

   ! output fields imported from the surface components
   call diag_surf(cam_in, state, pbuf)

   ! initialize the total tendency for the physics package
   call physics_tend_init(tend)

   call calc_te_and_aam_budgets(state, pbuf, 'pBF')

   ! Global mean total energy fixer
   call t_startf('energy_fixer')
   if (adiabatic .and. (.not. dycore_is('EUL'))) then
      call check_energy_fix(state, ptend, flx_heat)
      call physics_update(state, ptend, ztodt, tend)
      call check_energy_chng(state, tend, "chkengyfix", ztodt, zero, zero, zero, flx_heat)
   end if
   call t_stopf('energy_fixer')

   call calc_te_and_aam_budgets(state, pbuf, 'pBP')

   ! global energy diagnostics
   call check_energy_write(state, pbuf, ztodt)

   ! write out and save "before physics" state
   call diag_before_phys(state, pbuf)

   ! Compute package tendencies

   if (ideal_phys) then

      call held_suarez_tend(state, ptend, ztodt)
      call physics_update(state, ptend, ztodt, tend)

   else if (kessler_phys) then

      call kessler_tend(state, ptend, ztodt)
      call physics_update(state, ptend, ztodt, tend)

   else if (tj2016_phys) then

      ! Update surface, PBL and modified Held-Suarez forcings
      call thatcher_jablonowski_sfc_pbl_hs_tend(state, ptend, ztodt, cam_in)
      call physics_update(state, ptend, ztodt, tend)

      ! Compute the large-scale precipitation
      call thatcher_jablonowski_precip_tend(state, ptend, ztodt)
      call physics_update(state, ptend, ztodt, tend)

   end if

   ! Can't turn on conservation error messages unless the appropriate heat
   ! surface flux is computed and supplied as an argument to
   ! check_energy_chng to account for how the simplified physics forcings are
   ! changing the total exnergy.
   call check_energy_chng(state, tend, "physpkg", ztodt, zero, zero, zero, zero)

   if (chem_is_active()) then

      call chem_timestep_tend(state, ptend, cam_in, cam_out, ztodt, pbuf)
      call physics_update(state, ptend, ztodt, tend)
      call check_energy_chng(state, tend, "chempkg", ztodt, zero, zero, zero, zero)
      
   end if

   call diag_after_phys(state, pbuf, tend, ztodt)

   call cam_export(state, cam_out, pbuf)

   call calc_te_and_aam_budgets(state, pbuf, 'pAP')

   call check_energy_after_phys(state, pbuf)

 end subroutine physpkg_run

!=========================================================================================

end module physpkg
