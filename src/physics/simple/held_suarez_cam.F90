module held_suarez_cam

!-------------------------------------------------------------------------------
! 
!  Implement idealized Held-Suarez forcings:
!  Held, I. M., and M. J. Suarez, 1994: 'A proposal for the
!  intercomparison of the dynamical cores of atmospheric general
!  circulation models.'
!  Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
! 
!-------------------------------------------------------------------------------

use shr_kind_mod, only: r8 => shr_kind_r8
use ppgrid,       only: pcols, pver

implicit none
private
save

public :: held_suarez_init, held_suarez_tend

real(r8), parameter :: efoldf  =  1._r8  ! efolding time for wind dissipation
real(r8), parameter :: efolda  = 40._r8  ! efolding time for T dissipation
real(r8), parameter :: efolds  =  4._r8  ! efolding time for T dissipation
real(r8), parameter :: sigmab  =  0.7_r8 ! threshold sigma level
real(r8), parameter :: t00     = 200._r8 ! minimum reference temperature
real(r8), parameter :: kf      = 1._r8/(86400._r8*efoldf) ! 1./efolding_time for wind dissipation

real(r8), parameter :: onemsig = 1._r8 - sigmab ! 1. - sigma_reference

real(r8), parameter :: ka      = 1._r8/(86400._r8 * efolda) ! 1./efolding_time for temperature diss.
real(r8), parameter :: ks      = 1._r8/(86400._r8 * efolds)

!=========================================================================================
contains
!=========================================================================================

subroutine held_suarez_init()

   use cam_history,        only: addfld, add_default
   use physconst,          only: cappa, cpair
   use ref_pres,           only: pref_mid_norm, psurf_ref
   use held_suarez,        only: held_suarez_1994_init
   !----------------------------------------------------------------------------

   ! Set model constant values
   call held_suarez_1994_init(cappa, cpair, psurf_ref, pref_mid_norm)

   ! This field is added by radiation when full physics is used
   call addfld('QRS', (/ 'lev' /), 'A', 'K/s', &
        'Temperature tendency associated with the relaxation toward the equilibrium temperature profile')
   call add_default('QRS', 1, ' ')

end subroutine held_suarez_init

!=========================================================================================

subroutine held_suarez_tend(state, ptend, ztodt)

   use physconst,          only: cpairv
   use phys_grid,          only: get_rlat_all_p
   use physics_types,      only: physics_state, physics_ptend
   use physics_types,      only: physics_ptend_init
   use cam_history,        only: outfld
   use held_suarez,        only: held_suarez_1994

   ! arguments
   type(physics_state), intent(in)  :: state
   real(r8),            intent(in)  :: ztodt   ! physics package timestep
   type(physics_ptend), intent(out) :: ptend   ! Package tendencies

   ! Local workspace

   integer  :: lchnk            ! chunk identifier
   integer  :: ncol             ! number of atmospheric columns

   real(r8) :: clat(pcols)      ! latitudes(radians) for columns
   real(r8) :: qrs(pcols,pver)  ! heating rate
   integer  :: i, k             ! Longitude, level indices
   !----------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol

   call get_rlat_all_p(lchnk, ncol, clat)

   ! initialize individual parameterization tendencies
   call physics_ptend_init(ptend, state%psetcols, 'held_suarez', ls=.true., lu=.true., lv=.true.)

   call held_suarez_1994(pcols, ncol, clat, state%pmid, &
       state%u, state%v, state%t, ptend%u, ptend%v, ptend%s)

   qrs(:ncol,:) = ptend%s(:ncol,:)/cpairv(:ncol,:,lchnk)

   call outfld('QRS', qrs, pcols, lchnk)

end subroutine held_suarez_tend

end module held_suarez_cam
