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
module mo_aerosol_optics
  use mo_rrtmgp_kind,          only: wp
  use mo_gas_optics_specification,      only: ty_gas_optics_specification
  use mo_optical_props, only: ty_optical_props
  implicit none

  ! --- abstract base class: aerosol description
  type, abstract, public :: ty_aerosol_desc
  contains
    procedure(abstract_aer_optics), deferred :: aerosol_optics
    generic, public :: optics => aerosol_optics
  end type

  ! --- "Deferred" procedures. These need to be implmented in any derived class
  abstract interface
    subroutine abstract_aer_optics(this, k_dist, aer_props)
      import ty_aerosol_desc
      import ty_gas_optics_specification, ty_optical_props
      class(ty_aerosol_desc),  intent(in ) :: this
      type (ty_gas_optics_specification),      intent(in ) :: k_dist
      class(ty_optical_props), intent(out) :: aer_props
    end subroutine abstract_aer_optics
  end interface

end module mo_aerosol_optics
