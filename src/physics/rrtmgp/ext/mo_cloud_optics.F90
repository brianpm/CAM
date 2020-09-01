! Module: mo_cloud_optics

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

!
! This is the interface for routines that draw McICA samples of cloud physical properties
!    however provided and convert these into optical properties by band
!

module mo_cloud_optics
  use mo_rrtmgp_kind,          only: wp
  use mo_gas_optics_specification, & 
                              only: ty_gas_optics_specification
  use mo_rng,                 only: ty_rng
  use mo_optical_props,       only: ty_optical_props
  implicit none

  ! --- abstract base class: cloud description
  type, abstract, public :: ty_cloud_desc
  contains
    procedure(abstract_cloud_optics), deferred :: cloud_optics
    generic, public :: sample_and_optics => cloud_optics
  end type

! --- "Deferred" procedures. These need to be implmented in any derived class
! ------------------------------------------------------------------------------------------
  abstract interface
    function abstract_cloud_optics(this, rngs, spec_cfg, cloud_props) result(error_msg) 
      import ty_cloud_desc
      import ty_rng, ty_optical_props, ty_gas_optics_specification
      class(ty_cloud_desc), intent(in   ) :: this
      class(ty_rng), dimension(:), &
                            intent(inout) :: rngs
      type (ty_gas_optics_specification),   intent(in   ) :: spec_cfg
      class(ty_optical_props),     &
                            intent(  out) :: cloud_props
      character(len=128) :: error_msg                       
    end function abstract_cloud_optics
  end interface
! ------------------------------------------------------------------------------------------
end module mo_cloud_optics
