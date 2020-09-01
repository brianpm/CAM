! Module: mo_optical_props

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
! Description:  Sets up arrays needed for optical properties.

module mo_optical_props
  use mo_rrtmgp_kind,           only: wp
  use mo_optical_props_kernels, only: &
        increment_1scalar_by_1scalar, increment_1scalar_by_2stream, increment_1scalar_by_nstream, &
        increment_2stream_by_1scalar, increment_2stream_by_2stream, increment_2stream_by_nstream, &
        increment_nstream_by_1scalar, increment_nstream_by_2stream, increment_nstream_by_nstream, & 
        inc_1scalar_by_1scalar_bybnd, inc_1scalar_by_2stream_bybnd, inc_1scalar_by_nstream_bybnd, &
        inc_2stream_by_1scalar_bybnd, inc_2stream_by_2stream_bybnd, inc_2stream_by_nstream_bybnd, &
        inc_nstream_by_1scalar_bybnd, inc_nstream_by_2stream_bybnd, inc_nstream_by_nstream_bybnd
  implicit none

  ! --- 1 scalar ------------------------------------------------------------------------
  ! declaration of optical properties
  ! the 2 stream method
  type :: ty_optical_props
    real(wp), dimension(:,:,:), allocatable :: tau ! optical depth (ncol, nlay, ngpt)
  contains
    procedure, public :: init_1scalar
    procedure, public :: delta_scale_1scalar
    generic,   public :: delta_scale => delta_scale_1scalar
    procedure, private :: increment_gpt_by
    procedure, private :: increment_band_by
    generic,   public  :: increment_by => increment_gpt_by, increment_band_by
    procedure, private :: subset_1scalar_range
    procedure, public  :: get_subset => subset_1scalar_range
  end type

  ! --- 2 stream ------------------------------------------------------------------------
  ! declaration of optical properties
  ! the scalar method
  type, extends(ty_optical_props) :: ty_optical_props_2str
    real(wp), dimension(:,:,:), allocatable :: ssa ! single-scattering albedo (ncol, nlay, ngpt)
    real(wp), dimension(:,:,:), allocatable :: g ! asymmetry parameter (ncol, nlay, ngpt)
  contains
    procedure, public :: init_2stream
    procedure, public :: delta_scale_2stream
    generic,   public :: delta_scale => delta_scale_2stream
    procedure, private :: subset_2str_range
    procedure, public  :: get_subset => subset_2str_range
  end type

  ! --- n stream ------------------------------------------------------------------------
  ! declaration of optical properties
  ! the n stream method
  type, extends(ty_optical_props) :: ty_optical_props_nstr
    real(wp), dimension(:,:,:), allocatable :: ssa ! single-scattering albedo (ncol, nlay, ngpt)
    real(wp), dimension(:,:,:,:), allocatable :: p ! phase-function moments (nmom, ncol, nlay, ngpt)
  contains
    procedure, public :: init_nstream
    procedure, public :: delta_scale_nstream
    generic,   public :: delta_scale => delta_scale_nstream
    procedure, private :: subset_nstr_range
    procedure, public  :: get_subset => subset_nstr_range
  end type

contains

  ! --- 1 scalar ------------------------------------------------------------------------
  ! setup
  function init_1scalar(this, ncol, nlay, ngpt) result(err_message)
    class(ty_optical_props) :: this
    integer, intent(in)     :: ncol, nlay, ngpt
    character(len=128)      :: err_message

    err_message = "" 
    if(any([ncol, nlay, ngpt] <= 0)) then 
      err_message = "optical_props%init: must provide positive extents for ncol, nlay, ngpt"
    else
      !
      ! Assuming Fortran 2003 standard for allocation, meaning any current allocation is 
      !   deallocated 
      !
      allocate(this%tau(ncol,nlay,ngpt))
    end if 
  end function init_1scalar

  ! --- 2 stream ------------------------------------------------------------------------
  ! setup
  function init_2stream(this, ncol, nlay, ngpt) result(err_message)
    class(ty_optical_props_2str) :: this
    integer, intent(in)          :: ncol, nlay, ngpt
    character(len=128)           :: err_message

    err_message = this%ty_optical_props%init_1scalar(ncol, nlay, ngpt)
    if(err_message == "") allocate(this%ssa(ncol,nlay,ngpt), this%g(ncol,nlay,ngpt))
  end function init_2stream

  ! --- n stream ------------------------------------------------------------------------
  ! setup
  function init_nstream(this, nmom, ncol, nlay, ngpt) result(err_message)
    class(ty_optical_props_nstr) :: this
    integer, intent(in)          :: nmom ! number of moments
    integer, intent(in)          :: ncol, nlay, ngpt
    character(len=128)           :: err_message

    err_message = this%ty_optical_props%init_1scalar(ncol, nlay, ngpt)
    if(err_message == "") allocate(this%ssa(ncol,nlay,ngpt), this%p(nmom,ncol,nlay,ngpt))
  end function init_nstream
  ! ------------------------------------------------------------------------------------------
  ! --- delta scaling ------------------------------------------------------------------------
  ! ------------------------------------------------------------------------------------------
  function delta_scale_1scalar(this, f) result(err_message)
    class(ty_optical_props), intent(inout) :: this
    real, dimension(:,:,:),  intent(in   ) :: f
    character(128)                         :: err_message
    !
    ! Nothing to do
    !
    err_message = ""
  end function delta_scale_1scalar
  ! --------------------------------
  function delta_scale_2stream(this, for) result(err_message)
    class(ty_optical_props_2str), intent(inout) :: this
    real(wp), dimension(:,:,:), target, optional, &
                                  intent(in   ) :: for
    ! Forward scattering fraction; g**2 if not provided
    character(128)                              :: err_message

    integer :: ncol, nlay, ngpt
    integer :: icol, ilay, igpt
    real(wp), dimension(size(this%tau,dim=1), size(this%tau,dim=2), size(this%tau,dim=3)), target :: f_loc
    real(wp), dimension(:,:,:), pointer :: f
    real(wp) :: wf ! Temporary -- should it be vector?
    ! --------------------------------
    ncol = size(this%tau,dim=1)
    nlay = size(this%tau,dim=2)
    ngpt = size(this%tau,dim=3)
    err_message = "" 
    
    if(present(for)) then
      if(any([size(for, 1), size(for, 2), size(for, 3)] /= [ncol, nlay, ngpt])) then
        err_message = "delta_scale: dimension of 'for' don't match optical properties arrays" 
        return
      end if 
      if(any(for < 0._wp .or. for > 1._wp)) then 
        err_message = "delta_scale: values of 'for' out of bounds [0,1]" 
        return
      end if 
      f => for
    else
      f_loc(:,:,:) = this%g(:,:,:) * this%g(:,:,:)
      f => f_loc
    end if

    do igpt = 1, ngpt
      do ilay = 1, nlay
        do icol = 1, ncol
          wf = this%ssa(icol,ilay,igpt) * f(icol,ilay,igpt)
          this%tau(icol,ilay,igpt) = (1._wp - wf) * this%tau(icol,ilay,igpt)
          this%ssa(icol,ilay,igpt) = (this%ssa(icol,ilay,igpt) - wf) / (1.0_wp - wf)
          this%g  (icol,ilay,igpt) = (this%g  (icol,ilay,igpt) - f(icol,ilay,igpt)) / &
                                       (1._wp - f(icol,ilay,igpt))
        end do
      end do
    end do

  end function delta_scale_2stream
  ! --------------------------------

  function delta_scale_nstream(this, for) result(err_message)
    class(ty_optical_props_nstr), intent(inout) :: this
    real(wp), dimension(:,:,:), target, optional, &
                                 intent(in   ) :: for
    ! Forward scattering fraction; g**2 if not provided
    character(128)                             :: err_message

    err_message = 'delta_scale_nstream: Not yet implemented'
  end function delta_scale_nstream

  ! ------------------------------------------------------------------------------------------
  ! --- Return subsets of optical properties arrays along x (col) direction
  ! ------------------------------------------------------------------------------------------
  
  !
  ! Allocate class, then arrays; copy. Could probably be more efficient if 
  !   classes used pointers internally. 
  !
  
  ! ------------------------------------------------------------------------------------------
  ! This set takes start position and number as scalars
  ! ------------------------------------------------------------------------------------------
  function subset_1scalar_range(full, start, n, subset) result(err_message)
    class(ty_optical_props), intent(inout) :: full
    integer,                 intent(in   ) :: start, n
    class(ty_optical_props), intent(inout) :: subset 
    character(128)                         :: err_message

    integer :: ncol, nlay, ngpt

    ncol = size(full%tau,dim=1)
    nlay = size(full%tau,dim=2)
    ngpt = size(full%tau,dim=3)
    err_message = "" 
    if(start < 1 .or. start + n-1 > size(full%tau, 1)) & 
       err_message = "optical_props%subset: Asking for columns outside range" 
    if(err_message /= "") return 
    
    ! Seems like the deallocation statements should be needed under Fortran 2003    
    !   but Intel compiler doesn't run without them
    if(allocated(subset%tau)) deallocate(subset%tau)
    select type (subset)
      class is (ty_optical_props)  
        err_message = subset%init_1scalar(n, nlay, ngpt)
        if(err_message /= "") return 
      class is (ty_optical_props_2str)     
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%g  )) deallocate(subset%g  )
        err_message = subset%init_2stream(n, nlay, ngpt)
        if(err_message /= "") return 
        subset%ssa(1:n,:,:) = 0._wp
        subset%g  (1:n,:,:) = 0._wp
      class is (ty_optical_props_nstr) 
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%p  )) deallocate(subset%p  )
        err_message = subset%init_nstream(size(subset%p,1), n, nlay, ngpt)
        if(err_message /= "") return 
        subset%ssa(1:n,:,:) = 0._wp
        subset%p(:,1:n,:,:) = 0._wp
    end select 
    subset%tau(1:n,:,:) = full%tau(start:start+n-1,:,:)

  end function subset_1scalar_range
  ! ------------------------------------------------------------------------------------------
  function subset_2str_range(full, start, n, subset) result(err_message)
    class(ty_optical_props_2str), intent(inout) :: full
    integer,                      intent(in   ) :: start, n
    class(ty_optical_props),      intent(inout) :: subset 
    character(128)                              :: err_message

    integer :: ncol, nlay, ngpt

    ncol = size(full%tau,dim=1)
    nlay = size(full%tau,dim=2)
    ngpt = size(full%tau,dim=3)
    err_message = "" 
    if(start < 1 .or. start + n-1 > size(full%tau, 1)) & 
       err_message = "optical_props%subset: Asking for columns outside range" 
    if(err_message /= "") return 
    
    if(allocated(subset%tau)) deallocate(subset%tau)
    select type (subset)
      class is (ty_optical_props)  
        err_message = subset%init_1scalar(n, nlay, ngpt)
        if(err_message /= "") return 
        subset%tau(1:n,:,:) = full%tau(start:start+n-1,:,:) * & 
                     (1._wp - full%ssa(start:start+n-1,:,:))
      class is (ty_optical_props_2str)     
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%g  )) deallocate(subset%g  )
        err_message = subset%init_2stream(n, nlay, ngpt)
        if(err_message /= "") return 
        subset%tau(1:n,:,:) = full%tau(start:start+n-1,:,:) 
        subset%ssa(1:n,:,:) = full%ssa(start:start+n-1,:,:)
        subset%g  (1:n,:,:) = full%g  (start:start+n-1,:,:)
      class is (ty_optical_props_nstr) 
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%p  )) deallocate(subset%p  )
        err_message = subset%init_nstream(size(subset%p,1), n, nlay, ngpt)
        if(err_message /= "") return 
        subset%tau(1:n,:,:) = full%tau(start:start+n-1,:,:) 
        subset%ssa(1:n,:,:) = full%ssa(start:start+n-1,:,:)
        subset%p(1,1:n,:,:) = full%g  (start:start+n-1,:,:)
        subset%p(2:,:, :,:) = 0._wp
    end select 

  end function subset_2str_range
  ! ------------------------------------------------------------------------------------------
  function subset_nstr_range(full, start, n, subset) result(err_message)
    class(ty_optical_props_nstr), intent(inout) :: full
    integer,                      intent(in   ) :: start, n
    class(ty_optical_props),      intent(inout) :: subset 
    character(128)                              :: err_message

    integer :: ncol, nlay, ngpt, nmom

    ncol = size(full%tau,dim=1)
    nlay = size(full%tau,dim=2)
    ngpt = size(full%tau,dim=3)
    nmom = size(full%p  ,dim=1) 
    err_message = "" 
    if(start < 1 .or. start + n-1 > size(full%tau, 1)) & 
       err_message = "optical_props%subset: Asking for columns outside range" 
    if(err_message /= "") return 
    
    if(allocated(subset%tau)) deallocate(subset%tau)
    select type (subset)
      class is (ty_optical_props)  
        err_message = subset%init_1scalar(n, nlay, ngpt)
        if(err_message /= "") return 
        subset%tau(1:n,:,:) = full%tau(start:start+n-1,:,:) * & 
                     (1._wp - full%ssa(start:start+n-1,:,:))
      class is (ty_optical_props_2str)     
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%g  )) deallocate(subset%g  )
        err_message = subset%init_2stream(n, nlay, ngpt)
        if(err_message /= "") return 
        subset%tau(1:n,:,:) = full%tau(start:start+n-1,:,:) 
        subset%ssa(1:n,:,:) = full%ssa(start:start+n-1,:,:)
        subset%g  (1:n,:,:) = full%p(1,start:start+n-1,:,:)
      class is (ty_optical_props_nstr) 
        if(allocated(subset%ssa)) deallocate(subset%ssa)
        if(allocated(subset%p  )) deallocate(subset%p  )
        err_message = subset%init_nstream(nmom, n, nlay, ngpt)
        if(err_message /= "") return 
        subset%ssa(1:n,:,:) = full%ssa(start:start+n-1,:,:)
        subset%p(:,1:n,:,:) = full%p(:,start:start+n-1,:,:)
    end select 

  end function subset_nstr_range
  ! ------------------------------------------------------------------------------------------

  ! --- increment by ------------------------------------------------------------------------

  ! add 1scalar,2stream,nstream to 1scalar,2stream,nstream
  function increment_gpt_by(op1, op2)result(err_message)
    class(ty_optical_props), intent(inout) :: op1
    class(ty_optical_props), intent(in   ) :: op2
    character(128)                         :: err_message
    ! ----- 
    integer :: ncol, nlay, ngpt, nmom1
    ! ----- 
    err_message = "" 
    ncol = size(op1%tau,1) 
    nlay = size(op1%tau,2) 
    ngpt = size(op1%tau,3) 
    
    !
    ! Rudimentary error checking -- users are responsible for ensuring consistency of 
    !   array sizes within an object of ty_optical props 
    !
    if(any([size(op2%tau,1), size(op2%tau,2), size(op2%tau,3)] /= [ncol, nlay, ngpt])) then 
      err_message = "ty_optical_props%increment_by: optical properties objects are inconsistently sized" 
      return 
    end if 
    
    select type (op1)
      class is (ty_optical_props)
        select type (op2)
         class is (ty_optical_props)
           call increment_1scalar_by_1scalar(ncol, nlay, ngpt, &
                                             op1%tau,          &
                                             op2%tau)
         class is (ty_optical_props_2str)
           call increment_1scalar_by_2stream(ncol, nlay, ngpt, &
                                             op1%tau,          &
                                             op2%tau, op2%ssa)

         class is (ty_optical_props_nstr)
           ncol = size(op1%tau,1) 
           call increment_1scalar_by_nstream(ncol, nlay, ngpt, &
                                             op1%tau,          &
                                             op2%tau, op2%ssa)
        end select
      
    class is (ty_optical_props_2str)
      select type (op2)
        class is (ty_optical_props)
          call increment_2stream_by_1scalar(ncol, nlay, ngpt,   &
                                            op1%tau, op1%ssa,&
                                            op2%tau)
        class is (ty_optical_props_2str)
          call increment_2stream_by_2stream(ncol, nlay, ngpt,        &
                                            op1%tau, op1%ssa, op1%g, &
                                            op2%tau, op2%ssa, op2%g)
        class is (ty_optical_props_nstr)
          call increment_2stream_by_nstream(ncol, nlay, ngpt, size(op2%p, 1), &
                                            op1%tau, op1%ssa, op1%g, &
                                            op2%tau, op2%ssa, op2%p)
      end select
      
    class is (ty_optical_props_nstr)
      nmom1 = size(op1%p, 1) 
      select type (op2)
        class is (ty_optical_props)
          call increment_nstream_by_1scalar(ncol, nlay, ngpt, &
                                            op1%tau, op1%ssa, &
                                            op2%tau)
        class is (ty_optical_props_2str)
          call increment_nstream_by_2stream(ncol, nlay, ngpt, nmom1, &
                                            op1%tau, op1%ssa, op1%p, &
                                            op2%tau, op2%ssa, op2%g)
        class is (ty_optical_props_nstr)
          call increment_nstream_by_nstream(ncol, nlay, ngpt, nmom1, size(op2%p, 1), &
                                            op1%tau, op1%ssa, op1%p, &
                                            op2%tau, op2%ssa, op2%p)
      end select
    end select
  end function increment_gpt_by


  ! --- increment by (band) ------------------------------------------------------------------------

  ! add 1scalar,2stream,nstream to 1scalar,2stream,nstream with band dimension
  function increment_band_by(op1, op2, gpt_lims) result(err_message)
    class(ty_optical_props), intent(inout) :: op1
    class(ty_optical_props), intent(in   ) :: op2
    integer, dimension(:,:), intent(in   ) :: gpt_lims  ! (begin g-point, end g-point) = gpt_lims(2,band)
    character(128)                         :: err_message

    ! ----- 
    integer :: ncol, nlay, ngpt, nbnd, nmom1
    ! ----- 
    ncol = size(op1%tau,1) 
    nlay = size(op1%tau,2) 
    ngpt = size(op1%tau,3) 
    nbnd = size(gpt_lims, 2) 
    
    err_message = "" 
    !
    ! Rudimentary error checking -- users are responsible for ensuring consistency of 
    !   array sizes within an object of ty_optical props 
    !
    if(any([size(op2%tau,1), size(op2%tau,2)] /= [ncol, nlay])) & 
      err_message = "ty_optical_props%increment_by: optical properties objects are inconsistently sized" 
    if(size(op2%tau,3) /= nbnd)                          & 
      err_message = "ty_optical_props%increment_by: " // & 
                    "number of bands not consistent between optical properties objects, g-point limits" 
    if(minval(gpt_lims) < 1 .or. maxval(gpt_lims) > ngpt) & 
      err_message = "ty_optical_props%increment_by: " // & 
                    "band limits not consistent with number of gpoints" 
    if(err_message /= "") return 
        
    select type (op1)
      class is (ty_optical_props)
        select type (op2)
          class is (ty_optical_props)
            call inc_1scalar_by_1scalar_bybnd(ncol, nlay, ngpt, &
                                              op1%tau,          &
                                              op2%tau,          &
                                              nbnd, gpt_lims)
          class is (ty_optical_props_2str)
            call inc_1scalar_by_2stream_bybnd(ncol, nlay, ngpt, &
                                              op1%tau,          &
                                              op2%tau, op2%ssa, & 
                                              nbnd, gpt_lims)
          class is (ty_optical_props_nstr)
            call inc_1scalar_by_nstream_bybnd(ncol, nlay, ngpt, &
                                              op1%tau,          &
                                              op2%tau, op2%ssa, & 
                                              nbnd, gpt_lims)
        end select
      
      class is (ty_optical_props_2str)
        select type (op2)
          class is (ty_optical_props)
            call inc_2stream_by_1scalar_bybnd(ncol, nlay, ngpt, &
                                              op1%tau, op1%ssa, &
                                              op2%tau,          &
                                              nbnd, gpt_lims)
          class is (ty_optical_props_2str)
            call inc_2stream_by_2stream_bybnd(ncol, nlay, ngpt,        &
                                              op1%tau, op1%ssa, op1%g, &
                                              op2%tau, op2%ssa, op2%g, & 
                                              nbnd, gpt_lims)
          class is (ty_optical_props_nstr)
            call inc_2stream_by_nstream_bybnd(ncol, nlay, ngpt, size(op2%p, 1), &
                                              op1%tau, op1%ssa, op1%g, &
                                              op2%tau, op2%ssa, op2%p, & 
                                              nbnd, gpt_lims)
        end select
      
      class is (ty_optical_props_nstr)
        nmom1 = size(op1%p, 1) 
        select type (op2)
          class is (ty_optical_props)
            call inc_nstream_by_1scalar_bybnd(ncol, nlay, ngpt, &
                                              op1%tau, op1%ssa, &
                                              op2%tau,          &
                                              nbnd, gpt_lims)
          class is (ty_optical_props_2str)
            call inc_nstream_by_2stream_bybnd(ncol, nlay, ngpt, nmom1, &
                                              op1%tau, op1%ssa, op1%p, &
                                              op2%tau, op2%ssa, op2%g, & 
                                              nbnd, gpt_lims)
          class is (ty_optical_props_nstr)
            call inc_nstream_by_nstream_bybnd(ncol, nlay, ngpt, nmom1, size(op2%p, 1), &
                                              op1%tau, op1%ssa, op1%p, &
                                              op2%tau, op2%ssa, op2%p, & 
                                              nbnd, gpt_lims)
        end select
    end select
  end function increment_band_by
end module mo_optical_props
