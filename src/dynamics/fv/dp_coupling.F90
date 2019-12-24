module dp_coupling

!-------------------------------------------------------------------------------
! dynamics-physics coupling module
!-------------------------------------------------------------------------------

use shr_kind_mod,      only: r8 => shr_kind_r8
use shr_const_mod,     only: shr_const_rwv

use phys_control,      only: use_gw_front, use_gw_front_igw, waccmx_is
use constituents,      only: pcnst, cnst_get_ind, qmin
use physconst,         only: physconst_update, gravit, zvir, cpairv, rairv, &
                             physconst_calc_kappav

use pmgrid,            only: plev
use dynamics_vars,     only: T_FVDYCORE_GRID, t_fvdycore_state
use dyn_internal_state,only: get_dyn_state
use dyn_comp,          only: dyn_import_t, dyn_export_t, fv_print_dpcoup_warn, &
                             frontgf_idx, frontga_idx, uzm_idx, qini_idx
use d2a3dikj_mod,      only: d2a3dikj
use ctem,              only: ctem_diags, do_circulation_diags
use diag_module,       only: fv_diag_am_calc
use gravity_waves_sources, only: gws_src_fnct
use qbo,               only: qbo_use_forcing

use phys_grid,         only: get_ncols_p, get_lon_all_p, get_lat_all_p, &
                             block_to_chunk_send_pters, transpose_block_to_chunk,  &
                             block_to_chunk_recv_pters, chunk_to_block_send_pters, &
                             transpose_chunk_to_block, chunk_to_block_recv_pters

use ppgrid,            only: pcols, pver, begchunk, endchunk
use physics_types,     only: physics_state, physics_tend, set_dry_to_wet, &
                             set_state_pdry, set_wet_to_dry, physics_dme_adjust

use physics_buffer,    only: physics_buffer_desc, pbuf_get_chunk, &
                             pbuf_get_field

use geopotential,      only: geopotential_t
use check_energy,      only: check_energy_timestep_init
use qneg_module,       only: qneg3
use zonal_mean,        only: zonal_mean_3D

use cam_abortutils,    only: endrun
use cam_logfile,       only: iulog
use perf_mod,          only: t_startf, t_stopf, t_barrierf

#if ( defined OFFLINE_DYN )
use metdata,           only: get_met_fields
#endif
#if defined ( SPMD )
use spmd_dyn,          only: local_dp_map, block_buf_nrecs, chunk_buf_nrecs
#endif

implicit none
private
save

PUBLIC :: d_p_coupling, p_d_coupling

!=========================================================================================
CONTAINS
!=========================================================================================

subroutine d_p_coupling(grid, phys_state, phys_tend, pbuf2d, dyn_out)

   !------------------------------------------------------------------------------
   ! Coupler for converting dynamics output variables into physics input variables
   !------------------------------------------------------------------------------

   ! arguments
   type(t_fvdycore_grid),intent(in)    :: grid
   type(physics_state),  intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend ),  intent(inout) :: phys_tend(begchunk:endchunk)
   type(physics_buffer_desc), pointer  :: pbuf2d(:,:)
   type(dyn_export_t),   intent(in)    :: dyn_out  ! dynamics export

   ! LOCAL VARIABLES:

   type(t_fvdycore_state),    pointer :: dyn_state
   type(physics_buffer_desc), pointer :: pbuf_chnk(:)

   ! Variables from dynamics export container
   real(r8), pointer :: phisxy(:,:)              ! surface geopotential
   real(r8), pointer :: psxy (:,:)               ! surface pressure
   real(r8), pointer :: u3sxy(:,:,:)             ! u-wind on d-grid
   real(r8), pointer :: v3sxy(:,:,:)             ! v-wind on d-grid
   real(r8), pointer :: du3sxy(:,:,:)            ! u-wind increment on d-grid
   real(r8), pointer :: dv3sxy(:,:,:)            ! v-wind increment on d-grid
   real(r8), pointer :: dua3sxy(:,:,:)           ! u-wind adv. inc. on d-grid
   real(r8), pointer :: dva3sxy(:,:,:)           ! v-wind adv. inc. on d-grid
   real(r8), pointer :: duf3sxy(:,:,:)           ! u-wind fixer inc.on d-grid
   real(r8), pointer :: ptxy (:,:,:)             ! Virtual pot temp
   real(r8), pointer :: tracer(:,:,:,:)          ! constituents
   real(r8), pointer :: omgaxy(:,:,:)            ! vertical velocity
   real(r8), pointer :: pexy  (:,:,:)            ! edge pressure
   real(r8), pointer :: pelnxy(:,:,:)            ! log(pe)
   real(r8), pointer :: pkxy  (:,:,:)            ! pe**cappa
   real(r8), pointer :: pkzxy (:,:,:)            ! f-v mean of pk

   integer :: i,ib,j,k,m,lchnk      ! indices
   integer :: ncol                  ! number of columns in current chunk
   integer :: lats(pcols)           ! array of latitude indices
   integer :: lons(pcols)           ! array of longitude indices
   integer :: blksiz                ! number of columns in 2D block
   integer :: tsize                 ! amount of data per grid point passed to physics
   integer, allocatable, dimension(:,:) :: bpter
                                     ! offsets into block buffer for packing data
   integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for unpacking data

   real(r8) :: qmavl                ! available q at level pver-1
   real(r8) :: dqreq                ! q change at pver-1 required to remove q<qmin at pver
   real(r8) :: qbot                 ! bottom level q before change
   real(r8) :: qbotm1               ! bottom-1 level q before change
   real(r8) :: pic(pcols)           ! ps**cappa
   real(r8) :: fraction
   real(r8), allocatable :: u3(:, :, :)       ! u-wind on a-grid
   real(r8), allocatable :: v3(:, :, :)       ! v-wind on a-grid
   real(r8), allocatable :: du3 (:, :, :)     ! u-wind increment on a-grid
   real(r8), allocatable :: dv3 (:, :, :)     ! v-wind increment on a-grid
   real(r8), allocatable :: dua3(:, :, :)     ! u-wind adv. inc. on a-grid
   real(r8), allocatable :: dva3(:, :, :)     ! v-wind adv. inc. on a-grid
   real(r8), allocatable :: duf3(:, :, :)     ! u-wind fixer inc.on a-grid
   real(r8), allocatable :: dummy(:,:, :)     ! dummy work array ~ v3sxy*0

   real(r8), allocatable, dimension(:) :: bbuffer, cbuffer
                                     ! transpose buffers

   real(r8) :: zvirv(pcols,pver)    ! Local zvir array pointer

   real(r8) :: uzm(plev,grid%jfirstxy:grid%jlastxy) ! Zonal mean zonal wind

   integer  :: km, kmp1, iam
   integer  :: ifirstxy, ilastxy, jfirstxy, jlastxy
   integer  :: ic, jc
   integer  :: astat
   integer  :: boff

   ! frontogenesis function for gravity wave drag
   real(r8), allocatable :: frontgf(:,:,:)
   real(r8), pointer :: pbuf_frontgf(:,:)

   ! frontogenesis angle for gravity wave drag
   real(r8), allocatable :: frontga(:,:,:)
   real(r8), pointer :: pbuf_frontga(:,:)

   ! needed for qbo
   real(r8), pointer :: pbuf_uzm(:,:)

   ! Variables needed for WACCM-X
   integer  :: ixo, ixo2, ixh, ixh2, ixn  ! indices into state structure for O, O2, H, H2, and N
   real(r8) :: mmrSum_O_O2_H                ! Sum of mass mixing ratios for O, O2, and H
   real(r8), parameter :: mmrMin=1.e-20_r8  ! lower limit of o2, o, and h mixing ratios
   real(r8), parameter :: N2mmrMin=1.e-6_r8 ! lower limit of o2, o, and h mixing ratios

#if (! defined SPMD)
   integer :: block_buf_nrecs = 0
   integer :: chunk_buf_nrecs = 0
   logical :: local_dp_map=.true.
#endif
   !----------------------------------------------------------------------------

   dyn_state => get_dyn_state()

   if (use_gw_front .or. use_gw_front_igw) then

      allocate(frontgf(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy), stat=astat)
      if( astat /= 0 ) then
         write(iulog,*) 'd_p_coupling: failed to allocate frontgf; error = ',astat
         call endrun('d_p_coupling: failed to allocate frontgf')
      end if

      allocate(frontga(grid%ifirstxy:grid%ilastxy,plev,grid%jfirstxy:grid%jlastxy), stat=astat)
      if( astat /= 0 ) then
         write(iulog,*) 'd_p_coupling: failed to allocate frontga; error = ',astat
         call endrun('d_p_coupling: failed to allocate frontga')
      end if

   end if

   nullify(pbuf_chnk)
   nullify(pbuf_frontgf)
   nullify(pbuf_frontga)
   nullify(pbuf_uzm)

   fraction = 0.1_r8

   phisxy   => dyn_out%phis
   psxy     => dyn_out%ps
   u3sxy    => dyn_out%u3s
   v3sxy    => dyn_out%v3s
   ptxy     => dyn_out%pt
   tracer   => dyn_out%tracer

   omgaxy   => dyn_out%omga
   pexy     => dyn_out%pe
   pelnxy   => dyn_out%peln
   pkxy     => dyn_out%pk
   pkzxy    => dyn_out%pkz

   km       = grid%km
   kmp1     = km + 1

   ifirstxy = grid%ifirstxy
   ilastxy  = grid%ilastxy
   jfirstxy = grid%jfirstxy
   jlastxy  = grid%jlastxy

   iam      = grid%iam

   ! Transform dynamics staggered winds to physics grid (D=>A)
   call t_startf('d2a3dikj')
   allocate (u3(ifirstxy:ilastxy, km, jfirstxy:jlastxy))
   allocate (v3(ifirstxy:ilastxy, km, jfirstxy:jlastxy))
   if (iam < grid%npes_xy) then
      call d2a3dikj(grid, dyn_state%am_geom_crrct, u3sxy, v3sxy, u3, v3)
   end if
   call t_stopf('d2a3dikj')

   if (do_circulation_diags) then
      call t_startf('DP_CPLN_ctem')
      call ctem_diags(u3, v3, omgaxy, ptxy(:,jfirstxy:jlastxy,:), tracer(:,jfirstxy:jlastxy,:,1), &
                      psxy, pexy, grid)
      call t_stopf('DP_CPLN_ctem')
   end if

   if (dyn_state%am_diag) then
      du3sxy   => dyn_out%du3s
      dv3sxy   => dyn_out%dv3s
      dua3sxy  => dyn_out%dua3s
      dva3sxy  => dyn_out%dva3s
      duf3sxy  => dyn_out%duf3s
      allocate (du3 (ifirstxy:ilastxy, km, jfirstxy:jlastxy))
      allocate (dv3 (ifirstxy:ilastxy, km, jfirstxy:jlastxy))
      allocate (dua3(ifirstxy:ilastxy, km, jfirstxy:jlastxy))
      allocate (dva3(ifirstxy:ilastxy, km, jfirstxy:jlastxy))
      allocate (duf3(ifirstxy:ilastxy, km, jfirstxy:jlastxy))
      allocate (dummy(ifirstxy:ilastxy,jfirstxy:jlastxy, km))
      du3(:,:,:)   = 0._r8
      dv3(:,:,:)   = 0._r8
      dua3(:,:,:)  = 0._r8
      dva3(:,:,:)  = 0._r8
      duf3(:,:,:)  = 0._r8
      dummy(:,:,:) = 0._r8

      if (iam < grid%npes_xy) then
         ! (note dummy use of dva3 hence call order matters)
         call d2a3dikj(grid, dyn_state%am_geom_crrct,duf3sxy,   dummy, duf3 ,dva3)
         call d2a3dikj(grid, dyn_state%am_geom_crrct,dua3sxy, dva3sxy, dua3, dva3)
         call d2a3dikj(grid, dyn_state%am_geom_crrct, du3sxy,  dv3sxy, du3 , dv3 )
      end if

      call t_startf('DP_CPLN_fv_am')
      call fv_diag_am_calc(grid, psxy, pexy, du3, dv3, dua3, dva3, duf3)
      call t_stopf('DP_CPLN_fv_am')
   end if

   if ((iam .lt. grid%npes_xy) .and. (use_gw_front .or. use_gw_front_igw)) then
      call t_startf('DP_CPLN_gw_sources')
      call gws_src_fnct(grid, u3, v3, ptxy, tracer(:,jfirstxy:jlastxy,:,1), pexy, frontgf, frontga)
      call t_stopf('DP_CPLN_gw_sources')
   end if
   if (qbo_use_forcing) then
      call zonal_mean_3D(grid, plev, u3, uzm)
   end if

   !-----------------------------------------------------------------------
   ! Copy data from dynamics data structure to physics data structure
   !-----------------------------------------------------------------------
   if (local_dp_map) then

! This declaration is too long; this parallel section needs some stuff
! pulled out into routines.
!$omp parallel do private (lchnk, ncol, i, k, m, ic, jc, lons, lats, pic, pbuf_chnk, pbuf_uzm, pbuf_frontgf, pbuf_frontga)
      do lchnk = begchunk, endchunk
         ncol = phys_state(lchnk)%ncol
         call get_lon_all_p(lchnk, ncol, lons)
         call get_lat_all_p(lchnk, ncol, lats)

         pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)

         if (use_gw_front .or. use_gw_front_igw) then
            call pbuf_get_field(pbuf_chnk, frontgf_idx, pbuf_frontgf)
            call pbuf_get_field(pbuf_chnk, frontga_idx, pbuf_frontga)
         end if

         if (qbo_use_forcing) then
            call pbuf_get_field(pbuf_chnk, uzm_idx, pbuf_uzm)
         end if

         do i = 1, ncol
            ic = lons(i)
            jc = lats(i)
            phys_state(lchnk)%ps(i)   = psxy(ic,jc)
            phys_state(lchnk)%phis(i) = phisxy(ic,jc)
            pic(i) = pkxy(ic,jc,pver+1)
         end do
         do k = 1, km
            do i = 1, ncol
               ic = lons(i)
               jc = lats(i)
               phys_state(lchnk)%u    (i,k) = u3(ic,k,jc)
               phys_state(lchnk)%v    (i,k) = v3(ic,k,jc)
               phys_state(lchnk)%omega(i,k) = omgaxy(ic,k,jc)
               phys_state(lchnk)%t    (i,k) = ptxy(ic,jc,k) / (1.0_r8 + zvir*tracer(ic,jc,k,1))
               phys_state(lchnk)%exner(i,k) = pic(i) / pkzxy(ic,jc,k)

               if (use_gw_front .or. use_gw_front_igw) then
                  pbuf_frontgf(i,k) = frontgf(ic,k,jc)
                  pbuf_frontga(i,k) = frontga(ic,k,jc)
               end if

               if (qbo_use_forcing) then
                  pbuf_uzm(i,k) = uzm(k,jc)
               end if

            end do
         end do

         do k = 1, kmp1
            do i = 1, ncol

               ! edge-level pressure arrays: copy from the arrays computed by dynpkg
               ic = lons(i)
               jc = lats(i)
               phys_state(lchnk)%pint  (i,k) = pexy  (ic,k,jc)
               phys_state(lchnk)%lnpint(i,k) = pelnxy(ic,k,jc)
            end do
         end do

         ! Copy constituents
         ! Dry types converted from moist to dry m.r. at bottom of this routine
         do m = 1, pcnst
            do k = 1, km
               do i = 1, ncol
                  phys_state(lchnk)%q(i,k,m) = tracer(lons(i),lats(i),k,m)
               end do
            end do
         end do

      end do ! loop over chunks

   else ! map not local

      boff  = 6
      if (use_gw_front .or. use_gw_front_igw) boff = boff + 2
      if (qbo_use_forcing)  boff = boff + 1

      tsize = boff + 1 + pcnst

      blksiz = (jlastxy-jfirstxy+1)*(ilastxy-ifirstxy+1)
      allocate( bpter(blksiz,0:km),stat=astat )
      if( astat /= 0 ) then
         write(iulog,*) 'd_p_coupling: failed to allocate bpter; error = ',astat
         call endrun('d_p_coupling: failed to allocate bpter')
      end if
      allocate( bbuffer(tsize*block_buf_nrecs),stat=astat )
      if( astat /= 0 ) then
         write(iulog,*) 'd_p_coupling: failed to allocate bbuffer; error = ',astat
         call endrun('d_p_coupling: failed to allocate bbuffer')
      end if
      allocate( cbuffer(tsize*chunk_buf_nrecs),stat=astat )
      if( astat /= 0 ) then
         write(iulog,*) 'd_p_coupling: failed to allocate cbuffer; error = ',astat
         call endrun('d_p_coupling: failed to allocate cbuffer')
      end if

      if (iam .lt. grid%npes_xy) then
         call block_to_chunk_send_pters(iam+1,blksiz,kmp1,tsize,bpter)
      endif

!$omp parallel do private (j, i, ib, k, m)
      do j = jfirstxy, jlastxy
         do i = ifirstxy, ilastxy
            ib = (j-jfirstxy)*(ilastxy-ifirstxy+1) + (i-ifirstxy+1)

            bbuffer(bpter(ib,0)+4:bpter(ib,0)+boff+pcnst) = 0.0_r8

            bbuffer(bpter(ib,0))   = pexy(i,kmp1,j)
            bbuffer(bpter(ib,0)+1) = pelnxy(i,kmp1,j)
            bbuffer(bpter(ib,0)+2) = psxy(i,j)
            bbuffer(bpter(ib,0)+3) = phisxy(i,j)

            do k = 1, km

               bbuffer(bpter(ib,k))   = pexy(i,k,j)
               bbuffer(bpter(ib,k)+1) = pelnxy(i,k,j)
               bbuffer(bpter(ib,k)+2) = u3    (i,k,j)
               bbuffer(bpter(ib,k)+3) = v3    (i,k,j)
               bbuffer(bpter(ib,k)+4) = omgaxy(i,k,j)
               bbuffer(bpter(ib,k)+5) = ptxy(i,j,k) / (1.0_r8 + zvir*tracer(i,j,k,1))
               bbuffer(bpter(ib,k)+6) = pkxy(i,j,pver+1) / pkzxy(i,j,k)

               if (use_gw_front .or. use_gw_front_igw) then
                  bbuffer(bpter(ib,k)+7) = frontgf(i,k,j)
                  bbuffer(bpter(ib,k)+8) = frontga(i,k,j)
               end if

               if (qbo_use_forcing) then
                  bbuffer(bpter(ib,k)+9) = uzm(k,j)
               end if

               do m = 1, pcnst
                  bbuffer(bpter(ib,k)+boff+m) = tracer(i,j,k,m)
               end do

            end do
         end do
      end do

      call t_barrierf('sync_blk_to_chk', grid%commxy)
      call t_startf ('block_to_chunk')
      call transpose_block_to_chunk(tsize, bbuffer, cbuffer)
      call t_stopf  ('block_to_chunk')

!$omp parallel do private (lchnk, ncol, i, k, m, cpter, pbuf_chnk, pbuf_uzm, pbuf_frontgf, pbuf_frontga)
      do lchnk = begchunk, endchunk
         ncol = phys_state(lchnk)%ncol

         pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)

         if (use_gw_front .or. use_gw_front_igw) then
            call pbuf_get_field(pbuf_chnk, frontgf_idx, pbuf_frontgf)
            call pbuf_get_field(pbuf_chnk, frontga_idx, pbuf_frontga)
         end if

         if (qbo_use_forcing) then
            call pbuf_get_field(pbuf_chnk, uzm_idx, pbuf_uzm)
         end if

         call block_to_chunk_recv_pters(lchnk,pcols,pver+1,tsize,cpter)

         do i = 1, ncol

            phys_state(lchnk)%pint  (i,pver+1) = cbuffer(cpter(i,0))
            phys_state(lchnk)%lnpint(i,pver+1) = cbuffer(cpter(i,0)+1)
            phys_state(lchnk)%ps(i)            = cbuffer(cpter(i,0)+2)
            phys_state(lchnk)%phis(i)          = cbuffer(cpter(i,0)+3)

            do k = 1, km

               phys_state(lchnk)%pint  (i,k) = cbuffer(cpter(i,k))
               phys_state(lchnk)%lnpint(i,k) = cbuffer(cpter(i,k)+1)
               phys_state(lchnk)%u     (i,k) = cbuffer(cpter(i,k)+2)
               phys_state(lchnk)%v     (i,k) = cbuffer(cpter(i,k)+3)
               phys_state(lchnk)%omega (i,k) = cbuffer(cpter(i,k)+4)
               phys_state(lchnk)%t     (i,k) = cbuffer(cpter(i,k)+5)
               phys_state(lchnk)%exner (i,k) = cbuffer(cpter(i,k)+6)

               if (use_gw_front .or. use_gw_front_igw) then
                  pbuf_frontgf(i,k)  = cbuffer(cpter(i,k)+7)
                  pbuf_frontga(i,k)  = cbuffer(cpter(i,k)+8)
               end if

               if (qbo_use_forcing) then
                  pbuf_uzm(i,k) = cbuffer(cpter(i,k)+9)
               end if

               ! dry type constituents converted from moist to dry at bottom of routine
               do m = 1, pcnst
                  phys_state(lchnk)%q(i,k,m) = cbuffer(cpter(i,k)+boff+m)
               end do

            end do
         end do

      end do ! loop over chunks

      deallocate(bpter)
      deallocate(bbuffer)
      deallocate(cbuffer)

   end if ! (local_dp_map)

   ! Get indices to access O, O2, H, H2, and N species
   if (waccmx_is('ionosphere') .or. waccmx_is('neutral')) then
      call cnst_get_ind('O', ixo)
      call cnst_get_ind('O2', ixo2)
      call cnst_get_ind('H', ixh)
      call cnst_get_ind('H2', ixh2)
      call cnst_get_ind('N', ixn)
   end if

   ! Evaluate derived quantities
   call t_startf ('derived_fields')
!$omp parallel do private (lchnk, ncol, i, k, m, qmavl, dqreq, qbot, qbotm1, zvirv, pbuf_chnk, mmrSum_O_O2_H)
   do lchnk = begchunk, endchunk
      ncol = phys_state(lchnk)%ncol
      do k = 1, km
         do i = 1, ncol
            phys_state(lchnk)%pdel (i,k) = phys_state(lchnk)%pint(i,k+1) - phys_state(lchnk)%pint(i,k)
            phys_state(lchnk)%rpdel(i,k) = 1.0_r8/phys_state(lchnk)%pdel(i,k)
            phys_state(lchnk)%pmid (i,k) = 0.5_r8*(phys_state(lchnk)%pint(i,k) + phys_state(lchnk)%pint(i,k+1))
            phys_state(lchnk)%lnpmid(i,k) = log(phys_state(lchnk)%pmid(i,k))
         end do
      end do

      ! Attempt to remove negative constituents in bottom layer only by moving from next level
      ! This is a BAB kludge to avoid masses of warning messages for cloud water and ice, since
      ! the vertical remapping operator currently being used for cam is not strictly monotonic
      ! at the endpoints.
      do m = 1, pcnst
         do i = 1, ncol
            if (phys_state(lchnk)%q(i,pver,m) < qmin(m)) then
               ! available q in 2nd level
               qmavl = phys_state(lchnk)%q(i,pver-1,m) - qmin(m)
               ! required q change in bottom level rescaled to mass fraction in 2nd level
               dqreq = (qmin(m) - phys_state(lchnk)%q(i,pver,m))                         &
                      * phys_state(lchnk)%pdel(i,pver) / phys_state(lchnk)%pdel(i,pver-1)
               qbot   = phys_state(lchnk)%q(i,pver  ,m)
               qbotm1 = phys_state(lchnk)%q(i,pver-1,m)
               if (dqreq < qmavl) then
                  phys_state(lchnk)%q(i,pver  ,m) = qmin(m)
                  phys_state(lchnk)%q(i,pver-1,m) = phys_state(lchnk)%q(i,pver-1,m) - dqreq
                  ! these log messages are off by default.  use namelist to turn them on
                  if (dqreq>qmin(m) .and. dqreq>fraction*qbotm1 .and. (trim(fv_print_dpcoup_warn) == "full")) &
                      write(iulog,*) 'dpcoup dqreq', m, lchnk, i, qbot, qbotm1, dqreq
               else
                  if (dqreq>qmin(m) .and. (trim(fv_print_dpcoup_warn) == "full")) then
                     write(iulog,*) 'dpcoup cant adjust', m, lchnk, i, qbot, qbotm1, dqreq
                  end if
               end if
            end if
         end do
      end do

      ! Ensure O2 + O + H (N2) mmr greater than one.
      ! Check for unusually large H2 values and set to lower value
      if (waccmx_is('ionosphere') .or. waccmx_is('neutral')) then
         do i = 1, ncol
            do k = 1, pver

               if (phys_state(lchnk)%q(i,k,ixo) < mmrMin) phys_state(lchnk)%q(i,k,ixo) = mmrMin
               if (phys_state(lchnk)%q(i,k,ixo2) < mmrMin) phys_state(lchnk)%q(i,k,ixo2) = mmrMin

               mmrSum_O_O2_H = phys_state(lchnk)%q(i,k,ixo) + phys_state(lchnk)%q(i,k,ixo2) + &
                               phys_state(lchnk)%q(i,k,ixh)

               if ((1._r8-mmrMin-mmrSum_O_O2_H) < 0._r8) then

                  phys_state(lchnk)%q(i,k,ixo) = phys_state(lchnk)%q(i,k,ixo) * (1._r8 - N2mmrMin) / &
                                                 mmrSum_O_O2_H
                  phys_state(lchnk)%q(i,k,ixo2) = phys_state(lchnk)%q(i,k,ixo2) * (1._r8 - N2mmrMin) / &
                                                  mmrSum_O_O2_H
                  phys_state(lchnk)%q(i,k,ixh) = phys_state(lchnk)%q(i,k,ixh) * (1._r8 - N2mmrMin) / &
                                                  mmrSum_O_O2_H
               end if

               if (phys_state(lchnk)%q(i,k,ixh2) .gt. 6.e-5_r8) then
                  phys_state(lchnk)%q(i,k,ixh2) = 6.e-5_r8
               end if

            end do
         end do
      end if

      ! Call physconst_update to compute cpairv, rairv, mbarv, and cappav as constituent dependent variables
      ! and compute molecular viscosity (kmvis) and conductivity (kmcnd)
      if (waccmx_is('ionosphere') .or. waccmx_is('neutral')) then
         call physconst_update(phys_state(lchnk)%q, phys_state(lchnk)%t, lchnk, ncol)
      end if

      ! Fill local zvirv variable; calculated for WACCM-X
      if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
         zvirv(:,:) = shr_const_rwv / rairv(:,:,lchnk) -1._r8
      else
         zvirv(:,:) = zvir
      end if

      ! Compute initial geopotential heights
      call geopotential_t(phys_state(lchnk)%lnpint, phys_state(lchnk)%lnpmid  , phys_state(lchnk)%pint  , &
                          phys_state(lchnk)%pmid  , phys_state(lchnk)%pdel    , phys_state(lchnk)%rpdel , &
                          phys_state(lchnk)%t     , phys_state(lchnk)%q(:,:,1), rairv(:,:,lchnk), gravit, &
                          zvirv, phys_state(lchnk)%zi, phys_state(lchnk)%zm, ncol)

      ! Compute initial dry static energy, include surface geopotential
      do k = 1, pver
         do i = 1, ncol
            phys_state(lchnk)%s(i,k) = cpairv(i,k,lchnk)*phys_state(lchnk)%t(i,k) &
                                       + gravit*phys_state(lchnk)%zm(i,k) + phys_state(lchnk)%phis(i)
         end do
      end do

      ! Convert dry type constituents from moist to dry mixing ratio
      call set_state_pdry(phys_state(lchnk))    ! First get dry pressure to use for this timestep
      call set_wet_to_dry(phys_state(lchnk))    ! Dynamics had moist, physics wants dry.

      ! Ensure tracers are all positive
      call qneg3('D_P_COUPLING', lchnk, ncol, pcols, pver, 1, pcnst, qmin, phys_state(lchnk)%q)

      ! Compute energy and water integrals of input state

      pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
      call check_energy_timestep_init(phys_state(lchnk), pbuf_chnk)

   end do ! loop over chunks
   call t_stopf('derived_fields')

   deallocate(u3)
   deallocate(v3)
   if (dyn_state%am_diag) then
      deallocate(du3)
      deallocate(dv3)
      deallocate(dua3)
      deallocate(dva3)
      deallocate(duf3)
      deallocate(dummy)
   end if

end subroutine d_p_coupling

!=========================================================================================

subroutine p_d_coupling(grid, phys_state, phys_tend, pbuf2d, &
                        dyn_in, dtime, zvir, cappa, ptop)

   ! convert physics output to dynamics input
   ! Variables ending in xy are xy-decomposition instanciations.

   ! arguments
   type(T_FVDYCORE_GRID), intent(in)    :: grid ! FV Dynamics grid
   type(physics_state),   intent(inout) :: phys_state(begchunk:endchunk)
   type(physics_tend),    intent(inout) :: phys_tend(begchunk:endchunk)
   type(physics_buffer_desc), pointer   :: pbuf2d(:,:)
   type(dyn_import_t),    intent(inout) :: dyn_in

   real(r8), intent(in) :: dtime
   real(r8), intent(in) :: zvir
   real(r8), intent(in) :: cappa
   real(r8), intent(in) :: ptop

   ! LOCAL VARIABLES:

   type(physics_buffer_desc), pointer :: pbuf_chnk(:)
   real(r8),                  pointer :: qini(:,:)

   type(t_fvdycore_state),    pointer :: dyn_state

   ! Variables from the dynamics import container
   real(r8), pointer :: psxy(:,:)
   real(r8), pointer :: u3sxy(:,:,:)
   real(r8), pointer :: v3sxy(:,:,:)
   real(r8), pointer :: t3xy(:,:,:)                  !  Temperature
   real(r8), pointer :: ptxy(:,:,:)                  !  Virt. pot. temp.
   real(r8), pointer :: tracer(:,:,:,:)              !  Constituents

   real(r8), pointer :: pexy(:,:,:)
   real(r8), pointer :: delpxy(:,:,:)
   real(r8), pointer :: pkxy(:,:,:)
   real(r8), pointer :: pkzxy(:,:,:)

   real(r8):: dudtxy(grid%ifirstxy:grid%ilastxy,grid%km,grid%jfirstxy:grid%jlastxy)
   real(r8):: dvdtxy(grid%ifirstxy:grid%ilastxy,grid%km,grid%jfirstxy:grid%jlastxy)
   real(r8):: dummy_pelnxy(grid%ifirstxy:grid%ilastxy,grid%km+1,grid%jfirstxy:grid%jlastxy)

   integer :: i, ib, k, m, j, lchnk  ! indices
   integer :: ncol                   ! number of columns in current chunk
   integer :: lats(pcols)            ! array of latitude indices
   integer :: lons(pcols)            ! array of longitude indices
   integer :: blksiz                 ! number of columns in 2D block
   integer :: tsize                  ! amount of data per grid point passed to physics
   integer, allocatable, dimension(:,:) :: bpter
                                     ! offsets into block buffer for unpacking data
   integer :: cpter(pcols,0:pver)   ! offsets into chunk buffer for packing data

   real(r8) :: dt5
   real(r8), allocatable, dimension(:) :: bbuffer, cbuffer ! transpose buffers
#if (! defined SPMD)
   integer :: block_buf_nrecs = 0
   integer :: chunk_buf_nrecs = 0
   logical :: local_dp_map=.true.
#endif
   integer :: km, iam
   integer :: ifirstxy, ilastxy, jfirstxy, jlastxy

   real(r8) :: cappa3v(grid%ifirstxy:grid%ilastxy,grid%jfirstxy:grid%jlastxy,grid%km)
   !----------------------------------------------------------------------------

   ! Adjustments to the physics state before copying to dynamics state.
!$omp parallel do private(lchnk, pbuf_chnk, qini)
   do lchnk = begchunk, endchunk

      ! FV expects wet mmr
      call set_dry_to_wet(phys_state(lchnk))

      ! Adjust pressure and moist mmr to account for water vapor change due to physics.
      ! (assumes water vapor is only constituent that contributes to the atm mass)
      pbuf_chnk => pbuf_get_chunk(pbuf2d, lchnk)
      call pbuf_get_field(pbuf_chnk, qini_idx, qini)
      call physics_dme_adjust(phys_state(lchnk), qini)
   end do

   dyn_state => get_dyn_state()

   ! Pull the variables out of the dynamics export container
   psxy    => dyn_in%ps
   u3sxy   => dyn_in%u3s
   v3sxy   => dyn_in%v3s
   t3xy    => dyn_in%t3
   ptxy    => dyn_in%pt
   tracer  => dyn_in%tracer

   pexy    => dyn_in%pe
   delpxy  => dyn_in%delp
   pkxy    => dyn_in%pk
   pkzxy   => dyn_in%pkz

   km   = grid%km

   ifirstxy = grid%ifirstxy
   ilastxy  = grid%ilastxy
   jfirstxy = grid%jfirstxy
   jlastxy  = grid%jlastxy

   iam      = grid%iam

#if ( defined OFFLINE_DYN )
   ! set the dyn flds to offline meteorological data
   call get_met_fields(phys_state, phys_tend, dtime)
#endif

   ! Copy temperature, tendencies and constituents to dynamics data structures
   ! Copy onto xy decomposition, then transpose to yz decomposition

   if (local_dp_map) then

!$omp parallel do private(lchnk, i, k, ncol, m, lons, lats)
      do lchnk = begchunk,endchunk
         ncol = get_ncols_p(lchnk)
         call get_lon_all_p(lchnk, ncol, lons)
         call get_lat_all_p(lchnk, ncol, lats)

         do k = 1, km
            do i = 1, ncol
               dvdtxy(lons(i),k,lats(i)) = phys_tend(lchnk)%dvdt(i,k)
               dudtxy(lons(i),k,lats(i)) = phys_tend(lchnk)%dudt(i,k)
               ptxy  (lons(i),lats(i),k) = phys_state(lchnk)%t(i,k)
               delpxy(lons(i),lats(i),k) = phys_state(lchnk)%pdel(i,k)
            end do
         end do

         do m = 1, pcnst
            do k = 1, km
               do i = 1, ncol
                  tracer(lons(i),lats(i),k,m) = phys_state(lchnk)%q(i,k,m)
               end do
            end do
         end do

      end do

   else ! map not local

      tsize = 4 + pcnst

      blksiz = (jlastxy-jfirstxy+1)*(ilastxy-ifirstxy+1)
      allocate(bpter(blksiz,0:km))
      allocate(bbuffer(tsize*block_buf_nrecs))
      allocate(cbuffer(tsize*chunk_buf_nrecs))

!$omp parallel do private (lchnk, ncol, i, k, m, cpter)
      do lchnk = begchunk,endchunk
         ncol = get_ncols_p(lchnk)

         call chunk_to_block_send_pters(lchnk,pcols,km+1,tsize,cpter)

         do i=1,ncol
            cbuffer(cpter(i,0):cpter(i,0)+3+pcnst) = 0.0_r8
         end do

         do k=1,km
            do i=1,ncol

               cbuffer(cpter(i,k))   = phys_tend(lchnk)%dvdt(i,k)
               cbuffer(cpter(i,k)+1) = phys_tend(lchnk)%dudt(i,k)
               cbuffer(cpter(i,k)+2) = phys_state(lchnk)%t(i,k)
               cbuffer(cpter(i,k)+3) = phys_state(lchnk)%pdel(i,k)

               do m=1,pcnst
                  cbuffer(cpter(i,k)+3+m) = phys_state(lchnk)%q(i,k,m)
               end do

            end do

         end do

      end do

      call t_barrierf('sync_chk_to_blk', grid%commxy)
      call t_startf ('chunk_to_block')
      call transpose_chunk_to_block(tsize, cbuffer, bbuffer)
      call t_stopf  ('chunk_to_block')

      if (iam .lt. grid%npes_xy) then
         call chunk_to_block_recv_pters(iam+1,blksiz,km+1,tsize,bpter)
      endif

!$omp parallel do private (j, i, ib, k, m)
      do j = jfirstxy, jlastxy
         do k = 1, km
            do i = ifirstxy, ilastxy
               ib = (j-jfirstxy)*(ilastxy-ifirstxy+1) + (i-ifirstxy+1)

               dvdtxy(i,k,j) = bbuffer(bpter(ib,k))
               dudtxy(i,k,j) = bbuffer(bpter(ib,k)+1)
               ptxy  (i,j,k) = bbuffer(bpter(ib,k)+2)
               delpxy(i,j,k) = bbuffer(bpter(ib,k)+3)

               do m = 1, pcnst
                  tracer(i,j,k,m) = bbuffer(bpter(ib,k)+3+m)
               end do

            end do
         end do
      end do

      deallocate(bpter)
      deallocate(bbuffer)
      deallocate(cbuffer)

   end if ! local_dp_map

!$omp parallel do private(i,j,k)
   do k = 1, km
      do j = jfirstxy, jlastxy
         do i = ifirstxy, ilastxy
            t3xy(i,j,k) = ptxy(i,j,k)
         end do
      end do
   end do

   ! Update u3s and v3s from tendencies dudt and dvdt.

   dt5 = 0.5_r8*dtime

   call t_barrierf('sync_uv3s_update', grid%commxy)
   call t_startf('uv3s_update')
   if (iam .lt. grid%npes_xy) then
      call uv3s_update(grid, dudtxy, u3sxy, dvdtxy, v3sxy, dt5, &
                       dyn_state%am_geom_crrct)
   end if
   call t_stopf('uv3s_update')

   ! -------------------------------------------------------------------------
   ! Compute pt, q3, pe, delp, ps, peln, pkz and pk.
   ! For 2-D decomposition, delp is transposed to delpxy, pexy is computed
   ! from delpxy (and ptop), and pexy is transposed back to pe.
   ! Note that pt, q3, delp and pe are input parameters as well.
   ! -------------------------------------------------------------------------
   call t_barrierf('sync_p_d_adjust', grid%commxy)
   call t_startf('p_d_adjust')
   if (iam .lt. grid%npes_xy) then
      if (grid%high_alt) then
         call physconst_calc_kappav(ifirstxy, ilastxy, jfirstxy, jlastxy, 1, km, grid%ntotq, &
                                    tracer, cappa3v)
      else
         cappa3v = cappa
      endif
      call p_d_adjust(grid, tracer, dummy_pelnxy, pkxy, pkzxy, zvir, cappa3v, &
                      delpxy, ptxy, pexy, psxy, ptop)
   end if
   call t_stopf('p_d_adjust')

end subroutine p_d_coupling

!=========================================================================================

end module dp_coupling
