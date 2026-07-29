!
!  Copyright 2019-2020 SALMON developers
!
!  Licensed under the Apache License, Version 2.0 (the "License");
!  you may not use this file except in compliance with the License.
!  You may obtain a copy of the License at
!
!      http://www.apache.org/licenses/LICENSE-2.0
!
!  Unless required by applicable law or agreed to in writing, software
!  distributed under the License is distributed on an "AS IS" BASIS,
!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!  See the License for the specific language governing permissions and
!  limitations under the License.
!

#include "config.h"

module nonlocal_potential
  use nvtx_wrapper
  implicit none

! WARNING: We must not call these except for hpsi routine.

contains

subroutine dpseudo(tpsi,htpsi,info,nspin,ppg)
  use structures
  use communication, only: comm_summation
  use timer
  implicit none
  integer,intent(in) :: nspin
  type(s_parallel_info),intent(in) :: info
  type(s_pp_grid),intent(in) :: ppg
  type(s_orbital),intent(in) :: tpsi
  type(s_orbital) :: htpsi
  !
  integer :: ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e,norb
  integer :: ilma,ia,j,ix,iy,iz,Nlma,ilocal
  real(8) :: uVpsi,wrk
  real(8),allocatable :: uVpsibox (:,:,:,:,:)
  real(8),allocatable :: uVpsibox2(:,:,:,:,:)

  call timer_begin(LOG_UHPSI_PSEUDO)

  im_s = info%im_s
  im_e = info%im_e
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e
  norb = Nspin* info%numo * info%numk * info%numm

  Nlma = ppg%Nlma

  if(info%if_divide_rspace) then

    allocate(uVpsibox (Nlma,Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e))
    allocate(uVpsibox2(Nlma,Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e))

!$omp parallel do collapse(4) &
!$omp             private(im,ik,io,ispin,ilma,ia,uVpsi,j,ix,iy,iz)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin

      do ilma=1,Nlma
        ia = ppg%ia_tbl(ilma)
        uVpsi = 0.d0
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          uVpsi = uVpsi + ppg%uV(j,ilma) * tpsi%rwf(ix,iy,iz,ispin,io,ik,im)
        end do
        uVpsi = uVpsi * ppg%rinv_uvu(ilma)
        uVpsibox(ilma,ispin,io,ik,im) = uVpsi
      end do

    end do
    end do
    end do
    end do
!$omp end parallel do

    call timer_end(LOG_UHPSI_PSEUDO)

    call timer_begin(LOG_UHPSI_PSEUDO_COMM)
    call comm_summation(uVpsibox,uVpsibox2,Nlma*Norb,info%icomm_r)
    call timer_end(LOG_UHPSI_PSEUDO_COMM)

    call timer_begin(LOG_UHPSI_PSEUDO)

!$omp parallel do collapse(4) &
!$omp             private(im,ik,io,ispin,ilocal,ilma,ia,uVpsi,j,ix,iy,iz,wrk)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin

      do ilocal=1,ppg%ilocal_nlma
        ilma=ppg%ilocal_nlma2ilma(ilocal)
        ia  =ppg%ilocal_nlma2ia  (ilocal)
        uVpsi = uVpsibox2(ilma,ispin,io,ik,im)
!OCL norecurrence
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          wrk = uVpsi * ppg%uV(j,ilma)
          htpsi%rwf(ix,iy,iz,ispin,io,ik,im) = htpsi%rwf(ix,iy,iz,ispin,io,ik,im) + wrk
        end do
      end do

    end do
    end do
    end do
    end do
!$omp end parallel do

    deallocate(uVpsibox,uVpsibox2)

  else

#ifdef USE_OPENACC
!$acc parallel loop collapse(4) private(im,ik,io,ispin,ilma,ia,uVpsi,j,ix,iy,iz,wrk)
#else
!$omp parallel do collapse(4) &
!$omp             private(im,ik,io,ispin,ilma,ia,uVpsi,j,ix,iy,iz,wrk)
#endif
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin

      do ilma=1,Nlma
        ia = ppg%ia_tbl(ilma)
        uVpsi = 0.d0
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          uVpsi = uVpsi + ppg%uV(j,ilma) * tpsi%rwf(ix,iy,iz,ispin,io,ik,im)
        end do
        uVpsi = uVpsi * ppg%rinv_uvu(ilma)
!OCL norecurrence
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          wrk = uVpsi * ppg%uV(j,ilma)
          htpsi%rwf(ix,iy,iz,ispin,io,ik,im) = htpsi%rwf(ix,iy,iz,ispin,io,ik,im) + wrk
        end do
      end do

    end do
    end do
    end do
    end do
#ifdef USE_OPENACC
!$acc end parallel
#else
!$omp end parallel do
#endif

  end if

  call timer_end(LOG_UHPSI_PSEUDO)

  return
end subroutine dpseudo

!-----------------------------------------------------------------------------------------------------------------------------------

subroutine zpseudo(tpsi,htpsi,info,nspin,ppg)
  use structures
  use timer
  use iso_c_binding
#if defined(USE_OPENACC) && defined(USE_GEMM)
  use cublas
  use openacc
#endif
  implicit none
  intrinsic :: aimag
  integer,intent(in) :: nspin
  type(s_parallel_info),intent(in) :: info
  type(s_pp_grid),intent(in) :: ppg
  type(s_orbital),intent(in) :: tpsi
  type(s_orbital) :: htpsi
  !
  integer :: ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e
  integer :: ilma,ia,j,ix,iy,iz,Nlma,ilocal,vi,my_nlma,k
  complex(8) :: uVpsi,wrk
  complex(8),allocatable :: uVpsibox (:,:,:,:,:)
  complex(8),allocatable :: uVpsibox2(:,:,:,:,:)
#ifdef USE_OPENACC
  real(8),allocatable,save :: htpsi_zwf_r(:,:,:,:,:,:,:), htpsi_zwf_i(:,:,:,:,:,:,:)
  integer :: i1, i2, i3, i4, i5, i6, i7
  integer :: mps_max
  real(8) :: wrk_r
#endif
  complex(8) :: IMAGINARY_UNIT = (0, 1)
#if defined(USE_OPENACC) && defined(USE_CUDA)
  integer :: n,natom
  interface
    subroutine zpseudo_cuda(htpsi_zwf,n,im_s,im_e,ik_s,ik_e,io_s,io_e,Nspin,Nlma,ppg_nps,natom,mg_is_array_1,mg_ie_array_1,mg_is_array_2,mg_ie_array_2,mg_is_array_3,mg_ie_array_3,ppg_ia_tbl,ppg_mps,ppg_jxyz,ppg_zekr_uV,ppg_rinv_uvu,tpsi_zwf) bind(c)
      import
      ! Input integer
      integer(c_int), intent(in), value :: n
      integer(c_int), intent(in), value :: im_s
      integer(c_int), intent(in), value :: im_e
      integer(c_int), intent(in), value :: ik_s
      integer(c_int), intent(in), value :: ik_e
      integer(c_int), intent(in), value :: io_s
      integer(c_int), intent(in), value :: io_e
      integer(c_int), intent(in), value :: Nspin
      integer(c_int), intent(in), value :: Nlma
      integer(c_int), intent(in), value :: ppg_nps
      integer(c_int), intent(in), value :: natom
      integer(c_int), intent(in), value :: mg_is_array_1
      integer(c_int), intent(in), value :: mg_ie_array_1
      integer(c_int), intent(in), value :: mg_is_array_2
      integer(c_int), intent(in), value :: mg_ie_array_2
      integer(c_int), intent(in), value :: mg_is_array_3
      integer(c_int), intent(in), value :: mg_ie_array_3
      ! Output
      complex(c_double_complex), intent(inout) :: htpsi_zwf(mg_is_array_1:mg_ie_array_1, mg_is_array_2:mg_ie_array_2, mg_is_array_2:mg_ie_array_3, Nspin, io_e:io_s, ik_s:ik_e, im_s:im_e)
      ! Input
      integer(c_int)           , intent(in) :: ppg_ia_tbl(n * natom)
      integer(c_int)           , intent(in) :: ppg_mps(natom)
      integer(c_int)           , intent(in) :: ppg_jxyz(3, ppg_nps, natom)
      complex(c_double_complex), intent(in) :: ppg_zekr_uV(ppg_nps, ppg_Nlma, ik_s:ik_e)
      real(c_double)           , intent(in) :: ppg_rinv_uvu(n * natom)
      complex(c_double_complex), intent(in) :: tpsi_zwf(mg_is_array_1:mg_ie_array_1, mg_is_array_2:mg_ie_array_2,mg_is_array_2:mg_ie_array_3, Nspin, io_e:io_s, ik_s:ik_e, im_s:im_e)
    end subroutine zpseudo_cuda
  end interface
#endif
#if defined(USE_OPENACC) && defined(USE_GEMM)
  ! Batched-GEMM reformulation of phase 1 (projection), prototyped and
  ! validated in zpseudo-mini (see subwg2-benchmarks repo,
  ! salmon-gpu-optimization-ideas skill idea 2): a ~1.55x win over the
  ! vector(32) reduction and ~6.7x over the hand-written CUDA kernel at
  ! natom=486/nproj=15 (matching this real code's actual atom/projector
  ! counts). Phase 2 (back-projection) is completely unchanged, still the
  ! existing atomic-free inverse-gather-map algorithm below.
  !
  ! Real SALMON's own ppg%jxyz(3,ppg%nps,natom)/ppg%zekr_uV(ppg%nps,Nlma,ik)
  ! are allocated via plain `allocate` and only filled up to
  ! ppg%mps(ia)/the atom's actual projector count (confirmed by reading
  ! prep_pp.f90/hamiltonian.f90) -- rows beyond that are genuinely
  ! uninitialized, not zero. A batched GEMM with K=ppg%nps unconditionally
  ! reads every row for every batch element, so using those arrays directly
  ! would multiply real data against garbage. Fixed by building our OWN
  ! zero-padded, dense companion arrays ONCE (SAVE'd, matching this file's
  ! own htpsi_zwf_r/i idiom above) rather than trusting the original
  ! arrays' unused-region contents.
  integer, parameter :: gemm_io_block = 64
  integer,save :: gemm_max_nproj = -1
  integer,allocatable,save :: gemm_nproj_atom(:), gemm_l2g(:,:)
  complex(8),allocatable,save :: gemm_zekr_packed(:,:,:)
  real(8),allocatable,save :: gemm_rinv_packed(:,:)
  complex(8),allocatable,save :: gemm_wf_packed(:,:,:), gemm_out_packed(:,:,:)
  type(cublasHandle),save :: gemm_handle
  logical,save :: gemm_handle_created = .false.
  logical,save :: gemm_debug_printed = .false.
  integer,save :: gemm_call_count = 0
  integer :: gemm_ia, gemm_p, gemm_ilma, gemm_j, gemm_natom
  integer :: gemm_io_blk_s, gemm_this_block, gemm_io_local, gemm_stat
#endif

#ifdef FORTRAN_COMPILER_HAS_2MB_ALIGNED_ALLOCATION
!dir$ attributes align : 2097152 :: uVpsibox, uVpsibox2
#endif

  call nvtxStartRange('zpseudo', __LINE__)
  call timer_begin(LOG_UHPSI_PSEUDO)

  im_s = info%im_s
  im_e = info%im_e
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e

  Nlma = ppg%Nlma

  if(info%if_divide_rspace) then

    call timer_end(LOG_UHPSI_PSEUDO)

    call calc_uVpsi_rdivided(nspin,info,ppg,tpsi,uVpsibox,uVpsibox2)

    call timer_begin(LOG_UHPSI_PSEUDO)

#ifdef USE_OPENACC
! `acc atomic` does not support complex(8). Split it into the real and imag part.
    if (.not. allocated(htpsi_zwf_r)) then
      allocate( &
        htpsi_zwf_r( &
        lbound(htpsi%zwf,1):ubound(htpsi%zwf,1), &
        lbound(htpsi%zwf,2):ubound(htpsi%zwf,2), &
        lbound(htpsi%zwf,3):ubound(htpsi%zwf,3), &
        lbound(htpsi%zwf,4):ubound(htpsi%zwf,4), &
        lbound(htpsi%zwf,5):ubound(htpsi%zwf,5), &
        lbound(htpsi%zwf,6):ubound(htpsi%zwf,6), &
        lbound(htpsi%zwf,7):ubound(htpsi%zwf,7)  &
        ) &
        )
    end if
    if (.not. allocated(htpsi_zwf_i)) then
      allocate( &
        htpsi_zwf_i( &
        lbound(htpsi%zwf,1):ubound(htpsi%zwf,1), &
        lbound(htpsi%zwf,2):ubound(htpsi%zwf,2), &
        lbound(htpsi%zwf,3):ubound(htpsi%zwf,3), &
        lbound(htpsi%zwf,4):ubound(htpsi%zwf,4), &
        lbound(htpsi%zwf,5):ubound(htpsi%zwf,5), &
        lbound(htpsi%zwf,6):ubound(htpsi%zwf,6), &
        lbound(htpsi%zwf,7):ubound(htpsi%zwf,7)  &
        ) &
        )
    end if
    mps_max = size(ppg%jxyz, 2)
!$acc kernels
    htpsi_zwf_i = 0d0
    htpsi_zwf_r = 0d0

!$acc loop collapse(6) independent gang vector private(im,ik,io,ispin,ilocal,ilma,ia,uVpsi,j,ix,iy,iz,wrk,uVpsi)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin
      do ilocal=1,ppg%ilocal_nlma
!OCL norecurrence
        do j=1,mps_max
          ilma=ppg%ilocal_nlma2ilma(ilocal)
          ia  =ppg%ilocal_nlma2ia  (ilocal)
          uVpsi = uVpsibox2(ispin,io,ik,im,ilma)
          if (j <= ppg%mps(ia)) then
            ix = ppg%jxyz(1,j,ia)
            iy = ppg%jxyz(2,j,ia)
            iz = ppg%jxyz(3,j,ia)
            wrk = uVpsi * ppg%zekr_uV(j,ilma,ik)

            wrk_r = real(wrk)
            !$acc atomic update
            htpsi_zwf_r(ix,iy,iz,ispin,io,ik,im) = htpsi_zwf_r(ix,iy,iz,ispin,io,ik,im) + wrk_r
            !$acc end atomic

            wrk_r = aimag(wrk)
            !$acc atomic update
            htpsi_zwf_i(ix,iy,iz,ispin,io,ik,im) = htpsi_zwf_i(ix,iy,iz,ispin,io,ik,im) + wrk_r
            !$acc end atomic
          end if
        end do
      end do

    end do
    end do
    end do
    end do
!$acc loop collapse(7) independent gang vector
    do i1 = lbound(htpsi%zwf,1), ubound(htpsi%zwf,1)
    do i2 = lbound(htpsi%zwf,2), ubound(htpsi%zwf,2)
    do i3 = lbound(htpsi%zwf,3), ubound(htpsi%zwf,3)
    do i4 = lbound(htpsi%zwf,4), ubound(htpsi%zwf,4)
    do i5 = lbound(htpsi%zwf,5), ubound(htpsi%zwf,5)
    do i6 = lbound(htpsi%zwf,6), ubound(htpsi%zwf,6)
    do i7 = lbound(htpsi%zwf,7), ubound(htpsi%zwf,7)

      htpsi%zwf(i1,i2,i3,i4,i5,i6,i7) = htpsi%zwf(i1,i2,i3,i4,i5,i6,i7) + cmplx( &
          htpsi_zwf_r(i1,i2,i3,i4,i5,i6,i7), &
          htpsi_zwf_i(i1,i2,i3,i4,i5,i6,i7), &
          kind=8 &
      )

    end do
    end do
    end do
    end do
    end do
    end do
    end do
!$acc end kernels
#else
!$omp parallel do collapse(4) &
!$omp             private(im,ik,io,ispin,ilocal,ilma,ia,uVpsi,j,ix,iy,iz,wrk)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin

      do ilocal=1,ppg%ilocal_nlma
        ilma=ppg%ilocal_nlma2ilma(ilocal)
        ia  =ppg%ilocal_nlma2ia  (ilocal)
        uVpsi = uVpsibox2(ispin,io,ik,im,ilma)
!OCL norecurrence
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          wrk = uVpsi * ppg%zekr_uV(j,ilma,ik)

          htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) + wrk
        end do
      end do

    end do
    end do
    end do
    end do
!$omp end parallel do
#endif

    deallocate(uVpsibox,uVpsibox2)

  else
#ifdef USE_OPENACC
#ifdef USE_CUDA
    natom=size(ppg%mps)
    n=size(ppg%rinv_uvu)/natom
    call zpseudo_cuda(htpsi%zwf,&
        n,&
        im_s,im_e,&
        ik_s,ik_e,&
        io_s,io_e,&
        Nspin,&
        Nlma,&
        ppg%nps,&
        natom,&
        lbound(htpsi%zwf,1),ubound(htpsi%zwf,1),&
        lbound(htpsi%zwf,2),ubound(htpsi%zwf,2),&
        lbound(htpsi%zwf,3),ubound(htpsi%zwf,3),&
        ppg%ia_tbl,&
        ppg%mps,&
        ppg%jxyz,&
        ppg%zekr_uV,&
        ppg%rinv_uvu,&
        tpsi%zwf)
#elif defined(USE_GEMM)
    gemm_natom = size(ppg%mps)

    ! --- one-time setup: build zero-padded, dense companion arrays ---
    if (gemm_max_nproj < 0) then
      allocate(gemm_nproj_atom(gemm_natom))
      gemm_nproj_atom = 0
      do ilma=1,Nlma
        ia = ppg%ia_tbl(ilma)
        gemm_nproj_atom(ia) = gemm_nproj_atom(ia) + 1
      end do
      gemm_max_nproj = maxval(gemm_nproj_atom)

      allocate(gemm_l2g(gemm_max_nproj, gemm_natom))
      gemm_l2g = 0
      gemm_nproj_atom = 0   ! reused as a per-atom fill cursor below
      do ilma=1,Nlma
        ia = ppg%ia_tbl(ilma)
        gemm_nproj_atom(ia) = gemm_nproj_atom(ia) + 1
        gemm_l2g(gemm_nproj_atom(ia), ia) = ilma
      end do

      ! Zero-init THEN fill only the valid (j<=mps(ia), p<=nproj_atom(ia))
      ! region -- ppg%zekr_uV's own out-of-range rows are uninitialized
      ! (see header comment), so we deliberately never read them here.
      ! No ik dimension here (unlike ppg%zekr_uV itself): this prototype
      ! assumes a single k-point (ik_s=ik_e), matching every real-SALMON
      ! test this session -- keeps gemm_zekr_packed a plain 3D array
      ! passed whole (not sliced) to cublas, avoiding an array-section
      ! host_data/cublas interaction pattern that triggered an nvfortran
      ! internal compiler error when this was (nps,max_nproj,natom,ik).
      allocate(gemm_zekr_packed(ppg%nps, gemm_max_nproj, gemm_natom))
      allocate(gemm_rinv_packed(gemm_max_nproj, gemm_natom))
      gemm_zekr_packed = (0.d0,0.d0)
      gemm_rinv_packed = 0.d0
      do gemm_ia=1,gemm_natom
      do gemm_p=1,gemm_nproj_atom(gemm_ia)
        gemm_ilma = gemm_l2g(gemm_p, gemm_ia)
        do gemm_j=1,ppg%mps(gemm_ia)
          gemm_zekr_packed(gemm_j, gemm_p, gemm_ia) = ppg%zekr_uV(gemm_j, gemm_ilma, ik_s)
        end do
        gemm_rinv_packed(gemm_p, gemm_ia) = ppg%rinv_uvu(gemm_ilma)
      end do
      end do
!$acc enter data copyin(gemm_zekr_packed, gemm_rinv_packed, gemm_l2g, gemm_nproj_atom)

      ! gemm_wf_packed zero-initialized on the host THEN copied in once --
      ! its (j>mps(ia)) padding rows must reach the device as real zeros,
      ! not garbage, since every future call only ever overwrites the
      ! valid (j<=mps(ia)) sub-region, leaving the zero padding untouched
      ! (and correct) for the lifetime of the run.
      allocate(gemm_wf_packed(ppg%nps, gemm_io_block, gemm_natom))
      allocate(gemm_out_packed(gemm_max_nproj, gemm_io_block, gemm_natom))
      gemm_wf_packed = (0.d0, 0.d0)
!$acc enter data copyin(gemm_wf_packed) create(gemm_out_packed)

      gemm_stat = cublasCreate(gemm_handle)
      ! Without this, cublas launches on CUDA's default stream while the
      ! surrounding OpenACC gather/scatter kernels run on NVHPC's own
      ! accelerator queue -- nothing then guarantees the GEMM completes
      ! before the scatter kernel starts reading gemm_out_packed. Binding
      ! cublas to the same queue OpenACC's synchronous (non-async) regions
      ! use makes all three (gather, GEMM, scatter) serialize in-order,
      ! per NVIDIA's documented OpenACC/cuBLAS interop pattern.
      gemm_stat = cublasSetStream(gemm_handle, acc_get_cuda_stream(acc_async_sync))
      gemm_handle_created = .true.

      ! TEMPORARY DIAGNOSTIC -- remove before merging.
      write(*,'(A,I0,A,I0,A,I0,A,I0)') 'GEMM_DEBUG natom=',gemm_natom,' Nlma=',Nlma, &
          ' max_nproj=',gemm_max_nproj,' sum_nproj=',sum(gemm_nproj_atom)
      write(*,'(A,I0,A,I0)') 'GEMM_DEBUG mps(1)=',ppg%mps(1),' nps=',ppg%nps
      write(*,'(A,I0,A,I0)') 'GEMM_DEBUG l2g(1,1)=',gemm_l2g(1,1),' ia_tbl(l2g(1,1))=',ppg%ia_tbl(gemm_l2g(1,1))
      write(*,'(A,2ES16.8,A,2ES16.8)') 'GEMM_DEBUG zekr_packed(1,1,1)=',gemm_zekr_packed(1,1,1), &
          ' zekr_uV(1,l2g(1,1),ik_s)=',ppg%zekr_uV(1,gemm_l2g(1,1),ik_s)
      write(*,'(A,ES16.8,A,ES16.8)') 'GEMM_DEBUG rinv_packed(1,1)=',gemm_rinv_packed(1,1), &
          ' rinv_uvu(l2g(1,1))=',ppg%rinv_uvu(gemm_l2g(1,1))
      if (gemm_max_nproj > 1) then
        write(*,'(A,I0,A,I0)') 'GEMM_DEBUG l2g(2,1)=',gemm_l2g(2,1),' nproj_atom(1)=',gemm_nproj_atom(1)
      end if
    end if

    gemm_stat = cublasSetStream(gemm_handle, acc_get_cuda_stream(acc_async_sync))

    ! ROOT CAUSE of the electron-count/energy drift: ppg%zekr_uV is NOT
    ! constant during an RT run. structures.f90 documents it as
    ! zekr_uV = exp(-i(k+A/c)r)*uv, and time_evolution_step.f90 sets
    ! system%vec_Ac to the current step's vector potential and then calls
    ! update_kvector_nonlocalpt() -- which rewrites the whole array -- on
    ! EVERY time step. Packing it once at setup therefore freezes the t=0
    ! phase factors while the rest of the code advances with A(t), so the
    ! GEMM's projectors drift further out of date the longer the run goes.
    ! That matches every observation: exact at call 1, ~0.2% off by call
    ! 50, and a monotonically growing electron-count error.
    !
    ! Refresh the packed copy every call, on the device. Only zekr needs
    ! this: rinv_uvu is a genuine time-independent normalization (verified
    ! unchanged at call 50), and l2g/nproj_atom are pure topology. Padding
    ! rows stay untouched (and therefore still zero) because the guard
    ! only ever writes the valid (p,j) sub-region.
!$acc parallel loop collapse(3) present(ppg, gemm_zekr_packed, gemm_l2g, gemm_nproj_atom)
    do gemm_ia = 1, gemm_natom
    do gemm_p = 1, gemm_max_nproj
    do gemm_j = 1, ppg%nps
      if (gemm_p <= gemm_nproj_atom(gemm_ia) .and. gemm_j <= ppg%mps(gemm_ia)) then
        gemm_zekr_packed(gemm_j, gemm_p, gemm_ia) = &
            ppg%zekr_uV(gemm_j, gemm_l2g(gemm_p, gemm_ia), ik_s)
      end if
    end do
    end do
    end do

    ! --- per-call: batched GEMM projection, blocked over bands ---
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do ispin=1,Nspin
      gemm_io_blk_s = io_s
      do while (gemm_io_blk_s <= io_e)
        gemm_this_block = min(gemm_io_block, io_e - gemm_io_blk_s + 1)

        ! Gather: loop bound must be invariant (ppg%nps, the same for
        ! every atom) for collapse(3) to be legal -- OpenACC requires a
        ! rectangular iteration space, but ppg%mps(gemm_ia) is ragged
        ! (varies per atom). Guard with an if instead: only the valid
        ! (j<=mps(ia)) region is touched, matching the ONLY safe read
        ! range of ppg%jxyz (reading beyond mps(ia) would be an
        ! out-of-bounds array-index use, not just a wrong answer -- see
        ! header comment). Padding stays the zero it was copied in as
        ! above (the loop simply skips writing to it here).
!$acc parallel loop collapse(3) present(tpsi,ppg,gemm_wf_packed)
        do gemm_ia = 1, gemm_natom
        do gemm_io_local = 1, gemm_this_block
        do gemm_j = 1, ppg%nps
          if (gemm_j <= ppg%mps(gemm_ia)) then
            gemm_wf_packed(gemm_j, gemm_io_local, gemm_ia) = tpsi%zwf( &
                ppg%jxyz(1,gemm_j,gemm_ia), ppg%jxyz(2,gemm_j,gemm_ia), ppg%jxyz(3,gemm_j,gemm_ia), &
                ispin, gemm_io_blk_s+gemm_io_local-1, ik, im)
          end if
        end do
        end do
        end do

        ! Batched GEMM: (max_nproj x nps) @ (nps x this_block) per atom,
        ! natom batches. CUBLAS_OP_C applies conjg() as part of the GEMM.
!$acc host_data use_device(gemm_zekr_packed, gemm_wf_packed, gemm_out_packed)
        gemm_stat = cublasZgemmStridedBatched(gemm_handle, CUBLAS_OP_C, CUBLAS_OP_N, &
            gemm_max_nproj, gemm_this_block, ppg%nps, &
            (1.d0,0.d0), gemm_zekr_packed, ppg%nps, int(ppg%nps,8)*int(gemm_max_nproj,8), &
            gemm_wf_packed, ppg%nps, int(ppg%nps,8)*int(gemm_io_block,8), &
            (0.d0,0.d0), gemm_out_packed, gemm_max_nproj, int(gemm_max_nproj,8)*int(gemm_io_block,8), &
            gemm_natom)
!$acc end host_data

        ! Scale by rinv_uvu, scatter the GEMM's dense (p,io_local,ia)
        ! output back into ppg%uVpsibox's original ilma indexing (via
        ! gemm_l2g) -- padding projector rows (p>nproj_atom(ia)) are
        ! skipped, never written (their GEMM output is meaningless: a
        ! zero-projector-row dot product, harmless but not a real ilma).
!$acc parallel loop collapse(3) present(ppg,gemm_out_packed,gemm_rinv_packed,gemm_l2g,gemm_nproj_atom)
        do gemm_ia = 1, gemm_natom
        do gemm_io_local = 1, gemm_this_block
        do gemm_p = 1, gemm_max_nproj
          if (gemm_p <= gemm_nproj_atom(gemm_ia)) then
            ppg%uVpsibox(gemm_l2g(gemm_p,gemm_ia), ispin, gemm_io_blk_s+gemm_io_local-1, ik, im) = &
                gemm_out_packed(gemm_p, gemm_io_local, gemm_ia) * gemm_rinv_packed(gemm_p, gemm_ia)
          end if
        end do
        end do
        end do

        gemm_io_blk_s = gemm_io_blk_s + gemm_this_block
      end do
    end do
    end do
    end do

    ! TEMPORARY DIAGNOSTIC -- remove before merging. Compare the GEMM's
    ! actual uVpsibox(1,...) / uVpsibox(2,...) (atom 1's first two
    ! projectors, first orbital/spin/k/im) against a plain host-side
    ! reference computed directly from tpsi/zekr_uV/jxyz, mirroring the
    ! original reduction algorithm exactly.
    !
    ! Fires on call #GEMM_DEBUG_CALL (not call #1): every earlier check at
    ! call 1 matched exactly (including a full-Nlma-coverage sum), yet the
    ! electron-count drift only becomes visible by RT step 10 (~40-60 hpsi
    ! calls in) -- if this same check now fails at a later call, that
    ! points to state corruption building up across repeated calls rather
    ! than a one-shot logic bug in the per-call math itself.
    gemm_call_count = gemm_call_count + 1
    if (gemm_call_count == 50 .and. .not. gemm_debug_printed) then
      gemm_debug_printed = .true.
      uVpsi = 0.d0
      do j=1,ppg%mps(1)
        ix = ppg%jxyz(1,j,1); iy = ppg%jxyz(2,j,1); iz = ppg%jxyz(3,j,1)
        uVpsi = uVpsi + conjg(ppg%zekr_uV(j,1,ik_s)) * tpsi%zwf(ix,iy,iz,1,io_s,ik_s,im_s)
      end do
      uVpsi = uVpsi * ppg%rinv_uvu(1)
      write(*,'(A,2ES16.8,A,2ES16.8)') 'GEMM_DEBUG uVpsibox(ilma=1) gemm=',ppg%uVpsibox(1,1,io_s,ik_s,im_s), &
          ' ref=',uVpsi
      uVpsi = 0.d0
      do j=1,ppg%mps(1)
        ix = ppg%jxyz(1,j,1); iy = ppg%jxyz(2,j,1); iz = ppg%jxyz(3,j,1)
        uVpsi = uVpsi + conjg(ppg%zekr_uV(j,2,ik_s)) * tpsi%zwf(ix,iy,iz,1,io_s,ik_s,im_s)
      end do
      uVpsi = uVpsi * ppg%rinv_uvu(2)
      write(*,'(A,2ES16.8,A,2ES16.8)') 'GEMM_DEBUG uVpsibox(ilma=2) gemm=',ppg%uVpsibox(2,1,io_s,ik_s,im_s), &
          ' ref=',uVpsi

      ! Atom 1 has nproj_atom(1)=max_nproj (no padding) -- find an atom
      ! that actually needs padding (nproj_atom(ia) < max_nproj) and check
      ! both its first and its LAST valid projector, since that's the
      ! untested case (real SiO2 stoichiometry: Si atoms use all 18 slots,
      ! O atoms only use 13, so the majority of atoms are padded and were
      ! never actually exercised by the atom-1-only check above).
      gemm_ia = 0
      do gemm_j = 1, gemm_natom
        if (gemm_nproj_atom(gemm_j) < gemm_max_nproj) then
          gemm_ia = gemm_j
          exit
        end if
      end do
      write(*,'(A,I0,A,I0,A,I0)') 'GEMM_DEBUG padded_atom ia=',gemm_ia,' nproj=',gemm_nproj_atom(max(gemm_ia,1)), &
          ' mps=',ppg%mps(max(gemm_ia,1))
      if (gemm_ia > 0) then
        gemm_ilma = gemm_l2g(1, gemm_ia)
        uVpsi = 0.d0
        do j=1,ppg%mps(gemm_ia)
          ix = ppg%jxyz(1,j,gemm_ia); iy = ppg%jxyz(2,j,gemm_ia); iz = ppg%jxyz(3,j,gemm_ia)
          uVpsi = uVpsi + conjg(ppg%zekr_uV(j,gemm_ilma,ik_s)) * tpsi%zwf(ix,iy,iz,1,io_s,ik_s,im_s)
        end do
        uVpsi = uVpsi * ppg%rinv_uvu(gemm_ilma)
        write(*,'(A,I0,A,2ES16.8,A,2ES16.8)') 'GEMM_DEBUG padded_first ilma=',gemm_ilma,' gemm=', &
            ppg%uVpsibox(gemm_ilma,1,io_s,ik_s,im_s),' ref=',uVpsi

        gemm_ilma = gemm_l2g(gemm_nproj_atom(gemm_ia), gemm_ia)
        uVpsi = 0.d0
        do j=1,ppg%mps(gemm_ia)
          ix = ppg%jxyz(1,j,gemm_ia); iy = ppg%jxyz(2,j,gemm_ia); iz = ppg%jxyz(3,j,gemm_ia)
          uVpsi = uVpsi + conjg(ppg%zekr_uV(j,gemm_ilma,ik_s)) * tpsi%zwf(ix,iy,iz,1,io_s,ik_s,im_s)
        end do
        uVpsi = uVpsi * ppg%rinv_uvu(gemm_ilma)
        write(*,'(A,I0,A,2ES16.8,A,2ES16.8)') 'GEMM_DEBUG padded_last ilma=',gemm_ilma,' gemm=', &
            ppg%uVpsibox(gemm_ilma,1,io_s,ik_s,im_s),' ref=',uVpsi
      end if

      ! Everything checked so far is band io_s (block 1 of the do-while
      ! band-blocking loop). Check a band from the SECOND block
      ! (io_s+gemm_io_block) too, atom 1/ilma=1 -- this isolates whether
      ! the bug is specific to crossing a block boundary (e.g. an
      ! absolute-vs-relative band index mixup) rather than atom padding.
      if (io_e >= io_s + gemm_io_block) then
        gemm_j = io_s + gemm_io_block
        uVpsi = 0.d0
        do j=1,ppg%mps(1)
          ix = ppg%jxyz(1,j,1); iy = ppg%jxyz(2,j,1); iz = ppg%jxyz(3,j,1)
          uVpsi = uVpsi + conjg(ppg%zekr_uV(j,1,ik_s)) * tpsi%zwf(ix,iy,iz,1,gemm_j,ik_s,im_s)
        end do
        uVpsi = uVpsi * ppg%rinv_uvu(1)
        write(*,'(A,I0,A,2ES16.8,A,2ES16.8)') 'GEMM_DEBUG block2 io=',gemm_j,' gemm=', &
            ppg%uVpsibox(1,1,gemm_j,ik_s,im_s),' ref=',uVpsi
      end if

      ! Everything so far falls in a FULL-width band block (io_s is block 1,
      ! io_s+gemm_io_block is block 2, both size gemm_io_block). If
      ! (io_e-io_s+1) isn't a multiple of gemm_io_block, the LAST block is
      ! partial (this_block < gemm_io_block) -- untested until now. Check
      ! the very last band (io_e), atom 1/ilma=1.
      write(*,'(A,I0,A,I0)') 'GEMM_DEBUG nstate=',io_e-io_s+1,' io_block=',gemm_io_block
      uVpsi = 0.d0
      do j=1,ppg%mps(1)
        ix = ppg%jxyz(1,j,1); iy = ppg%jxyz(2,j,1); iz = ppg%jxyz(3,j,1)
        uVpsi = uVpsi + conjg(ppg%zekr_uV(j,1,ik_s)) * tpsi%zwf(ix,iy,iz,1,io_e,ik_s,im_s)
      end do
      uVpsi = uVpsi * ppg%rinv_uvu(1)
      write(*,'(A,I0,A,2ES16.8,A,2ES16.8)') 'GEMM_DEBUG lastblock io=',io_e,' gemm=', &
          ppg%uVpsibox(1,1,io_e,ik_s,im_s),' ref=',uVpsi

      ! Full-coverage cross-check: every individual spot-check so far
      ! (unpadded/padded atoms, band 1/second-block/last-partial-block) has
      ! matched exactly, yet the bug still reproduces -- sum ALL Nlma
      ! projectors' uVpsibox (band io_s) via GEMM vs the ORIGINAL per-ilma
      ! reduction algorithm, covering the entire ilma range at once instead
      ! of a handful of examples.
      block
        complex(8) :: gemm_total, ref_total, ref_val
        integer :: gemm_check_ilma, gemm_check_ia
        gemm_total = (0.d0,0.d0)
        ref_total = (0.d0,0.d0)
        do gemm_check_ilma = 1, Nlma
          gemm_check_ia = ppg%ia_tbl(gemm_check_ilma)
          ref_val = (0.d0,0.d0)
          do j=1,ppg%mps(gemm_check_ia)
            ix = ppg%jxyz(1,j,gemm_check_ia); iy = ppg%jxyz(2,j,gemm_check_ia); iz = ppg%jxyz(3,j,gemm_check_ia)
            ref_val = ref_val + conjg(ppg%zekr_uV(j,gemm_check_ilma,ik_s)) * tpsi%zwf(ix,iy,iz,1,io_s,ik_s,im_s)
          end do
          ref_val = ref_val * ppg%rinv_uvu(gemm_check_ilma)
          gemm_total = gemm_total + ppg%uVpsibox(gemm_check_ilma,1,io_s,ik_s,im_s)
          ref_total = ref_total + ref_val
        end do
        write(*,'(A,2ES16.8,A,2ES16.8)') 'GEMM_DEBUG fullcheck io=io_s gemm_total=',gemm_total,' ref_total=',ref_total
      end block

      ! gemm_wf_packed's (j>mps(ia)) padding rows are only ever zeroed ONCE
      ! at setup and never rewritten by any gather step -- the whole design
      ! assumes they stay zero for the run's lifetime. If something (most
      ! plausibly cublas's own internal workspace churn reusing that
      ! memory) corrupts them over repeated calls, that would explain
      ! exactly this pattern: correct at call 1, wrong by call 50. Check a
      ! definitely-padding cell for atom 1 directly (mps(1)=960 < nps=968).
      write(*,'(A,2ES16.8)') 'GEMM_DEBUG padding_check gemm_wf_packed(965,1,1) (should be 0,0) =', &
          gemm_wf_packed(965,1,1)

      ! gemm_zekr_packed/gemm_rinv_packed are also built ONCE at setup and
      ! never rewritten -- re-verify they still match their source arrays
      ! at call 50 (same check the setup diagnostic did at call 1), to
      ! rule out these "constant" arrays drifting/getting corrupted too.
      write(*,'(A,2ES16.8,A,2ES16.8)') 'GEMM_DEBUG zekr_recheck gemm=',gemm_zekr_packed(1,1,1), &
          ' src=',ppg%zekr_uV(1,gemm_l2g(1,1),ik_s)
      write(*,'(A,ES16.8,A,ES16.8)') 'GEMM_DEBUG rinv_recheck gemm=',gemm_rinv_packed(1,1), &
          ' src=',ppg%rinv_uvu(gemm_l2g(1,1))
    end if

    ! Phase 2 (back-projection): completely unchanged from the plain-acc
    ! path below -- same atomic-free inverse-gather-map algorithm, own
    ! acc kernels region since phase 1 no longer shares one with it here.
!$acc kernels present(ppg,tpsi,htpsi)
!$acc loop collapse(5) independent private(ilocal,ilma,ia,uVpsi,vi,my_nlma,k,j,ix,iy,iz,wrk)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin
      do vi=0,ppg%max_vi-1
        my_nlma = ppg%v2nlma(vi)
        if (my_nlma < 1) cycle

        wrk = 0d0
!$acc loop seq
        do k=1,my_nlma
          ilma = ppg%k2ilma(vi,k)
          j    = ppg%k2j(vi,k)
          wrk  = wrk + ppg%uVpsibox(ilma,ispin,io,ik,im) * ppg%zekr_uV(j,ilma,ik)
        end do

        ix = ppg%v2j(1,vi)
        iy = ppg%v2j(2,vi)
        iz = ppg%v2j(3,vi)
        htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) + wrk
      end do
    end do
    end do
    end do
    end do
!$acc end kernels
#else
!$acc kernels present(ppg,tpsi,htpsi)
!$acc loop collapse(5) independent gang private(ilocal,ilma,ia,uVpsi,vi,my_nlma,k,j,ix,iy,iz,wrk)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin
      do ilma=1,Nlma
        ia = ppg%ia_tbl(ilma)
        uVpsi = 0.d0
!$acc loop vector(32) reduction(+:uVpsi)
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          uVpsi = uVpsi + conjg(ppg%zekr_uV(j,ilma,ik)) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
        end do
        ppg%uVpsibox(ilma,ispin,io,ik,im) = uVpsi * ppg%rinv_uvu(ilma)
      end do
#ifdef USE_OPENACC
    end do
    end do
    end do
    end do
!$acc loop collapse(5) independent private(ilocal,ilma,ia,uVpsi,vi,my_nlma,k,j,ix,iy,iz,wrk)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin
#endif
      do vi=0,ppg%max_vi-1
        my_nlma = ppg%v2nlma(vi)
        if (my_nlma < 1) cycle

        wrk = 0d0
!$acc loop seq
        do k=1,my_nlma
          ilma = ppg%k2ilma(vi,k)
          j    = ppg%k2j(vi,k)
          wrk  = wrk + ppg%uVpsibox(ilma,ispin,io,ik,im) * ppg%zekr_uV(j,ilma,ik)
        end do

        ix = ppg%v2j(1,vi)
        iy = ppg%v2j(2,vi)
        iz = ppg%v2j(3,vi)
        htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) + wrk
      end do
    end do
    end do
    end do
    end do
!$acc end kernels
#endif
#else
!$omp parallel do collapse(4) &
!$omp             private(im,ik,io,ispin,ilma,ia,uVpsi,j,ix,iy,iz,wrk)
    do im=im_s,im_e
    do ik=ik_s,ik_e
    do io=io_s,io_e
    do ispin=1,Nspin

      do ilma=1,Nlma
        ia = ppg%ia_tbl(ilma)
        uVpsi = 0.d0
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          uVpsi = uVpsi + conjg(ppg%zekr_uV(j,ilma,ik)) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
        end do
        uVpsi = uVpsi * ppg%rinv_uvu(ilma)
!OCL norecurrence
        do j=1,ppg%mps(ia)
          ix = ppg%jxyz(1,j,ia)
          iy = ppg%jxyz(2,j,ia)
          iz = ppg%jxyz(3,j,ia)
          wrk = uVpsi * ppg%zekr_uV(j,ilma,ik)
          htpsi%zwf(ix,iy,iz,ispin,io,ik,im) = htpsi%zwf(ix,iy,iz,ispin,io,ik,im) + wrk
        end do
      end do

    end do
    end do
    end do
    end do
!$omp end parallel do
#endif

  end if

  call nvtxEndRange()
  call timer_end(LOG_UHPSI_PSEUDO)

  return
end subroutine zpseudo

subroutine calc_uVpsi(nspin,info,ppg,tpsi,uVpsibox)
  use structures
  use timer
  implicit none
  integer        ,intent(in) :: nspin
  type(s_parallel_info),intent(in) :: info
  type(s_pp_grid),intent(in) :: ppg
  type(s_orbital),intent(in) :: tpsi
  complex(8)    ,allocatable :: uVpsibox (:,:,:,:,:)
  integer :: ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e,norb
  integer :: ilma,ia,j,ix,iy,iz,Nlma
  complex(8) :: uVpsi
  integer :: ilocal

  call timer_begin(LOG_UHPSI_PSEUDO)

  im_s = info%im_s
  im_e = info%im_e
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e
  norb = Nspin* info%numo * info%numk * info%numm

  Nlma = ppg%Nlma

  allocate(uVpsibox(Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e,Nlma))

#ifdef USE_OPENACC
!$acc parallel loop collapse(4) private(im,ik,io,ispin,ilocal,ilma,ia,uVpsi,j,ix,iy,iz)
#else
!$omp parallel do collapse(4) &
!$omp             private(im,ik,io,ispin,ilocal,ilma,ia,uVpsi,j,ix,iy,iz)
#endif
  do im=im_s,im_e
  do ik=ik_s,ik_e
  do io=io_s,io_e
  do ispin=1,Nspin

    do ilocal=1,ppg%ilocal_nlma
      ilma=ppg%ilocal_nlma2ilma(ilocal)
      ia  =ppg%ilocal_nlma2ia  (ilocal)
      uVpsi = 0.d0
      do j=1,ppg%mps(ia)
        ix = ppg%jxyz(1,j,ia)
        iy = ppg%jxyz(2,j,ia)
        iz = ppg%jxyz(3,j,ia)
        uVpsi = uVpsi + conjg(ppg%zekr_uV(j,ilma,ik)) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
      end do
      uVpsi = uVpsi * ppg%rinv_uvu(ilma)
      uVpsibox(ispin,io,ik,im,ilma) = uVpsi
    end do

  end do
  end do
  end do
  end do
#ifdef USE_OPENACC
!$acc end parallel
#else
!$omp end parallel do
#endif

  call timer_end(LOG_UHPSI_PSEUDO)

  return
end subroutine calc_uVpsi

subroutine calc_uVpsi_rdivided(nspin,info,ppg,tpsi,uVpsibox,uVpsibox2)
  use structures
  use timer
#ifdef FORTRAN_COMPILER_HAS_MPI_VERSION3
  use salmon_global, only: natom
  use communication, only: comm_wait_all,comm_show_error
  use mpi, only: MPI_SUM,MPI_DOUBLE_COMPLEX
#else
  use communication, only: comm_summation
#endif
  implicit none
  integer        ,intent(in) :: nspin
  type(s_parallel_info),intent(in) :: info
  type(s_pp_grid),intent(in) :: ppg
  type(s_orbital),intent(in) :: tpsi
  complex(8)    ,allocatable :: uVpsibox (:,:,:,:,:)
  complex(8)    ,allocatable :: uVpsibox2(:,:,:,:,:)
#ifdef USE_OPENACC
  complex(8)    ,allocatable, save, device :: dev_uVpsibox (:,:,:,:,:)
  complex(8)    ,allocatable, save, device :: dev_uVpsibox2(:,:,:,:,:)
#endif
  integer :: ispin,io,ik,im,im_s,im_e,ik_s,ik_e,io_s,io_e,norb
  integer :: ilma,ia,j,ix,iy,iz,Nlma,ilocal
  complex(8) :: uVpsi
#ifdef FORTRAN_COMPILER_HAS_MPI_VERSION3
  integer :: ireqs(natom),nreq,ierr,is,ie,ns
#endif

  call nvtxStartRange('calc_uVpsi_rdivided', __LINE__)
  call timer_begin(LOG_UHPSI_PSEUDO)

  im_s = info%im_s
  im_e = info%im_e
  ik_s = info%ik_s
  ik_e = info%ik_e
  io_s = info%io_s
  io_e = info%io_e
  norb = Nspin* info%numo * info%numk * info%numm

  Nlma = ppg%Nlma

  allocate(uVpsibox (Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e,Nlma))
  allocate(uVpsibox2(Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e,Nlma))

#ifdef USE_OPENACC
  if (.not. allocated(dev_uVpsibox)) then
    allocate(dev_uVpsibox (Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e,Nlma))
  end if

  if (.not. allocated(dev_uVpsibox2)) then
    allocate(dev_uVpsibox2(Nspin,io_s:io_e,ik_s:ik_e,im_s:im_e,Nlma))
  end if

!$acc kernels
  dev_uVpsibox = 0d0
  dev_uVpsibox2 = 0d0
!$acc end kernels
#else
#ifdef FORTRAN_COMPILER_HAS_MPI_VERSION3
  uVpsibox2 = 0d0
#else
  uVpsibox  = 0d0
#endif
#endif

#ifdef USE_OPENACC
!$acc kernels
#else
!$omp parallel do collapse(4) &
!$omp             private(im,ik,io,ispin,ilocal,ilma,ia,uVpsi,j,ix,iy,iz)
#endif
  do im=im_s,im_e
  do ik=ik_s,ik_e
  do io=io_s,io_e
  do ispin=1,Nspin

    do ilocal=1,ppg%ilocal_nlma
      ilma=ppg%ilocal_nlma2ilma(ilocal)
      ia  =ppg%ilocal_nlma2ia  (ilocal)
      uVpsi = 0.d0
      do j=1,ppg%mps(ia)
        ix = ppg%jxyz(1,j,ia)
        iy = ppg%jxyz(2,j,ia)
        iz = ppg%jxyz(3,j,ia)
        uVpsi = uVpsi + conjg(ppg%zekr_uV(j,ilma,ik)) * tpsi%zwf(ix,iy,iz,ispin,io,ik,im)
      end do
      uVpsi = uVpsi * ppg%rinv_uvu(ilma)
#ifdef USE_OPENACC
      dev_uVpsibox(ispin,io,ik,im,ilma) = uVpsi
#else
      uVpsibox(ispin,io,ik,im,ilma) = uVpsi
#endif
    end do

  end do
  end do
  end do
  end do
#ifdef USE_OPENACC
!$acc end kernels
#else
!$omp end parallel do
#endif

  call timer_end(LOG_UHPSI_PSEUDO)

  call timer_begin(LOG_UHPSI_PSEUDO_COMM)
#ifdef FORTRAN_COMPILER_HAS_MPI_VERSION3
! FIXME: This subroutine uses MPI functions directly...
  nreq = 0
  do ia=1,natom
    is = ppg%irange_atom(1,ia)
    ie = ppg%irange_atom(2,ia)
    ns = ie - is + 1

    if (ppg%ireferred_atom(ia)) then
      nreq = nreq + 1
#ifdef USE_OPENACC
      call MPI_Iallreduce( dev_uVpsibox (1,io_s,ik_s,im_s,is) &
                         , dev_uVpsibox2(1,io_s,ik_s,im_s,is) &
                         , ns*norb, MPI_DOUBLE_COMPLEX, MPI_SUM, ppg%icomm_atom(ia) &
                         , ireqs(nreq), ierr )
#else
      call MPI_Iallreduce( uVpsibox (1,io_s,ik_s,im_s,is) &
                         , uVpsibox2(1,io_s,ik_s,im_s,is) &
                         , ns*norb, MPI_DOUBLE_COMPLEX, MPI_SUM, ppg%icomm_atom(ia) &
                         , ireqs(nreq), ierr )
#endif
      call comm_show_error(ierr)
    !else
      ! uvpsibox2(:,:,:,:,ppg%irange_ia(1:2,ia)) does not use in this process...
      ! We can skip self copy, but zero clear required
    end if
  end do
  call comm_wait_all(ireqs(1:nreq))
#ifdef USE_OPENACC
!$acc kernels
uVpsibox2 = dev_uVpsibox2
uVpsibox = dev_uVpsibox
!$acc end kernels
#endif
#else
  call comm_summation(uVpsibox,uVpsibox2,Nlma*Norb,info%icomm_r)
#endif
  call timer_end(LOG_UHPSI_PSEUDO_COMM)

  call nvtxEndRange()
  return
end subroutine calc_uVpsi_rdivided

end module nonlocal_potential
