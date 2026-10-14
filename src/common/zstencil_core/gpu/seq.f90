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

#if defined(USE_OPENACC) && defined(USE_CUDA)
module zstencil_cuda
  use cudafor
  implicit none

contains

  ! One thread per real(8) component instead of per complex element, so
  ! consecutive threads read consecutive 8-byte words: 8 sectors per request
  ! rather than the 16 a 16-byte-strided scalar load costs. lapt/nabt are
  ! real, so the neighbour accumulation is componentwise and only the -zI*w
  ! term crosses components -- one shfl_xor between the paired lanes.
  ! Expression order matches the OpenACC version, so results are bit-exact.
  attributes(global) subroutine zstencil_x_cuda(tr, hr, vf, ix_map, &
      nwork, nx, ny, nz, igs1, igs2, igs3, ixlb, &
      lb1, lb2, lb3, n1, n12, n123, n123s, &
      vlb1, vlb2, vlb3, m1, m12, ispin, lap0, &
      lp1, lp2, lp3, lp4, nb1, nb2, nb3, nb4)
    implicit none
    real(8) :: tr(*), hr(*), vf(*)
    integer :: ix_map(*)
    integer, value :: nwork, nx, ny, nz, igs1, igs2, igs3, ixlb
    integer, value :: lb1, lb2, lb3, n1, n12, n123, n123s
    integer, value :: vlb1, vlb2, vlb3, m1, m12, ispin
    real(8), value :: lap0, lp1, lp2, lp3, lp4, nb1, nb2, nb3, nb4
    integer :: i, c, q, ix, iy, iz, iol, obase, cen_w, jp, jm
    real(8) :: vc, wc, wp, t0, t1, cen, vlocal, outv
    logical :: live

    i    = (blockidx%x - 1) * blockdim%x + threadidx%x
    live = (i <= nwork)
    vc   = 0d0
    wc   = 0d0
    obase = 0
    cen_w = 1

    if (live) then
      c   = mod(i - 1, 2)
      q   = (i - 1) / 2
      ix  = igs1 + mod(q, nx)
      iy  = igs2 + mod(q / nx, ny)
      iz  = igs3 + mod(q / (nx * ny), nz)
      iol = q / (nx * ny * nz)

      ! complex-element offset of (*, iy, iz, ispin, io); ix added per access
      obase = (iy - lb2) * n1 + (iz - lb3) * n12 &
            + (ispin - 1) * n123 + iol * n123s

      jp = ix_map(ix + 1 - ixlb + 1); jm = ix_map(ix - 1 - ixlb + 1)
      t0 = tr(2 * ((jp - lb1) + obase) + c + 1)
      t1 = tr(2 * ((jm - lb1) + obase) + c + 1)
      vc = lp1 * (t0 + t1) + vc
      wc = nb1 * (t0 - t1) + wc

      jp = ix_map(ix + 2 - ixlb + 1); jm = ix_map(ix - 2 - ixlb + 1)
      t0 = tr(2 * ((jp - lb1) + obase) + c + 1)
      t1 = tr(2 * ((jm - lb1) + obase) + c + 1)
      vc = lp2 * (t0 + t1) + vc
      wc = nb2 * (t0 - t1) + wc

      jp = ix_map(ix + 3 - ixlb + 1); jm = ix_map(ix - 3 - ixlb + 1)
      t0 = tr(2 * ((jp - lb1) + obase) + c + 1)
      t1 = tr(2 * ((jm - lb1) + obase) + c + 1)
      vc = lp3 * (t0 + t1) + vc
      wc = nb3 * (t0 - t1) + wc

      jp = ix_map(ix + 4 - ixlb + 1); jm = ix_map(ix - 4 - ixlb + 1)
      t0 = tr(2 * ((jp - lb1) + obase) + c + 1)
      t1 = tr(2 * ((jm - lb1) + obase) + c + 1)
      vc = lp4 * (t0 + t1) + vc
      wc = nb4 * (t0 - t1) + wc

      cen_w  = 2 * ((ix - lb1) + obase) + c + 1
      cen    = tr(cen_w)
      vlocal = vf((ix - vlb1) + (iy - vlb2) * m1 + (iz - vlb3) * m12 + 1)
    end if

    ! every lane must reach the shuffle, so this sits outside the guard
    wp = __shfl_xor(wc, 1)

    if (live) then
      if (c == 0) then
        outv = vlocal * cen + lap0 * cen - 0.5d0 * vc + wp
      else
        outv = vlocal * cen + lap0 * cen - 0.5d0 * vc - wp
      end if
      hr(cen_w) = outv
    end if
  end subroutine zstencil_x_cuda

end module zstencil_cuda
#endif

subroutine zstencil_typical_gpu(io_s,io_e,Nspin,is_array,ie_array,is,ie,idx,idy,idz,igs,ige &
                               ,tpsi,htpsi,V_local,lap0,lapt,nabt &
                               )
  use structures
#ifdef USE_OPENACC
  use cudafor, only: cudaGetDevice, cudaGetDeviceProperties, cudaDeviceProp
#endif
#if defined(USE_OPENACC) && defined(USE_CUDA)
  use zstencil_cuda, only: zstencil_x_cuda
  use openacc, only: acc_get_cuda_stream, acc_async_sync
  use iso_c_binding, only: c_loc, c_f_pointer
#endif
  implicit none


  integer,intent(in) :: io_s, io_e, Nspin

  integer,intent(in) :: is_array(3),ie_array(3),is(3),ie(3)
  integer,intent(in),target :: idx(is(1)-4:ie(1)+4)
  integer,intent(in) :: idy(is(2)-4:ie(2)+4),idz(is(3)-4:ie(3)+4)
  integer,intent(in) :: igs(3),ige(3)

  complex(8),intent(in) ,target :: tpsi   (is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3),1:Nspin,io_s:io_e)
  complex(8),intent(out),target :: htpsi  (is_array(1):ie_array(1),is_array(2):ie_array(2),is_array(3):ie_array(3),1:Nspin,io_s:io_e)
  ! type(s_scalar) ,intent(in) :: V_local(is(1):ie(1),is(2):ie(2),is(3):ie(3), 1:Nspin)
  type(s_scalar) ,intent(in) :: V_local(1:Nspin)
  real(8),   intent(in)  :: lap0
  real(8),   intent(in)  :: lapt(12), nabt(12)

  complex(8), parameter :: zI=(0.d0,1.d0)
  complex(8), parameter :: z0=(0.d0,0.d0)

  integer    :: ix,iy,iz,io,ispin
  complex(8) :: v,w
  !complex(8) :: t(8)
  complex(8) :: t_0,t_1
  integer    :: i

#ifdef USE_OPENACC
  logical,save :: probed = .false.
  logical,save :: io_restructure = .false.
  type(cudaDeviceProp) :: prop
  integer :: cuda_stat, cuda_dev
#endif
#if defined(USE_OPENACC) && defined(USE_CUDA)
  real(8),pointer  :: tpsi_r(:), htpsi_r(:)
  integer,pointer  :: idx_p(:)
  integer :: nw1,nw2,nw3,nwspin,nwio,xg_nwork,xg_grid
  integer :: xn1,xn12,xn123,xn123s,xm1,xm12
#endif

#ifdef __INTEL_COMPILER
#if defined(__KNC__) || defined(__AVX512F__)
#   define MEM_ALIGN   64
#   define VECTOR_SIZE 4
# else
#   define MEM_ALIGN   32
#   define VECTOR_SIZE 2
# endif

! !dir$ assume_aligned V_local:MEM_ALIGN
! !dir$ assume_aligned tpsi   :MEM_ALIGN
! !dir$ assume_aligned htpsi  :MEM_ALIGN
#endif

#define DX(dt) idx(ix+(dt)),iy,iz
#define DY(dt) ix,idy(iy+(dt)),iz
#define DZ(dt) ix,iy,idz(iz+(dt))
#ifdef USE_OPENACC
  ! Pulling io out of the collapse is ~1.4x faster on Blackwell, ~1.6x slower on
  ! Hopper; the binary is multi-architecture, so dispatch at runtime.
  if (.not. probed) then
    cuda_stat = cudaGetDevice(cuda_dev)
    cuda_stat = cudaGetDeviceProperties(prop, cuda_dev)
    io_restructure = (cuda_stat == 0 .and. prop%major >= 10)
    probed = .true.
  end if

  if (io_restructure) then

#if defined(USE_CUDA)
  ! x direction: coalesced component-per-thread kernel (see zstencil_cuda).
  nw1 = ige(1)-igs(1)+1 ; nw2 = ige(2)-igs(2)+1 ; nw3 = ige(3)-igs(3)+1
  nwio = io_e-io_s+1
  xn1    = ie_array(1)-is_array(1)+1
  xn12   = xn1*(ie_array(2)-is_array(2)+1)
  xn123  = xn12*(ie_array(3)-is_array(3)+1)
  xn123s = xn123*Nspin
  xm1  = ie(1)-is(1)+1
  xm12 = xm1*(ie(2)-is(2)+1)
  call c_f_pointer(c_loc(tpsi),  tpsi_r,  [2*xn123s*nwio])
  call c_f_pointer(c_loc(htpsi), htpsi_r, [2*xn123s*nwio])
  call c_f_pointer(c_loc(idx),   idx_p,   [ie(1)-is(1)+9])
  xg_nwork = 2*nw1*nw2*nw3*nwio
  xg_grid  = (xg_nwork + 127)/128
  do ispin=1,Nspin
!$acc host_data use_device(tpsi_r, htpsi_r, idx_p, V_local(ispin)%f)
    call zstencil_x_cuda<<<xg_grid,128,0,acc_get_cuda_stream(acc_async_sync)>>>( &
        tpsi_r, htpsi_r, V_local(ispin)%f, idx_p, &
        xg_nwork, nw1, nw2, nw3, igs(1), igs(2), igs(3), is(1)-4, &
        is_array(1), is_array(2), is_array(3), xn1, xn12, xn123, xn123s, &
        is(1), is(2), is(3), xm1, xm12, ispin, lap0, &
        lapt(1), lapt(2), lapt(3), lapt(4), nabt(1), nabt(2), nabt(3), nabt(4))
!$acc end host_data
  end do

  ! y+z direction: unchanged OpenACC, now its own region
!$acc kernels copyin(V_local, tpsi) copy(htpsi)
!$acc loop collapse(4)
  do ispin=1,Nspin
  do iz=igs(3),ige(3)
  do iy=igs(2),ige(2)
  do ix=igs(1),ige(1)
  !$acc loop seq
#else
!$acc kernels copyin(V_local, tpsi) copy(htpsi)
!$acc loop collapse(4)
  do ispin=1,Nspin
  do iz=igs(3),ige(3)
  do iy=igs(2),ige(2)
  do ix=igs(1),ige(1)
  !$acc loop seq
  do io=io_s,io_e
    v = z0
    w = z0
    !$acc loop seq
    do i = 1,4
      t_0 = tpsi(DX( i), ispin, io)
      t_1 = tpsi(DX(-i), ispin, io)
      v=lapt(i)*(t_0+t_1) + v
      w=nabt(i)*(t_0-t_1) + w
    end do
    htpsi(ix,iy,iz,ispin,io) = V_local(ispin)%f(ix,iy,iz)*tpsi(ix,iy,iz,ispin,io) &
                    + lap0*tpsi(ix,iy,iz,ispin,io) &
                    - 0.5d0 * v - zI * w
  end do

  !$acc loop seq
#endif
  do io=io_s,io_e
    v = z0
    w = z0
    !$acc loop seq
    do i = 1,4
      t_0 = tpsi(DY( i), ispin, io)
      t_1 = tpsi(DY(-i), ispin, io)
      v=lapt(i+4)*(t_0+t_1) + v
      w=nabt(i+4)*(t_0-t_1) + w
    end do
    !$acc loop seq
    do i = 1,4
      t_0 = tpsi(DZ( i), ispin, io)
      t_1 = tpsi(DZ(-i), ispin, io)
      v=lapt(i+8)*(t_0+t_1) + v
      w=nabt(i+8)*(t_0-t_1) + w
    end do
    htpsi(ix,iy,iz,ispin,io) = htpsi(ix,iy,iz,ispin,io) &
                    - 0.5d0 * v - zI * w
  end do
  end do
  end do
  end do
  end do
!$acc end kernels

  else

!$acc kernels copyin(V_local, tpsi) copy(htpsi)
!$acc loop collapse(5)
  do io=io_s,io_e
  do ispin=1,Nspin
  do iz=igs(3),ige(3)
  do iy=igs(2),ige(2)
  do ix=igs(1),ige(1)
    v = z0
    w = z0
    !$acc loop seq
    do i = 1,4
      t_0 = tpsi(DX( i), ispin, io)
      t_1 = tpsi(DX(-i), ispin, io)
      v=lapt(i)*(t_0+t_1) + v
      w=nabt(i)*(t_0-t_1) + w
    end do
    htpsi(ix,iy,iz,ispin,io) = V_local(ispin)%f(ix,iy,iz)*tpsi(ix,iy,iz,ispin,io) &
                    + lap0*tpsi(ix,iy,iz,ispin,io) &
                    - 0.5d0 * v - zI * w
  end do
  end do
  end do
  end do
  end do

!$acc loop collapse(5)
  do io=io_s,io_e
  do ispin=1,Nspin
  do iz=igs(3),ige(3)
  do iy=igs(2),ige(2)
  do ix=igs(1),ige(1)
    v = z0
    w = z0
    !$acc loop seq
    do i = 1,4
      t_0 = tpsi(DY( i), ispin, io)
      t_1 = tpsi(DY(-i), ispin, io)
      v=lapt(i+4)*(t_0+t_1) + v
      w=nabt(i+4)*(t_0-t_1) + w
    end do
    !$acc loop seq
    do i = 1,4
      t_0 = tpsi(DZ( i), ispin, io)
      t_1 = tpsi(DZ(-i), ispin, io)
      v=lapt(i+8)*(t_0+t_1) + v
      w=nabt(i+8)*(t_0-t_1) + w
    end do
    htpsi(ix,iy,iz,ispin,io) = htpsi(ix,iy,iz,ispin,io) &
                    - 0.5d0 * v - zI * w
  end do
  end do
  end do
  end do
  end do
!$acc end kernels

  end if
#else
  do io=io_s,io_e
  do ispin=1,Nspin
  do iz=igs(3),ige(3)
  do iy=igs(2),ige(2)
! !dir$ assume_aligned V_local(is(1),iy,iz):MEM_ALIGN
! !dir$ assume_aligned tpsi(is_array(1),iy,iz)   :MEM_ALIGN
! !dir$ assume_aligned htpsi(is_array(1),iy,iz)  :MEM_ALIGN
  do ix=igs(1),ige(1)
    v = z0
    w = z0
    do i = 1,4
      t_0 = tpsi(DX( i), ispin, io)
      t_1 = tpsi(DX(-i), ispin, io)
      v=lapt(i)*(t_0+t_1) + v
      w=nabt(i)*(t_0-t_1) + w
    end do
    do i = 1,4
      t_0 = tpsi(DY( i), ispin, io)
      t_1 = tpsi(DY(-i), ispin, io)
      v=lapt(i+4)*(t_0+t_1) + v
      w=nabt(i+4)*(t_0-t_1) + w
    end do
    do i = 1,4
      t_0 = tpsi(DZ( i), ispin, io)
      t_1 = tpsi(DZ(-i), ispin, io)
      v=lapt(i+8)*(t_0+t_1) + v
      w=nabt(i+8)*(t_0-t_1) + w
    end do
    htpsi(ix,iy,iz,ispin,io) = V_local(ispin)%f(ix,iy,iz)*tpsi(ix,iy,iz,ispin,io) &
                    + lap0*tpsi(ix,iy,iz,ispin,io) &
                    - 0.5d0 * v - zI * w
  end do
  end do
  end do
  end do
  end do
#endif
end subroutine
