program dos_get
  implicit none
  double precision, allocatable :: ampl(:), l_ampl_c(:), l_ampl_o(:)
  integer :: KT, ND, k, it, n, size
  double precision :: time_max, mag, ampr, ampi, Emax, Emin, time_interval
  double precision :: vol, time, pi, t, ampr2, ampi2, fac, tact, w, eps, p, ECM
  open(unit=90,file='evolvew2.dat')
  open(unit=80,file='local_evolvew2.dat')
  pi=4.0d0*datan(1.0d0)
  ND = 1728
  read(90,*) KT, time_max, Emax, Emin, vol, p
  time_interval = time_max/dble(KT)
  KT = KT ! must be power of two
  size = 2*KT ! for zero padding
  ECM = (Emax + Emin)/2.0d0
  print *, "p = ",p
  print *, "KT =",KT

  time_interval=time_interval*4.135667d0/((Emax-Emin)*pi)
!  time_interval = 2*pi/(Emax-Emin)
!  time_max = time_interval*KT

  allocate(ampl(2*KT))
  allocate(l_ampl_c(2*KT))
  allocate(l_ampl_o(2*KT))

! Read in amplitudes(And analytically continue)
  eps=0.0035d0
  do it = 1,2*KT,2
    read(90,*) t, ampl(it), ampl(it+1), mag
    read(80,*) t, l_ampl_c(it), l_ampl_c(it+1), l_ampl_o(it), l_ampl_o(it+1)
    ampl(it) = ampl(it)*exp(-t*eps)
    ampl(it+1) = ampl(it+1)*exp(-t*eps)
    l_ampl_c(it) = l_ampl_c(it)*exp(-t*eps)
    l_ampl_c(it+1) = l_ampl_c(it+1)*exp(-t*eps)
    l_ampl_o(it) = l_ampl_o(it)*exp(-t*eps)
    l_ampl_o(it+1) = l_ampl_o(it+1)*exp(-t*eps)
  enddo

! zero padding
!  do n = size,KT-1 ! need to start size on odd
!    if(n == size) print *, "Doing zero padding"
!    ampl(2*n+1)=0
!    ampl(2*n+2)=0
!  end do

! Fourier Transform
  call four1(ampl,KT,1)
  call four1(l_ampl_c,KT,1)
  call four1(l_ampl_o,KT,1)

! Multiply by time interval
  do k=1,2*KT
    ampl(k)=ampl(k)*time_interval
    l_ampl_c(k)=l_ampl_c(k)*time_interval
    l_ampl_o(k)=l_ampl_o(k)*time_interval
  end do

  open(unit=50,file='DOS.dat')
  open(unit=60,file='l_DOS.dat')
! Write DOS's to file
  do n = KT/2, KT-1 ! do negative frequencies first
!    w = -dble(n-KT)*pi*(Emax-Emin)/(KT*time_interval)
    w = dble(n-KT)*4.135667d0/(KT*time_interval) + ECM
    mag = dsqrt(ampl(2*n+1)**2+ampl(2*n+2)**2)
    write(50,*) w,ampl(2*n+1),ampl(2*n+2), mag
    write(60,*) w,l_ampl_c(2*n+1),l_ampl_c(2*n+2),l_ampl_o(2*n+1),l_ampl_o(2*n+2)
  end do

  do n = 0, KT/2 ! positive frequencies
!    w = -dble(n)*pi*(Emax-Emin)/(dble(KT)*time_interval)
    w = dble(n)*4.135667d0/(dble(KT)*time_interval) + ECM
    mag = dsqrt(ampl(2*n+1)**2+ampl(2*n+2)**2)
    write(50,*) w,ampl(2*n+1),ampl(2*n+2), mag
    write(60,*) w,l_ampl_c(2*n+1),l_ampl_c(2*n+2),l_ampl_o(2*n+1),l_ampl_o(2*n+2)
  end do

contains

! Fast Fourier Transform
! data = array length 2*nn; alternates between real and imaginary parts
! nn = number of data points (only powers of two)
! isign = +1, fourier transform; -1, inverse fourier transform
subroutine four1(data,nn,isign)
  implicit none
  integer :: isign,nn
  double precision,dimension(2*nn) :: data
  integer :: i,istep,j,m,mmax,n
  double precision :: tempi,tempr
  double precision :: theta,wi,wpi,wpr,wr,wtemp
  n=2*nn
  j=1
! bit reversal
  do i=1,n,2
    if(j .gt. i) then
      tempr=data(j)
      tempi=data(j+1)
      data(j)=data(i)
      data(j+1)=data(i+1)
      data(i)=tempr
      data(i+1)=tempi
    end if
    m=nn
1   if ((m .ge. 2) .and. (j .gt. m)) then
      j=j-m
      m=m/2
      goto 1
    end if
    j = j+m
  end do
! Danielson-Lanczos Lemma
  mmax=2
2  if (n .gt. mmax) then
    istep = 2*mmax
    theta=6.28318530717959d0/dble(isign*mmax)
    wpr=-2.d0*dsin(0.5d0*theta)**2
    wpi=dsin(theta)
    wr=1.d0
    wi=0.d0
    do m=1,mmax,2
      do i=m,n,istep
        j=i+mmax
        tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
        tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
        data(j)=data(i)-tempr
        data(j+1)=data(i+1)-tempi
        data(i)=data(i)+tempr
        data(i+1)=data(i+1)+tempi
      end do
      wtemp=wr
      wr=wr*wpr-wi*wpi+wr
      wi=wi*wpr+wtemp*wpi+wi
    end do
    mmax=istep
  goto 2
  end if
  return
end subroutine four1

! Factor to obtain results in Ohm**-1 cm**1
!  fac=(2.0d0*pi*1.602d0)**2/(vol*4.135667**2)
! Take into account time in fs for integration, factor of 2 for
! positive/negative time integration
!  fac=2.0d0*fac*tact
! Account for normalization of Hamiltonian for velocity operators
! and the fact that the delta function defined in terms of the normalized
! energy and Hamiltonian
! Final result to be integrated units 10**-15 Ohm**-1 cm**-1
! fac=fac*((Emax-Emin)/2.0d0)

! Diagonalize Hamiltonian
!  hmatrix = 0.0d0
!  do n = 1,ND
!    do j = 1,neighb(n,0)
!      hmatrix(n,neighb(n,j)) = hamilt(n,j,1,1)
!     hmatrix(neighb(n,j),n) = hmatrix(n,neighb(n,j))
!    end do
!  end do
!  print *, "before diag call"
!  l = ND*(3+ND/2)
!  allocate(work(l))
!  call dsyev('N','U',ND,hmatrix,ND,eig,work,l,info)
!  print *,"Info=",info
!  open(unit=14,file='eigen.dat')
!  do n=1,ND
!    write(14,*) n, eig(n)
!  end do
!  close(14)
!  print *, "after diag call"
!
! * - save the Kernel coefficients
!
!  open(unit=17,file='g.d',status='unknown')
!  do n = 0, NT-1
!     write(17,"(i4,1x,e10.3)")n,g(n)
!  end do
!  close(unit=17)
!
! Calculate Correction Terms
!    hub_cor = 0.0d0
!    dc_cor = 0.0d0
!    do j=1,ND
!      if (ntype(j).eq.1) then
!        hub_cor = hub_cor + U * nocc(j,1) * nocc(j,2)
!        dc_cor = dc_cor + U/2d0*(nocc(j,1)+nocc(j,2))*(nocc(j,1)+nocc(j,2) - 1)
!      end if
!    end do
!    print *, "DC Correction: ", dc_cor
!    print *, "Hubbard Correction: ", hub_cor

end program dos_get
