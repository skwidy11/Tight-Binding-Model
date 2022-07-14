subroutine update_band_AF(kx,ky,matrix,t,eps,nocc,site,U,numstates)
    complex, dimension(12,12), intent(out) :: matrix
    double precision,intent(in) :: kx,ky,U
    double precision, dimension(4,4,9,9) :: t
    double precision, dimension(4,9) :: eps
    double precision, dimension(numstates,2) :: nocc
    integer :: site,i,j,k,n
    double precision :: tpd,tpp
    tpd = Abs(t(1,2,1,1))
    tpp = Abs(t(2,2,1,1))
    matrix = 0.0d0

    ! H_A matrix set
    do i=0,3 ! Iterate over both CuO2 units and spins
      if(i.eq.0.or.i.eq.3) then
        matrix(1+i*3,1+i*3) = U*(0.5d0-nocc(site,2))
      else
        matrix(1+i*3,1+i*3) = U*(0.5d0-nocc(site,1))
      end if
      matrix(2+i*3,2+i*3) = eps(2,1)
      matrix(3+i*3,3+i*3) = eps(2,1)
      matrix(1+i*3,2+i*3) = tpd*exp((0d0,1d0)*kx) ! Copper to oxygen_x
      matrix(2+i*3,1+i*3) = conjg(matrix(1+i*3,2+i*3))
      matrix(1+i*3,3+i*3) = tpd*exp((0d0,1d0)*ky) ! Copper to oxygen_y
      matrix(3+i*3,1+i*3) = conjg(matrix(1+i*3,3+i*3))
      matrix(2+i*3,3+i*3) = 2d0*tpp*dcos(kx-ky) ! Oxygen_x to Oxygen_y
      matrix(3+i*3,2+i*3) = conjg(matrix(2+i*3,3+i*3))
    enddo

    ! H_M matrix set
    do j=0,1 ! Handles spin
    do i=0,1 ! Handles which off-diagonal we filling out.
      k = 3*(1-i)+j*6
      n = i*3+j*6
      matrix(1+k,2+n) = -tpd*exp((0d0,-1d0)*kx) ! Copper to oxygen_x
      matrix(2+n,1+k) = conjg(matrix(1+k,2+n))
      matrix(1+k,3+n) = -tpd*exp((0,-1d0)*ky) ! Copper to oxygen_y
      matrix(3+n,1+k) = conjg(matrix(1+k,3+n))
      matrix(2+k,3+n) = -2d0*tpp*dcos(kx+ky) ! Oxygen to Oxygen
      matrix(3+n,2+k) = conjg(matrix(2+k,3+n))
    end do
    end do


end subroutine
