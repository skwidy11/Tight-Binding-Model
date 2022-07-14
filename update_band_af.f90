subroutine update_band_AF(kx,ky,matrix,t,eps,nocc,site,U,numstates)
    complex, dimension(20,20), intent(out) :: matrix
    double precision,intent(in) :: kx,ky,U
    double precision, dimension(4,4,9,9) :: t
    double precision, dimension(4,9) :: eps
    double precision, dimension(numstates,2) :: nocc
    integer :: site,i,j,k,n
    double precision :: tpd,tpp,tpp_p
    tpd = Abs(t(1,2,1,1))
    tpp = Abs(t(2,2,2,1))
    tpp_p_11 = Abs(t(2,2,1,1))
    tpp_p_22 = Abs(t(2,2,2,2))
    matrix = 0.0d0

    ! H_A matrix set
    do i=0,3
      if(i.eq.0.or.i.eq.3) then
        matrix(1+i*5,1+i*5) = U*(0.5d0-nocc(site,2))
      else
        matrix(1+i*5,1+i*5) = U*(0.5d0-nocc(site,1))
      end if
      matrix(2+i*5,2+i*5) = eps(2,1)
      matrix(3+i*5,3+i*5) = eps(2,2)
      matrix(4+i*5,4+i*5) = eps(2,2)
      matrix(5+i*5,5+i*5) = eps(2,1)
      matrix(1+i*5,2+i*5) = tpd*exp((0d0,1d0)*kx)
      matrix(2+i*5,1+i*5) = conjg(matrix(1+i*5,2+i*5))
      matrix(1+i*5,5+i*5) = tpd*exp((0d0,1d0)*ky) ! Changed sign out front
      matrix(5+i*5,1+i*5) = conjg(matrix(1+i*5,5+i*5))
      matrix(2+i*5,4+i*5) = -2d0*tpp*dcos(kx-ky)
      matrix(4+i*5,2+i*5) = conjg(matrix(2+i*5,4+i*5))
      matrix(3+i*5,5+i*5) = -2d0*tpp*dcos(kx-ky)
      matrix(5+i*5,3+i*5) = conjg(matrix(3+i*5,5+i*5))
      matrix(2+i*5,5+i*5) = 2d0*tpp_p_11*dcos(kx-ky) ! Changed sign out front
      matrix(5+i*5,2+i*5) = conjg(matrix(2+i*5,5+i*5))
      matrix(3+i*5,4+i*5) = 2d0*tpp_p_22*dcos(kx-ky) ! Changed sign out front
      matrix(4+i*5,3+i*5) = conjg(matrix(3+i*5,4+i*5))
    enddo

    ! H_M matrix set
    do j=0,1
    do i=0,1
      k = 5*(1-i)+j*10
      n = i*5+j*10
      matrix(1+k,2+n) = -tpd*exp((0d0,-1d0)*kx) ! Changed power sign
      matrix(2+n,1+k) = conjg(matrix(1+k,2+n))
      matrix(1+k,5+n) = -tpd*exp((0,-1d0)*ky) ! Changed sign out front and prower
      matrix(5+n,1+k) = conjg(matrix(1+k,5+n))
      matrix(2+k,4+n) = -2d0*tpp*dcos(kx+ky)
      matrix(4+n,2+k) = conjg(matrix(2+k,4+n))
      matrix(2+k,5+n) =  -2d0*tpp_p_11*dcos(kx+ky) ! Changed sign out front
      matrix(5+n,2+k) = conjg(matrix(2+k,5+n))
      matrix(3+k,4+n) =  -2d0*tpp_p_22*dcos(kx+ky) ! Changed sign out front
      matrix(4+n,3+k) = conjg(matrix(3+k,4+n))
      matrix(3+k,5+n) = -2d0*tpp*dcos(kx+ky)
      matrix(5+n,3+k) = conjg(matrix(3+k,5+n))
    end do
    end do


end subroutine
