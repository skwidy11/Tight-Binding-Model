subroutine update_band(kx,ky,matrix,t,eps,nocc,site,U)
    complex, dimension(10,10), intent(out) :: matrix
    double precision,intent(in) :: kx,ky,U
    double precision, dimension(3,3,9,9) :: t
    double precision, dimension(3,9) :: eps
    double precision, dimension(2880,2) :: nocc
    integer :: site,i
    double precision :: tpd,tpp,tpp_p
    tpd = Abs(t(1,2,1,1))
    tpp = Abs(t(2,2,2,1))
    tpp_p = Abs(t(2,2,1,1))
    matrix = 0.0d0
    do i=0,1
      if(i.eq.0) then
        matrix(1,1) = U*(0.5-nocc(site,1))
      else
        matrix(6,6) = U*(0.5-nocc(site,2))
      end if
      matrix(2+i*5,2+i*5) = eps(2,1)
      matrix(3+i*5,3+i*5) = eps(2,2)
      matrix(4+i*5,4+i*5) = eps(2,2)
      matrix(5+i*5,5+i*5) = eps(2,1)
      matrix(1+i*5,2+i*5) = (0,-1d0)*2d0*tpd*dsin(kx)!*exp((0,1d0)*kx)
      matrix(2+i*5,1+i*5) = conjg(matrix(1+i*5,2+i*5))
      matrix(1+i*5,5+i*5) = (0,1d0)*2d0*tpd*dsin(ky)!*exp((0,1d0)*ky)
      matrix(5+i*5,1+i*5) = conjg(matrix(1+i*5,5+i*5))
      matrix(2+i*5,4+i*5) = -4d0*tpp*dcos(kx)*dcos(ky)!*exp((0,1d0)*(ky-kx))
      matrix(4+i*5,2+i*5) = conjg(matrix(2+i*5,4+i*5))
      matrix(3+i*5,5+i*5) = -4d0*tpp*dcos(kx)*dcos(ky)!*exp((0,1d0)*(ky-kx))
      matrix(5+i*5,3+i*5) = conjg(matrix(3+i*5,5+i*5))
      matrix(2+i*5,5+i*5) = -4d0*tpp_p*dsin(kx)*dsin(ky)!*exp((0,1d0)*(ky-kx))
      matrix(5+i*5,2+i*5) = conjg(matrix(2+i*5,5+i*5))
      matrix(3+i*5,4+i*5) = -4d0*tpp_p*dsin(kx)*dsin(ky)!*exp((0,1d0)*(ky-kx))
      matrix(4+i*5,3+i*5) = conjg(matrix(3+i*5,4+i*5))
    enddo
end subroutine
