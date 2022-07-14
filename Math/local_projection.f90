! Subroutine to calculate local projection of oxygen/copper atoms separately.

subroutine local_projection ( ND, KT, it, norb, psi0, psi, ntype, indx, ampl )
  implicit none
  integer, intent(in) :: ND, norb, indx, KT, it
  double complex, intent(in), dimension(ND,norb) :: psi0, psi
  integer, intent(in), dimension(ND) :: ntype
  double complex, intent(out) :: ampl(0:KT-1,2)
  integer :: k, mu

  do k = 1, ND
    if (ntype(k) == 2) then
      ampl(it,2) = ampl(it,2) + dconjg(psi0(k,indx))*psi(k,indx)
    else
      ampl(it,1) = ampl(it,1) + dconjg(psi0(k,indx))*psi(k,indx)
    endif
  enddo

!  print *, "Ox=",ampl(it,2)," Cu=",ampl(it,1)

end subroutine local_projection
