! Subroutine to initialize the random number generator.
!
! E. Mucciolo, UCF (Apr 01, 2014)
!
  subroutine initrand ( iseed )
!
  implicit none
!
  integer, intent(in) :: iseed
!
  integer :: k, n
  integer, dimension(:), allocatable:: seed 
!
! initialize the random sequence
!
  call random_seed (size=n)
  allocate (seed(n))
  do k = 1, n
     seed(k) = iseed + n - 1
  end do
!
! generate random numbers
!
  call random_seed (put=seed)
  deallocate (seed)
!
  end subroutine initrand
