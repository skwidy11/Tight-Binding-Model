! Subroutine to convert an integer number into a character (4 digits).
!
  subroutine numchar4 ( ireal, ext )
!
  implicit none
!
  integer, intent(in) :: ireal
  character (len=4), intent(out) :: ext
!
  integer :: ir1, ir2, ir3, ir4
!
  ir1 = ireal/1000
  ir2 = (ireal - ir1*1000)/100
  ir3 = (ireal - ir1*1000 - ir2*100)/10
  ir4 = ireal - ir1*1000 - ir2*100 - ir3*10
  ext = achar(ir1-10*(ir1/10)+48)//achar(ir2-10*(ir2/10)+48) &
         //achar(ir3-10*(ir3/10)+48)//achar(ir4-10*(ir4/10)+48)
!
  end subroutine numchar4
