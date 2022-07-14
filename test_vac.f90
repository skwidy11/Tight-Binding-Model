program test_vac
  integer :: total, i, n, typ
  double precision :: rx,ry,rz
  total = 0
  open(unit=1, file="structure")

  read (1,*)
  read (1,*)
  read (1,*)
  read (1,*)
  read (1,*)

  do i=1,24*24*3
    read (1,*) n, rx, ry, rz, typ
    !print *, typ
    if (typ==4) then
      total = total + 1
    endif
  end do
  print *, total

end program
