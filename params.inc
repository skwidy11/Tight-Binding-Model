        
        implicit none 

 	integer :: ierr, mynod, nproc
        INTEGER, PARAMETER :: Prec14=SELECTED_REAL_KIND(14)
        INTEGER, PARAMETER :: NZ=30
	REAL(Kind=Prec14) :: rc,a

        COMMON /params1/ rc,a
