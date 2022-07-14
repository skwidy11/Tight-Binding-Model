! Invert a 3x3 matrix, obtain inverted matrix and also the determinant of original matrix
subroutine matinv(a,b,c)
        INTEGER, PARAMETER :: Prec14=SELECTED_REAL_KIND(14)
        REAL(Kind=Prec14), DIMENSION(3,3) :: a,b
        REAL(Kind=Prec14) :: d11,d12,d13
        REAL(Kind=Prec14) :: d21,d22,d23
        REAL(Kind=Prec14) :: d31,d32,d33
        REAL(Kind=Prec14) :: c

        d11=a(2,2)*a(3,3)-a(2,3)*a(3,2)
        d22=a(3,3)*a(1,1)-a(3,1)*a(1,3)
        d33=a(1,1)*a(2,2)-a(1,2)*a(2,1)
        d12=a(2,3)*a(3,1)-a(2,1)*a(3,3)
        d23=a(3,1)*a(1,2)-a(3,2)*a(1,1)
        d31=a(1,2)*a(2,3)-a(1,3)*a(2,2)
        d13=a(2,1)*a(3,2)-a(3,1)*a(2,2)
        d21=a(3,2)*a(1,3)-a(1,2)*a(3,3)
        d32=a(1,3)*a(2,1)-a(2,3)*a(1,1)
        c=a(1,1)*d11+a(1,2)*d12+a(1,3)*d13
        b(1,1)=d11/c
        b(2,2)=d22/c
        b(3,3)=d33/c
        b(1,2)=d21/c
        b(2,3)=d32/c
        b(3,1)=d13/c
        b(2,1)=d12/c
        b(3,2)=d23/c
        b(1,3)=d31/c

        return
        end
