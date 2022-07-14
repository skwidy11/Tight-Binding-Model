! Subroutine to generate a Chebyshev series using recursion relation (order n > 1).
!
! Input: Energy
!
! Output: T_n(E) in array tn, Chebyshev expansion
!
!
! Patrick Schelling, UCF (August 8, 2018)
!
  subroutine cheb(NT, nesteps, epsmin, epsmax,  g, tn )
!
!
  include 'params.inc'
  double precision E, epsmin, epsmax
  double precision, dimension(1:nesteps+1,0:NT-1) :: tn
  double precision, dimension(0:NT-1) :: g
  integer  n, NT, nesteps, ne

  do ne=1, nesteps+1
        E=epsmin+dfloat(ne-1)*(epsmax-epsmin)/dfloat(nesteps)
  
        tn(ne,0)=1.0d0; tn(ne,1)=E

        do n=2,NT-1
                tn(ne,n)=2.0d0*E*tn(ne,n-1)-tn(ne,n-2)
        enddo
  enddo

  do ne=1, nesteps+1
  tn(ne,0)=tn(ne,0)*g(0)      
  do n=1,NT-1  
        tn(ne,n)=tn(ne,n)*2.0d0*g(n)      
  enddo
  enddo
 
!
  end subroutine cheb

