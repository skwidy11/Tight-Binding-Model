double precision function bessj1(x)
integer, parameter :: double = KIND(0.d0)
real(kind=double) x
! Returns the Bessel function J1(x) for any real x.
real(kind=double) ax,xx,z
real(kind=double) p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6,y
SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,s5,s6
data r1,r2,r3,r4,r5,r6/72362614232.d0,-7895059235.d0,242396853.1d0, -2972611.439d0,15704.48260d0,-30.16036606d0/
data s1,s2,s3,s4,s5,s6/144725228442.d0,2300535178.d0, 18583304.74d0,99447.43394d0,376.9991397d0,1.d0/
data p1,p2,p3,p4,p5/1.d0,.183105d-2,-.3516396496d-4,.2457520174d-5, -.240337019d-6/
data q1,q2,q3,q4,q5/.04687499995d0,-.2002690873d-3, .8449199096d-5,-.88228987d-6,.105787412d-6/

if(dabs(x).lt.8.0d0) then
 y=x**2
 bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))))) /(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
else
 ax=dabs(x)
 z=8.0d0/ax
 y=z**2
 xx=ax-2.356194491d0
 bessj1=dsqrt(.636619772d0/ax)*(dcos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))-z*dsin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*sign(1.0d0,x)
endif
return
end
