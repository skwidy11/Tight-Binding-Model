double precision function bessjn(n,x)
integer, parameter :: double = KIND(0.d0)
integer n,IACC
real(kind=double) bessj,x,BIGNO,BIGNI
parameter (IACC=40,BIGNO=1.d10,BIGNI=1.d-10)
! uses bessj0,bessj1
! Returns the Bessel function Jn(x) for any real x and n > 2.
integer j,jsum,m
real(kind=double) ax,bj,bjm,bjp,sum,tox,bessj0,bessj1
! if(n.lt.2) pause
 ax=dabs(x) 
 if(ax.eq.0.d0) then
  bessjn=0.0d0
 else if(ax.gt.dfloat(n))then
  tox=2.0d0/ax
  bjm=bessj0(ax)
  bj=bessj1(ax)
  do j=1,n-1
    bjp=j*tox*bj-bjm
    bjm=bj
    bj=bjp
  enddo
  bessjn=bj
 else
  tox=2.0d0/ax
  m=2*((n+int(sqrt(float(IACC*n))))/2)
  bessj=0.0d0
  jsum=0
  sum=0.0d0
  bjp=0.0d0
  bj=1.0d0
  do j=m,1,-1
   bjm=j*tox*bj-bjp
   bjp=bj
   bj=bjm
   if(dabs(bj).gt.BIGNO) then
    bj=bj*BIGNI
    bjp=bjp*BIGNI
    bessjn=bessjn*BIGNI
    sum=sum*BIGNI
   endif
   if(jsum.ne.0)sum=sum+bj
   jsum=1-jsum
   if(j.eq.n) bessjn=bjp
  enddo
  sum=2.0d0*sum-bj
  bessjn=bessjn/sum
 endif

 if(x.lt.0..and.mod(n,2).eq.1) bessjn=-bessjn

return
end
