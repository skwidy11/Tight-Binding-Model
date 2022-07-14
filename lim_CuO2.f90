
! Generate square, planar CuO2 lattice, for use with electron transport codes
! P. Schelling June 8, 2020

IMPLICIT NONE
INTEGER, PARAMETER :: Prec14=SELECTED_REAL_KIND(14)
INTEGER, PARAMETER :: n1=12,n2=12,n3=1
INTEGER, PARAMETER :: natoms=3*n1*n2*n3
REAL(KIND=Prec14), PARAMETER :: a=1.0d0 ! Lattice parameter in Angstroms
REAL(KIND=Prec14), DIMENSION(3,3) :: r
REAL(KIND=Prec14), DIMENSION(3,3) :: h
REAL(KIND=Prec14) :: rx,ry,rz,amp,rc,p,vac,percVac
CHARACTER*1 :: label
CHARACTER*3 :: sname
INTEGER :: i,j,k,n,nt,na,nb,nc,nd,norb,typ,maxvac,numvac,numOx
INTEGER, DIMENSION(3) :: ntype
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: neighbor

! Max vac thing doesn't work. Maybe I should fix someday?

! What do the following variables do?
! a, amp, rc, label, sname, na, nb, nc, nd, norb
! Is the variable nt necessary?
amp=0.0d0

! p - Probability for an oxygen to be a vacancy
p=0.00d0
! maxvac - maximum number of vacancies adjacent to a Cu atom
maxvac=3

call random_seed()

norb=2

open(unit=1, file="structure")

write(1,*) natoms, norb, p
! write(1,*) n1, n2, n3

rc=1.1d0


h=0.0d0

h(1,1)=dfloat(n1)
h(2,2)=dfloat(n2)
h(3,3)=100.0d0

write(1,100) h(1,1),h(1,2),h(1,3)
write(1,100) h(2,1),h(2,2),h(2,3)
write(1,100) h(3,1),h(3,2),h(3,3)

100 format(3f12.6)

write(1,*) a,rc

! Basis vectors for atoms within cell
r(1,1)=0.0d0
r(1,2)=0.0d0
r(1,3)=0.0d0

r(2,1)=0.5d0
r(2,2)=0.0d0
r(2,3)=0.0d0

r(3,1)=0.0d0
r(3,2)=0.5d0
r(3,3)=0.0d0

! typ physical significance:
! 1 - Cu not adjacent to vacancy
! 3 - Cu adjacent to vacancy
! 2 - O (non-vacancy)
! 2 - O (vacancy)
ntype(1)=1
ntype(2)=2
ntype(3)=2

nt=0; na=0; nb=0; nc=0

allocate(neighbor(n1+1,n2+1,n3+1))
neighbor=0
numvac=0
do i=0,n1-1
do j=0,n2-1
do k=0,n3-1
do n=1,3
rx=r(n,1)+(dfloat(i)/dfloat(n1))*h(1,1)+(dfloat(j)/dfloat(n2))*h(2,1)+(dfloat(k)/dfloat(n3))*h(3,1)
ry=r(n,2)+(dfloat(i)/dfloat(n1))*h(1,2)+(dfloat(j)/dfloat(n2))*h(2,2)+(dfloat(k)/dfloat(n3))*h(3,2)
rz=r(n,3)+(dfloat(i)/dfloat(n1))*h(1,3)+(dfloat(j)/dfloat(n2))*h(2,3)+(dfloat(k)/dfloat(n3))*h(3,3)
rx=rx-h(1,1)/2.0d0; ry=ry-h(2,2)/2.0d0; rz=1.0d0

typ=ntype(n)
call random_number(vac)
if(vac.lt.p) then ! Adjust here to make new site
  if(n.eq.2 .and. neighbor(i+1,j+1,k+1).lt.maxvac .and. neighbor(i+2,j+1,k+1).lt.maxvac) then
    typ=4
    neighbor(i+1,j+1,k+1)=neighbor(i+1,j+1,k+1)+1
    neighbor(i+2,j+1,k+1)=neighbor(i+2,j+1,k+1)+1
    numvac=numvac+1
  elseif(n.eq.3 .and. neighbor(i+1,j+1,k+1).lt.maxvac .and. neighbor(i+1,j+2,k+1).lt.maxvac) then
    typ=4
    neighbor(i+1,j+1,k+1)=neighbor(i+1,j+1,k+1)+1
    neighbor(i+1,j+2,k+1)=neighbor(i+1,j+2,k+1)+1
    numvac=numvac+1
  endif
endif

if(n.ne.1) then
  nt=nt+1
  write(1,200) nt,rx,ry,rz,typ
  numOx=numOx+1
endif
enddo
enddo
enddo
enddo


! Assigning copper's neighboring vacancies
do i=0,n1-1
do j=0,n2-1
do k=0,n3-1
  typ=1
  rx=r(1,1)+(dfloat(i)/dfloat(n1))*h(1,1)+(dfloat(j)/dfloat(n2))*h(2,1)+(dfloat(k)/dfloat(n3))*h(3,1)
  ry=r(1,2)+(dfloat(i)/dfloat(n1))*h(1,2)+(dfloat(j)/dfloat(n2))*h(2,2)+(dfloat(k)/dfloat(n3))*h(3,2)
  rz=r(1,3)+(dfloat(i)/dfloat(n1))*h(1,3)+(dfloat(j)/dfloat(n2))*h(2,3)+(dfloat(k)/dfloat(n3))*h(3,3)
  rx=rx-h(1,1)/2.0d0; ry=ry-h(2,2)/2.0d0; rz=1.0d0
  if(neighbor(i+1,j+1,k+1).ge.1) typ=3
  nt=nt+1
  write(1,200) nt,rx,ry,rz,typ
enddo
enddo
enddo
close(1)

percVac=dfloat(numvac)/dfloat(numOx) * 100.0
print *, numvac, numOx
write(6,300) 'Percent of Oxygens with vacancies:', percVac, '%'

200 format(i6,3f16.9,i4)
300 format(A,f6.2,A)
stop
end
