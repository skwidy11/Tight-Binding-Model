c here I only look at spectral weight at ef
c see Science 345, 187 (2014)
c 3/17/15 This is now running but it is very slow because
c it reads the wavefunction and eigenvalue files again
c for each kx,ky.  The first run,   which didn't quite finish
c, at p=0.7 was disappointing. The largest density was at
c  0,pi not in the pi,pi direction as seen in experiments.
c  However to see it more clearly I should look at a range
c  of p's . Also I think the energy width of the range
c  of eignevalue energies sampled may be too big.
c  Finally, I need to look a different fillings.
c  Program: 1) rewrite to speed up by readiing the energy and
c wave function files just once and putting the results for the
c function around the fermi energy in new files 2) run
c for a lower filling and for a range of p's.
c 4/3/15 starting the revisions noted above.
c 8/1/16 starting revision for 2 band model.  Code called
c spectralweight2b.f
c 8/3/16 seems to be running correctly.  Gets killed if run from the
c keyboard, probably because it's too long. Try running batch.
      double precision energy,delta,ef
      double precision fkx,fky,pi,phasearg
      double precision costrans,sintrans
      double precision,dimension(:,:,:),allocatable :: wavenearfs
c         first index: wflabel, then ix,iy
      double precision, dimension(:), allocatable :: energynearfs
      double precision, dimension(:,:), allocatable :: weight
      integer lx,ly,states,ifill,hx,hy
c          eigenenergies for selected wfs
      open(2,file='eigenvalues2bandp1b')
      open(3,file='eigenvectors2bandp.7b')
      open(7,file='weightatfsp2bandp.7fill85')
c set filling
      lx = 70
      ly = 70
      states = lx*ly
      hx = lx/2
      hy = ly/2

      ifill=3120
      icount=0
      pi=4.0*datan(1.0d0)
      allocate(wavenearfs(states,-hx+1:hx,-hy+1:hy))
      allocate(weight(-hx:hx-1,-hy:hy-1),energynearfs(states))


      do ie=1,states ! sets fermi energy I am pretty sure
        ef=energy
        read(2,*) energy
        if(icount.eq.ifill)then
           goto 11
        endif
        icount=icount+1
      enddo
11    continue

      rewind(2)
      delta =4.0d0/70.0d0
      write(*,*) 'delta',delta
      write(*,*) 'ef,icount',ef,icount
c     read (*,13) dummy
13    format(F12.4)
c 8/1/16 fermi energy checks at p=1
c look at all the states at each kx,ky; removed rewinds
        icount=0
        do ie=1,states
          read(2,*) energy
          if(energy.gt.ef+delta.or.energy.le.ef-delta)then
c read past the wave function if not near fs
            do il=1,states
            read(3,*)
            enddo !il
            goto 12
          endif
          icount=icount+1
          do il=1,states
            ix=mod(il,lx)-hx ! 8/1/16check this
            iy=int(il/ly)-hy+1
            read(3,*) wave
            wavenearfs(icount,ix,iy)=wave
            energynearfs(icount)=energy
          enddo ! il
12      continue
        enddo ! ie

c compute overlap
        write(*,*) 'number near fs',icount
c       read(*,13) dummy

      do ikx=-hx,hx-1 ! the BZ doesn't change for 3 band case
c but need finer grained fkx
      do iky=-hy,hy-1
        weight(ikx,iky)=0.0d0
        fkx=(dfloat(ikx)*pi)/dble(hx)
        fky=(dfloat(iky)*pi)/dble(hy)

         do icc=1,icount ! icc for eigenvectors
        costrans=0.0d0
        sintrans=0.0d0

          do iy=-hx+1,hx !  changed to include  4
c              atoms/cell
          do ix=-hy+1,hy
           phasearg=(dfloat(ix)*fkx+dfloat(iy)*fky)
           wave=wavenearfs(icc,ix,iy)
            costrans=costrans+wave*dcos(phasearg)
            sintrans=sintrans+wave*dsin(phasearg)

c           write(*,*) 'weight for this fkx,fky'
            overlap=(costrans**2+sintrans**2)
            weight(ikx,iky)=weight(ikx,iky)+overlap
          enddo ! ix
          enddo ! iy
       enddo ! icc

c           write(*,*) fkx,fky,weight(ikx,iky)
            write(7,*)  fkx,fky,weight(ikx,iky)
c           pause
       enddo ! iky
       enddo ! ikx
c     save weight
c      do ikx=-50,0
c      do iky=-50,0
c          fkx=dfloat(ikx)*2.*pi/dfloat(100)
c          fky=dfloat(iky)*2.*pi/dfloat(100)
c          write(7,*) fkx,fky,weight(ikx,iky)
c      enddo ! iky
c      enddo ! ikx
       rewind(2)
       rewind(3)
       close(2)
       close(3)
       end
