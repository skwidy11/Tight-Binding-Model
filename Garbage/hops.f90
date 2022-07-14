!    s-s hops

            e=ess_sigma; f=fss_sigma; g=gss_sigma 
            tm=hop(e,f,g,rij)
            mu=1; nu=1  
            t=tm
            hamilt(i,nn,mu,nu)=t

!    s-p hops

            e=esp_sigma; f=fsp_sigma; g=gsp_sigma 
            tm=hop(e,f,g,rij)

            mu=1; nu=2 
            t=xl*tm
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=-t

            mu=1; nu=3  
            t=yl*tm
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=-t

            mu=1; nu=4  
            t=zl*tm
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=-t

!    s-d hops

            e=esd_sigma; f=fsd_sigma; g=gsd_sigma 
            tm=hop(e,f,g,rij)

            mu=1; nu=5  
            t=xl*yl*tm
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=t

            mu=1; nu=6
            t=yl*zl*tm
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=t

            mu=1; nu=7 
            t=xl*zl*tm
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=t

            mu=1; nu=8  
            t=(xl**2-yl**2)*tm
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=t

            mu=1; nu=9  
            t=(zl**2-0.50d0*(xl**2+yl**2))*tm
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=t

!    p-p hops

	    e=epp_sigma; f=fpp_sigma; g=gpp_sigma
	    tm1=hop(e,f,g,rij)
	    e=epp_pi; f=fpp_pi; g=gpp_pi
	    tm2=hop(e,f,g,rij)
	    
            mu=2; nu=2
            t=xl**2*tm1+(1.0d0-xl**2)*tm2
            hamilt(i,nn,mu,nu)=t

            mu=2; nu=3
            t=xl*yl*(tm1-tm2)
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=t

            mu=2; nu=4
            t=xl*zl*(tm1-tm2)
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=t

            mu=3; nu=3
            t=yl**2*tm1+(1.0d0-yl**2)*tm2
            hamilt(i,nn,mu,nu)=t

            mu=3; nu=4
            t=yl*zl*(tm1-tm2)
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=t

            mu=4; nu=4
            t=zl**2*tm1+(1.0d0-zl**2)*tm2
            hamilt(i,nn,mu,nu)=t

!    p-d hops

	    e=epd_sigma; f=fpd_sigma; g=gpd_sigma
	    tm1=hop(e,f,g,rij)
	    e=epd_pi; f=fpd_pi; g=gpd_pi
	    tm2=hop(e,f,g,rij)
	    mu=2; nu=5
            t=dsqrt(3.0d0)*xl**2*yl*tm1+yl*(1.0d0-2.0d0*xl**2)*tm2
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=-t
	    
	    mu=2; nu=6
            t=dsqrt(3.0d0)*xl*yl*zl*tm1-2.0d0*xl*yl*zl*tm2
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=-t

	    mu=2; nu=7
	    t=dsqrt(3.0d0)*xl**2*zl*tm1+zl*(1.0d0-2.0d0*xl**2)*tm2
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=-t

	    mu=2; nu=8
	    t=0.50d0*dsqrt(3.0d0)*xl*(xl**2-yl**2)*tm1+xl*(1.0d0-xl**2+yl**2)*tm2
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=-t

	    mu=2; nu=9
	    t=xl*(zl**2-0.50d0*(xl**2+yl**2))*tm1+dsqrt(3.0d0)*xl*zl**2*tm2
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=-t

	    mu=3; nu=5
	    t=dsqrt(3.0d0)*xl*yl**2*tm1+xl*(1.0d0-2.0d0*yl**2)*tm2
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=-t

	    mu=3; nu=6
	    t=dsqrt(3.0d0)*yl**2*zl*tm1+zl*(1.0d0-2.0d0*yl**2)*tm2
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=-t

	    mu=3; nu=7
	    t=dsqrt(3.0d0)*xl*yl*zl*tm1-2.0d0*xl*yl*zl*tm2
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=-t

	    mu=3; nu=8
	    t=0.50d0*dsqrt(3.0d0)*yl*(xl**2-yl**2)*tm1-yl*(1.0d0+xl**2-yl**2)*tm2
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=-t

	    mu=3; nu=9
    	    t=yl*(zl**2-0.50d0*(xl**2+yl**2))*tm1-dsqrt(3.0d0)*yl*zl**2*tm2
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=-t

	    mu=4; nu=5
	    t=dsqrt(3.0d0)*xl*yl*zl*tm1-2.0d0*xl*yl*zl*tm2
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=-t

	    mu=4; nu=6
	    t=dsqrt(3.0d0)*zl**2*yl*tm1+yl*(1.0d0-2.0d0*zl**2)*tm2
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=-t

	    mu=4; nu=7
	    t=dsqrt(3.0d0)*zl**2*xl*tm1+xl*(1.0d0-2.0d0*zl**2)*tm2
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=-t

	    mu=4; nu=8
	    t=0.50d0*dsqrt(3.0d0)*zl*(xl**2-yl**2)*tm1-zl*(xl**2-yl**2)*tm2
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=-t

	    mu=4; nu=9
	    t=zl*(zl**2-0.50d0*(xl**2+yl**2))*tm1+dsqrt(3.0d0)*zl*(xl**2+yl**2)*tm2
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=-t

!    d-d hops

	    e=edd_sigma; f=fdd_sigma; g=gdd_sigma
	    tm1=hop(e,f,g,rij)
	    e=edd_pi; f=fdd_pi; g=gdd_pi
	    tm2=hop(e,f,g,rij)
	    e=edd_delta; f=fdd_delta; g=gdd_delta
	    tm3=hop(e,f,g,rij)

	    mu=5; nu=5
	    t=3.0d0*xl**2*yl**2*tm1+(xl**2+yl**2-4.0d0*xl**2*yl**2)*tm2+(zl**2+xl**2*yl**2)*tm3
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=t

	    mu=5; nu=6
   	    t=3.0d0*xl*yl**2*zl*tm1+xl*zl*(1.0d0-4.0d0*yl**2)*tm2+xl*zl*(yl**2-1.0d0)*tm3
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=t

	    mu=5; nu=7
 	    t=3.0d0*xl**2*yl*zl*tm1+yl*zl*(1.0d0-4.0d0*xl**2)*tm2+yl*zl*(xl**2-1.0d0)*tm3
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=t

	    mu=5; nu=8
	    t=1.5d0*xl*yl*(xl**2-yl**2)*tm1+2.0d0*xl*yl*(yl**2-xl**2)*tm2+0.50d0*xl*yl*(xl**2-yl**2)*tm3
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=t

	    mu=5; nu=9
            t=dsqrt(3.0d0)*xl*yl*(zl**2-0.50d0*(xl**2+yl**2))*tm1-2.0d0*dsqrt(3.0d0)*xl*yl*zl**2*tm2+0.50d0*dsqrt(3.0d0)*xl*yl*(1.0d0+zl**2)*tm3
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=t

	    mu=6; nu=6
	    t=3.0d0*yl**2*zl**2*tm1+(yl**2+zl**2-4.0d0*yl**2*zl**2)*tm2+(xl**2+yl**2*zl**2)*tm3
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=t

	    mu=6; nu=7
	    t=3.0d0*xl*yl*zl**2*tm1+xl*yl*(1.0d0-4.0d0*zl**2)*tm2+xl*yl*(zl**2-1.0d0)*tm3
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=t

	    mu=6; nu=8
	    t=1.50d0*yl*zl*(xl**2-yl**2)*tm1-yl*zl*(1.0d0+2.0d0*(xl**2-yl**2))*tm2+yl*zl*(1.0d0+0.50d0*(xl**2-yl**2))*tm3
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=t

	    mu=6; nu=9
	    t=dsqrt(3.0d0)*yl*zl*(zl**2-0.50d0*(xl**2+yl**2))*tm1+dsqrt(3.0d0)*yl*zl*(xl**2+yl**2-zl**2)*tm2+0.50d0*dsqrt(3.0d0)*yl*zl*(xl**2+yl**2)*tm3
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=t

	    mu=7; nu=7
	    t=3.0d0*xl**2*zl**2*tm1+(zl**2+xl**2-4.0*xl**2*zl**2)*tm2+(yl**2+xl**2*zl**2)*tm3
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=t

	    mu=7; nu=8
	    t=1.5d0*xl*zl*(xl**2-yl**2)*tm1+xl*zl*(1.0d0-2.0d0*(xl**2-yl**2))*tm2-xl*zl*(1.0d0-0.50d0*(xl**2-yl**2))*tm3
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=t

	    mu=7; nu=9
	    t=dsqrt(3.0d0)*xl*zl*(zl**2-0.50d0*(xl**2+yl**2))*tm1+dsqrt(3.0d0)*xl*zl*(xl**2+yl**2-zl**2)*tm2-0.50d0*dsqrt(3.0d0)*xl*zl*(xl**2+yl**2)*tm3
            
 	    hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=t

	    mu=8; nu=8
	    t=0.75d0*((xl**2-yl**2)**2)*tm1+(xl**2+yl**2-((xl**2-yl**2)**2))*tm2+(zl**2+0.25*((xl**2-yl**2)**2))*tm3

            hamilt(i,nn,mu,nu)=t

	    mu=8; nu=9
	    t=0.50d0*dsqrt(3.0d0)*(xl**2-yl**2)*(zl**2-0.5d0*(xl**2+yl**2))*tm1+dsqrt(3.0d0)*zl**2*(yl**2-xl**2)*tm2+0.25d0*dsqrt(3.0d0)*(1.0d0+zl**2)*(xl**2-yl**2)*tm3
            hamilt(i,nn,mu,nu)=t
            hamilt(i,nn,nu,mu)=t

	    mu=9; nu=9
	    t=((zl**2-0.50d0*(xl**2+yl**2))**2)*tm1+3.0d0*zl**2*(xl**2+yl**2)*tm2+0.75d0*((xl**2+yl**2)**2)*tm3
            hamilt(i,nn,mu,nu)=t
