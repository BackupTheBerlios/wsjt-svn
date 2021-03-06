      real function fchisq(cx,npts,fsample,nflip,a,ccfmax,dtmax)

      parameter (NMAX=60*96000)          !Samples per 60 s
      complex cx(npts)
      real a(5)
      complex w,wstep,za,zb,z
      real ss(2600)
      complex csx(0:NMAX/64)
      data twopi/6.283185307/a1,a2,a3/99.,99.,99./
      save

      baud=11025.0/4096.0
      if(a(1).ne.a1 .or. a(2).ne.a2 .or. a(3).ne.a3) then
         a1=a(1)
         a2=a(2)
         a3=a(3)

C  Mix and integrate the complex X and Y signals
         csx(0)=0.
         w=1.0
         x0=0.5*(npts+1)
         s=2.0/npts
         do i=1,npts
            x=s*(i-x0)
            if(mod(i,100).eq.1) then
               p2=1.5*x*x - 0.5
!               p3=2.5*(x**3) - 1.5*x
!               p4=4.375*(x**4) - 3.75*(x**2) + 0.375
               dphi=(a(1) + x*a(2) + p2*a(3)) * (twopi/fsample)
               wstep=cmplx(cos(dphi),sin(dphi))
            endif
            w=w*wstep
            csx(i)=csx(i-1) + w*cx(i)
         enddo
      endif

C  Compute 1/2-symbol powers at 1/16-symbol steps.
      fac=1.e-4
      nsps=nint(fsample/baud)                  !Samples per symbol
      nsph=nsps/2                              !Samples per half-symbol

      ndiv=16                                  !Output ss() steps per symbol
      nout=ndiv*npts/nsps
      dtstep=1.0/(ndiv*baud)                   !Time per output step

      do i=1,nout
         j=i*nsps/ndiv
         k=j-nsph
         ss(i)=0.
         if(k.ge.1) then
            za=csx(j)-csx(k)
            z=za
            ss(i)=fac*(real(z)**2 + aimag(z)**2)
         endif
      enddo

      ccfmax=0.
      call ccf2(ss,nout,nflip,ccf,lagpk)
      if(ccf.gt.ccfmax) then
         ccfmax=ccf
         dtmax=lagpk*dtstep
      endif
      fchisq=-ccfmax

      return
      end
