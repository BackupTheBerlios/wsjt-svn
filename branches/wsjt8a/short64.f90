subroutine short64(dat,jz,NFreeze,MouseDF,DFTolerance,                 &
     mode64,nspecialbest,nstest,dfsh,iderrbest,idriftbest,              &
     snrdb,ss1a,ss2a,nwsh,idfsh)

! Check to see if this might be a shorthand message.
! This is done before zapping, downsampling, or normal decoding.

  parameter (NP2=60*12000)               !Size of dat array
  parameter (NFFT=16384)                 !FFT length
  parameter (NH=NFFT/2)                  !Step size
  parameter (MAXSTEPS=90)       !Max # of steps

  real dat(jz)
  integer DFTolerance
  real s2(NH,MAXSTEPS)                   !2d spectrum
  real ss(NH,4)                          !Save spectra in four phase bins
  real psavg(NH)
  real sigmax(4)                         !Peak of spectrum at each phase
  real ss1a(-224:224)                    !Lower magenta curve
  real ss2a(-224:224)                    !Upper magenta curve
  real ss1(-473:1784)                    !Lower magenta curve (temp)
  real ss2(-473:1784)                    !Upper magenta curve (temp)
!  real ssavg(-11:11)
  real ssavg(-18:18)
  integer ipk(4)                         !Peak bin at each phase
  save

  nspecialbest=0                         !Default return value
  nstest=0
  df=12000.0/NFFT
  f0=1270.46

! Do 16 k FFTs, stepped by 8k.  (*** Maybe should step by 4k? ***)
  call zero(psavg,NH)
  nsteps=(jz-NH)/(4*NH)
  nsteps=4*nsteps                        !Number of steps
  do j=1,nsteps
     k=(j-1)*NH + 1
     call ps(dat(k),NFFT,s2(1,j))       !Get power spectra
     if(mode64.eq.4) then
        call smooth(s2(1,j),NH)
        call smooth(s2(1,j),NH)
     endif
     call add(psavg,s2(1,j),psavg,NH)
  enddo
  
  call flat1(psavg,s2,NH,nsteps,NH,MAXSTEPS)

  xfac=(16384./6480.0) * 21 * mode64
  dtstep=0.5/df
  fac=dtstep/(60.0*df)

! Define range of frequencies to be searched
  fa=max(200.0,f0+MouseDF-600.0)
  fb=min(4800.0,f0+MouseDF+600.0)
  ia=fa/df
  ib=fb/df + 3.1*xfac  !Tone for '73' is above sync tone by 3*xfac*df Hz
  if(NFreeze.eq.1) then
     fa=max(200.0,f0+MouseDF-DFTolerance)
     fb=min(4800.0,f0+MouseDF+DFTolerance)
  endif
  ia2=fa/df
  ib2=fb/df + 3.1*xfac  !Upper tone is above sync tone by 3*xfac*df Hz
  if(ib2.gt.NH) ib2=NH

! Find strongest line in each of the 4 phases, repeating for each drift rate.
  sbest=0.
  snrbest=0.
  idz=6.0/df                !Is this the right drift range?
  do idrift=-idz,idz
     drift=idrift*df*60.0/49.04
     call zero(ss,4*NH)     !Clear the accumulating array
     do j=1,nsteps
        n=mod(j-1,4)+1
        k=nint((j-nsteps/2)*drift*fac) + ia
        call add(ss(ia,n),s2(k,j),ss(ia,n),ib-ia+1)
     enddo

     do n=1,4
        sigmax(n)=0.
        do i=ia2,ib2
           sig=ss(i,n)
           if(sig.ge.sigmax(n)) then
              ipk(n)=i
              sigmax(n)=sig
              if(sig.ge.sbest) then
                 sbest=sig
                 nbest=n
              endif
           endif
        enddo
     enddo
     n2best=nbest+2
     if(n2best.gt.4) n2best=nbest-2
     xdf=min(ipk(nbest),ipk(n2best))*df - f0
     if(NFreeze.eq.1 .and. abs(xdf-mousedf).gt.DFTolerance) goto 10

     idiff=abs(ipk(nbest)-ipk(n2best))
     xk=float(idiff)/xfac
     k=nint(xk)
     iderr=nint((xk-k)*xfac)
     nspecial=0
     maxerr=nint(0.008*abs(idiff) + 0.51)
     if(abs(iderr).le.maxerr .and. k.ge.1 .and. k.le.3) nspecial=k
     if(nspecial.gt.0) then
        call getsnr(ss(ia2,nbest),ib2-ia2+1,snr1)
        call getsnr(ss(ia2,n2best),ib2-ia2+1,snr2)
        snr=0.5*(snr1+snr2)
        if(snr.gt.snrbest) then
           snrbest=snr
           nspecialbest=nspecial
           nstest=snr/2.0 - 2.0             !Threshold set here
           if(nstest.lt.0) nstest=0
           if(nstest.gt.10) nstest=10
           dfsh=nint(xdf)
           iderrbest=iderr
           idriftbest=idrift
           snrdb=db(snr) - db(2500.0/df) - db(sqrt(nsteps/4.0))+1.8
           n1=nbest
           n2=n2best
           ipk1=ipk(n1)
           ipk2=ipk(n2)
        endif
     endif
     if(nstest.eq.0) nspecial=0
10   continue
  enddo

  if(nstest.eq.0) nspecialbest=0
  df4=(16384.0/6480.0)*df
  if(nstest.gt.0) then

     if(ipk1.gt.ipk2) then
        ntmp=n1
        n1=n2
        n2=ntmp
        ntmp=ipk1
        ipk1=ipk2
        ipk2=ntmp
     endif

     call zero(ss1,2258)
     call zero(ss2,2258)

     ka=(ia2*df-f0)/df4
     kb=(ib2*df-f0)/df4
     do k=ka,kb
        i=nint((f0+k*df4)/df)
        ss1(k)=0.5 * (ss(i-1,n1) + ss(i,n1) + ss(i+1,n1))
        ss2(k)=0.5 * (ss(i-1,n2) + ss(i,n2) + ss(i+1,n2))
     enddo

     kpk1=nint(0.25*ipk1-472.0)
     kpk2=kpk1 + nspecial*mode64*10
     ssmax=0.
     do i=-10,10
        ssavg(i)=ss1(kpk1+i) + ss2(kpk2+i)
        if(ssavg(i).gt.ssmax) then
           ssmax=ssavg(i)
           itop=i
        endif
     enddo
     base=0.25*(ssavg(-10)+ssavg(-9)+ssavg(9)+ssavg(10))
     shalf=0.5*(ssmax+base)
     do k=1,8
        if(ssavg(itop-k).lt.shalf) go to 110
     enddo
     k=8
110  x=(ssavg(itop-(k-1))-shalf)/(ssavg(itop-(k-1))-ssavg(itop-k))
     do k=1,8
        if(ssavg(itop+k).lt.shalf) go to 120
     enddo
     k=8
120  x=x+(ssavg(itop+(k-1))-shalf)/(ssavg(itop+(k-1))-ssavg(itop+k))
     nwsh=nint(x*df4)
  endif

! See if orange/magenta curves need to be shifted:
  idfsh=0
  if(mousedf.lt.-600) idfsh=-670
  if(mousedf.gt.600) idfsh=1000
  if(mousedf.gt.1600) idfsh=2000
  if(mousedf.gt.2600) idfsh=3000
  i0=nint(idfsh/df4)

  do i=-224,224
     ss1a(i)=ss1(i+i0)
     ss2a(i)=ss2(i+i0)
  enddo

  return
end subroutine short64
