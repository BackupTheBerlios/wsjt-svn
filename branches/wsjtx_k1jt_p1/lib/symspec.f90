subroutine symspec(k,ntrperiod,nsps,ingain,pxdb,s,red,df3,ihsym,npts8)

! Input:
!  k         pointer to the most recent new data
!  ntrperiod T/R sequence length, minutes
!  nsps      samples per symbol, at 12000 Hz
!  ndiskdat  0/1 to indicate if data from disk
!  nb        0/1 status of noise blanker (off/on)
!  nbslider  NB setting, 0-100

! Output:
!  pxdb      power (0-60 dB)
!  s()       spectrum for waterfall display
!  red()     first cut at JT9 sync amplitude
!  ihsym     index number of this half-symbol (1-184)

  parameter (NTMAX=120)
  parameter (NMAX=NTMAX*12000)        !Total sample intervals per 30 minutes
  parameter (NDMAX=NTMAX*1500)        !Sample intervals at 1500 Hz rate
  parameter (NSMAX=1365)              !Max length of saved spectra
  parameter (NFFT1=1024)
  parameter (MAXFFT3=16384)
  real*4 s(NSMAX),w3(MAXFFT3)
  real*4 x1(NFFT1)
  real*4 x2(NFFT1+105)
  real*4 ssum(NSMAX)
  real*4 red(NSMAX)
  real*4 xc(0:MAXFFT3-1)
  complex cx(0:MAXFFT3/2)
  integer*2 id2
  complex c0
  common/jt9com/ss(184,NSMAX),savg(NSMAX),c0(NDMAX),id2(NMAX),nutc,ndiskdat, &
       ntr,mousefqso,newdat,nfa,nfb,ntol,kin,nzhsym,nsynced,ndecoded
  data rms/999.0/,k0/99999999/,ntrperiod0/0/,nfft3z/0/
  equivalence (xc,cx)
  save

  nfft3=16384                                  !df=12000.0/16384 = 0.732422
!  nfft3=2048
  jstep=nsps/2                                 !Step size = half-symbol in id2()
  jstep8=nsps/16                               !Step size = half-symbol in c0()
  if(k.gt.NMAX) go to 999
!  if(k.lt.nfft3) then
  if(k.lt.2048) then
     ihsym=0
     go to 999                                 !Wait for enough samples to start
  endif
  if(nfft3.ne.nfft3z) then                     !New nfft3, compute window
     pi=4.0*atan(1.0)
     if(ntrperiod.eq.1) then                 !Compute window for nfft3 spectrun
        do i=1,nfft3
           xx=float(i-1)/(nfft3-1)
           w3(i)=0.40897 -0.5*cos(2.0*pi*xx) + 0.09103*cos(4.0*pi*xx)
!           w3(i)=0.355768 - 0.487306*cos(2.0*pi*xx) + 0.144232*cos(4.0*pi*xx) &
!                - 0.012604*cos(6.0*pi*xx)
        enddo
     else
        do i=1,nfft3
           w3(i)=2.0*(sin(i*pi/nfft3))**2         !Window for nfft3 spectrum
        enddo
     endif
     nfft3z=nfft3
  endif

  if(k.lt.k0) then                             !Start a new data block
     ja=0
     ssum=0.
     ihsym=0
     k1=0
     k8=0
     x2=0.
     if(ndiskdat.eq.0) then
        id2(k+1:)=0
        c0=0.          !This is necessary to prevent "ghosts".  Not sure why.
     endif
  endif
  k0=k
  kstep1=NFFT1
  fac=2.0/NFFT1
  nblks=(k-k1)/kstep1
  gain=10.0**(0.05*ingain)
  sq=0.
  do nblk=1,nblks
     do i=1,NFFT1
        x1(i)=gain*id2(k1+i)
     enddo
     sq=sq + dot_product(x1,x1)
! Mix at 1500 Hz, lowpass at +/-750 Hz, and downsample to 1500 Hz complex.
     x2(106:105+kstep1)=x1(1:kstep1)
     call fil3(x2,kstep1+105,c0(k8+1),n2)
     x2(1:105)=x1(kstep1-104:kstep1)   !Save 105 trailing samples
     k1=k1+kstep1
     k8=k8+kstep1/8
  enddo

  npts8=k8
  ja=ja+jstep8                         !Index of first sample
  rms=sqrt(sq/(nblks*NFFT1))
  pxdb=0.
  if(rms.gt.0.0) pxdb=20.0*log10(rms)
  if(pxdb.gt.60.0) pxdb=60.0

  fac0=0.1
  do i=0,nfft3-1                      !Copy data into cx
     j=8*ja+i-(nfft3-1)
     xc(i)=0.
!     if(j.ge.1 .and. j.le.NDMAX) cx(i)=c0(j)
     if(j.ge.1) xc(i)=fac0*id2(j)
  enddo

  if(ihsym.lt.184) ihsym=ihsym+1

!  write(69,3001) nsps,jstep8,jstep,kstep1,ja,8*ja,ihsym
!3001 format(7i9)
!  call flush(69)

!  cx(0:nfft3-1)=w3(1:nfft3)*cx(0:nfft3-1)  !Apply window w3
!  call four2a(cx,nfft3,1,1,1)              !Third FFT (forward)
  xc(0:nfft3-1)=w3(1:nfft3)*xc(0:nfft3-1)   !Apply window w3
  call four2a(xc,nfft3,1,-1,0)               !Real-to-complex FFT

  n=min(184,ihsym)
!   df3=1500.0/nfft3                    !JT9-1: 0.732 Hz = 0.42 * tone spacing
  df3=12000.0/nfft3                   !JT9-1: 0.732 Hz = 0.42 * tone spacing
!  i0=nint(-500.0/df3)
  i0=nint(1000.0/df3)
  iz=min(NSMAX,nint(1000.0/df3))
  fac=(1.0/nfft3)**2
  do i=1,iz
     j=i0+i-1
     if(j.lt.0) j=j+nfft3
     sx=fac*(real(cx(j))**2 + aimag(cx(j))**2)
     ss(n,i)=sx
     ssum(i)=ssum(i) + sx
     s(i)=sx
  enddo

999 continue

  fac00=0.35
  npct=20
  call pctile(s,iz,npct,xmed0)
  fac0=fac00/max(xmed0,0.006)
  s(1:iz)=fac0*s(1:iz)
  call pctile(ssum,iz,npct,xmed1)
  fac1=fac00/max(xmed1,0.006*ihsym)
  savg(1:iz)=fac1*ssum(1:iz)
!  savg(iz+1:iz+20)=savg(iz)
  call redsync(ss,ntrperiod,ihsym,iz,red)

  return
end subroutine symspec