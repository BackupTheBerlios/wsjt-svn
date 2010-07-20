subroutine decodems(dat,npts,cfile6,t2,mswidth,ndb,nrpt,Nfreeze,       &
     DFTolerance,MouseDF)

! Decode a JTMS ping

  parameter (NZ=30*11025)
  real dat(npts)                        !Raw data
  integer DFTolerance
  character*6 cfile6
  complex cdat(NZ)                      !Analytic form of signal
  real s(NZ)
  real fs(56)
  real sm(0:63)
  complex c(NZ)
  complex cw(56,0:63)                   !Complex waveforms for codewords
  complex cwb(56)                       !Complex waveform for 'space'
  complex z
  logical first
  character msg*400,msg28*28
  character*90 line
  character cc*64
!                    1         2         3         4         5         6
!          0123456789012345678901234567890123456789012345678901234567890123
  data cc/'0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ./?-                 _     @'/
  data first/.true./
  save first,smax,cw,cwb           !Why is this needed for save?  But it is!
  common/ccom/nline,tping(100),line(100)

  if(first) call setupms(cw,cwb)        !Calculate waveforms for codewords
  first=.false.

  nsps=8
  f0=1155.46875
  xn=log(float(npts))/log(2.0)
  n=xn
  if(xn-n .gt.0.001) n=n+1
  nfft1=2**n
  nh1=nfft1/2
  nq1=nfft1/4
  df1=11025.0/nfft1
  ia=nint(11025.0/(3.0*nsps*df1))
  ib=nint(dftolerance/df1)
  i0=nint(f0/df1)

  call analytic(dat,npts,nfft1,s,cdat)        !Convert to analytic signal

  fac=1.0/(nfft1**2)
  do i=1,npts
     c(i)=fac*cdat(i)**2
  enddo
  c(npts+1:nfft1)=0.
  call four2a(c,nfft1,1,-1,1)

! In the "doubled-frequencies" spectrum of squared cdat:
  fa=2.0*(f0-400)
  fb=2.0*(f0+400)
  if(NFreeze.gt.0) then
     fa=2.0*(f0+MouseDF-DFtolerance)
     fb=2.0*(f0+MouseDF+DFtolerance)
  endif  
  ja=nint(fa/df1)
  jb=nint(fb/df1)
  jd=nint(1378.125/df1)

  smax=0.
  do j=ja,jb
     ss=real(c(j))**2 + aimag(c(j))**2 + real(c(j+jd))**2 + aimag(c(j+jd))**2
     if(ss.gt.smax) then
        smax=ss
        fpk=(j-1)*df1
     endif
  enddo
  dfx=0.5*fpk-f0

! Should have a test here to reject non-JTMS signals

! We know DF, now find character sync.
  fs=0.
  do i1=1,npts-55
     z=0.
     do i=1,56
        z=z + cdat(i+i1-1)*conjg(cwb(i))
     enddo
     z=0.001*z
     j=mod(i1-1,56)+1
     ss=abs(z)
     fs(j)=fs(j)+ss
  enddo

  smax=0.
  do j=1,56
     if(fs(j).gt.smax) then
        smax=fs(j)
        i1=j
     endif
  enddo

!  rewind 51

  msg=' '
  nchar=(npts-55-i1)/56
  if(nchar.gt.400) nchar=400
  do j=1,nchar
     ia=i1 + (j-1)*56
     smax=0.
     do k=0,63
        z=0.
        do i=1,56
           z=z + cdat(ia+i)*conjg(cw(i,k))
        enddo
        ss=abs(z)
        sm(k)=ss
        if(ss.gt.smax) then
           smax=ss
           phapk=atan2(aimag(z),real(z))
           kpk=k
        endif
     enddo
     sm(kpk)=0.
     smax2=0.
     do k=0,63
        smax2=max(smax2,sm(k))
     enddo
     if(kpk.lt.1) then
        kpk=64
     endif
     msg(j:j)=cc(kpk:kpk)
     if(kpk.eq.58) msg(j:j)=' '
     if(smax/smax2.lt.1.05) msg(j:j)=' '               !Threshold test
!     write(51,3007) j,smax,phapk,phapk+6.283185307
!3007 format(i5,3f12.3)
  enddo
!  call flush(51)

  ia=max(1,nchar/3)
  ib=min(ia+27,nchar)
  msg28=msg(ia:ib)
  ndf=nint(dfx)

  if(nline.le.99) nline=nline+1
  tping(nline)=t2
  write(line(nline),1110) cfile6,t2,mswidth,ndb,nrpt,ndf,msg28
1110 format(a6,f5.1,i5,i3,1x,i2.2,i5,5x,a28,f8.1,f6.2,i3,a2)

  return
end subroutine decodems
