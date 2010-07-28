subroutine jtms(dat,npts,cfile6,t2,mswidth,ndb,nrpt,Nfreeze,       &
     DFTolerance,MouseDF,mycall,hiscall)

! Decode a JTMS ping

  parameter (NZ=30*11025)
  real dat(npts)                        !Raw data
  complex cdat(NZ)                      !Analytic form of signal
  character*6 cfile6                    !FileID
  integer DFTolerance
  character*12 mycall,hiscall
  real s(NZ)                            !Power spectrum
  real sm(0:63)
  real s2(0:63,400)
  real r(60000)
  complex cw(56,0:63)                   !Complex waveforms for codewords
  complex cwb(56)                       !Complex waveform for 'space'
  logical first
  character msg*400,msg29*29,frag*29
  character*90 line
  common/ccom/nline,tping(100),line(100)
  data first/.true./
  data nsum/0/,nrec/0/
  save nsum,nrec
  save first,cw,cwb
  save cdat

  nrec=nrec+1
  if(first) call setupms(cw,cwb)        !Calculate waveforms for codewords
  first=.false.

  nsps=8                                !Samples per symbol
  f0=1155.46875                         !Nominal frequency for bit=0
  xn=log(float(npts))/log(2.0)
  n=xn
  if(xn-n .gt.0.001) n=n+1
  nfft1=2**n
  df1=11025.0/nfft1

  call analytic(dat,npts,nfft1,s,cdat)        !Convert to analytic signal

  call msdf(cdat,npts,nfft1,f0,mousedf,dftolerance,dfx,ferr)    !Get DF

  if(abs(ferr).gt.0.002) go to 900           !Reject non-JTMS signals
  call tweak1(cdat,npts,-dfx,cdat)           !Mix to standard frequency

! DF is known, now establish character sync.

  call syncms(cdat,npts,cwb,r,i1)

  call lenms(r,npts,msglen)

  s2=0.
  nchar=(npts-55-i1)/56
  if(nchar.gt.400) nchar=400

  frag=' '//mycall
  call searchms(cdat(i1),npts-i1,frag,ndi1,rmax1)
  frag=' '//hiscall
  call searchms(cdat(i1),npts-i1,frag,ndi2,rmax2)
  frag=' CQ'
  call searchms(cdat(i1),npts-i1,frag,ndi3,rmax3)

  ndi=99
  if(max(rmax1,rmax2,rmax3).ge.0.6) then
     if(max(rmax1,rmax2,rmax3).eq.rmax1 .and. abs(ndi1).le.5) ndi=ndi1
     if(max(rmax1,rmax2,rmax3).eq.rmax2 .and. abs(ndi2).le.5) ndi=ndi2
     if(max(rmax1,rmax2,rmax3).eq.rmax3 .and. abs(ndi3).le.5) ndi=ndi3
     if(abs(ndi).le.5) i1=i1+ndi
  endif

  call decodems(cdat,npts,cw,i1,nchar,s2,sm,msg)

  ia=max(1,nchar/3)
  ib=min(ia+27,nchar)
  msg29=msg(ia:ib)
  ndf=nint(dfx)

  if(msglen.eq.0 .or. nchar.lt.max(20,2*msglen)) then

     if(nline.le.99) nline=nline+1
     tping(nline)=t2
     call cs_lock('decodems')
     if(abs(ndi).le.5) then
        write(line(nline),1110) cfile6,t2,mswidth,ndb,nrpt,ndf,msg29,ndi
1110    format(a6,f5.1,i5,i3,1x,i2.2,i5,5x,a29,12x,i3)
     else
        write(line(nline),1110) cfile6,t2,mswidth,ndb,nrpt,ndf,msg29
     endif
     call cs_unlock
  else if(msglen.gt.0) then

     call foldms(s2,msglen,nchar,mycall,msg,msg29)

     if(nline.le.99) nline=nline+1
     tping(nline)=t2
     call cs_lock('decodems')
     if(abs(ndi).le.5) then
        write(line(nline),1120) cfile6,t2,mswidth,ndb,nrpt,ndf,msg29,msglen,ndi
1120    format(a6,f5.1,i5,i3,1x,i2.2,i5,5x,a29,8x,i3,'*',i3)
     else
        write(line(nline),1120) cfile6,t2,mswidth,ndb,nrpt,ndf,msg29,msglen
     endif
!     write(*,1130) nrec,cfile6,t2,mswidth,ndb,nrpt,ndf,msg29,msglen
!1130 format(i3,1x,a6,f5.1,i5,i3,1x,i2.2,i5,5x,a29,10x,i5'*')
     call cs_unlock
   endif

900 continue
  call flush(52)

  return
end subroutine jtms

