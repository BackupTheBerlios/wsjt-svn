subroutine diana(dat,npts,cfile6,MinSigdB,DFTolerance,NFreeze,       &
     MouseDF,ccfblue,ccfred)

! Decode an ISCAT_2 signal

  parameter (NMAX=512*1024)
  parameter (NSZ=646)
  real dat(NMAX)                          !Raw signal, 30 s at 11025 sps
  character cfile6*6                      !File time
  character c42*42
  character msg*28
  real s0(1024,NSZ)
  real fs0(1024,108)                       !108 = 96 + 3*4
  real ccfblue(-5:540)
  real psavg(1024)         !Average spectrum of the whole file
  integer dftolerance
  data nsps/2048/,nsync/4/,nlen/2/,ndat/18/

! Define some constants
  nsym=npts/nsps                      !Total symbol intervals in file
  nblk=nsync+nlen+ndat                !Frame size
  nfft=4096
  nq=nfft/4
  df=11025.0/nfft
  kstep=nsps/4

! Compute spectra at 1/4-symbol steps, s0, and folded spectra, fs0
  call specdiana(dat,npts,s0,jsym,fs0)

  call syncdiana(fs0,kstep,nfreeze,mousedf,dftolerance,xsync,     &
     ipk,jpk,dfx,dtx,ccfblue)

  call lendiana(fs0,ipk,jpk,msglen)
  
  msg=' '
  if(msglen.gt.0) call decdiana(s0,jsym,ipk,jpk,msglen,msg,snrx)
  nsnr=nint(snrx)
  if(nsnr.lt.-28) then
     nsnr=-28
     msg=' '
  endif
  jsync=xsync
  jdf=nint(dfx)
  nwidth=0

  call cs_lock('iscat')
!  write(*,1020) cfile6,jsync,nsnr,dtx,jdf,nwidth,msg
  write(11,1020) cfile6,jsync,nsnr,dtx,jdf,nwidth,msg
  write(21,1020) cfile6,jsync,nsnr,dtx,jdf,nwidth,msg
1020 format(a6,i3,i5,f5.1,i5,i3,7x,a28)
  call cs_unlock

  return
end subroutine diana
