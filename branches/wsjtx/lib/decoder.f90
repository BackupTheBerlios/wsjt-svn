subroutine decoder(ntrSeconds,ndepth,nRxLog,c00)

! Decoder for JT9.

  parameter (NMAX=1800*12000)        !Total sample intervals per 30 minutes
  parameter (NDMAX=1800*1500)        !Sample intervals at 1500 Hz rate
  parameter (NSMAX=22000)            !Max length of saved spectra
  character*22 msg
  character*33 line
  character*80 fmt
  real*4 ccfred(NSMAX)
  integer*1 i1SoftSymbols(207)
  integer*2 id2
  complex c0(NDMAX),c00(NDMAX)
  common/jt9com/ss0(184,NSMAX),savg(NSMAX),id2(NMAX),nutc0,ndiskdat,    &
       ntr,nfqso,newdat,npts80,nfb,ntol,kin,nsynced,ndecoded
  common/jt9comB/ss(184,NSMAX),c0
  logical first
  data first/.true./
  save

  if(newdat.ne.0) then
     ss=ss0
     c0=c00
     nutc=nutc0
     npts8=npts80
  endif

  ntrMinutes=ntrSeconds/60
  newdat=1
  nsynced=0
  ndecoded=0
  limit=1000
  if(ndepth.ge.2) limit=20000
  if(ndepth.ge.3) limit=100000

  nsps=0
  if(ntrMinutes.eq.1) then
     nsps=6912
     df3=1500.0/2048.0
     fmt='(i4.4,i4,i5,f6.1,f8.0,f6.2,3x,a22)'
  else if(ntrMinutes.eq.2) then
     nsps=15360
     df3=1500.0/2048.0
     fmt='(i4.4,i4,i5,f6.1,f8.1,f6.2,3x,a22)'
  else if(ntrMinutes.eq.5) then
     nsps=40960
     df3=1500.0/6144.0
     fmt='(i4.4,i4,i5,f6.1,f8.1,f6.2,3x,a22)'
  else if(ntrMinutes.eq.10) then
     nsps=82944
     df3=1500.0/12288.0
     fmt='(i4.4,i4,i5,f6.1,f8.2,f6.2,3x,a22)'
  else if(ntrMinutes.eq.30) then
     nsps=252000
     df3=1500.0/32768.0
     fmt='(i4.4,i4,i5,f6.1,f8.2,f6.2,3x,a22)'
  endif
  if(nsps.eq.0) stop 'Error: bad TRperiod'    !Better: return an error code###

  kstep=nsps/2
  tstep=kstep/12000.0

  call sync9(ss,tstep,df3,ntol,nfqso,ccfred,ia,ib,ipk)  ! Get sync, approx freq

  open(13,file='decoded.txt',status='unknown')
  rewind 13
  if(first) then
     open(14,file='wsjtx_rx.log',status='unknown',position='append')
     first=.false.
  endif
  if(iand(nRxLog,2).ne.0) rewind 14
  if(iand(nRxLog,1).ne.0) then
! Write date and time to lu 14     
  endif

  fgood=0.
  df8=1500.0/(nsps/8)
  sbest=0.
  do i=ia,ib
     f=(i-1)*df3
     if((i.eq.ipk .or. ccfred(i).ge.3.0) .and. f.gt.fgood+10.0*df8) then
        call spec9(c0,npts8,nsps,f,fpk,xdt,snr,i1SoftSymbols)
        call decode9(i1SoftSymbols,limit,nlim,msg)
        sync=ccfred(i) - 2.0
        if(sync.lt.0.0) sync=0.0
        nsync=sync
        if(nsync.gt.10) nsync=10
        nsnr=nint(snr)
        drift=0.0

        if(ccfred(i).gt.sbest .and. fgood.eq.0.0) then
           sbest=ccfred(i)
           write(line,fmt) nutc,nsync,nsnr,xdt,1000.0+fpk,drift
           if(nsync.gt.0) nsynced=1
        endif

        if(msg.ne.'                      ') then
           write(13,fmt) nutc,nsync,nsnr,xdt,1000.0+fpk,width,msg
           write(14,fmt) nutc,nsync,nsnr,xdt,1000.0+fpk,width,msg
           fgood=f
           nsynced=1
           ndecoded=1
        endif
     endif
  enddo

  if(fgood.eq.0.0) then
     write(13,1020) line
     write(14,1020) line
1020 format(a33)
  endif

  call flush(13)
  call flush(14)
  close(13)

  return
end subroutine decoder
