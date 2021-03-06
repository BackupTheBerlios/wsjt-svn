subroutine fivehz

!  Called at interrupt level from the PortAudio callback routine.
!  For nspb=2048 the callback rate is nfsample/nspb = 5.38 Hz.
!  Thus, we should be able to control the timing of T/R sequence events
!  here to within about 0.2 s.

!  Do not do anything very time consuming in this routine!!
!  Disk I/O is a bad idea.  Writing to stdout (for diagnostic purposes)
!  seems to be OK.

#ifdef Win32
  use dflib
  use dfport
#endif

  real*8 tstart,tstop,t60
  logical first,txtime,debug
  integer ptt
  integer TxOKz
  real*8 fs,fsample,tt,tt0,u
  include 'gcom1.f90'
  include 'gcom2.f90'
  data first/.true./,nc0/1/,nc1/1/
  save

  n1=time()
  n2=mod(n1,86400)
  tt=n1-n2+tsec-0.1d0*ndsec

  if(first) then
     rxdelay=0.2
     txdelay=0.2
     tlatency=1.0
     first=.false.
     iptt=0
     ntr0=-99
     debug=.false.
     rxdone=.false.
     ibuf00=-99
     ncall=-1
     tt0=tt
     u=0.05d0
     fsample=11025.d0
     maxms=0
     mfsample=110250
  endif

  if(txdelay.lt.0.2d0) txdelay=0.2d0

! Measure average sampling frequency over a recent interval

  ncall=ncall+1
  if(ncall.gt.0) then
     fs=ncall*2048.d0/(tt-tt0)
     fsample=u*fs + (1.d0-u)*fsample
     mfsample=nint(10.d0*fsample)
  endif

  if(trperiod.le.0) trperiod=30
  tx1=0.0                              !Time to start a TX sequence
  tx2=trperiod-(tlatency+txdelay)      !Time to turn TX off
  if(mode(1:4).eq.'JT65') then
     if(nwave.lt.126*4096) nwave=126*4096
     tx2=nwave/11025.0
  endif

  if(TxFirst.eq.0) then
     tx1=tx1+trperiod
     tx2=tx2+trperiod
  endif

  t=mod(Tsec,2.d0*trperiod)
  txtime = t.ge.tx1 .and. t.lt.tx2

! If we're transmitting, freeze the input buffer pointers where they were.
  receiving=1
  if(((txtime .and. (lauto.eq.1)) .or. TxOK.eq.1 .or. transmitting.eq.1) & 
       .and. (mute.eq.0)) then
     receiving=0
     ibuf=ibuf000
     iwrite=iwrite000
  endif
  ibuf000=ibuf
  iwrite000=iwrite

  nsec=Tsec

  ntr=mod(nsec/trperiod,2)             !ntr=0 in 1st sequence, 1 in 2nd

  if(ntr.ne.ntr0) then
     ibuf0=ibuf                        !Start of new sequence, save ibuf
     ibuf0=ibuf0-3
     if(ibuf0.lt.1) ibuf0=ibuf0+1024
!     if(mode(1:4).ne.'JT65') then
!        ibuf0=ibuf0+3                  !So we don't copy our own Tx
!        if(ibuf0.gt.1024) ibuf0=ibuf0-1024
!     endif
     ntime=time()                      !Save start time
     if(mantx.eq.1 .and. iptt.eq.1) then
        mantx=0
        TxOK=0
     endif
  endif

! Switch PTT line and TxOK appropriately
  if(lauto.eq.1) then
     if(txtime .and. iptt.eq.0 .and.          &
          mute.eq.0) i1=ptt(nport,1,iptt)                !Raise PTT
     if(.not.txtime .or. mute.eq.1) TxOK=0               !Lower TxOK
  else
     if(mantx.eq.1 .and. iptt.eq.0 .and.      &
          mute.eq.0) i2=ptt(nport,1,iptt)                !Raise PTT
     if(mantx.eq.0 .or. mute.eq.1) TxOK=0                !Lower TxOK
  endif

! Calculate Tx waveform as needed
  if((iptt.eq.1 .and. iptt0.eq.0) .or. nrestart.eq.1) then
     call wsjtgen
     nrestart=0
  endif

! If PTT was just raised, start a countdown for raising TxOK:
  nc1a=txdelay/0.18576
  if(nc1a.lt.2) nc1a=2
  if(mode(1:4).eq.'JT65') nc1a=2                    !No extra delay for JT65
  if(iptt.eq.1 .and. iptt0.eq.0) nc1=-nc1a
  if(nc1.le.0) nc1=nc1+1
  if(nc1.eq.0) TxOK=1                               ! We are transmitting

! If TxOK was just lowered, start a countdown for lowering PTT:
  nc0a=txdelay/0.18576
  if(nc0a.lt.4) nc0a=4
  if(TxOK.eq.0 .and. TxOKz.eq.1 .and. iptt.eq.1) nc0=-nc0a
  if(nc0.le.0) nc0=nc0+1
  if(nc0.eq.0) i3=ptt(nport,0,iptt)

  if(iptt.eq.0 .and.TxOK.eq.0) then
     sending="                      "
     sendingsh=0
  endif

  nbufs=ibuf-ibuf0
  if(nbufs.lt.0) nbufs=nbufs+1024
  tdata=nbufs*2048.0/11025.0
  if(mode(1:4).eq.'JT65' .and. monitoring.eq.1 .and. tdata.gt.53.0    &
       .and. ibuf0.ne.ibuf00) then
     rxdone=.true.
     ibuf00=ibuf0
  endif

!  if(ndebug.ne.0) then
!     t60=mod(tsec,60.d0)
!     if(iptt.ne.iptt0) then
!        if(iptt.eq.1) tstart=tsec
!        if(iptt.eq.0) write(*,1101) tsec-tstop,t60
!1101    format('Delay from TxOFF to PTT was',f6.2,' s at t=',f6.2)
!     endif
!     if(TxOK.ne.TxOKz) then
!        if(TxOK.eq.0) tstop=tsec
!        if(TxOK.eq.1) write(*,1102) tsec-tstart,t60
!1102    format('Delay from PTT to TxON was ',f6.2,' s at t=',f6.2)
!     endif
!  endif

  iptt0=iptt
  TxOKz=TxOK
  ntr0=ntr

  return
end subroutine fivehz

subroutine fivehztx

!  Called at interrupt level from the PortAudio output callback.

#ifdef Win32
  use dflib
  use dfport
#endif

  logical first
  real*8 fs,fsample,tt,tt0,u
  include 'gcom1.f90'
  data first/.true./
  save

  n1=time()
  n2=mod(n1,86400)
  tt=n1-n2+tsec-0.1d0*ndsec

  if(first) then
     first=.false.
     ncall=-1
     tt0=tt
     fsample=11025.d0
     nsec0=-999
     u=0.05d0
     mfsample2=110250
  endif

  ncall=ncall+1
  if(ncall.gt.0) then
     fs=ncall*2048.d0/(tt-tt0)
     fsample=u*fs + (1.d0-u)*fsample
     mfsample2=nint(10.d0*fsample)
  endif

  return
end subroutine fivehztx
