subroutine fivehz

!  Called at interrupt level from the PortAudio callback routine.
!  For nspb=2048 the callback rate is nfsample/nspb = 5.38 Hz.
!  Thus, we should be able to control the timing of T/R sequence events
!  here to within about 0.2 s.

!  Do not do anything very time consuming in this routine!!
!  Disk I/O is a bad idea.  Writing to stdout (for diagnostic purposes)
!  seems to be OK.

  integer time

  parameter (NTRING=64)
  real*8 tt1(0:NTRING-1)
  logical first,txtime,filled
  integer ptt
  integer TxOKz
  real*8 fs,fsample,tt,u
  include 'gcom1.f90'
  include 'gcom2.f90'
  data first/.true./,nc0/1/,nc1/1/,nsec0/-99/
  save

  n1=time()
  n2=mod(n1,86400)
  tt=n1-n2+tsec-0.1d0*ndsec

  if(first) then
     rxdelay=0.2
     txdelay=0.4
     tlatency=1.0
     first=.false.
     iptt=0
     ntr0=-99
     rxdone=.false.
     ibuf00=-99
     ncall=-1
     u=0.05d0
     fsample=11025.d0
     mfsample=110250
     filled=.false.
     techo=2.5
  endif

  if(txdelay.lt.0.2d0) txdelay=0.2d0

! Measure average sampling frequency over a recent interval
  ncall=ncall+1
  if(ncall.eq.9) then
     ntt0=0
     ntt1=0
     tt1(ntt1)=tt
  endif
!  if(ncall.ge.10 .and. mod(ncall,2).eq.1) then
  if(ncall.ge.10) then
     ntt1=iand(ntt1+1,NTRING-1)
     tt1(ntt1)=tt
     if(ntt1.eq.NTRING-1) filled=.true.
     if(filled) ntt0=iand(ntt1+1,NTRING-1)
     if(mod(ncall,2).eq.1) then
        nd=ntt1-ntt0
        if(nd.lt.0) nd=nd+NTRING
        fs=nd*2048.d0/(tt1(ntt1)-tt1(ntt0))
        fsample=u*fs + (1.d0-u)*fsample
        mfsample=nint(10.d0*fsample)
     endif
  endif

! Echo mode
  if(mode(1:4).ne.'Echo') n3=0
  if(mode(1:4).eq.'Echo') then
     s6=mod(tsec,6.d0)
     if(lauto.eq.0 .and. (n3.le.0 .or. n3.ge.6)) go to 10

! When s6 has wrapped back to zero, start a new cycle
     if(s6.lt.s6z) then
        TxOK=0                              !Lower TxOK
        i1=ptt(nport,pttport,1,ndtr,iptt)   !Raise PTT
        call wsjtgen                        !Generate the waveform
        t1a=s6
        n3=1
!        if(ndebug.gt.0) write(*,3001) n3,s6,0.0,' Raise PTT'
!3001    format(i1,2f7.2,10x,a)
        go to 10
     endif

     if(n3.eq.1 .and. s6.gt.s6z) then
        TxOK=1
        t2a=s6                              !Save start time of Tx audio
        n3=2
!        if(ndebug.gt.0) write(*,3001) n3,s6,s6-t1a,' Start Tx audio'
        go to 10
     endif

     if(n3.eq.2 .and. s6.gt.2.6) then
        TxOK=0                              !Stop Tx audio
        t3a=s6
        n3=3
!        if(ndebug.gt.0) write(*,3001) n3,s6,s6-t2a,' Stop Tx audio'
        go to 10
     endif

!     if(n3.eq.3 .and. s6.gt.s6z) then
     if(n3.eq.3 .and. s6.gt.3.2) then
        i1=ptt(nport,pttport,0,ndtr,iptt)        !Lower PTT
        t4a=s6
        n3=4
!        if(ndebug.gt.0) write(*,3001) n3,s6,s6-t3a,' Lower PTT'
        go to 10
     endif

     if(n3.eq.4 .and. s6.ge.t2a+techo) then
        t2az=t2a
        f1z=f1
        ibuf0=ibuf
        t5a=s6
        n3=5
!        if(ndebug.gt.0) write(*,3002) n3,s6,s6-t4a,ibuf0,' Start Rx'
!3002    format(i1,2f7.2,i8,2x,a)
        go to 10
     endif

     if(n3.eq.5 .and. s6.gt.t5a+2.1) then
        ndecoding=1
        t6a=s6
        n3=6
!        if(ndebug.gt.0) write(*,3002) n3,s6,s6-t5a,ibuf,' Stop Rx'
!        if(ndebug.gt.0) write(*,*)
        go to 10
     endif

10   s6z=s6
     go to 900
  endif
!End of Echo mode

  if(trperiod.le.0) trperiod=30
  tx1=0.0                              !Time to start a TX sequence
  tx2=trperiod-(tlatency+txdelay)      !Time to turn TX off
  if(mode(1:4).eq.'JT65' .or. mode(1:3).eq.'JT4') then
     if(nwave.lt.126*4096) nwave=126*4096
     tx2=txdelay + nwave/11025.0
     if(tx2.gt.(trperiod-2.0)) tx2=trperiod-tlatency-1.0
  endif

  if(TxFirst.eq.0 .and. ntxboth.eq.0) then
     tx1=tx1+trperiod
     tx2=tx2+trperiod
  endif

  t=mod(Tsec,2.d0*trperiod)
  if(ntxboth.eq.0) then
     txtime = t.ge.tx1 .and. t.lt.tx2
  else
     txtime = (t.ge.tx1 .and. t.lt.tx2) .or.                     &
          (t.ge.tx1+trperiod .and. t.lt.tx2+trperiod)
  endif

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
          mute.eq.0) i1=ptt(nport,pttport,1,ndtr,iptt)        !Raise PTT
     if(.not.txtime .or. mute.eq.1) TxOK=0                    !Lower TxOK
  else
     if(mantx.eq.1 .and. iptt.eq.0 .and.      &
          mute.eq.0) i2=ptt(nport,pttport,1,ndtr,iptt)        !Raise PTT
     if(mantx.eq.0 .or. mute.eq.1) TxOK=0                     !Lower TxOK
  endif

! Calculate Tx waveform as needed
  if((iptt.eq.1 .and. iptt0.eq.0) .or. nrestart.eq.1) then
     call wsjtgen
     nrestart=0
  endif

! If PTT was just raised, start a countdown for raising TxOK:
  nc1a=txdelay/0.18576
  if(nc1a.lt.2) nc1a=2
  if(iptt.eq.1 .and. iptt0.eq.0) nc1=-nc1a-1
  if(nc1.le.0) nc1=nc1+1
  if(nc1.eq.0) TxOK=1                               ! We are transmitting

! If TxOK was just lowered, start a countdown for lowering PTT:
  nc0a=(tlatency+txdelay)/0.18576
  if(nc0a.lt.5) nc0a=5
  if(TxOK.eq.0 .and. TxOKz.eq.1 .and. iptt.eq.1) nc0=-nc0a-1
  if(nc0.le.0) nc0=nc0+1
  if(nc0.eq.0) i3=ptt(nport,pttport,0,ndtr,iptt)

  if(iptt.eq.0 .and.TxOK.eq.0) then
     sending="                      "
     sendingsh=0
  endif

  nbufs=i1+i2+i3                             !Silence g95 warning
  nbufs=ibuf-ibuf0
  if(nbufs.lt.0) nbufs=nbufs+1024
  tdata=nbufs*2048.0/11025.0

  if((mode(1:4).eq.'JT65' .or. mode(1:3).eq.'JT4'                      &
       .or. mode(1:2).eq.'CW') .and. monitoring.eq.1                   &
       .and. tdata.gt.float(ntdecode) .and. ibuf0.ne.ibuf00) then
     rxdone=.true.
     ibuf00=ibuf0
  endif

  iptt0=iptt
  TxOKz=TxOK
  ntr0=ntr

900 return
end subroutine fivehz

subroutine fivehztx

!  Called at interrupt level from the PortAudio output callback.

  parameter (NTRING=64)
  real*8 tt1(0:NTRING-1)
  logical first,filled
  real*8 fs,fsample,tt,u
  include 'gcom1.f90'
  data first/.true./
  save

  n1=time()
  n2=mod(n1,86400)
  tt=n1-n2+tsec-0.1d0*ndsec

  if(first) then
     first=.false.
     ncall=-1
     fsample=11025.d0
     u=0.05d0
     mfsample2=110250
     filled=.false.
  endif

! Measure average sampling frequency over a recent interval
  ncall=ncall+1
  if(ncall.eq.9) then
     ntt0=0
     ntt1=0
     tt1(ntt1)=tt
  endif
  if(ncall.ge.10) then
     ntt1=iand(ntt1+1,NTRING-1)
     tt1(ntt1)=tt
     if(ntt1.eq.NTRING-1) filled=.true.
     if(filled) ntt0=iand(ntt1+1,NTRING-1)
     if(mod(ncall,2).eq.1) then
        nd=ntt1-ntt0
        if(nd.lt.0) nd=nd+NTRING
        fs=nd*2048.d0/(tt1(ntt1)-tt1(ntt0))
        fsample=u*fs + (1.d0-u)*fsample
        mfsample2=nint(10.d0*fsample)
     endif
  endif

  return
end subroutine fivehztx

subroutine addnoise(n)
  integer*2 n
  real*8 txsnrdb0
  include 'gcom1.f90'
  data idum/0/
  save

  if(txsnrdb.gt.40.0) return
  if(txsnrdb.ne.txsnrdb0) then
     snr=10.0**(0.05*(txsnrdb-1))
     fac=3000.0
     if(snr.gt.1.0) fac=3000.0/snr
     txsnrdb0=txsnrdb
  endif
  i=fac*(gran(idum) + n*snr/32768.0)
  if(i>32767) i=32767;
  if(i<-32767) i=-32767;
  n=i

  return
end subroutine addnoise
