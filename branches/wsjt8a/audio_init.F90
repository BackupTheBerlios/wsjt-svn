subroutine audio_init(ndin,ndout)
!f2py threadsafe

#ifdef CVF
  use dfmt
  integer Thread1,Thread2,ABOVE_NORMAL_PRIORITY_CLASS
  parameter (ABOVE_NORMAL_PRIORITY_CLASS=32768)
  external a2d,decode1
#endif
  integer start_threads

  include 'gcom1.f90'
  include 'gcom2.f90'

  nmode=1
  if(mode(1:4).eq.'JT64') then
     nmode=2
     if(mode(5:5).eq.'A') mode64=1
     if(mode(5:5).eq.'B') mode64=2
     if(mode(5:5).eq.'C') mode64=4
  endif
  if(mode(1:4).eq.'Echo') nmode=3
  if(mode(1:5).eq.'ISCAT') nmode=4
  if(mode(1:3).eq.'JT8') nmode=5

  ndevin=ndin
  ndevout=ndout
  TxOK=0
  Transmitting=0
  nfsample=12000
  nspb=1024
  nbufs=2048
  nmax=nbufs*nspb
  nwave=60*nfsample
  ngo=1
  f0=800.0
  do i=1,nwave
     iwave(i)=nint(32767.0*sin(6.283185307*i*f0/nfsample))
  enddo

#ifdef CVF
!  Priority classes (for processes):
!     IDLE_PRIORITY_CLASS               64
!     NORMAL_PRIORITY_CLASS             32
!     HIGH_PRIORITY_CLASS              128

!  Priority definitions (for threads):
!     THREAD_PRIORITY_IDLE             -15
!     THREAD_PRIORITY_LOWEST            -2
!     THREAD_PRIORITY_BELOW_NORMAL      -1
!     THREAD_PRIORITY_NORMAL             0
!     THREAD_PRIORITY_ABOVE_NORMAL       1
!     THREAD_PRIORITY_HIGHEST            2
!     THREAD_PRIORITY_TIME_CRITICAL     15

  npri=NORMAL_PRIORITY_CLASS
  if(nhighpri.ne.0) npri=ABOVE_NORMAL_PRIORITY_CLASS
  m0=SetPriorityClass(GetCurrentProcess(),npri)

! Start a thread for doing A/D and D/A with sound card.
  Thread1=CreateThread(0,0,a2d,0,CREATE_SUSPENDED,id1)
  m1=SetThreadPriority(Thread1,THREAD_PRIORITY_ABOVE_NORMAL)
  m2=ResumeThread(Thread1)

! Start a thread for background decoding.
  Thread2=CreateThread(0,0,decode1,0,CREATE_SUSPENDED,id2)
  m3=SetThreadPriority(Thread2,THREAD_PRIORITY_BELOW_NORMAL)
  m4=ResumeThread(Thread2)
#else
  ierr=start_threads()
#endif

  return
end subroutine audio_init
