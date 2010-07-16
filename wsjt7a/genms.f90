subroutine genms(msg,txsnrdb,iwave,nwave)

! Generate a JTMS wavefile.

  parameter (NMAX=30*11025)     !Max length of wave file
  character*28 msg              !Message to be generated
  character cc*64
  integer imsg(28)
  integer sent(168)
  real*8 txsnrdb,t
  real*8 dt,phi,f,f0,dfgen,dphi,twopi
  integer*2 iwave(NMAX)         !Generated wave file
  complex cw
  common/mscom/cw(48,0:63)
  data idum/-1/
!           1234567890123456789012345678901234567890123456789012345678901234
  data cc/' 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ./?-                      @'/
  save idum

  sig=sqrt(10.0**(0.1*txsnrdb))
  do i=28,1,-1                                 !Find user's message length
     if(msg(i:i).ne.' ') go to 1
  enddo
1 msglen=i

! Convert message to a bit sequence, 6 bits per character
! Start with two blank characters (all 0's)
  k=12
  sent=0
  do j=1,msglen
     do i=1,64
        if(msg(j:j).eq.cc(i:i)) then
           imsg(j)=i-1
           go to 5
        endif
     enddo
5    do n=5,0,1                            !Each has 6 bits, 6*nsps samples
        k=k+1
        sent(k)=iand(1,ishft(i,-n))
     enddo
  enddo
  nsym=k
  print*,nsym

 ! Set up necessary constants
  nsps=8
  dt=1.d0/11025.d0
  f0=1500.d0
  dfgen=750.d0
  t=0.d0
  k=0
  phi=0.d0
  nrpt=30.0*12000.0/(nsym*nsps)
  do irpt=1,nrpt
     do j=1,nsym
        if(sent(j).eq.1) then
           f=f0 + 0.5d0*dfgen
        else
           f=f0 - 0.5d0*dfgen
        endif
        dphi=twopi*f*dt
        do i=1,nsps
           k=k+1
           phi=phi+dphi
           iwave(k)=nint(32767.0*sin(phi))
        enddo
     enddo
  enddo

900 iwave(k+1:)=0
  nwave=k

  if(txsnrdb.lt.40.d0) then
! ###  Make some pings (for tests only) ###
     do i=1,nwave
        iping=i/(3*12000)
        if(iping.ne.iping0) then
           ip=mod(iping,3)
           w=0.05*(ip+1)
           ig=(iping-1)/3
           amp=sqrt((3.0-ig)/3.0)
           t0=dt*(iping+0.5)*(3*12000)
           iping0=iping
        endif
        t=(i*dt-t0)/w
        if(t.lt.0.d0 .and. t.lt.10.d0) then
           fac=0.
        else
           fac=2.718*t*dexp(-t)
        endif
        iwave(i)=nint(fac*amp*iwave(i))
     enddo
  endif

  return
end subroutine genms
