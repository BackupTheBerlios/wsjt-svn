subroutine genms(msg,iwave,nwave)

! Generate a JTMS wavefile.

  parameter (NMAX=30*11025)     !Max length of wave file
  character*28 msg              !Message to be generated
  character cc*64
  integer sent(168)
  real*8 dt,phi,f,f0,dfgen,dphi,twopi
  integer*2 iwave(NMAX)         !Generated wave file
  complex cw
  common/mscom/cw(48,0:63)
  data idum/-1/
!           1234567890123456789012345678901234567890123456789012345678901234
  data cc/' 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ./?-                      @'/
  save idum

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
        if(msg(j:j).eq.cc(i:i)) go to 5
     enddo
5    do n=5,0,-1                            !Each character gets 6 bits
        k=k+1
        sent(k)=iand(1,ishft(i-1,-n))
     enddo
  enddo
  nsym=k

 ! Set up necessary constants
  twopi=8.d0*atan(1.d0)
  nsps=8
  dt=1.d0/11025.d0
  f0=11025.d0/nsps                               ! 1575.0 Hz
  dfgen=11025.d0/(2*nsps)                        !  787.5 Hz
  t=0.d0
  k=0
  phi=0.d0
  nrpt=NMAX/(nsym*nsps)

!  write(*,3001) (sent(k),k=1,nsym)
!3001 format(10(1x,6i1))

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

  return
end subroutine genms
