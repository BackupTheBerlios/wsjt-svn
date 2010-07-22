subroutine genms(msg,iwave,nwave)

! Generate a JTMS wavefile.

  parameter (NMAX=30*11025)     !Max length of wave file
  character*28 msg              !Message to be generated
  character cc*64
  integer sent(196)
  real*8 dt,phi,f,f0,dfgen,dphi,twopi,foffset
  integer*2 iwave(NMAX)         !Generated wave file
!                   1         2         3         4         5         6
!          0123456789012345678901234567890123456789012345678901234567890123
  data cc/'0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ./?-                 _     @'/

  do i=28,1,-1                                 !Find user's message length
     if(msg(i:i).ne.' ') go to 1
  enddo
1 msglen=i+1                                   !Add one for space at EOM
  if(msglen.gt.28) msglen=28

! Convert message to a bit sequence, 7 bits per character (6 + even parity)
  sent=0
  k=0
  do j=1,msglen
     if(msg(j:j).eq.' ') then
        i=58
        go to 5
     else
        do i=1,64
           if(msg(j:j).eq.cc(i:i)) go to 5
        enddo
     endif
5    m=0
     do n=5,0,-1                            !Each character gets 6 bits
        k=k+1
        sent(k)=iand(1,ishft(i-1,-n))
        m=m+sent(k)
     enddo
     k=k+1
     sent(k)=iand(m,1)                      !Insert parity bit
  enddo
  nsym=k

 ! Set up necessary constants
  twopi=8.d0*atan(1.d0)
  nsps=8
  dt=1.d0/11025.d0
  f0=11025.d0/nsps                               ! 1575.0 Hz
  dfgen=0.5d0*f0                                 !  787.5 Hz
  foffset=1500.d0 - f0
  t=0.d0
  k=0
  phi=0.d0
  nrpt=NMAX/(nsym*nsps)

  do irpt=1,nrpt
     do j=1,nsym
        if(sent(j).eq.1) then
           f=f0 + 0.5d0*dfgen + foffset
        else
           f=f0 - 0.5d0*dfgen + foffset
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
