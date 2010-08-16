subroutine gendiana(msg,msglen,samfac,iwave,nwave,msgsent,sendingsh)

! Generate waveform for Diana mode.

  parameter (NMAX=30*11025,NSZ=126,NSPS=2048)
  character msg*28,msgsent*28
  integer*2 iwave(NMAX)
  integer imsg(28)
  integer itone(NSZ)
  character c*42
  real*8 twopi,dt,f0,f,df,pha,dpha,samfac
  integer isync(4)                              !Sync pattern
  integer irpt(31)
  integer sendingsh
  data isync/8,16,32,24/
  data nsync/4/,nlen/2/,ndat/18/
  data c/'0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ /.?+-'/
  data irpt/1,2,3,4,6,7,8,9,11,12,13,14,15,16,17,18,19,21,22,23,24,     &
       25,26,27,28,29,31,32,33,34,35/

  twopi=8.d0*atan(1.d0)
  df=11025.d0/NSPS
  dt=1.0/(samfac*11025.0)
  f0=236*df
  nsym=126

  nblk=nsync+nlen+ndat
  k=0
  kk=1
  do i=1,msglen                        !Define tone sequence for user message
     imsg(i)=36
     do j=1,42
        if(msg(i:i).eq.c(j:j)) imsg(i)=j-1
     enddo
  enddo

  do i=1,nsym                          !Total symbols in whole transmission
     j=mod(i-1,nblk)+1
     if(j.le.nsync) then
        itone(i)=isync(j)
     else if(j.gt.nsync .and. j.le.nsync+nlen) then
        itone(i)=msglen
        if(j.ge.nsync+2) then
           n=msglen + 5*(j-nsync-1)
           if(n.gt.41) n=n-42
           itone(i)=n
        endif
     else
        k=k+1
        kk=mod(k-1,msglen)+1
        itone(i)=imsg(kk)
     endif
  enddo
  msgsent=msg

  k=0
  pha=0.
  do m=1,nsym                                    !Generate iwave
     f=f0 + itone(m)*df
     dpha=twopi*f*dt
     do i=1,NSPS
        k=k+1
        pha=pha+dpha
        iwave(k)=nint(32767.0*sin(pha))
     enddo
  enddo
  nwave=k
  print*,'Diana ',f0,nwave
  write(*,3001) itone
3001 format(24i3)

  return
end subroutine gendiana
