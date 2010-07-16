subroutine geniscat(msg,nmsg,shok,iwave,nwave,sendingsh,msgsent)

  parameter (NMAX=30*11025,NSZ=1291,NSPS=256)
  character msg*28,msgsent*22
  integer*2 iwave(NMAX)
  integer sendingsh
  logical first
  integer shok
  integer imsg(30)
  integer itone(NSZ)
  character c*42
  real*8 twopi,dt,f0,f,df,pha,dpha
  integer icos(4)
  integer irpt(31)
  data icos/0,1,3,2/
  data nsync/4/,nlen/2/,ndat/18/,jz/645/
  data c/'0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ /.?@-'/
  data irpt/1,2,3,4,6,7,8,9,11,12,13,14,15,16,17,18,19,21,22,23,24,     &
       25,26,27,28,29,31,32,33,34,35/

  twopi=8.d0*atan(1.d0)
  df=11025.d0/NSPS
  dt=1.d0/11025.d0
  f0=13*df
  nsym=NMAX/NSPS
  sendingsh=0

! Check for shorthand message
  if(shok.eq.1 .and. nmsg.le.4 .and.                                   &
       (msg(1:1).eq.'R' .or. msg(1:1).eq.'7')) then
     n=0
     m=0
     if(nmsg.eq.2 .and. msg(1:3).eq.'RO') n=5
     if(nmsg.eq.3 .and. msg(1:3).eq.'R26') n=10
     if(nmsg.eq.3 .and. msg(1:3).eq.'R27') n=20
     if(nmsg.eq.3 .and. msg(1:3).eq.'RRR') n=30
     if(nmsg.eq.2 .and. msg(1:2).eq.'73') n=40
     if(n.eq.0 .and. msg(1:1).eq.'R') then
        read(msg(2:4),*,err=10) m
        if(m.lt.-20) m=-20
        if(m.gt.10) m=10
        write(msgsent,1002) m
1002    format('R',i3)
        if(msgsent(2:2).eq.' ') msgsent=msgsent(1:1)//msgsent(3:)
        if(msgsent(2:2).eq.' ') msgsent=msgsent(1:1)//msgsent(3:)
        n=irpt(m+21)
     endif

     if(n.ne.0) then
        do i=1,nsym-1,2
           itone(i)=0
           itone(i+1)=n
        enddo
        sendingsh=1
     endif
  else
10   nblk=nsync+nlen+ndat
     msglen=nmsg+1
     k=0
     kk=1
     imsg(1)=40
     do i=1,nmsg
        imsg(i+1)=36
        do j=1,42
           if(msg(i:i).eq.c(j:j)) imsg(i+1)=j-1
        enddo
     enddo

     do i=1,nsym
        j=mod(i-1,nblk)+1
        if(j.le.nsync) then
           itone(i)=icos(j)
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
  endif

! Generate iwave
  k=0
  pha=0.
  do m=1,nsym
     f=f0 + itone(m)*df
     dpha=twopi*f*dt
     do i=1,NSPS
        k=k+1
        pha=pha+dpha
        iwave(k)=nint(32767.0*sin(pha))
     enddo
  enddo
  nwave=k

  return
end subroutine geniscat