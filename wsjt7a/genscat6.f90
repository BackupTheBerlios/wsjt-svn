subroutine genscat6(msg,itone)

  character msg*28,msg1*29
  character c*43
  integer imsg(29)
  integer itone(645)
  integer icos6(6)
  data icos6/0,1,4,3,5,2/
  data c/'0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ .,?/#$'/  !1-43
  data nsync/7/,ndat/14/

  do i=28,1,-1
     if(msg(i:i).ne.' ') go to 10
  enddo
10 msglen=i+1                      !Message length, including leading '$'
  msg1='$'//msg

  do i=1,msglen
     i1=index(c,msg1(i:i))
     imsg(i)=i1-1
  enddo

  ntot=ndat+nsync
  k=0
  kk=1
  do i=1,645
     j=mod(i-1,ntot)+1
     if(j.lt.nsync) then
        itone(i)=icos6(j)
     else if(j.eq.nsync) then
        itone(i)=msglen
     else
        k=k+1
        kk=mod(k-1,msglen)+1
        itone(i)=imsg(kk)
     endif
  enddo

  return

end subroutine genscat6
