subroutine genscat6(msg,itone)

  character*28 msg
  integer itone(645)
  integer icos6(6)
  integer nt(6)
  character c*43
  data icos6/1,2,3,4,5,6/
  data nt/13,17,19,23,29,31/
  data c/'0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ .,?/#$'/

  do i=24,1,-1
     if(msg(i:i).ne.' ') go to 10
  enddo
10 msglen=i

  do i=1,6
     ntot=nt(i)
     if(ntot.ge.msglen+7) go to 20
  enddo

20 ndext=ntot-7 
  itone(1:6)=icos6
  itone(7)=ntot
  do i=1,msglen
     i1=index(c,msg(i:i))
     itone(7+i)=i1-1
  enddo
  if(ntot.gt.msglen+7) itone(msglen+8:ntot)=ntot

  write(*,3001) msglen,ndext,ntot,(itone(i),i=1,ntot)
3001 format(3i3,2x,7i3/(20i3))

  j=ntot
  do nrpt=1,999
     do i=1,ntot
        j=j+1
        if(j.gt.645) go to 900
        itone(j)=itone(i)
     enddo
  enddo

900 return
end subroutine genscat6
