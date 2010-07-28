subroutine foldms(s2,msglen,nchar,mycall,msg,msg28)

  real s2(0:63,400)
  real fs2(0:63,29)
  integer nfs2(29)
  real sm(0:63)
  character mycall*12
  character msg*400,msg28*28
  character cc*64
!                    1         2         3         4         5         6
!          0123456789012345678901234567890123456789012345678901234567890123
  data cc/'0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ./?-                 _     @'/

  fs2=0.
  nfs2=0
  do j=1,nchar                           !Fold s2 into fs2, modulo msglen
     jj=mod(j-1,msglen)+1
     nfs2(jj)=nfs2(jj)+1
     do i=0,40
        fs2(i,jj)=fs2(i,jj) + s2(i,j)
     enddo
  enddo

  msg=' '
  do j=1,msglen
     smax=0.
     do k=0,40
        if(fs2(k,j).gt.smax) then
           smax=fs2(k,j)
           kpk=k
        endif
     enddo
     sm(kpk)=0.
     smax2=0.
     do k=0,40
        smax2=max(smax2,sm(k))
     enddo
     if(kpk.eq.40) kpk=57
     msg(j:j)=cc(kpk+1:kpk+1)
     if(kpk.eq.57) msg(j:j)=' '
  enddo
  msg28=msg(1:msglen)
  call match(mycall,msg28(1:msglen),nstart,nmatch)
  call match(' CQ ',msg28(1:msglen),nstart2,nmatch2)
  if(nmatch.ge.3 .and.nstart.gt.1) then
     msg28=msg(nstart:msglen)//msg(1:nstart-1)
  else if(nmatch.ge.3 .and.nstart.eq.1) then
     msg28=msg(1:msglen)
  else if(nmatch2.ge.3 .and.nstart2.gt.1) then
     msg28=msg(nstart2:msglen)//msg(1:nstart2-1)
  else if(nmatch2.ge.3 .and.nstart2.eq.1) then
     msg28=msg(2:msglen)
  else
     i3=index(msg,'  ')
     if(i3.gt.0 .and. i3.le.msglen-2) then
        msg28=msg(i3+2:msglen)//msg(1:msglen)
     else
        i3=index(msg,' ')
        if(i3.gt.0 .and. i3.lt.msglen) msg28=msg(i3:msglen)//msg(1:msglen)
     endif
  endif
  if(msg28(1:1).eq.' ') msg28=msg28(2:)

  return
end subroutine foldms
