subroutine decdiana(s0,jsym,ipk,jpk,msglen,msg,avg)

  parameter (NSZ=646)
  real s0(1024,NSZ)
  real fs1(0:41,30)
  character msg*28
  character c42*42
  data c42/'0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ /.?+-'/

  nblk=24
  fs1=0.
  jb=(jsym-4*nblk+1)/4
  jb=4*jb
  k=0
  n=0
  do j=jpk,jpk+4*125,4                    !Fold s0 data symbols into fs1
     k=k+1                                ! ... modulo msglen
     if(mod(k-1,nblk)+1.gt.6) then
        n=n+1
        m=mod(n-1,msglen)+1
        iblk=(j-jpk)/(4*nblk)             !iblk runs from 0 to 4
        ioffset=7*iblk
        do i=0,41
           ii=i+ioffset
           if(ii.ge.42) ii=ii-42
           fs1(i,m)=fs1(i,m) + s0(ipk+2*ii,j)
        enddo
     endif
  enddo

! Read out the message:
  msg='                            '
  worst=9999.
  sum=0.
  do m=1,msglen
     smax=0.
     smax2=0.
     do i=0,41
        if(fs1(i,m).gt.smax) then
           smax=fs1(i,m)
           ipk3=i
        endif
     enddo
     do i=0,41
        if(fs1(i,m).gt.smax2 .and. i.ne.ipk3) smax2=fs1(i,m)
     enddo
     rr=smax/smax2
     sum=sum + rr
     if(rr.lt.worst) worst=rr
     msg(m:m)=c42(ipk3+1:ipk3+1)
  enddo

  avg=sum/msglen

  return
end subroutine decdiana
