subroutine decodems(dat,npts,nchar,imsg)

  parameter (NZ=30*11025)
  real dat(NZ)
  character msg*28
  complex zsum
  integer imsg(5250)
  real ss(0:63)
  real*8 t
  complex cw
  common/mscom/cw(63,0:63)

  nchar=5250
  fac=1.0/300000.0
  nerr=0
  do k=1,nchar
     i0=(k-1)*63
     smax1=0.
     do j=0,63                               !Test 64 possible characters
        zsum=0.
        do i=1,63                            !Sum over 63 samples
           zsum=zsum + conjg(cw(i,j))*dat(i0+i)
        enddo
        sum=abs(zsum)
        if(sum.gt.smax1) then
           smax1=sum
           jpk=j
        endif
        ss(j)=sum
     enddo

     ss(jpk)=0.
     smax2=0.
     do j=0,63                               !Test 64 possible characters
        smax2=max(ss(j),smax2)
     enddo
     rr=smax1/smax2
     if(jpk.ne.imsg(k)) nerr=nerr+1
!     write(*,3001) k,imsg(k),jpk,rr
!3001 format(3i6,f12.3)
  enddo

  ferr=float(nerr)/nchar
  write(*,1110) nerr,ferr
1110 format('Errors:',i5,f10.4)

  return
999 end subroutine decodems
