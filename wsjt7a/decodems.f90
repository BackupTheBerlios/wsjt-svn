subroutine decodems(dat,npts,cfile6,t2,mswidth,ndb,nrpt,Nfreeze,DFTolerance)

  parameter (NZ=30*11025)
  real dat(NZ)
  integer DFTolerance
  character*6 cfile6

  ndf=0
  write(*,3001) cfile6,t2,mswidth,ndb,nrpt,ndf,npts,NFreeze,DFTolerance
3001 format(a6,f5.1,i5,i3,i3.2,i5,2x,4i6)
! 1050    format(a6,f5.1,i5,i3,1x,i2.2,i5,1x,a3,1x,a40)

  if(t2.lt.3.0) then
     do i=1,npts
        write(51,3002) i,dat(i)
3002    format(i8,f10.1)
     enddo
  endif

  return
999 end subroutine decodems
