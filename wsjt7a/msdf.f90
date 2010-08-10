subroutine msdf(cdat,npts,t2,nfft1,f0,nfreeze,mousedf,dftolerance,dfx,ferr)

! Determine DF for a JTMS signal.  Also find ferr, a measure of
! frequency differerence between 1st and 2nd harmonic.  
! (Should be 0.000)

  parameter (NZ=32768)
  complex cdat(npts)
  integer dftolerance
  real sq(NZ)
  real tmp(NZ)
  complex c(NZ)
  data nsps/8/
  save c

  rewind 50
  rewind 51

  df1=11025.0/nfft1
  nh=nfft1/2
  fac=1.0/(nfft1**2)

  do i=1,npts
     c(i)=fac*cdat(i)**2
  enddo
  c(npts+1:nfft1)=0.
  call four2a(c,nfft1,1,-1,1)

! In the "doubled-frequencies" spectrum of squared cdat:
  fa=2.0*(f0-400)
  fb=2.0*(f0+400)
  if(NFreeze.gt.0) then
     fa=2.0*(f0+MouseDF-DFtolerance)
     fb=2.0*(f0+MouseDF+DFtolerance)
  endif  
  ja=nint(fa/df1)
  jb=nint(fb/df1)
  jd=nfft1/nsps

  do j=1,nh+1
     sq(j)=real(c(j))**2 + aimag(c(j))**2
  enddo

!  call smo(sq,nh+1,tmp,9)
!  call smo(sq,nh+1,tmp,9)

  do j=1,nh+1
     sq2=0.
     if(j+jd.le.nfft1/2+1) sq2=sq(j+jd)
     f=0.5*(j-1)*df1 - f0
     write(50,3001) f,sq(j),sq2,j
3001 format(f10.2,2f12.3,i8)
  enddo

  smax=0.
  smax1=0.
  smax2=0.
  do j=ja,jb
     ss=sq(j)+sq(j+jd)
     f=0.5*(j-1)*df1 - f0
     write(51,3001) f,ss
     if(ss.gt.smax) then
        smax=sq(j)+sq(j+jd)
        jpk=j
     endif
     if(sq(j).gt.smax1) then
        smax1=sq(j)
        jpk1=j
     endif
     if(sq(j+jd).gt.smax2) then
        smax2=sq(j+jd)
        jpk2=j+jd
     endif
  enddo

  fpk=(jpk-1)*df1  
  dfx=0.5*fpk-f0

  fpk1=(jpk1-1)*df1
  fpk2=(jpk2-1)*df1
  ferr=(fpk2-fpk1)/1378.125 - 1.0
  
  write(*,3501) t2,dfx,ferr,fpk1,fpk2
3501 format(2f8.1,f12.6,2f10.2)
  call flush(50)
  call flush(51)

  return
end subroutine msdf
