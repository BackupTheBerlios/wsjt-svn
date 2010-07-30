subroutine msdf(cdat,npts,nfft1,f0,nfreeze,mousedf,dftolerance,dfx,ferr)

! Determine DF for a JTMS signal.  Also find ferr, a measure of
! frequency differerence between 1st and 2nd harmonic.  
! (Should be 0.000)

  parameter (NZ=32768)
  complex cdat(npts)
  integer dftolerance
  real sq(NZ)
  complex c(NZ)
  data nsps/8/
  save c

  df1=11025.0/nfft1
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

  do j=1,nfft1/2+1
     sq(j)=real(c(j))**2 + aimag(c(j))**2
  enddo

  smax=0.
  smax1=0.
  do j=ja,jb
     if(sq(j)+sq(j+jd).gt.smax) then
        smax=sq(j)+sq(j+jd)
        jpk=j
     endif
     if(sq(j).gt.smax1) then
        smax1=sq(j)
        jpk1=j
     endif
  enddo

  smax2=0.
  do j=jpk1+jd,jb+jd
     if(sq(j).gt.smax2) then
        smax2=sq(j)
        jpk2=j
     endif
  enddo

  fpk=(jpk-1)*df1  
  dfx=0.5*fpk-f0

  fpk1=(jpk1-1)*df1
  fpk2=(jpk2-1)*df1
  ferr=(fpk2-fpk1)/1378.125 - 1.0

  return
end subroutine msdf
