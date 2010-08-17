subroutine specdiana(dat,npts,s0,jsym,fs0)

  parameter (NSZ=646)
  real dat(npts)                          !Raw signal, 30 s at 11025 sps
  real x(4096),x2(4096)
  complex c(0:4096)
  real s0(1024,NSZ)
  real fs0(1024,108)                       !108 = 96 + 3*4
  real savg(1024)
  real b(1024)
  equivalence (x,c)

  nsps=2048
  nblk=24
  nsym=npts/nsps                      !Total symbol intervals in file
  nfft=4096                           !Do FFTs at twice the symbol length
  kstep=nsps/4                        !Step by 1/4 symbol
  nq=nfft/4
  fac=1.0/1000.0                      !Somewhat arbitrary
  savg=0.

  ia=1-kstep
  do j=1,4*nsym                       !Compute symbol spectra
     ia=ia+kstep
     ib=ia+nsps-1
     if(ib.gt.npts) go to 10
     x(1:nsps)=fac*dat(ia:ib)
     x(nsps+1:nfft)=0.
     call four2a(x,nfft,1,-1,0)
     do i=1,nq
        s0(i,j)=real(c(i))**2 + aimag(c(i))**2
        savg(i)=savg(i) + s0(i,j)     !Accumulate average spectrum
     enddo
  enddo

10 jsym=j-1

  savg=savg/jsym
  do i=1,nq                           !Find baseline
     x(1:jsym)=s0(i,1:jsym)
     call pctile(x,x2,jsym,30,b(i))
  enddo
  b(1:10)=b(11)
  nadd=51
  call smo(b,nq,x2,nadd)              !Smooth the baseline

!  do i=nadd/1+1,nq-nadd/2
!     write(51,3001) i*df,savg(i),db(savg(i)),b(i),db(b(i))
!3001 format(5f12.3)
!  enddo

  do j=1,jsym                         !Normalize the spectra
     s0(1:nq,j)=s0(1:nq,j)/b(1:nq)
  enddo

  fs0=0.
  jb=(jsym-4*nblk+1)/4
  jb=4*jb
  do j=1,jb                           !Fold s0 into fs0, modulo 4*nblk
     k=mod(j-1,4*nblk)+1
     fs0(1:nq,k)=fs0(1:nq,k) + s0(1:nq,j)
  enddo

  do j=1,12                           !Replicate first 12 spectra at end
     fs0(1:nq,96+j)=fs0(1:nq,j)
  enddo

  return
end subroutine specdiana
