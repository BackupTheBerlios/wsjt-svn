subroutine symspec65(dd,npts,ss,nhsym,savg)

  parameter (NFFT=8192)
  parameter (NSZ=3413)                         !NFFT*5000/12000
  parameter (MAXHSYM=322)
  real*8 hstep
  real*4 dd(npts)
  real*4 ss(MAXHSYM,NSZ)
  real*4 savg(NSZ)
  real*4 ref(NSZ)
  real*4 x(NFFT)
  complex c(0:NFFT/2)
  equivalence (x,c)

  hstep=2048.d0*12000.d0/11025.d0              !half-symbol = 2229.116 samples
  nsps=nint(2*hstep)
  df=12000.0/NFFT
  nhsym=npts/hstep - 1.0
  savg=0.
  fac1=1.e-3

  do j=1,nhsym
     i0=(j-1)*hstep
     x(1:nsps)=fac1*dd(i0+1:i0+nsps)
     x(nsps+1:)=0.
     call four2a(c,NFFT,1,-1,0)                !r2c forward FFT
     do i=1,NSZ
        s=real(c(i))**2 + aimag(c(i))**2
        ss(j,i)=s
        savg(i)=savg(i)+s
     enddo
  enddo
  savg=savg/nhsym

  call flat65(ss,nhsym,MAXHSYM,NSZ,ref)

!  do i=1,NSZ
!     write(71,3001) i*df,savg(i),db(savg(i)),ref(i),db(ref(i))
!3001 format(5f15.3)
!  enddo

  savg=savg/ref
  do j=1,nhsym
     ss(j,1:NSZ)=ss(j,1:NSZ)/ref
  enddo

!  do i=1,NSZ
!     write(72,3001) i*df,savg(i),db(savg(i))
!  enddo

  return
end subroutine symspec65
