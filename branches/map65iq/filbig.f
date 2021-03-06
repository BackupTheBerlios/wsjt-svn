      subroutine filbig(dd,nmax,f0,newdat,c4a,n4)

C  Filter and downsample the complex data for one polarization,
C  stored in array dd(2,nmax).  Output is downsampled from 96000 Hz
C  to 1378.125 Hz.  

!      parameter (NFFT1=5376000,NFFT2=77175)
      parameter (NFFT1=5120000,NFFT2=74088)
      real*4 dd(2,nmax)                          !Input data
      complex c4a(NFFT2)                         !Output data
      complex ca(NFFT1)                          !FFT of input
      real*8 df
C Impulse response of filter (one side)
      real halfpulse(8)
      complex cfilt(NFFT2)
                       !Filter (complex; imag = 0)
      real rfilt(NFFT2)                          !Filter (real)
      integer plan1,plan3,plan5
      logical first
      include 'fftw3.f'
      equivalence (rfilt,cfilt)
      data first/.true./
      data halfpulse/114.97547150,36.57879257,-20.93789101,
     +  5.89886379,1.59355187,-2.49138308,0.60910773,-0.04248129/
      save

      if(nmax.lt.0) go to 900
      if(first) then
         nflag=FFTW_ESTIMATE_PATIENT
!         nflag=FFTW_MEASURE
C  Plan the FFTs just once
         call sfftw_plan_dft_1d_(plan1,NFFT1,ca,ca,
     +        FFTW_BACKWARD,nflag)
         call sfftw_plan_dft_1d_(plan3,NFFT2,c4a,c4a,
     +        FFTW_FORWARD,nflag)
         call sfftw_plan_dft_1d_(plan5,NFFT2,cfilt,cfilt,
     +        FFTW_BACKWARD,nflag)

C  Convert impulse response to filter function
         do i=1,NFFT2
            cfilt(i)=0.
         enddo
         fac=0.00625/NFFT1
         cfilt(1)=fac*halfpulse(1)
         do i=2,8
            cfilt(i)=fac*halfpulse(i)
            cfilt(NFFT2+2-i)=fac*halfpulse(i)
         enddo
         call sfftw_execute_(plan5)

         base=cfilt(NFFT2/2+1)
         do i=1,NFFT2
            rfilt(i)=real(cfilt(i))-base
         enddo

!         df=96000.d0/NFFT1
         df=95238.1d0/NFFT1
         first=.false.
      endif

C  When new data comes along, we need to compute a new "big FFT"
C  If we just have a new f0, continue with the existing ca and cb.

      if(newdat.ne.0) then
         nz=min(nmax,NFFT1)
         do i=1,nz
            ca(i)=cmplx(dd(1,i),dd(2,i))
         enddo

         if(nmax.lt.NFFT1) then
            do i=nmax+1,NFFT1
               ca(i)=0.
            enddo
         endif
         call sfftw_execute_(plan1)
         newdat=0
      endif

C  NB: f0 is the frequency at which we want our filter centered.
C      i0 is the bin number in ca and cb closest to f0.

      i0=nint(f0/df) + 1
      nh=NFFT2/2
      do i=1,nh                                !Copy data into c4a
         j=i0+i-1                              !and apply the filter function
         c4a(i)=rfilt(i)*ca(j)
      enddo
      do i=nh+1,NFFT2
         j=i0+i-1-NFFT2
         if(j.lt.1) j=j+NFFT1                  !NFFT1 was NFFT2
         c4a(i)=rfilt(i)*ca(j)
      enddo

C  Do the short reverse transform, to go back to time domain.
      call sfftw_execute_(plan3)

      n4=min(nmax/64,NFFT2)
      go to 999

 900  call sfftw_destroy_plan_(plan1)
      call sfftw_destroy_plan_(plan3)
      call sfftw_destroy_plan_(plan5)

 999  return
      end
