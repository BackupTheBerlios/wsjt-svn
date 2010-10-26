program fmtave

  implicit real*8 (a-h,o-z)
  character infile*80
  character*8 cutc,cutc1

  nargs=iargc()
  if(nargs.ne.1) then
     print*,'Usage: fmtave <infile>'
     go to 999
  endif
  call getarg(1,infile)

  open(10,file=infile,status='old')
  read(10,*) cutc
  read(10,*) cutc
  read(10,*) cutc

  nkhz0=0
  sum=0.d0
  sumsq=0.d0
  n=0
  do i=1,99999
     read(10,*,end=10) cutc,nkHz,noffset,faudio,snr,yave,yrms
     if(nkHz.ne.nkHz0 .and. i.ne.1) then
        ave=sum/n
        rms=0.d0
        if(n.gt.1) then
           rms=sqrt(abs(sumsq - sum*sum/n)/(n-1.d0))
        endif
        err=rms/sqrt(n-1.d0)
        if(err.lt.0.1) err=0.1
        fMHz=0.001d0*nkHz0
        write(*,1010) fMHz,ave,rms,err,ave/fMHz,cutc1
1010    format(f8.3,4f8.2,2x,a8)
        sum=0.d0
        sumsq=0.d0
        n=0
     endif
     dial_error=faudio-noffset
     sum=sum + dial_error
     sumsq=sumsq + dial_error**2
     n=n+1
     if(n.eq.1) cutc1=cutc
     nkHz0=nkHz
  enddo

10 ave=sum/n
  rms=0.d0
  if(n.gt.0) then
     rms=sqrt((sumsq - sum*sum/n)/(n-1.d0))
  endif
  err=rms/sqrt(n-1.d0)
  if(err.lt.0.1) err=0.1
  fMHz=0.001d0*nkHz
  write(*,1010) fMHz,ave,rms,err,ave/fMHz,cutc1

999 end program fmtave
