program genmet

  character*12 arg
  real*4 r(0:255)
  integer hist(2,-128:128)
  data hist/514*0/,idum/-1/

  lim(x)=min(127,max(-128,nint(scale*x)))

  nargs=iargc()
  if(nargs.ne.4) then
     print*,'Usage: genmet nadd m0 snr iters'
     go to 999
  endif
  call getarg(1,arg)
  read(arg,*) nadd
  call getarg(2,arg)
  read(arg,*) m0
  call getarg(3,arg)
  read(arg,*) snr
  call getarg(4,arg)
  read(arg,*) iters

  ntones=2**m0
  scale=5.0
  fac=sqrt(1.0/nadd)
  s=sqrt(10.0**(0.1*snr))

  do i=-128,128
     hist(1,i)=0
     hist(2,i)=0
  enddo
  nerr=1
  do iter=1,iters
     do i=0,ntones-1
        r(i)=0.
        do n=1,nadd
           x1=0.707*gasdev(idum)
           y1=0.707*gasdev(idum)
           if(i.eq.0) x1=x1+s
           r(i)=r(i) + x1*x1 + y1*y1
        enddo
        r(i)=fac*r(i)
     enddo
     do m=0,m0-1
        n=2**m
        r1=0.
        r2=0.
        do i=0,ntones-1
           if(iand(i,n).ne.0) r1=max(r1,r(i))
           if(iand(i,n).eq.0) r2=max(r2,r(i))
        enddo
        don=r2-r1
        doff=r1-r2
        if(don.lt.0.0) nerr=nerr+1
        j1=lim(doff)
        hist(1,j1)=hist(1,j1)+1
        j2=lim(don)
        hist(2,j2)=hist(2,j2)+1
     enddo
  enddo

  do i=-128,127
     write(13,1010) i/scale,hist(1,i)/(6.0*iters),hist(2,i)/(6.0*iters)
1010 format(f8.3,2f12.9)
  enddo

  ber=nerr/(6.0*iters)
  write(*,1020) nadd,m0,snr,ber
1020 format('nadd:',i3,'   m0:',i2,'   snr: 'f5.1,'   BER:',f8.3)
      
  xln2=log(2.0)
  do i=-128,127
     p1=hist(2,i)/(6.0*iters)
     p0=hist(1,i)/(6.0*iters)
     if(p0+p1.eq.0.0 .and. i.lt.0) p0=1.e-6
     if(p0+p1.eq.0.0 .and. i.gt.0) p1=1.e-6
     xlhd0=log(max(0.001,2.0*p0/(p0+p1)))/xln2
     xlhd1=log(max(0.001,2.0*p1/(p0+p1)))/xln2
     write(24,1012) i/scale,xlhd0,xlhd1,p0/(p0+p1),p1/(p0+p1)
1012 format(f7.1,2f8.3,2f9.6)
  enddo
      
999 end program genmet

