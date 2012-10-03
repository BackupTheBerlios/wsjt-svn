subroutine spec9(c0,npts8,nsps,f0a,lagpk,fpk)

  parameter (MAXFFT=31500)
  complex c0(0:npts8-1)
  real s(0:MAXFFT-1)
  real ssym(0:8,184)
  real ssymg(0:8,184)
  complex c(0:MAXFFT-1)
  integer*1 i1SoftSymbolsScrambled(207)
  integer*1 i1SoftSymbols(207)
  integer ibit(207)

  integer*1 t1(13)              !72 bits and zero tail as 8-bit bytes
  integer*4 t4(69)              !Symbols from t5, values 0-7
  integer*4 mettab(0:255,0:1)
  integer*1 tmp(72)
  character*22 msg

  integer isync(85)
  integer ii(16)                       !Locations of sync symbols
  data ii/1,6,11,16,21,26,31,39,45,51,57,63,69,75,81,85/
  integer ig(0:7)
  data ig/0,1,3,2,7,6,4,5/
!  data ig/0,1,3,2,6,7,5,4/

  isync=0
  do i=1,16
     isync(ii(i))=1
  enddo

  idt=-400
  idf=0.
  fshift=fpk-f0a + 0.1*idf
  twopi=8.0*atan(1.0)
  dphi=twopi*fshift/1500.0
  nsps8=nsps/8
  nfft=nsps8
  df=1500.0/nfft
  s=0.
  istart=lagpk*nsps8 + idt
  nsym=min((npts8-istart)/nsps8,85)

  do j=0,nsym-1
     ia=j*nsps8 + istart
     ib=ia+nsps8-1
     c(0:nfft-1)=c0(ia:ib)

     phi=0.
     do i=0,nfft-1
        phi=phi + dphi
        c(i)=c(i) * cmplx(cos(phi),-sin(phi))
     enddo

     call four2a(c,nfft,1,-1,1)
     do i=0,nfft-1
        sx=real(c(i))**2 + aimag(c(i))**2
        if(i.le.8) ssym(i,1+j)=sx
        s(i)=s(i) + sx
     enddo
  enddo

  ssymg=ssym
  do j=1,nsym
     ssym(0,j)=ssymg(0,j)
     do i=0,7
        ssym(ig(i)+1,j)=ssymg(i+1,j)
     enddo
  enddo

  m0=3
  ntones=8
  k=0
  do j=1,nsym
     if(isync(j).eq.1) cycle
     do m=m0-1,0,-1                   !Get bit-wise soft symbols
        n=2**m
        r1=0.
        r2=0.
        do i=0,ntones-1
           if(iand(i,n).ne.0) then
              r1=max(r1,ssym(i+1,j))
           else
              r2=max(r2,ssym(i+1,j))
           endif
        enddo
        k=k+1
        i1SoftSymbolsScrambled(k)=min(127,max(-127,nint(10.0*(r1-r2)))) + 128
     enddo
  enddo

  ibit=0
  do i=1,207
     if(i1SoftSymbolsScrambled(i).lt.0) ibit(i)=1
  enddo

! Get the metric table
  bias=0.37                          !To be optimized, in decoder program
  scale=10                           !  ... ditto ...
  open(19,file='met8.21',status='old')

  do i=0,255
     read(19,*) x00,x0,x1
     mettab(i,0)=nint(scale*(x0-bias))
     mettab(i,1)=nint(scale*(x1-bias))    !### Check range, etc.  ###
  enddo
  close(19)
  nbits=72
  ndelta=17
  limit=1000

  call interleave9(i1SoftSymbolsScrambled,-1,i1SoftSymbols)
  call fano232(i1SoftSymbols,nbits+31,mettab,ndelta,limit,t1,ncycles,    &
       metric,ierr,maxmetric,maxnp)

  nbytes=(nbits+7)/8
  do i=1,nbytes
     n=t1(i)
     t4(i)=iand(n,255)
  enddo
  call unpackbits(t4,nbytes,8,tmp)
  call packbits(tmp,12,6,t4)
  do i=1,12
     if(t4(i).lt.128) t1(i)=t4(i)
     if(t4(i).ge.128) t1(i)=t4(i)-256
  enddo
  do i=1,12
     t4(i)=t1(i)
  enddo
  call unpackmsg(t4,msg)         !Unpack decoded msg
  print*,msg

  return
end subroutine spec9
