subroutine sync41(dat,npts,cfile6)
  parameter (NMAX=512*1024)
  parameter (NSZ=4*1292)
  real dat(NMAX)                          !Raw signal, 30 s at 11025 sps
  character cfile6*6                      !File time
  real x(NSZ),x2(NSZ)
  complex c(128)
  real s0(128,NSZ)
  real fs0(128,96)
  real fs1(0:40,30)
  real savg(128)
  real ccfred(-10:10)
  real ccfblue(0:95)
  integer dftolerance
  integer icos(4)
  equivalence (x,c)
  data icos/0,1,3,2/
  
  nsps=256
  nsym=npts/nsps
  kstep=nsps/4
  nfft=512
  nh=nfft/2
  nq=nfft/4
  df=11025.0/nfft
  fac=1.0/32768.0
  savg=0.

  ia=1-kstep
  do j=1,4*nsym
     ia=ia+kstep
     ib=ia+nsps-1
     if(ib.gt.npts) go to 10
     x(1:nsps)=fac*dat(ia:ib)
     x(nsps+1:nfft)=0.
     call four2a(x,nfft,1,-1,0)
     do i=1,nq
        s0(i,j)=real(c(i))**2 + aimag(c(i))**2
        savg(i)=savg(i) + s0(i,j)
     enddo
  enddo

10 jsym=j
  print*,'A',npts,jsym

!  do i=1,nq
!     do j=1,jsym
!        x(j)=s0(i,j)
!     enddo
!     call pctile(x,x2,jsym,20,base)
!     savg(i)=savg(i)/(jsym*base)               !May be a problem for shorthands
!  enddo

  do i=3,nq
     write(13,3001) i*df,savg(i)
3001 format(2f10.3)
  enddo

  fs0=0.
  nblk=4+2+18
  jb=(jsym-4*nblk+1)/4
  jb=4*jb

  do j=1,jb
     k=mod(j-1,4*nblk)+1
     fs0(1:nq,k)=fs0(1:nq,k) + s0(1:nq,j)
  enddo

  i0=2*16
  smax=0.
  ipk=9999
  jpk=9999
  do j=0,4*nblk-1
     do i=-10,10
        ss=0.
        do n=1,4
           k=j+4*n-3
           if(k.gt.4*nblk) k=k-4*nblk
           ss=ss + fs0(i0+i+2*icos(n),k)
        enddo
        if(ss.gt.smax) then
           smax=ss
           ipk=i0+i
           jpk=j+1
        endif
     enddo
  enddo

  print*,'b',ipk,jpk,smax

  smax=0.
  do i=ipk,ipk+40,2
     ss=fs0(i,jpk+16) + fs0(i+10,jpk+20)
     if(ss.gt.smax) then
        smax=ss
        ipk2=i
     endif
     write(19,3003) i,ss
3003 format(i5,e12.3)
  enddo
  msglen=(ipk2-i0)/2
  print*,'d',msglen,smax

  fs1=0.
  jb=(jsym-4*nblk+1)/4
  jb=4*jb

! Fold the information symbols
  k=0
  n=0
  do j=jpk,jsym,4
     k=k+1
     if(mod(k-1,nblk)+1.gt.6) then
        n=n+1
        m=mod(n-1,msglen)+1
        do i=0,40
           fs1(i,m)=fs1(i,m) + s0(ipk+2*i,k)
        enddo
     endif
  enddo

! Read out the message:

  do m=1,msglen
     smax=0.
     do i=0,40
        if(fs1(i,m).gt.smax) then
           smax=fs1(i,m)
           ipk3=i
        endif
     enddo
     print*,'e',ipk3,smax
  enddo

!  j=jpk
!  do i=-10,10
!     ss=0.
!     do n=1,4
!        k=j+4*n-3
!        if(k.gt.4*nblk) k=k-4*nblk
!        ss=ss + fs0(i0+i+2*icos(n),k)
!     enddo
!     ccfred(i)=ss
!     write(17,3003) i,ccfred(i)
!  enddo

!  i=ipk
!  do j=0,4*nblk-1
!     ss=0.
!     do n=1,4
!        k=j+4*n-3
!        if(k.gt.4*nblk) k=k-4*nblk
!        ss=ss + fs0(i0+i+2*icos(n),k)
!     enddo
!     ccfblue(j)=ss
!     write(18,3003) j,ccfblue(j)
!  enddo

  return
end subroutine sync41
