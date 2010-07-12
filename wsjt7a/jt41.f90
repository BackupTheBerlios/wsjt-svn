subroutine jt41(dat,npts,cfile6)

! Find sync 
  parameter (NMAX=512*1024)
  parameter (NSZ=4*1292)
  real dat(NMAX)                          !Raw signal, 30 s at 11025 sps
  character cfile6*6                      !File time
  character c41*41
  character msg*28,msg1*28
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
  data nsps/256/,nsync/4/,nlen/2/,ndat/18/
  data c41/'0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ /.?@'/
  
  nsym=npts/nsps
  nblk=nsync+nlen+ndat
  nfft=512                             !FFTs at twice the symbol length
  kstep=nsps/4                         !Step by 1/4 symbol
  nh=nfft/2
  nq=nfft/4
  df=11025.0/nfft
  fac=1.0/32768.0                      !Somewhat arbitrary
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

!  do i=1,nq
!     do j=1,jsym
!        x(j)=s0(i,j)
!     enddo
!     call pctile(x,x2,jsym,20,base)
!     savg(i)=savg(i)/(jsym*base)               !May be a problem for shorthands
!  enddo

!  do i=3,nq
!     write(13,3001) i*df,savg(i)
!3001 format(2f10.3)
!  enddo

  fs0=0.
  jb=(jsym-4*nblk+1)/4
  jb=4*jb
  do j=1,jb                                  !Fold s0 modulo 4*nblk into fs0
     k=mod(j-1,4*nblk)+1
     fs0(1:nq,k)=fs0(1:nq,k) + s0(1:nq,j)
  enddo

  i0=2*16
  smax=0.
  ipk=9999
  jpk=9999
  do j=0,4*nblk-1                            !Find the sync pattern
     do i=-10,10
        ss=0.
        do n=1,4
           k=j+4*n-3
           if(k.gt.4*nblk) k=k-4*nblk
           ss=ss + fs0(i0+i+2*icos(n),k)
        enddo
        if(ss.gt.smax) then
           smax=ss
           ipk=i0+i                          !Frequency offset, DF
           jpk=j+1                           !Time offset, DT
        endif
     enddo
  enddo

  if(ipk.gt.100 .or. jpk.gt.96) then
     print*,'ipk:',ipk,'   jpk:',jpk
     go to 900
  endif
  smax=0.
  do i=ipk,ipk+40,2                         !Find User's message length
     ss=fs0(i,jpk+16) + fs0(i+10,jpk+20)
     if(ss.gt.smax) then
        smax=ss
        ipk2=i
     endif
!     write(19,3003) i,ss
!3003 format(i5,e12.3)
  enddo
  msglen=(ipk2-i0)/2

  fs1=0.
  jb=(jsym-4*nblk+1)/4
  jb=4*jb
  k=0
  n=0
  do j=jpk,jsym,4                         !Fold information symbols into fs1
     k=k+1
     if(mod(k-1,nblk)+1.gt.6) then
        n=n+1
        m=mod(n-1,msglen)+1
!        write(15,4001) j,k,m
!4001    format(3i6)
        do i=0,40
           fs1(i,m)=fs1(i,m) + s0(ipk+2*i,j)
        enddo
     endif
  enddo

!  do i=0,40
!     write(16,4002) i,(nint(10*fs1(i,j)),j=1,10)
!4002 format(i2,10i7)
!  enddo

! Read out the message:

  msg1='                            '
  mpk=0
  do m=1,msglen
     smax=0.
     do i=0,40
        if(fs1(i,m).gt.smax) then
           smax=fs1(i,m)
           ipk3=i+1
        endif
     enddo
     if(ipk3.eq.41) mpk=m
     msg1(m:m)=c41(ipk3:ipk3)
  enddo

  if(mpk.eq.1) then
     msg=msg1(2:)
  else if(mpk.lt.msglen) then
     msg=msg1(mpk+1:msglen)//msg1(1:mpk-1)
  else
     msg=msg1(1:msglen-1)
  endif

  tping=jpk*kstep/11025.0
  width=0.0
  nsig=nint(db(smax))
  ndf0=nint((ipk-i0) * 11025.0/nfft)
  write(*,1010) cfile6,tping,width,nsig,ndf0,msg
  write(11,1010) cfile6,tping,width,nsig,ndf0,msg
  write(21,1010) cfile6,tping,width,nsig,ndf0,msg
1010 format(a6,2f5.1,i4,i5,6x,a28)

900 return
end subroutine jt41
