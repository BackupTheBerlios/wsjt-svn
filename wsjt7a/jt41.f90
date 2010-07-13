subroutine jt41(dat,npts,cfile6,MinSigdB,DFTolerance,NFreeze,MouseDF,ccf,psavg)

! Decoder for JT41 

  parameter (NMAX=512*1024)
  parameter (NSZ=4*1292)
  real dat(NMAX)                          !Raw signal, 30 s at 11025 sps
  character cfile6*6                      !File time
  character c41*41
  character msg*28,msg1*28
  real x(NSZ),x2(NSZ)
  complex c(0:512)
  real s0(128,NSZ)
  real fs0(128,108)                       !108 = 96 + 3*4
  real fs1(0:40,30)
  real savg(128)
  real savg2(128)
  real b(128)
  real ccfred(-10:10)
  real ccfblue(0:95)
  real ccf(-5:540)
  real psavg(450)         !Average spectrum of the whole file
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
  fac=1.0/1000.0                      !Somewhat arbitrary
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

10 jsym=j-1

  savg=savg/jsym
  do i=1,nq
     x(1:jsym)=s0(i,1:jsym)
     call pctile(x,x2,jsym,30,b(i))
  enddo
  b(1:10)=b(11)
  do i=1,nq/2
     psavg(i)=2*db(savg(2*i)) + 10.0
  enddo

!  rewind 53
!  do i=1,nq
!     savg2(i)=savg(i)/b(i)
!     write(53,3001) i*df,savg(i),b(i),savg2(i)
!3001 format(4f10.3)
!  enddo
!  call flush(53)

  do j=1,jsym
     s0(1:nq,j)=s0(1:nq,j)/b(1:nq)
  enddo

  fs0=0.
  jb=(jsym-4*nblk+1)/4
  jb=4*jb
  do j=1,jb                                  !Fold s0 modulo 4*nblk into fs0
     k=mod(j-1,4*nblk)+1
     fs0(1:nq,k)=fs0(1:nq,k) + s0(1:nq,j)
  enddo

  do j=1,12
     fs0(1:nq,96+j)=fs0(1:nq,j)
  enddo

  i0=2*13
  smax=0.
  ipk=9999
  jpk=9999
  ia=-10
  ib=10
  if(nfreeze.eq.1) then
     ia=(mousedf-dftolerance)/df
     ib=(mousedf+dftolerance)/df
  endif

  do j=0,4*nblk-1                            !Find sync pattern, lags 0-95
     do i=ia,ib
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

  ref=fs0(ipk+2,jpk) + fs0(ipk+4,jpk) + fs0(ipk+6,jpk)  +        &
      fs0(ipk,jpk+4) + fs0(ipk+4,jpk+4) + fs0(ipk+6,jpk+4) +     &
      fs0(ipk,jpk+8) + fs0(ipk+2,jpk+8) + fs0(ipk+4,jpk+8) +     &
      fs0(ipk,jpk+12) + fs0(ipk+2,jpk+12) + fs0(ipk+6,jpk+12)
  ref=ref/3.0

  kk=0
!  rewind 54
  do j=0,4*nblk-1
     ss=0.
     do n=1,4
        k=j+4*n-3
        if(k.gt.4*nblk) k=k-4*nblk
        ss=ss + fs0(ipk+2*icos(n),k)
     enddo
     kk=kk+1
     ccf(kk)=ss/ref
!     write(54,3101) kk,ss/ref
!3101 format(i5,f12.3)
  enddo
!  call flush(54)

  tping=jpk*kstep/11025.0
  xsync=smax/ref
  nsig=nint(db(smax/ref - 1.0) -15.0)
  if(nsig.lt.-20) nsig=-20
  ndf0=nint((ipk-i0) * 11025.0/nfft)
  if(nsig.lt.MinSigdB) go to 800

  if(ipk.gt.100 .or. jpk.gt.96) then
     print*,'ipk:',ipk,'   jpk:',jpk
     go to 900
  endif
  smax=0.
  ja=jpk+16
  if(ja.gt.4*nblk) ja=ja-4*nblk
  jb=jpk+20
  if(jb.gt.4*nblk) jb=jb-4*nblk
  do i=ipk,ipk+60,2                         !Find User's message length
     ss=fs0(i,ja) + fs0(i+10,jb)
     if(ss.gt.smax) then
        smax=ss
        ipk2=i
     endif
  enddo

  msglen=(ipk2-ipk)/2
  if(msglen.lt.1 .or. msglen.gt.28) msglen=2         !### tests only ###
!  if(msglen.lt.1 .or. msglen.gt.28) then
!     print*,'msglen:',msglen
!     go to 900
!  endif

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
        do i=0,40
           fs1(i,m)=fs1(i,m) + s0(ipk+2*i,j)
        enddo
     endif
  enddo

! Read out the message:
  msg1='                            '
  mpk=0
  worst=9999.
  sum=0.
  do m=1,msglen
     smax=0.
     smax2=0.
     do i=0,40
        if(fs1(i,m).gt.smax) then
           smax=fs1(i,m)
           ipk3=i
        endif
     enddo
     do i=0,40
        if(fs1(i,m).gt.smax2 .and. i.ne.ipk3) smax2=fs1(i,m)
     enddo
     rr=smax/smax2
     sum=sum + rr
     if(rr.lt.worst) worst=rr
     if(ipk3.eq.40) mpk=m
     msg1(m:m)=c41(ipk3+1:ipk3+1)
  enddo

  avg=sum/msglen
  if(mpk.eq.1) then
     msg=msg1(2:)
  else if(mpk.lt.msglen) then
     msg=msg1(mpk+1:msglen)//msg1(1:mpk-1)
  else
     msg=msg1(1:msglen-1)
  endif

800 continue
  if(nsig.lt.MinSigdB) then
     msglen=0
     worst=1.
     avg=1.
  endif
  nworst=10.0*(worst-1.0)
  navg=10.0*(avg-1.0)
  if(nworst.gt.10) nworst=10
  if(navg.gt.10) navg=10
  xsync=xsync-0.3
  isync=xsync
  if(navg.le.0) msg=' '

!  write(*,1020)  cfile6,isync,width,nsig,ndf0,msg,msglen,nworst,navg
  write(11,1020) cfile6,nsig,ndf0,msg,msglen,nworst,navg
  write(21,1020) cfile6,nsig,ndf0,msg,msglen,nworst,navg
1020 format(a6,i5,i5,6x,a28,i4,2i3)
  call flush(11)
  call flush(21)

900 return
end subroutine jt41
