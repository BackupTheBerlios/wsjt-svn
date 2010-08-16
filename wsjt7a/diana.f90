subroutine diana(dat,npts,cfile6,MinSigdB,DFTolerance,NFreeze,       &
     MouseDF,ccfblue,ccfred)

! Decode an ISCAT_2 signal

  parameter (NMAX=512*1024)
  parameter (NSZ=646)
  real dat(NMAX)                          !Raw signal, 30 s at 11025 sps
  character cfile6*6                      !File time
  character c42*42
  character msg*28
  real x(4096),x2(4096)
  complex c(0:4096)
  real s0(1024,NSZ)
  real fs0(1024,108)                       !108 = 96 + 3*4
  real fs1(0:41,30)
  real savg(1024)
  real b(1024)
  real ccfblue(-5:540)
  real psavg(1024)         !Average spectrum of the whole file
  integer dftolerance
  integer isync(4)
  equivalence (x,c)
  data isync/8,16,32,24/
  data nsps/2048/,nsync/4/,nlen/2/,ndat/18/
  data c42/'0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ /.?@-'/

! Define some constants
  nsym=npts/nsps                      !Total symbol intervals in file
  nblk=nsync+nlen+ndat                !Frame size
  nfft=4096                           !Do FFTs at twice the symbol length
  kstep=nsps/4                        !Step by 1/4 symbol
  nh=nfft/2
  nq=nfft/4
  df=11025.0/nfft                     !Bin spacing in frequency
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

  i0=2*236
  smax=0.
  ipk=9999
  jpk=9999
  ia=-10
  ib=10
  if(nfreeze.eq.1) then
     ia=(mousedf-dftolerance)/df
     ib=(mousedf+dftolerance)/df
  endif

  do j=0,4*nblk-1                     !Find sync pattern, lags 0-95
     do i=ia,ib                       !Search over DF range
        ss=0.
        do n=1,4                      !Sum the four sync tones
           k=j+4*n-3
           if(k.gt.4*nblk) k=k-4*nblk
           ss=ss + fs0(i0+i+2*isync(n),k)
        enddo
        if(ss.gt.smax) then
           smax=ss
           ipk=i0+i                   !Frequency offset, DF
           jpk=j+1                    !Time offset, DT
        endif
     enddo
  enddo

  dfx=(ipk-i0)*df
  dtx=jpk*kstep/11025.0

  ref=fs0(ipk+2,jpk) + fs0(ipk+4,jpk) + fs0(ipk+6,jpk)  +        &
      fs0(ipk,jpk+4) + fs0(ipk+4,jpk+4) + fs0(ipk+6,jpk+4) +     &
      fs0(ipk,jpk+8) + fs0(ipk+2,jpk+8) + fs0(ipk+4,jpk+8) +     &
      fs0(ipk,jpk+12) + fs0(ipk+2,jpk+12) + fs0(ipk+6,jpk+12)
  ref=ref/3.0                         !Reference level near (DF,DT)

  kk=0
  do j=0,4*nblk-1                     !Compute ccfblue
     ss=0.
     do n=1,4
        k=j+4*n-3
        if(k.gt.4*nblk) k=k-4*nblk
        ss=ss + fs0(ipk+2*isync(n),k)
     enddo
     ccfblue(j+1)=ss/ref
!     write(52,3002) j+1,ccfblue(j+1)
!3002 format(i8,f12.3)
  enddo

  xsync=smax/ref
  nsig=nint(db(smax/ref - 1.0) -15.0)
  if(nsig.lt.-20) nsig=-20
  ndf0=nint((ipk-i0) * 11025.0/nfft)
!  if(nsig.lt.MinSigdB) go to 800

  
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
  fs1=0.
  jb=(jsym-4*nblk+1)/4
  jb=4*jb
  k=0
  n=0
  do j=jpk,jpk+4*125,4                    !Fold s0 data symbols into fs1
     k=k+1                                ! ... modulo msglen
     if(mod(k-1,nblk)+1.gt.6) then
        n=n+1
        m=mod(n-1,msglen)+1
        iblk=(j-jpk)/(4*nblk)             !iblk runs from 0 to 4
        ioffset=7*iblk
        do i=0,41
           ii=i+ioffset
           if(ii.lt.0) ii=ii+42
           fs1(i,m)=fs1(i,m) + s0(ipk+2*ii,j)
        enddo
     endif
  enddo

! Read out the message:
  msg='                            '
  worst=9999.
  sum=0.
  do m=1,msglen
     smax=0.
     smax2=0.
     do i=0,41
        if(fs1(i,m).gt.smax) then
           smax=fs1(i,m)
           ipk3=i
        endif
     enddo
     do i=0,41
        if(fs1(i,m).gt.smax2 .and. i.ne.ipk3) smax2=fs1(i,m)
     enddo
     rr=smax/smax2
     sum=sum + rr
     if(rr.lt.worst) worst=rr
     msg(m:m)=c42(ipk3+1:ipk3+1)
  enddo

  avg=sum/msglen

800 continue
  if(nsig.lt.MinSigdB) then
     msglen=0
     worst=1.
     avg=1.
  endif
  nsnr=10.0*(avg-1.0)
  if(nsnr.gt.10) nsnr=10
  xsync=xsync-0.3
  jsync=xsync
  if(nsnr.le.-99) msg=' '
  jdf=nint(dfx)
  nwidth=0

  call cs_lock('iscat')
!  write(*,1020) cfile6,jsync,nsnr,dtx,jdf,nwidth,msg
  write(11,1020) cfile6,jsync,nsnr,dtx,jdf,nwidth,msg
  write(21,1020) cfile6,jsync,nsnr,dtx,jdf,nwidth,msg
1020 format(a6,i3,i5,f5.1,i5,i3,7x,a28)
  call cs_unlock

900 return
end subroutine diana
