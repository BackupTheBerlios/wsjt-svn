subroutine decodems(dat,npts,cfile6,t2,mswidth,ndb,nrpt,Nfreeze,       &
     DFTolerance,MouseDF,mycall)

! Decode a JTMS ping

  parameter (NZ=30*11025)
  real dat(npts)                        !Raw data
  complex cdat(NZ)                      !Analytic form of signal
  character*6 cfile6                    !FileID
  integer DFTolerance
  character*12 mycall
  real s(NZ)                            !Power spectrum
  real sq(NZ)                           !Double-frequency spectrum
  real sm(0:63)
  real s2(0:63,400)
  real fs2(0:63,28)
  integer nfs2(28)
  real r(60000)
  real fr(56)
  integer nfr(56)
  real acf(1624)
  complex c(NZ)
  complex cw(56,0:63)                   !Complex waveforms for codewords
  complex cwb(56)                       !Complex waveform for 'space'
  complex z,zmax,zmax0
  logical first
  character msg*400,msg28*28
  character cc*64
  integer np(8)
  character*90 line
  common/ccom/nline,tping(100),line(100)
  data np/5,7,11,13,17,19,23,29/        !Permissible message lengths
!                    1         2         3         4         5         6
!          0123456789012345678901234567890123456789012345678901234567890123
  data cc/'0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ./?-                 _     @'/
  data first/.true./
  data nsum/0/,nrec/0/
  save nsum,nrec
  save first,smax,cw,cwb              !Why is this needed for save?  But it is!
  save c,cdat

  nrec=nrec+1
  if(first) call setupms(cw,cwb)        !Calculate waveforms for codewords
  first=.false.

  nsps=8                                !Samples per symbol
  f0=1155.46875                         !Nominal frequency for bit=0
  xn=log(float(npts))/log(2.0)
  n=xn
  if(xn-n .gt.0.001) n=n+1
  nfft1=2**n
  df1=11025.0/nfft1

  call analytic(dat,npts,nfft1,s,cdat)        !Convert to analytic signal

  fac=1.0/(nfft1**2)
  do i=1,npts
     c(i)=fac*cdat(i)**2
  enddo
  c(npts+1:nfft1)=0.
  call four2a(c,nfft1,1,-1,1)

  rewind 52
  do i=1,nfft1/2+1
     ss=real(c(i))**2 + aimag(c(i))**2
     write(52,5001) df1*(i-1),ss
5001 format(f10.3,f12.3)
  enddo

! In the "doubled-frequencies" spectrum of squared cdat:
  fa=2.0*(f0-400)
  fb=2.0*(f0+400)
  if(NFreeze.gt.0) then
     fa=2.0*(f0+MouseDF-DFtolerance)
     fb=2.0*(f0+MouseDF+DFtolerance)
  endif  
  ja=nint(fa/df1)
  jb=nint(fb/df1)
  jd=nfft1/nsps

  smax=0.
  do j=ja,jb
     ss=real(c(j))**2 + aimag(c(j))**2 + real(c(j+jd))**2 + aimag(c(j+jd))**2
     sq(j)=ss
     if(ss.gt.smax) then
        smax=ss
        jpk=j
        fpk=(j-1)*df1
     endif
  enddo
  call pctile (sq(ja),r,jb-ja+1,50,base2)
  ss1=real(c(jpk))**2 + aimag(c(jpk))**2
  ss2=real(c(jpk+jd))**2 + aimag(c(jpk+jd))**2
  if(smax/base2 .lt. 6.0) go to 900                   !Reject non-JTMS signals
  if(ss1.lt.0.1*smax .or. ss2.lt.0.1*smax) go to 900
  if(ss1/base2.lt.1.0 .or. ss2/base2.lt.1.0) go to 900
  dfx=0.5*fpk-f0

! DF is known, now find character sync.
  r=0.
  fr=0.
  nfr=0
  rmax=0.
  do i1=1,npts-55
     z=0.
     ss=0.
     do i=1,56
        ss=ss + abs(cdat(i+i1-1))
        z=z + cdat(i+i1-1)*conjg(cwb(i))
     enddo
!     r(i1)=abs(z)/ss
     r(i1)=abs(z)
     k=mod(i1-1,56)+1
     fr(k)=fr(k)+r(i1)
     nfr(k)=nfr(k)+1
!     write(52,5001) i1,abs(z),ss,r(i1)
!5001 format(i6,3f12.3)
     if(r(i1).gt.rmax) then
        rmax=r(i1)
        jpk=i1
     endif
  enddo

  rmax=0.
  do i=1,56
     fr(i)=fr(i)/nfr(i)
     if(fr(i).gt.rmax) then
        rmax=fr(i)
        ipk=i
     endif
!     write(52,5002) i,fr(i),nfr(i)
!5002 format(i6,f12.3,i6)
  enddo

!  i1=mod(jpk-1,56)+1-3                     !### Better solution needed?  ###
!  if(i1.lt.1) i1=i1+56
!  print*,'A',jpk,i1,ipk
  i1=ipk-3
  if(i1.lt.1) i1=i1+56

  msglen=0                                 !Use ACF to find msg length
  if(npts.ge.8*56) then
     r=r-sum(r(1:npts))/npts
     acfmax=0.
     acf0=dot_product(r(1:npts),r(1:npts))
     kz=min(npts/2,28*56)
     do k=8,kz
        fac=float(npts)/(npts-k)
           acf(k)=fac*dot_product(r(1:npts),r(1+k:npts+k))/acf0
        if(acf(k).gt.acfmax) then
           acfmax=acf(k)
           kpk=k
        endif
     enddo

     amax2=0
     do i=1,8
        k=56*np(i)
        if(acf(k).gt.0.45 .and. acf(k).gt.amax2) then
           amax2=acf(k)
           msglen=np(i)
        endif
     enddo
  endif

  msg=' '                                !Decode the message
  zmax0=1.0
  s2=0.
  nchar=(npts-55-i1)/56
  if(nchar.gt.400) nchar=400
!  call searchms(cdat(i1),npts,nchar,dfx)
  do j=1,nchar
     ia=i1 + (j-1)*56
     smax=0.
     do k=0,40
        kk=k
        if(k.eq.40) kk=57
        z=0.
        do i=1,56
           z=z + cdat(ia+i)*conjg(cw(i,kk))
        enddo
        ss=abs(z)
        s2(k,j)=ss
        sm(k)=ss
        if(ss.gt.smax) then
           smax=ss
           zmax=z
           kpk=kk
        endif
     enddo
     sm(kpk)=0.
     smax2=0.
     do k=0,40
        smax2=max(smax2,sm(k))
     enddo
     msg(j:j)=cc(kpk+1:kpk+1)
     if(kpk.eq.57) msg(j:j)=' '
!     if(smax/smax2.lt.1.025) msg(j:j)=' '               !Threshold test
     phapk=atan2(aimag(zmax),real(zmax))
     dphapk=atan2(aimag(zmax/zmax0),real(zmax/zmax0))
!     write(51,3007) j,smax,phapk,dphapk
!3007 format(i5,3f12.3)
     zmax0=zmax
  enddo

  ia=max(1,nchar/3)
  ib=min(ia+27,nchar)
  msg28=msg(ia:ib)
  ndf=nint(dfx)

  if(msglen.eq.0 .or. nchar.lt.max(20,2*msglen)) then

     if(nline.le.99) nline=nline+1
     tping(nline)=t2
     call cs_lock('decodems')
!     write(*,1110) cfile6,t2,mswidth,ndb,nrpt,ndf,msg28,nmatch,nsum
     write(line(nline),1110) cfile6,t2,mswidth,ndb,nrpt,ndf,msg28
!     write(21,1110) cfile6,t2,mswidth,ndb,nrpt,ndf,msg28
1110 format(a6,f5.1,i5,i3,1x,i2.2,i5,5x,a28,2i5)
     call cs_unlock
  else if(msglen.gt.0) then
!   if(msglen.gt.0) then
     fs2=0.
     nfs2=0
     do j=1,nchar                           !Fold s2 into fs2, modulo msglen
        jj=mod(j-1,msglen)+1
        nfs2(jj)=nfs2(jj)+1
        do i=0,40
           fs2(i,jj)=fs2(i,jj) + s2(i,j)
        enddo
     enddo

     msg=' '
     do j=1,msglen
        smax=0.
        do k=0,40
           if(fs2(k,j).gt.smax) then
              smax=fs2(k,j)
              kpk=k
           endif
        enddo
        sm(kpk)=0.
        smax2=0.
        do k=0,40
           smax2=max(smax2,sm(k))
        enddo
        if(kpk.eq.40) kpk=57
        msg(j:j)=cc(kpk+1:kpk+1)
        if(kpk.eq.57) msg(j:j)=' '
     enddo
     msg28=msg(1:msglen)
     call match(mycall,msg28(1:msglen),nstart,nmatch)
     if(nmatch.ge.3 .and.nstart.gt.1) then
        msg28=msg(nstart:msglen)//msg(1:nstart-1)
     else if(nmatch.ge.3 .and.nstart.eq.1) then
        msg28=msg(1:msglen)
     else
        i3=index(msg,'  ')
        if(i3.gt.0 .and. i3.le.msglen-2) then
           msg28=msg(i3+2:msglen)//msg(1:msglen)
        else
           i3=index(msg,' ')
           if(i3.gt.0 .and. i3.lt.msglen) msg28=msg(i3:msglen)//msg(1:msglen)
        endif
     endif
     if(nline.le.99) nline=nline+1
     tping(nline)=t2
     call cs_lock('decodems')
     write(line(nline),1120) cfile6,t2,mswidth,ndb,nrpt,ndf,msg28,msglen
!     write(21,1120) cfile6,t2,mswidth,ndb,nrpt,ndf,msg28,msglen
1120 format(a6,f5.1,i5,i3,1x,i2.2,i5,5x,a28,10x,i5,'*')
!     write(*,1130) nrec,cfile6,t2,mswidth,ndb,nrpt,ndf,msg28,msglen
1130 format(i3,1x,a6,f5.1,i5,i3,1x,i2.2,i5,5x,a28,10x,i5'*')
     call cs_unlock
   endif

900 continue
  call flush(52)

  return
end subroutine decodems
