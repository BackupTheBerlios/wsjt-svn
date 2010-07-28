subroutine decodems(dat,npts,cfile6,t2,mswidth,ndb,nrpt,Nfreeze,       &
     DFTolerance,MouseDF,mycall,hiscall)

! Decode a JTMS ping

  parameter (NZ=30*11025)
  real dat(npts)                        !Raw data
  complex cdat(NZ)                      !Analytic form of signal
  character*6 cfile6                    !FileID
  integer DFTolerance
  character*12 mycall,hiscall
  real s(NZ)                            !Power spectrum
  real sm(0:63)
  real s2(0:63,400)
  real fs2(0:63,29)
  integer nfs2(29)
  real r(60000)
  complex cw(56,0:63)                   !Complex waveforms for codewords
  complex cwb(56)                       !Complex waveform for 'space'
  complex z,zmax,zmax0
  logical first
  character msg*400,msg28*28,frag*28
  character cc*64
  character*90 line
  common/ccom/nline,tping(100),line(100)
!                    1         2         3         4         5         6
!          0123456789012345678901234567890123456789012345678901234567890123
  data cc/'0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ./?-                 _     @'/
  data first/.true./
  data nsum/0/,nrec/0/
  save nsum,nrec
  save first,smax,cw,cwb
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

  call msdf(cdat,npts,nfft1,f0,mousedf,dftolerance,dfx,ferr)    !Get DF

!  write(*,2001) t2,fpk1,fpk2,ferr
!2001 format(f6.1,2f10.1,f10.3)
  if(abs(ferr).gt.0.002) go to 900           !Reject non-JTMS signals
  call tweak1(cdat,npts,-dfx,cdat)           !Mix to standard frequency

! DF is known, now establish character sync.

  call syncms(cdat,npts,cwb,r,i1)

  call lenms(r,npts,msglen)

  msg=' '                                !Decode the message
  zmax0=1.0
  s2=0.
  nchar=(npts-55-i1)/56
  if(nchar.gt.400) nchar=400

  frag=' '//mycall
  call searchms(cdat(i1),npts-i1,frag,ndi1,rmax1)
  frag=' '//hiscall
  call searchms(cdat(i1),npts-i1,frag,ndi2,rmax2)
  frag=' CQ'
  call searchms(cdat(i1),npts-i1,frag,ndi3,rmax3)

!  write(*,2002) t2,ndi1,ndi2,ndi3,rmax1,rmax2,rmax3
!2002 format(f7.1,3i8,3f10.2)
  ndi=99
  if(max(rmax1,rmax2,rmax3).ge.0.6) then
     if(max(rmax1,rmax2,rmax3).eq.rmax1 .and. abs(ndi1).le.5) ndi=ndi1
     if(max(rmax1,rmax2,rmax3).eq.rmax2 .and. abs(ndi2).le.5) ndi=ndi2
     if(max(rmax1,rmax2,rmax3).eq.rmax3 .and. abs(ndi3).le.5) ndi=ndi3
     if(abs(ndi).le.5) i1=i1+ndi
  endif

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
!     write(*,1110) cfile6,t2,mswidth,ndb,nrpt,ndf,msg28
     if(abs(ndi).le.5) then
        write(line(nline),1110) cfile6,t2,mswidth,ndb,nrpt,ndf,msg28,ndi
1110    format(a6,f5.1,i5,i3,1x,i2.2,i5,5x,a28,12x,i3)
     else
        write(line(nline),1110) cfile6,t2,mswidth,ndb,nrpt,ndf,msg28
     endif
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
     call match(' CQ ',msg28(1:msglen),nstart2,nmatch2)
     if(nmatch.ge.3 .and.nstart.gt.1) then
        msg28=msg(nstart:msglen)//msg(1:nstart-1)
     else if(nmatch.ge.3 .and.nstart.eq.1) then
        msg28=msg(1:msglen)
     else if(nmatch2.ge.3 .and.nstart2.gt.1) then
        msg28=msg(nstart2:msglen)//msg(1:nstart2-1)
     else if(nmatch2.ge.3 .and.nstart2.eq.1) then
        msg28=msg(2:msglen)
     else
        i3=index(msg,'  ')
        if(i3.gt.0 .and. i3.le.msglen-2) then
           msg28=msg(i3+2:msglen)//msg(1:msglen)
        else
           i3=index(msg,' ')
           if(i3.gt.0 .and. i3.lt.msglen) msg28=msg(i3:msglen)//msg(1:msglen)
        endif
     endif
     if(msg28(1:1).eq.' ') msg28=msg28(2:)
     if(nline.le.99) nline=nline+1
     tping(nline)=t2
     call cs_lock('decodems')
     if(abs(ndi).le.5) then
        write(line(nline),1120) cfile6,t2,mswidth,ndb,nrpt,ndf,msg28,msglen,ndi
1120    format(a6,f5.1,i5,i3,1x,i2.2,i5,5x,a28,8x,i3,'*',i3)
     else
        write(line(nline),1120) cfile6,t2,mswidth,ndb,nrpt,ndf,msg28,msglen
     endif
!     write(*,1130) nrec,cfile6,t2,mswidth,ndb,nrpt,ndf,msg28,msglen
!1130 format(i3,1x,a6,f5.1,i5,i3,1x,i2.2,i5,5x,a28,10x,i5'*')
     call cs_unlock
   endif

900 continue
  call flush(52)

  return
end subroutine decodems

