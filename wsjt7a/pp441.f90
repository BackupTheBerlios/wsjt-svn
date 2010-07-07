subroutine pp441(dat,jz,cfile6,tstart,t2,width,npeak,nrpt,     &
     dftolerance,frag,ncon)

! FSK441++ decoder

  parameter (NMAX=512*1024)
  parameter (MAXFFT=8192)
  real dat(NMAX)                          !Raw signal, 30 s at 11025 sps
  complex cdat(NMAX)                      !Analytic form of signal
  character cfile6*6                      !File time
  character frag*28,frag0*28              !Message fragment to be matched
  character msg*40
  character*28 msg0,msg1,msg2,msg3
  character c1*1,tok1*4,tok2*4
  complex cfrag(2100)                     !Complex waveform of message fragment
  complex ct0(25)
  complex ct1(25)
  complex ct2(25)
  complex ct3(25)
  complex z
  complex zz(0:3)
  real r(0:3)
  real rf(0:3)
  real spec(0:3)
  integer nspec(0:3)
  integer itone(84)                       !Generated tones for msg fragment
  real s(NMAX)
  real ccf(-6000:6000)
  integer dftolerance
  integer dit(-3000:3000)
  real y(0:3,-3000:3000)
  complex za(0:3,-3000:3000)
  real yf(0:3,0:86)
  integer nf(0:86)
  integer ditf(0:86)
  character c*48
  character*90 line
  common/ccom/nline,tping(100),line(100)
  common/scratch/work(NMAX)
  data c/' 123456789.,?/# $ABCD FGHIJKLMNOPQRSTUVWXY 0EZ*!'/
  data frag0/'xxxxx'/
  save frag0,cfrag,ct0,ct1,ct2,ct3,ndits

  if(frag.ne.frag0) then
! Generate waveform for message fragment
     do i=28,1,-1                          !Get length of fragment
        if(frag(i:i).ne.' ') go to 10
     enddo
10   nfrag=1000+i
     call abc441(frag,nfrag,itone,ndits)
     call gen441(itone,ndits,cfrag)

! Generate symbol (single-dit) waveforms
     call gen441(1,1,ct0)
     call gen441(2,1,ct1)
     call gen441(3,1,ct2)
     call gen441(4,1,ct3)
     frag0=frag
  endif

! Initialize variables
  nsps=25
  nsam=nsps*ndits
  mswidth=10*nint(100.0*width)
  dt=1.0/11025.0
  i0=(tstart-0.02)/dt
  if(i0.lt.1) i0=1
  npts=nint((width+0.02)/dt)+1
  npts=min(npts,jz+1-i0)
  xn=log(float(npts))/log(2.0)
  n=xn
  if(xn-n .gt.0.001) n=n+1
  nfft1=2**n
  df1=11025.0/nfft1

  call analytic(dat(i0),npts,nfft1,s,cdat)    !Convert to analytic signal

! Get DF by looking for the four FSK441 tones
  ia=dftolerance/df1
  ccfmax=0.
  do i=-ia,ia
     ccf(i)=s(i+nint(882.0/df1)) + s(i+nint(1323.0/df1)) +           &
          s(i+nint(1764.0/df1)) + s(i+nint(2205.0/df1))
  enddo
  ccf(:-ia-1)=0.
  ccf(ia+1:)=0.
  nadd=2*nint(5.0/df1)+1
  call smo(ccf(-ia),2*ia+1,work,nadd)         !Smooth CCF by nadd

  do i=-ia,ia                                 !Fin max of smoothed CCF
     if(ccf(i).gt.ccfmax) then
        ccfmax=ccf(i)
        ipk=i
        dfx=i*df1
     endif
  enddo
  ib=min(nint(220.5/df1),ia)                  !Search range +/- 220.5 Hz
  call pctile(ccf(ipk-ib),work,2*ib+1,50,base)
  ccfmax=ccfmax/base
  if(ccfmax.lt.4.0) go to 800                 !Is CCF search successful?

! We seem to have an FSK441 ping, and we know DF; now find DT.
  call tweak1(cdat,npts,-dfx,cdat)            !Mix to standard frequency

! Look for best match to "frag", find its DT
  sbest=0.
  fac=1.e-6/base
  do i=1,npts-nsam
     z=0.
!     sq=0.
     do j=1,nsam
        z=z + cdat(j+i-1)*cfrag(j)
!        sq=sq + real(cdat(j+i-1))**2 + aimag(cdat(j+i-1))**2
     enddo
!     ss=(real(z)**2 + aimag(z)**2)/sq         !??? Is this right ???
     ss=(real(z)**2 + aimag(z)**2)*fac
     if(ss.gt.sbest) then
        sbest=ss
        ibest=i
        tbest=(i+i0-1)*dt
     endif
  enddo

  if(sbest.lt.ccfmax) go to 800         !Skip if not FSK441 data

! We know DF and DT; now demodulate and decode.
  spec=0.
  nspec=0
  n=ibest/nsps - 1
  i1a=ibest-n*nsps
  n=(npts-nsps+1)/nsps - 1
  i1b=i1a+n*nsps

! Full range of potentially useful symbols is is1 to is2:
  is1=(1-ibest)/nsps
  is2=(npts-nsps+1-ibest)/nsps

! Demodulate the symbols
  do i1=i1a,i1b,nsps                      
     is=(i1-ibest)/nsps
     sq=dot_product(cdat(i1:i1+nsps-1),cdat(i1:i1+nsps-1))
     rms=sqrt(sq)
     zz(0)=dot_product(cdat(i1:i1+nsps-1),conjg(ct0))/rms
     zz(1)=dot_product(cdat(i1:i1+nsps-1),conjg(ct1))/rms
     zz(2)=dot_product(cdat(i1:i1+nsps-1),conjg(ct2))/rms
     zz(3)=dot_product(cdat(i1:i1+nsps-1),conjg(ct3))/rms
     
     rmax=0.
     do i=0,3
        r(i)=abs(zz(i))
        za(i,is)=zz(i)
        if(r(i).gt.rmax) then
           rmax=r(i)                        !Non-coherent demodulation
           ipk=i
        endif
     enddo
     do i=0,3
        if(i.ne.ipk) then
           spec(i)=spec(i)+r(i)             !Accumulate an avg 4-pt spectrum
           nspec(i)=nspec(i)+1
        endif
     enddo
  enddo
  
  do i=0,3                                  !Normalize the 4-pt spectrum
     if(nspec(i).gt.0) then
        spec(i)=spec(i)/nspec(i)
     else
        spec(i)=1.0
     endif
  enddo
  
  do i1=i1a,i1b,nsps                     !Get the dit values
     is=(i1-ibest)/nsps
     rmax=0.
     do i=0,3
        za(i,is)=za(i,is)/spec(i)
        zz(i)=za(i,is)
        r(i)=abs(zz(i))                  !Non-coherent amplitude
        y(i,is)=r(i)
        if(r(i).gt.rmax) then
           rmax=r(i)
           ipk=i
        endif
     enddo
     dit(is)=ipk
  enddo

  nc=16*dit(6) + 4*dit(7) +dit(8)
  c1=' '
  if(nc.le.47 .and. nc.ge.0) c1=c(nc+1:nc+1)
  call token(c1,n1,tok1,n2)                     !Get length encoded in msg, n2
  nc=16*dit(9) + 4*dit(10) +dit(11)
  c1=' '
  if(nc.le.47 .and. nc.ge.0) c1=c(nc+1:nc+1)
  call token(c1,n3,tok2,n4)                     !Get encoded token, if any

  msg='                                        '
  msglen=min((is2-is1+1)/3,40)         !Legth of potentially decodable text
  j=is1/3                              !Set starting location
  j=3*j

  do i=1,msglen                        !Read off the hard-decision message
     j=j+3
     nc=16*dit(j) + 4*dit(j+1) +dit(j+2)
     msg(i:i)=' '
     if(nc.le.47 .and. nc.ge.0) msg(i:i)=c(nc+1:nc+1)
  enddo

  call cs_lock('pp441')
! Probably shouldn't write multiple times, just use the best one:
  do i=1,msglen-n2-1
     if(msg(i:i+1).eq.'$!') then
        msg1=msg(i:i+n2-1)
        call dec441(msg1,msg2)
        i3=index(msg2,'$')
        if(i3.gt.1) msg2=msg2(:i3-1)
        i4=index(msg2,'!')
        if(i4.gt.1) msg2=msg2(:i4-2)
        if(msg2.ne.msg0) then
           if(ncon.ne.0) write(*,1110) cfile6,tbest,mswidth,npeak,     &
                nrpt,nint(dfx),msg2,'B',n2
1110       format(a6,f5.1,i5,i3,1x,i2.2,i5,5x,a28,10x,5x,a1,i3)
           if(nline.le.99) nline=nline+1
           tping(nline)=t2
           write(line(nline),1110) cfile6,t2,mswidth,npeak,        &
                nrpt,nint(dfx),msg2,'B',n2
           msg0=msg2
        endif
     endif
  enddo
  call cs_unlock

  if(n2.ge.4) then
! Fold the y() array to get average message
     yf=0.
     nf=0
     do j=is1,is2                          !Use weighted averages ?
        k=mod(j+300*n2,3*n2)
        yf(0,k)=yf(0,k) + y(0,j)
        yf(1,k)=yf(1,k) + y(1,j)
        yf(2,k)=yf(2,k) + y(2,j)
        yf(3,k)=yf(3,k) + y(3,j)
        nf(k)=nf(k)+1
     enddo

     do k=0,3*n2-1                     !Get dit values for averaged spectrum
        if(nf(k).gt.0) then
           rmax=0.
           do i=0,3
              rf(i)=yf(i,k)/nf(k)
              if(rf(i).gt.rmax) then
                 rmax=rf(i)
                 ipk=i
              endif
           enddo
        endif
        ditf(k)=ipk
     enddo

     j=-3
     msg='                                        '
     do i=1,n2                                !Read off the averaged message
        j=j+3
        nc=16*ditf(j) + 4*ditf(j+1) +ditf(j+2)
        msg(i:i)=' '
        if(nc.le.47 .and. nc.ge.0) msg(i:i)=c(nc+1:nc+1)
     enddo
     call token(msg(4:4),n3,tok2,n4)              !Get encoded token, if any
     msg3=msg(5:n2)
     if(n3.ge.1 .and. n3.le.3) msg3=tok2//msg3
     if(n3.ge.4 .and. n3.le.11) then
        do i=28,1,-1
           if(msg3(i:i).ne.' ') go to 100
        enddo
100     msg3=msg3(1:i)//tok2
     endif

     call cs_lock('pp441')
     if(ncon.ne.0) write(*,1110) cfile6,tbest,mswidth,npeak,nrpt,     &
          nint(dfx),msg3,'C',n2
        if(nline.le.99) nline=nline+1
        tping(nline)=t2
        write(line(nline),1110) cfile6,t2,mswidth,npeak,              &
             nrpt,nint(dfx),msg3,'C',n2
     call cs_unlock
  endif

800 continue

  return
end subroutine pp441
