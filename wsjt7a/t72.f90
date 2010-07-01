program t72

! Tests experimental FSK441 decoder

  parameter (NMAX=512*1024)
  parameter (MAXFFT=8192)
  real dat0(NMAX)                          !Raw signal, 30 s at 11025 sps
  real dat(NMAX)                          !Raw signal, 30 s at 11025 sps
  complex cdat(NMAX)                      !Analytic form of signal
  character arg*12                        !Command-line argument
  character cfile6*6                      !File time
  character frag*28                       !Message fragment to be matched
  character msg*40
  complex cfrag(2100)                     !Complex waveform of message fragment
  complex ct0(25)
  complex ct1(25)
  complex ct2(25)
  complex ct3(25)
  complex z
  complex zz(0:3)
  real r(0:3)
  real rf(0:3)
  real th(0:3)
  real s(NMAX)
  real ccf(-2400:2400)
  integer itone(84)                       !Generated tones for msg fragment
  integer dftolerance
  integer dit(-3000:3000)
  real y(0:3,-3000:3000)
  real yf(0:3,0:86)
  integer nf(0:86)
  integer ditf(0:86)
  logical pick
  real pingdat(3,100)                     !Detected pings
  character c*48
  common/scratch/work(10000)
  data c/' 123456789.,?/# $ABCD FGHIJKLMNOPQRSTUVWXY 0EZ  '/

  nargs=iargc()
  if(nargs.ne.2) then
     print*,'Usage: t72 nfile frag'
     go to 999
  endif
  call getarg(1,arg)
  read(arg,*) nfile
  call getarg(2,frag)
  open(72,file='dat.72',form='unformatted',status='old')

  do i=28,1,-1                          !Get length of fragment
     if(frag(i:i).ne.' ') go to 10
  enddo
10 nfrag=i

! Generate waveform for message fragment
  call abc441(frag,nfrag,itone,ndits)
  call gen441(itone,ndits,cfrag)

! Generate symbol (single-dit) waveforms
  call gen441(1,1,ct0)
  call gen441(2,1,ct1)
  call gen441(3,1,ct2)
  call gen441(4,1,ct3)

! Initialize variables
  jpk=0
  dt=1.0/11025.0
  minsigdb=2
  minwidth=40
  dftolerance=400
  pick=.false.
  nsps=25
  nsam=nsps*ndits
  xn=log(float(nsam))/log(2.0)
  n=xn
  if(xn-n .gt.0.001) n=n+1
  nfft=2**n
  nh=nfft/2
  df=11025.0/nfft
  ja=dftolerance/df
  ibest=0
  dfx=0

  do ifile=1,nfile
     read(72,end=999) jz,nz,cfile6,(dat0(j),j=1,jz)
     if(ifile.ne.nfile .and. nfile.ne.999) go to 900

! If necessary, correct sample-rate errors
     if(ifile.eq.3 .or. ifile.eq.4 .or. ifile.eq.5) then
        j=0
        do i=1,jz
           j=j+1
           if(mod(i,147).eq.0) then
              j=j-1
           endif
           dat(i)=dat0(j)
        enddo
        dat(j:jz)=0.
     else
        dat(1:jz)=dat0(1:jz)
     endif

     call ping441(dat,jz,nz,MinSigdB,MinWidth,pick,pingdat,nping)   !Find pings

     do iping=1,nping                        !Process each ping
        tstart=pingdat(1,iping)
        width=pingdat(2,iping)
        peak=pingdat(3,iping)
        mswidth=10*nint(100.0*width)
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
        ia=dftolerance/df1

        ccfmax=0.
        do i=-ia,ia
           ccf(i)=s(i+nint(882.0/df1)) + s(i+nint(1323.0/df1)) +           &
                s(i+nint(1764.0/df1)) + s(i+nint(2205.0/df1))
        enddo
        ccf(:-ia-1)=0.
        ccf(ia+1:)=0.
        nadd=2*nint(5.0/df1)+1
        call smo(ccf(-ia),2*ia+1,nadd)

        do i=-ia,ia
!           write(14,3001) i*df1,ccf(i)
!3001       format(2f10.3)
           if(ccf(i).gt.ccfmax) then
              ccfmax=ccf(i)
              ipk=i
              dfx=i*df1
           endif
        enddo
        ib=min(nint(220.5/df1),ia)
        call pctile(ccf(ipk-ib),work,2*ib+1,50,base)
        ccfmax=ccfmax/base
        if(ccfmax.lt.4.0) go to 800

! We seem to have an FSK441 ping, and we know DF; now find DT.
        call tweak1(cdat,npts,-dfx,cdat)
        sbest=0.
        do i=1,npts-nsam                       !Look for matches to frag
           z=0.
           sq=0.
           do j=1,nsam
              z=z + cdat(j+i-1)*cfrag(j)
              sq=sq + real(cdat(j+i-1))**2 + aimag(cdat(j+i-1))**2
           enddo
           ss=(real(z)**2 + aimag(z)**2)/sq
!           write(15,3002) i,(i+i0-1)*dt,ss
!3002       format(i8,f12.6,f12.3)
           if(ss.gt.sbest) then
              sbest=ss
              ibest=i
              tbest=(i+i0-1)*dt
           endif
        enddo

! We know DF and DT; now decode the message.
        k1=0
        do i1=ibest,npts-nsps+1,nsps                  !Positive is
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
              y(i,is)=r(i)
              th(i)=atan2(aimag(zz(i)),real(zz(i)))
              if(r(i).gt.rmax) then
                 rmax=r(i)
                 ipk=i
              endif
           enddo
           k1=k1+1
           dit(is)=ipk
           write(17,3004) is,zz
3004       format(i3,4(f8.2,f6.2))
           write(18,3005) is,ipk,rmax,th(ipk),is/3
3005       format(2i3,2f8.2,i6)
           if(ipk.eq.0) write(20,3005) is,ipk,rmax,th(ipk),is/3
           if(ipk.eq.1) write(21,3005) is,ipk,rmax,th(ipk),is/3
           if(ipk.eq.2) write(22,3005) is,ipk,rmax,th(ipk),is/3
           if(ipk.eq.3) write(23,3005) is,ipk,rmax,th(ipk),is/3
        enddo

        k2=0
        do i1=ibest-nsps,1,-nsps                             !Negative is
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
              y(i,is)=r(i)
              th(i)=atan2(aimag(zz(i)),real(zz(i)))
              if(r(i).gt.rmax) then
                 rmax=r(i)
                 ipk=i
              endif
           enddo
           k2=k2+1
           dit(is)=ipk
           write(17,3004) is,zz
           write(18,3005) is,ipk,rmax,th(ipk),(is-2)/3
           if(ipk.eq.0) write(20,3005) is,ipk,rmax,th(ipk),(is-2)/3
           if(ipk.eq.1) write(21,3005) is,ipk,rmax,th(ipk),(is-2)/3
           if(ipk.eq.2) write(22,3005) is,ipk,rmax,th(ipk),(is-2)/3
           if(ipk.eq.3) write(23,3005) is,ipk,rmax,th(ipk),(is-2)/3
        enddo

! Get message length
        smax=0.
        sum0=1.
        do lag=0,28
           sum=0.
           do j=is,k1-3*28
              sum=sum + y(0,j)*y(0,j+3*lag) + y(1,j)*y(1,j+3*lag) +      &
                        y(2,j)*y(2,j+3*lag) + y(3,j)*y(3,j+3*lag) 
           enddo
           if(lag.eq.0) then
              sum0=sum
              sum=1.0
           else
              sum=sum/sum0
              if(sum.gt.smax) then
                 smax=sum
                 len=lag
              endif
!              if(smax.gt.0.85) go to 20
           endif
           write(19,3006) lag,sum
3006       format(i3,f12.3)
        enddo
!20      continue

        msg='                                        '
        msglen=min((k1+k2)/3,40)
        j=is/3
        j=3*j

        do i=1,msglen
           j=j+3
           nc=16*dit(j) + 4*dit(j+1) +dit(j+2)
           msg(i:i)=' '
           if(nc.le.47 .and. nc.ge.0) msg(i:i)=c(nc+1:nc+1)
        enddo

        write(*,1110) cfile6,tbest,mswidth,nint(dfx),     &
             ccfmax,sbest,len,msg(:msglen)
1110    format(a6,f5.1,i5,i5,2f6.1,i3,2x,a)

! Fold the y() array
        yf=0.
        nf=0
        do j=is,k1
           k=mod(j+300*len,3*len)
           yf(0,k)=yf(0,k) + y(0,j)
           yf(1,k)=yf(1,k) + y(1,j)
           yf(2,k)=yf(2,k) + y(2,j)
           yf(3,k)=yf(3,k) + y(3,j)
           nf(k)=nf(k)+1
        enddo

        do k=0,3*len-1
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
!        write(*,4001) (ditf(k),k=0,3*len-1)
!4001    format(10(i2,2i1))

        msglen=len
        j=-3
        msg='                                        '
        do i=1,msglen
           j=j+3
           nc=16*ditf(j) + 4*ditf(j+1) +ditf(j+2)
           msg(i:i)=' '
           if(nc.le.47 .and. nc.ge.0) msg(i:i)=c(nc+1:nc+1)
        enddo

        write(*,1120) msg(:msglen)
1120    format(37x,a)


800     continue
     enddo
900  continue
  enddo

999 end program t72
