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
  real spec(0:3)
  integer nspec(0:3)
  real s(NMAX)
  real ccf(-6000:6000)
  integer itone(84)                       !Generated tones for msg fragment
  integer dftolerance
  integer dit(-3000:3000)
  real y(0:3,-3000:3000)
  complex za(0:3,-3000:3000)
  real yf(0:3,0:86)
  integer nf(0:86)
  integer ditf(0:86)
  logical pick
  real pingdat(3,100)                     !Detected pings
  character c*48
  common/scratch/work(NMAX)
  data c/' 123456789.,?/# $ABCD FGHIJKLMNOPQRSTUVWXY 0EZ  '/

  nargs=iargc()
  if(nargs.ne.2) then
     print*,'Usage: t72 nfile frag'
     go to 999
  endif
  call getarg(1,arg)
  read(arg,*) nfile
  call getarg(2,frag)
  open(72,file='dat.72a',form='unformatted',status='old')

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
  ipk=0
  twopi=8.0*atan(1.0)
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
        call smo(ccf(-ia),2*ia+1,work,nadd)

        do i=-ia,ia
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
           if(ss.gt.sbest) then
              sbest=ss
              ibest=i
              tbest=(i+i0-1)*dt
           endif
        enddo

        if(sbest.lt.ccfmax) go to 800         !Skip if not FSK441 data

! We know DF and DT; now decode the message.
        spec=0.
        nspec=0
        n=ibest/nsps - 1
        i1a=ibest-n*nsps
        n=(npts-nsps+1)/nsps - 1
        i1b=i1a+n*nsps
        is1=(1-ibest)/nsps
        is2=(npts-nsps+1-ibest)/nsps

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
                 rmax=r(i)
                 ipk=i
              endif
           enddo
           do i=0,3
              if(i.ne.ipk) then
                 spec(i)=spec(i)+r(i)
                 nspec(i)=nspec(i)+1
              endif
           enddo
        enddo

        do i=0,3
           if(nspec(i).gt.0) then
              spec(i)=spec(i)/nspec(i)
           else
              spec(i)=1.0
           endif
        enddo

        do i1=i1a,i1b,nsps
           is=(i1-ibest)/nsps
           rmax=0.
           do i=0,3
              za(i,is)=za(i,is)/spec(i)
              zz(i)=za(i,is)
              r(i)=abs(zz(i))
              y(i,is)=r(i)
              if(r(i).gt.rmax) then
                 rmax=r(i)
                 ipk=i
              endif
           enddo
           dit(is)=ipk
        enddo

! Get message length
        smax=0.
        sum0=1.
        lenavg=1
        do lag=0,28
           sum=0.
           do j=is1,is2-3*28
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
                 lenavg=lag
              endif
           endif
        enddo
        if(lenavg.eq.1) go to 800

        msg='                                        '
        msglen=min((is2-is1+1)/3,40)
        j=is1/3
        j=3*j

        do i=1,msglen
           j=j+3
           nc=16*dit(j) + 4*dit(j+1) +dit(j+2)
           msg(i:i)=' '
           if(nc.le.47 .and. nc.ge.0) msg(i:i)=c(nc+1:nc+1)
        enddo

        write(*,1110) cfile6,tbest,mswidth,nint(dfx),     &
             ccfmax,sbest,lenavg,msg(:msglen)
1110    format(a6,f5.1,i5,i5,2f6.1,i3,2x,a)

! Fold the y() array
        yf=0.
        nf=0
        do j=is1,is2
           k=mod(j+300*lenavg,3*lenavg)
           yf(0,k)=yf(0,k) + y(0,j)
           yf(1,k)=yf(1,k) + y(1,j)
           yf(2,k)=yf(2,k) + y(2,j)
           yf(3,k)=yf(3,k) + y(3,j)
           nf(k)=nf(k)+1
        enddo

        do k=0,3*lenavg-1
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
!        write(*,4001) (ditf(k),k=0,3*lenavg-1)
!4001    format(10(i2,2i1))

        j=-3
        msg='                                        '
        do i=1,lenavg
           j=j+3
           nc=16*ditf(j) + 4*ditf(j+1) +ditf(j+2)
           msg(i:i)=' '
           if(nc.le.47 .and. nc.ge.0) msg(i:i)=c(nc+1:nc+1)
        enddo

        write(*,1120) msg(:lenavg)
1120    format(37x,a)


800     continue
     enddo
900  continue
  enddo

999 end program t72
