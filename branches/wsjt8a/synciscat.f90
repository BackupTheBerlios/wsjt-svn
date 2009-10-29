subroutine synciscat(dat,jz,DFTolerance,NFreeze,MouseDF,dtx,dfx,      &
     snrx,snrsync,isbest,ccfblue,ccfred,s2)

! Synchronizes ISCAT data, finding the best-fit DT and DF.  

  parameter (NFFTMAX=1024)         !Max length of FFTs
  parameter (NHMAX=NFFTMAX/2)      !Max length of power spectra
  parameter (NSMAX=292)            !Max number of quarter-symbol steps
  integer DFTolerance              !Range of DF search
  real dat(jz)                     !Raw data, downsampled to 6 kHz
  real xs1(NHMAX)
  real s1(NHMAX,NSMAX)             !2d spectrum, stepped by half-symbols
  real s2(64,63)                   !2d spectrum, synced data symbols only
  real x(NFFTMAX)                  !Temp array for computing FFTs
  real ccfblue(-5:540)             !CCF with pseudorandom sequence
  real ccfred(-224:224)            !Peak of ccfblue, as function of freq
  integer isync(10,3)
  integer ic10(10)
  data ic10/0,1,3,7,4,9,8,6,2,5/     !10x10 Costas array

! Set up the ISCAT sync pattern
  nsync=10
  do i=1,10
     isync(i,1)=ic10(i)
     isync(i,2)=ic10(11-i)
     isync(i,3)=9-ic10(i)
  enddo
  nsym=nsync+63

! Do FFTs of twice symbol length, stepped by quarter symbols.  
  nfft=1024
  nh=nfft/2
  nsteps=4*(jz-NH)/nh
  nsteps=nsteps/292
  nsteps=nsteps*292
  kstep=nh/4
  df=12000.0/nfft

! Compute power spectrum for each quarter-symbol step
  s1=0.
  do j=1,nsteps
     k=(j-1)*kstep + 1
     jj=mod(j-1,292)+1
     do i=1,nh
        x(i)=dat(k+i-1)
        x(i+nh)=0.
     enddo
     call ps(x,nfft,xs1)
     call add(s1(1,jj),xs1,s1(1,jj),nh)
  enddo

! Determine the search range in frequency
  famin=3.
  fbmax=1200.
  f0=700.0
  fa=famin
  fb=fbmax
  if(NFreeze.eq.1) then
     fa=max(famin,f0+MouseDF-DFTolerance)
     fb=min(fbmax,f0+MouseDF+DFTolerance)
  else
     fa=max(famin,f0+MouseDF-600)
     fb=min(fbmax,f0+MouseDF+600)
  endif
  ia=fa/df
  ib=fb/df
  i0=nint(f0/df)

! Find best frequency bin and best sync pattern
  syncbest=-1.e30
  ss=0.
  nss=0
  do i=ia,ib
     smax=-1.e30
     do lag=0,291
        sum1=0.
        sum2=0.
        sum3=0.
        do j=1,nsync
           j0=4*j - 3 + lag
           jj0=mod(j0-1,292)+1
           sum1=sum1 + s1(i+2*isync(j,1),jj0)
           sum2=sum2 + s1(i+2*isync(j,2),jj0)
           sum3=sum3 + s1(i+2*isync(j,3),jj0)
        enddo
        ccf1=sum1/nsync
        ccf2=sum2/nsync
        ccf3=sum3/nsync
        if(ccf1.gt.smax) then
           smax=ccf1
           ispk=1
        endif
        if(ccf2.gt.smax) then
           smax=ccf2
           ispk=2
        endif
        if(ccf3.gt.smax) then
           smax=ccf3
           ispk=3
        endif
     enddo

     j=i-i0
     if(abs(j).le.224) then
        ccfred(i-i0)=smax
        ss=ss+smax
        nss=nss+1
     endif
     if(smax.gt.syncbest) then
        syncbest=smax
        ipk=i
        isbest=ispk
     endif
  enddo

  avered=ss/nss

! Once more, using best frequency and best sync pattern:
  ccfblue=0.
  syncbest=-1.e30
  do lag=0,291
     sum=0.
     do j=1,nsync
        j0=4*j - 3 + lag
        jj0=mod(j0-1,292)+1
        sum=sum + s1(ipk+2*isync(j,isbest),jj0)
     enddo
     ccfblue(lag)=sum/nsync
     if(ccfblue(lag).gt.syncbest) then
        lagpk=lag
        syncbest=ccfblue(lag)
     endif
  enddo

  sum=0.
  nsum=0
  do j=0,291
     if(abs(j-lagpk).gt.2) then
        sum=sum + ccfblue(j)
        nsum=nsum + 1
     endif
  enddo
  ave=sum/nsum
  do j=0,291
     ccfblue(j)=ccfblue(j)-ave
  enddo

  snrsync=syncbest/ave - 1.0
  snrx=-31.
  if(syncbest.gt.1.0) snrx=db(snrsync) - 20.0
  dtstep=kstep/12000.d0
  dtx=dtstep*lagpk
  dfx=(ipk-i0)*df

  ja=nint(dftolerance/df)
  do j=-ja,ja
     ccfred(j)=0.5*(ccfred(j)-avered)
  enddo
  ccfred(-224:-ja)=0.
  ccfred(ja:224)=0.

! Copy synchronized data symbols from s1 into s2
  do j=1,63
     j0=4*j - 3 + lagpk + 40
     jj0=mod(j0-1,292)+1
     do i=1,64
        s2(i,j)=s1(ipk+2*(i-1),jj0)
     enddo
  enddo

  return
end subroutine synciscat