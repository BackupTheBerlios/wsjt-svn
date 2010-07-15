subroutine genms(msg,txsnrdb,iwave,nwave)

! Generate a JTMS_2 wave file

  parameter (NZ=30*11025)
  character msg*28
  character cc*64
  integer imsg(5250)
  integer*2 iwave(NZ)
  integer*2 jwave(NZ)
  integer icw(0:63)
  logical first
  real*8 txsnrdb,t
  complex cw
  common/mscom/cw(63,0:63)
  data idum/-1/,first/.true./
  data icw/                                                 &
    z'000000',z'1C8AF6',z'12CF8D',z'0E457B',z'15ED30',z'0967C6',  &
    z'0722BD',z'1BA84B',z'167C6E',z'0AF698',z'04B3E3',z'183915',  &
    z'03915E',z'1F1BA8',z'115ED3',z'0DD425',z'0B3E37',z'17B4C1',  &
    z'19F1BA',z'057B4C',z'1ED307',z'0259F1',z'0C1C8A',z'10967C',  &
    z'1D4259',z'01C8AF',z'0F8DD4',z'130722',z'08AF69',z'14259F',  &
    z'1A60E4',z'06EA12',z'1915ED',z'059F1B',z'0BDA60',z'175096',  &
    z'0CF8DD',z'10722B',z'1E3750',z'02BDA6',z'0F6983',z'13E375',  &
    z'1DA60E',z'012CF8',z'1A84B3',z'060E45',z'084B3E',z'14C1C8',  &
    z'122BDA',z'0EA12C',z'00E457',z'1C6EA1',z'07C6EA',z'1B4C1C',  &
    z'150967',z'098391',z'0457B4',z'18DD42',z'169839',z'0A12CF',  &
    z'11BA84',z'0D3072',z'037509',z'1FFFFF'/

!           1234567890123456789012345678901234567890123456789012345678901234
  data cc/' 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ./?-                      @'/
  save first,idum

  if(first) then
! Generate 64 single-character waveforms.
     twopi=8.0*atan(1.0)
     dt=1.0/11025.0
     dpha0=twopi*dt*800.0
     dpha1=twopi*dt*2000.0
     do j=0,63
        pha=0.
        n=icw(j)
        do i=1,63
           i3=(i+2)/3
           if(iand(n,1).eq.0) then
              pha=pha+dpha0
           else
              pha=pha+dpha1
           endif
           cw(i,j)=cmplx(cos(pha),sin(pha))
           if(i.eq.3*i3) n=n/2
        enddo
     enddo
     first=.false.
  endif

  sig=sqrt(10.0**(0.1*txsnrdb))
  do i=28,1,-1                                 !Find user's message length
     if(msg(i:i).ne.' ') go to 1
  enddo
1 msglen=i+1                                   !Add 1 for EOM character

! Convert message to an array of waveform numbers, imsg(1:msglen).
  do j=1,msglen
     imsg(j)=0
     do i=1,64
        if(msg(j:j).eq.cc(i:i)) then
           imsg(j)=i-1
           go to 5
        endif
     enddo
5    continue
  enddo

! Replicate imsg contents sufficiently for a 30-second transmission
  nrpt=5250/msglen
  j=msglen
  do n=2,nrpt
     do i=1,msglen
        j=j+1
        imsg(j)=imsg(i)
     enddo
  enddo
  nchar=j                             !Total # characters to be transmitted

! Generate the wave file
  k=0
  do n=1,nchar
     do i=1,63
        k=k+1
        iwave(k)=32767.0*real(cw(i,imsg(n)))
     enddo
  enddo
  nwave=k

  if(txsnrdb.lt.40.0) then
! ###  Make some pings (for tests only) ###
     do i=1,nwave
        iping=i/(3*12000)
        if(iping.ne.iping0) then
           ip=mod(iping,3)
           w=0.05*(ip+1)
           ig=(iping-1)/3
           amp=sqrt((3.0-ig)/3.0)
           t0=dt*(iping+0.5)*(3*12000)
           iping0=iping
        endif
        t=(i*dt-t0)/w
        if(t.lt.0.d0 .and. t.lt.10.d0) then
           fac=0.
        else
           fac=2.718*t*dexp(-t)
        endif
!        x=sig*fac*amp*iwave(i)/32767.0 + gran(idum)
!        jwave(i)=nint(1000.0*x)
        iwave(i)=nint(fac*amp*iwave(i))
     enddo
  endif

  return
end subroutine genms
