subroutine syncdiana(fs0,kstep,nfreeze,mousedf,dftolerance,syncx,     &
     ipk,jpk,dfx,dtx,msglen,ccfblue)

  real fs0(1024,96)                       !Folded-for-sync spectra
  real ccfblue(-5:540)
  integer dftolerance
  integer isync(4)
  data isync/8,16,32,24/

  df=11025.0/4096.0
  nblk=24
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
  enddo
  syncx=smax/ref - 1.0

  smax=0.
  ja=jpk+16
  if(ja.gt.4*nblk) ja=ja-4*nblk
  jb=jpk+20
  if(jb.gt.4*nblk) jb=jb-4*nblk
  do i=ipk+2,ipk+56,2                         !Find User's message length
     ss=fs0(i,ja) + fs0(i+10,jb)
     if(ss.gt.smax) then
        smax=ss
        ipk2=i
     endif
  enddo
  msglen=(ipk2-ipk)/2

  return
end subroutine syncdiana
