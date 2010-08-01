subroutine len441(r,npts,msglen)

! Determine length of the user message in a JTMS ping.

  real r(60000)
  real rf(2250)
  integer nrf(2250)
  real acf(2250)
  save acf                              !Why necessary?  (But don't remove!)

  msglen=0                              !Use ACF to find msg length
  if(npts.ge.8*75) then
     r=r-sum(r(1:npts))/npts
     acfmax=0.
     acf0=dot_product(r(1:npts),r(1:npts))
     kz=min(npts/2,30*75)
     do k=8,kz
        fac=float(npts)/(npts-k)
        acf(k)=fac*dot_product(r(1:npts),r(1+k:npts+k))/acf0
     enddo
     call hipass(acf(8),kz-7,50)

     do k=8,kz                          !Find acfmax, kpk
        if(acf(k).gt.acfmax) then
           acfmax=acf(k)
           kpk=k
        endif
     enddo

     sumsq=0.
     n=0
     do k=8,kz                          !Find rms, skipping around kpk
        if(abs(k-kpk).gt.10) then
           sumsq=sumsq+acf(k)**2
           n=n+1
        endif
     enddo
     rms=sqrt(sumsq/n)
     acf=acf/rms                        !Normalize the acf

     rewind 55
     do k=8,kz                          !Find acfmax, kpk
        write(55,3001) k,k/75.0,acf(k)
3001    format(i8,2f12.4)
     enddo

     sbest=0.
     npz=min(28,npts/150)
     do np=5,npz
        npp=75*np
        rf=0.
        nrf=0
        do i=1,npts
           j=mod(i-1,npp)+1
           rf(j)=rf(j) + r(i)
           nrf(j)=nrf(j)+1
        enddo
        smax=0.
        do j=1,npp
           rf(j)=rf(j)/nrf(j)
           if(rf(j).gt.smax) smax=rf(j)
        enddo
        if(smax.gt.sbest) then
           npbest=np
           sbest=smax
        endif
        print*,np,smax,npbest,sbest
     enddo

     call flush(55)

     iz=min(28,kz/75)
     amax2=-1.e30
     do k=375,kz
        if(acf(k).lt.amax2 .and. k.ge.kpk+10) go to 100
        if(acf(k).gt.3.8) then  
           amax2=acf(k)                 !Save best value >3.8 sigma
           kpk=k                        !Save message length
        endif
     enddo
  endif

100 msglen=nint(kpk/75.0)
  chk=msglen*75.0/kpk
  if(chk.lt.0.992 .or. chk.gt.1.008) msglen=0
  
  return
end subroutine len441
