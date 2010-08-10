subroutine lenms(r,npts,msglen)

! Determine length of the user message in a JTMS ping.

  real r(60000)
  real acf(1624)
  integer np(8)
  data np/5,7,11,13,17,19,23,29/        !Permissible message lengths
  save acf                              !Why necessary?  (But don't remove!)

  msglen=0                              !Use ACF to find msg length
  if(npts.ge.8*56) then
     r=r-sum(r(1:npts))/npts
     acfmax=0.
     acf0=dot_product(r(1:npts),r(1:npts))
     kz=min(npts/2,29*56)
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

     amax2=0.
     do i=1,8
        k=56*np(i)                      !Check only the permitted lengths
        if(k.gt.kz) go to 10
        if(acf(k).gt.3.5 .and. acf(k).gt.amax2) then  
           amax2=acf(k)                 !Save best value >3.5 sigma
           msglen=np(i)                 !Save message length
           kpk2=k
        endif
     enddo
10   continue
  endif

  return
end subroutine lenms
