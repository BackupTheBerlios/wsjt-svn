subroutine syncms(cdat,npts,cwb,r,i1)

  complex cdat(npts)                    !Analytic signal
  complex cwb(56)                       !Complex waveform for 'space'
  real r(60000)
  real fr(56)
  integer nfr(56)
  complex z

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
     r(i1)=abs(z)
     k=mod(i1-1,56)+1
     fr(k)=fr(k)+r(i1)
     nfr(k)=nfr(k)+1
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
  enddo

  i1=ipk
  if(i1.lt.1) i1=i1+56

  return
end subroutine syncms
