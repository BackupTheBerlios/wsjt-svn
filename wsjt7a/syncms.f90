subroutine syncms(cdat,npts,cwb,r,i1)

! Establish character sync within a JTMS ping.

  complex cdat(npts)                    !Analytic signal
  complex cwb(56)                       !Complex waveform for 'space'
  real r(60000)
  integer hist(56),hmax(1)
  complex z

  r=0.
  hist=0
  rlim=0.9
  do j=1,npts-55
     z=0.
     ss=0.
     do i=1,56
        ss=ss + abs(cdat(i+j-1))          !Total power
        z=z + cdat(i+j-1)*conjg(cwb(i))   !Signal matching <space>
     enddo
     r(j)=abs(z)/ss                       !Goodness-of-fit to <space>
     k=mod(j-1,56)+1
     if(r(j).gt.rlim) hist(k)=hist(k)+1
  enddo

  hmax=maxloc(hist)
  i1=hmax(1)

  return
end subroutine syncms
