subroutine syncms(cdat,npts,cwb,r,i1)

  complex cdat(npts)                    !Analytic signal
  complex cwb(56)                       !Complex waveform for 'space'
  real r(60000)
  complex z

  r=0.
  rmax=0.
  do j=1,npts-55
     z=0.
     ss=0.
     do i=1,56
        ss=ss + abs(cdat(i+j-1))          !Accumulate total power
        z=z + cdat(i+j-1)*conjg(cwb(i))   !Accumulate signal matching <space>
     enddo
     r(j)=abs(z)/ss                       !Goodness-of-fit to <space>
     if(r(j).gt.rmax) then
        rmax=r(j)
        jpk=j                             !Location of best match
     endif
  enddo

  i1=mod(jpk-1,56)+1

  return
end subroutine syncms
