subroutine lendiana(fs0,ipk,jpk,msglen)
  real fs0(1024,108)                       !108 = 96 + 12 extra

  nblk=24
  smax=0.
  ja=jpk+16
  if(ja.gt.4*nblk) ja=ja-4*nblk
  jb=jpk+20
  if(jb.gt.4*nblk) jb=jb-4*nblk
  do i=ipk,ipk+60,2                         !Find User's message length
     ss=fs0(i,ja) + fs0(i+10,jb)
     if(ss.gt.smax) then
        smax=ss
        ipk2=i
     endif
  enddo

  msglen=(ipk2-ipk)/2

  return
end subroutine lendiana
