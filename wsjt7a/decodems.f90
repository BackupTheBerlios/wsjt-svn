subroutine decodems(cdat,npts,cw,i1,nchar,s2,msg)

  complex cdat(npts)
  complex cw(56,0:63)                   !Complex waveforms for codewords
  real s2(0:63,400)
  character msg*400
  complex z,zmax
!  complex zmax0
  character cc*64
!                    1         2         3         4         5         6
!          0123456789012345678901234567890123456789012345678901234567890123
  data cc/'0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ./?-                 _     @'/

  msg=' '
!  zmax0=1.0
  do j=1,nchar
     ia=i1 + (j-1)*56
     smax=0.
     do k=0,40
        kk=k
        if(k.eq.40) kk=57
        z=0.
        do i=1,56
           z=z + cdat(ia+i)*conjg(cw(i,kk))
        enddo
        ss=abs(z)
        s2(k,j)=ss
        if(ss.gt.smax) then
           smax=ss
           zmax=z
           kpk=kk
        endif
     enddo
     msg(j:j)=cc(kpk+1:kpk+1)
     if(kpk.eq.57) msg(j:j)=' '
!     phapk=atan2(aimag(zmax),real(zmax))
!     dphapk=atan2(aimag(zmax/zmax0),real(zmax/zmax0))
!     zmax0=zmax
  enddo

  return
end subroutine decodems
