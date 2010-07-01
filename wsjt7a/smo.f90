subroutine smo(x,npts,nadd)

  real x(npts)
  common/scratch/y(10000)

  nh=nadd/2
  do i=1+nh,npts-nh
     sum=0.
     do j=-nh,nh
        sum=sum + x(i+j)
     enddo
     y(i)=sum
  enddo
  y(:nh)=0.
  y(npts-nh+1:)=0.

  fac=1.0/nadd
  do i=1,npts
     x(i)=fac*y(i)
  enddo

  return
end subroutine smo
