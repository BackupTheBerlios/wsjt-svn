      subroutine twkfreq(c4aa,n5,a)

      complex c4aa(n5)
      real a(5)
      complex w,wstep
      data twopi/6.283185307/

C  Apply AFC corrections to the c4aa data
      w=1.0
      x0=0.5*(n5+1)
      s=2.0/n5
      do i=1,n5
         x=s*(i-x0)
         if(mod(i,1000).eq.1) then
            p2=1.5*x*x - 0.5
!            p3=2.5*(x**3) - 1.5*x
!            p4=4.375*(x**4) - 3.75*(x**2) + 0.375
            dphi=(a(1) + x*a(2) + p2*a(3)) * (twopi/1378.125)
            wstep=cmplx(cos(dphi),sin(dphi))
         endif
         w=w*wstep
         c4aa(i)=w*c4aa(i)
      enddo

      return
      end
