subroutine gen441(msg,nmsg,itone)

  character msg*28
  integer itone(84)
  integer lookup(0:91)
  character cc*43
  data cc/' 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ.,?/#$'/
  data lookup/13, 15, 17, 46, 47, 45, 44, 12, 11, 14, &
               1,  2,  3,  4,  5,  6,  7,  8,  9, 10, &
              16, 48, 18, 19, 20, 21, 22, 23, 24, 25, &
              26, 27, 15, 29, 30, 14, 16, 42, 46, 35, &
              36, 37, 21,  0, 11, 41, 10, 13, 43,  1, &
               2,  3,  4,  5,  6,  7,  8,  9, 49, 56, &
              52, 55, 54, 12, 63, 17, 18, 19, 20, 44, &
              22, 23, 24, 25, 26, 27, 28, 29, 30, 31, &
              32, 33, 34, 35, 36, 37, 38, 39, 40, 41, &
              45, 63/
  save

  do i=28,1,-1
     if(msg(i:i).ne.' ') go to 1
  enddo
1 nmsg=i
  if(nmsg.eq.0) nmsg=1
     
  do i=1,nmsg
     n=ichar(msg(i:i))
     if(n.eq.95) n=32
     if(n.ge.97 .and. n.le.122) n=n-32   !Promote lower case to upper
     if(n.lt.0 .or. n.gt.91) n=32        !Replace illegal char with blank 
     n=lookup(n)
     itone(3*i-2)=n/16
     itone(3*i-1)=mod(n/4,4)
     itone(3*i)=mod(n,4)
  enddo

!  do i=1,43
!     n=ichar(cc(i:i))
!     if(n.eq.95) n=32
!     if(n.ge.97 .and. n.le.122) n=n-32   !Promote lower case to upper
!     if(n.lt.0 .or. n.gt.91) n=32        !Replace illegal char with blank 
!     n=lookup(n)
!     m1=n/16
!     m2=mod(n/4,4)
!     m3=mod(n,4)
!     mm=100*m1+10*m2+m3
!     write(*,3001) i,cc(i:i),mm
!3001 format(i2,2x,a1,i5.3)
!  enddo

  return
end subroutine gen441
