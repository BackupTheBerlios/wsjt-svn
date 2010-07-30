subroutine abc441(msg,itone,ndits)

  character msg*28
  integer itone(84)
  integer lookup(0:91)
  character cc*43
  data cc/' 0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ.,?/#$'/
  data lookup/13, 15, 17, 46, 47, 45, 44, 12, 11, 14, &
               1,  2,  3,  4,  5,  6,  7,  8,  9, 10, &
              16, 48, 18, 19, 20, 21, 22, 23, 24, 25, &
              26, 27, 15, 47, 30, 14, 16, 42, 46, 35, &
              36, 37, 21,  0, 11, 41, 10, 13, 43,  1, &
               2,  3,  4,  5,  6,  7,  8,  9, 49, 56, &
              52, 55, 54, 12, 63, 17, 18, 19, 20, 44, &
              22, 23, 24, 25, 26, 27, 28, 29, 30, 31, &
              32, 33, 34, 35, 36, 37, 38, 39, 40, 41, &
              45, 63/

  do i=1,nmsg
     j=ichar(msg(i:i))
     if(j.lt.0 .or. j.gt.91) j=32 !Replace illegal char with blank 
     n=lookup(j)
     itone(3*i-2)=n/16 + 1
     itone(3*i-1)=mod(n/4,4) + 1
     itone(3*i)=mod(n,4) + 1
  enddo
  ndits=3*nmsg
  return
end subroutine abc441
