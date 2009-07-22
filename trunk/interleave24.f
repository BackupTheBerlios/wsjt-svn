      subroutine interleave24(id,ndir)
      integer*1 id(0:205),itmp(0:205)
      integer j0(0:205)
      character*8 c0,c1
      logical first
      data first/.true./
      save first

      if(first) then
         call cs_lock('interleave24')
         k=-1
         do i=0,255
            write(c0,1001) i
 1001       format(b8.8)
            c1=c0(8:8)//c0(7:7)//c0(6:6)//c0(5:5)//c0(4:4)//c0(3:3)//
     +         c0(2:2)//c0(1:1)
            read(c1,1001) n
            if(n.le.205) then
               k=k+1
               j0(k)=n
            endif
         enddo
         first=.false.
         call cs_unlock
      endif

      if(ndir.eq.1) then
         do i=0,205
            itmp(j0(i))=id(i)
         enddo
      else
         do i=0,205
            itmp(i)=id(j0(i))
         enddo
      endif

      do i=0,205
         id(i)=itmp(i)
      enddo

      return
      end
