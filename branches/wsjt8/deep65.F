      subroutine deep65(s3,mode65,neme,flip,mycall,hiscall,
     +  hisgrid,decoded,qual)

      parameter (MAXCALLS=7000,MAXRPT=63)
      real s3(64,63)
      character callsign*12,grid*4,message*22,hisgrid*6,c*1,ceme*3
      character*12 mycall,hiscall
      character*22 decoded
      character*22 testmsg(2*MAXCALLS + 2 + MAXRPT)
      character*15 callgrid(MAXCALLS)
      character*180 line
      character*4 rpt(MAXRPT)
      integer ncode(63,2*MAXCALLS + 2 + MAXRPT)
!      character*1 cflag
      character*36 cc
      integer dgen(12)
      common/tmp9/ mrs(63),mrs2(63)
      data cc/'0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
#ifdef CVF
C  This prevents some optimizations that break this subroutine.
      volatile p1,p2,bias
#endif

      data neme0/-99/
      data rpt/'-01','-02','-03','-04','-05',
     +         '-06','-07','-08','-09','-10',
     +         '-11','-12','-13','-14','-15',
     +         '-16','-17','-18','-19','-20',
     +         '-21','-22','-23','-24','-25',
     +         '-26','-27','-28','-29','-30',
     +         'R-01','R-02','R-03','R-04','R-05',
     +         'R-06','R-07','R-08','R-09','R-10',
     +         'R-11','R-12','R-13','R-14','R-15',
     +         'R-16','R-17','R-18','R-19','R-20',
     +         'R-21','R-22','R-23','R-24','R-25',
     +         'R-26','R-27','R-28','R-29','R-30',
     +         'RO','RRR','73'/
      data iseed/1000000001/,dgen(1)/999/
      save

      rewind 23
      k=0
      icall=0
      do n=1,MAXCALLS
         if(n.eq.1) then
            callsign=hiscall
            do i=4,12
               if(ichar(callsign(i:i)).eq.0) callsign(i:i)=' '
            enddo
            grid=hisgrid(1:4)
            if(ichar(grid(3:3)).eq.0) grid(3:3)=' '
            if(ichar(grid(4:4)).eq.0) grid(4:4)=' '
         else
            read(23,1002,end=20) line
 1002       format (A80)
            if(line(1:4).eq.'ZZZZ') go to 20
            if(line(1:2).eq.'//') go to 10
            i1=index(line,',')
            if(i1.lt.4) go to 10
            i2=index(line(i1+1:),',')
            if(i2.lt.5) go to 10
            i2=i2+i1
            i3=index(line(i2+1:),',')
            if(i3.lt.1) i3=index(line(i2+1:),' ')
            i3=i2+i3
            callsign=line(1:i1-1)
            grid=line(i1+1:i2-1)
            ceme=line(i2+1:i3-1)
            if(callsign.eq.hiscall .and. grid.eq.hisgrid(1:4)) go to 10
            if(neme.eq.1 .and. ceme.ne.'EME') go to 10
         endif

         icall=icall+1
         j1=index(mycall,' ') - 1
         if(j1.le.-1) j1=12
         if(j1.lt.3) j1=6
         j2=index(callsign,' ') - 1
         if(j2.le.-1) j2=12
         if(j2.lt.3) j2=6
         j3=index(mycall,'/')
         j4=index(callsign,'/')
         callgrid(icall)=callsign(1:j2)

         mz=1
         if(n.eq.1 .and. j3.lt.1 .and. j4.lt.1 .and.
     +      flip.gt.0.0 .and. callsign(1:6).ne.'      ') mz=MAXRPT+1
C  Test for messages with MyCall + HisCall + report
         do m=1,mz
            if(m.gt.1) grid=rpt(m-1)
            if(j3.lt.1 .and.j4.lt.1) 
     +         callgrid(icall)=callsign(1:j2)//' '//grid
            message=mycall(1:j1)//' '//callgrid(icall)
            k=k+1
            testmsg(k)=message
            call encode65(message,ncode(1,k))
C  Insert CQ message unless sync=OOO (flip=-1).
            if(m.eq.1 .and. flip.gt.0.0) then
               message='CQ '//callgrid(icall)
               k=k+1
               testmsg(k)=message
               call encode65(message,ncode(1,k))
            endif
         enddo
 10      continue
      enddo
 20   ntot=k
      neme0=neme

      ref0=0.
      do j=1,63
         ref0=ref0 + s3(mrs(j),j)
      enddo

      p1=0.
      p2=0.
      do k=1,ntot
         sum=0.
         ref=ref0
         do j=1,63
            i=ncode(j,k)+1
            sum=sum + s3(i,j)
            if(i.eq.mrs(j)) then
               ref=ref - s3(i,j) + s3(mrs2(j),j)
            endif
         enddo
         p=sum/ref
         if(p.gt.p2) then
            if(p.gt.p1) then
               p2=p1
!               ip2=ip1
               p1=p
               ip1=k
            else
               p2=p
!               ip2=k
            endif
         endif
      enddo

      if(mode65.eq.1) bias=max(1.12*p2,0.335)
      if(mode65.eq.2) bias=max(1.08*p2,0.405)
      if(mode65.ge.4) bias=max(1.04*p2,0.505)
      qual=100.0*(p1-bias)
      if(qual.lt.0.0) qual=0.0
      decoded='                      '
      c=' '

      if(qual.gt.1.0) then
         if(ip1.le.ntot) then
            if(qual.lt.6.0) c='?'
            decoded=testmsg(ip1)
         else
            i=ip1-ntot-1
            i1=i/(36**3)
            i2=(i - i1*36**3)/(36**2)
            i3=(i - i1*36**3 - i2*36**2)/36
            i4=mod(i,36)
            decoded=cc(i1+1:i1+1)//cc(i2+1:i2+1)//cc(i3+1:i3+1)//
     +           cc(i4+1:i4+1)//'                  '
         endif
      endif
      decoded(22:22)=c

C  Neutralize the SM2CEW foolishness
      call packmsg(decoded,dgen)
      call unpackmsg(dgen,decoded)

      return
      end
