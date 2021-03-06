      subroutine deep65(s3,mode65,neme,nsked,flip,mycall,hiscall,
     +  hisgrid,decoded,qual)

      parameter (MAXCALLS=7000,MAXRPT=63)
      real s3(64,63)
      character callsign*12,grid*4,message*22,hisgrid*6,c*1,ceme*3
      character*12 mycall,hiscall
      character*22 decoded,deepmsg
      character*22 testmsg(2*MAXCALLS + 2 + MAXRPT)
      character*15 callgrid(MAXCALLS)
      character*80 line
      character*4 rpt(MAXRPT)
      logical first
      integer ncode(63,2*MAXCALLS)
      common/tmp8/ p(64,63)

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
 1002       format(a80)
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
            if(neme.eq.1 .and. ceme.ne.'EME') go to 10
         endif

 5       icall=icall+1
         j1=index(mycall,' ') - 1
         if(j1.lt.3) j1=6
         j2=index(callsign,' ') - 1
         if(j2.lt.3) j2=6
         j3=index(mycall,'/')
         j4=index(callsign,'/')
         callgrid(icall)=callsign(1:j2)

         mz=1
         if(n.eq.1) mz=MAXRPT+1
         do m=1,mz
            if(m.gt.1) grid=rpt(m-1)
            if(j3.lt.1 .and.j4.lt.1) 
     +         callgrid(icall)=callsign(1:j2)//' '//grid
            message=mycall(1:j1)//' '//callgrid(icall)
            k=k+1
            testmsg(k)=message
            call encode65(message,ncode(1,k))
            if(m.eq.1) then
               message='CQ '//callgrid(icall)
               k=k+1
               testmsg(k)=message
               call encode65(message,ncode(1,k))
            endif
         enddo
         if(nsked.eq.1) go to 20
 10   enddo
 20   ntot=k
      neme0=neme

      sum0=0.
      do j=1,63
         smax=-1.e30
         do i=1,64
            smax=max(smax,s3(i,j))
         enddo
         sum0=sum0+smax
      enddo

      p1=-1.e30
      ip1=0
      p2=-1.e30
      ip2=0
      do k=1,ntot
C  If sync=OOO, no CQ messages
         if(flip.lt.0.0 .and. testmsg(k)(1:3).eq.'CQ ') go to 30
         sum=0.
         sum2=0.
         do j=1,63
            i=ncode(j,k)+1
            sum=sum + s3(i,j)
            sum2=sum2 + p(i,j)
         enddo
         if(sum.gt.p1) then
            p1=sum
            ip1=k
         endif
         if(sum2.gt.p2) then
            p2=sum2
            ip2=k
         endif
 30   enddo

      p1=p1/sum0
      qual=100.0*(p1-0.40)
      if(mode65.eq.1) qual=100.0*(p1-0.33)
      if(mode65.eq.4) qual=100.0*(p1-0.50)
      if(qual.lt.0.) qual=0.
      if(qual.gt.10.) qual=10.
      decoded='                      '
      c=' '
      if(qual.gt.0.0) then
         if(qual.lt.4.0) c='?'
         decoded=testmsg(ip1)
      endif
      decoded(22:22)=c
      deepmsg=decoded

      q2=0.27*p2 + 81.3
!      if(mode65.eq.1) qual=100.0*(p1-0.33)
!      if(mode65.eq.4) qual=100.0*(p1-0.50)
      if(q2.lt.0.) q2=0.
      if(q2.gt.10.) q2=10.
      decoded='                      '
      c=' '
      if(q2.gt.0.0) then
         if(q2.lt.4.0) c='?'
         decoded=testmsg(ip2)
      endif
      decoded(22:22)=c

!      qual=q2
      decoded=deepmsg

      return
      end
