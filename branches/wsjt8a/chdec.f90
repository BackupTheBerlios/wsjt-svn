subroutine chdec(cmode,nbit,gsym,iu)

! Decode channel symbols to recover source-encoded user message

  character*5 cmode
  integer gsym(372)
  integer gsym2(372)
  integer iu(3)
  integer era(63)
  integer dat4(13)
  integer*1 dbits(96)
  integer*1 symbols(372)
  integer*1 ddec(10)
  integer mettab(0:255,0:1)
  logical first
  integer igray0(0:7)
  data igray0/0,1,3,2,7,6,4,5/    !Use this to remove the gray code
  data first/.true./
  save first,mettab

  if(cmode.eq.'JTMS' .or. cmode.eq.'JT8') then

     if(first) then
! Get the metric table
        bias=0.0                        !Metric bias: viterbi=0, seq=rate
        scale=10                        !Optimize?
        if(cmode.eq.'JTMS') open(19,file='met2.21',status='old')
        if(cmode.eq.'JT8')  open(19,file='met8.21',status='old')

        do i=0,255
           read(19,*) xjunk,d0,d1
           mettab(i,0)=nint(scale*(d0-bias))
           mettab(i,1)=nint(scale*(d1-bias))    !### Check range, etc.  ###
        enddo
        first=.false.
     endif

     if(cmode.eq.'JTMS') then
        nhdata=nbit+12
        nsym=2*nhdata                          !Number of binary symbols
        gsym2(1:2*nhdata)=gsym(1:2*nhdata)
        do i=1,nhdata                          !Remove the interleaving
           gsym(2*i-1)=gsym2(i)
           gsym(2*i)=gsym2(nhdata+i)
        enddo        
     else if(cmode.eq.'JT8') then
        nsym=4*(nbit+15)                       !Number of binary symbols
! Remove gray code and convert 3-bit symbols to bit-level soft symbols
        gsym2=gsym
        do i=1,124
           n=igray0(gsym2(i))
           gsym(3*i-2)=n/4
           gsym(3*i-1)=iand(n,2)/2
           gsym(3*i)=iand(n,1)
        enddo
        do i1=0,30                    !Remove interleaving
           do i2=0,11
              i=31*i2+i1
              j=12*i1+i2
              gsym2(j+1)=gsym(i+1)
           enddo
        enddo
        gsym=gsym2
     endif

     do i=1,nsym
        n=127
        if(gsym(i).eq.1) n=-127
        symbols(i)=n
     enddo

     if(cmode.eq.'JTMS') call vit213(symbols,nbit,mettab,ddec,metric)
     if(cmode.eq.'JT8')  call vit416(symbols,nbit,mettab,ddec,metric)
     iz=(nbit+7)/8
     ddec(iz+1:)=0
     n1=ddec(1)
     n2=ddec(2)
     n3=ddec(3)
     n4=ddec(4)
     iu(1)=ishft(iand(n1,255),24) + ishft(iand(n2,255),16) +         &
           ishft(iand(n3,255),8) + iand(n4,255)
     n1=ddec(5)
     n2=ddec(6)
     n3=ddec(7)
     n4=ddec(8)
     iu(2)=ishft(iand(n1,255),24) + ishft(iand(n2,255),16) +         &
           ishft(iand(n3,255),8) + iand(n4,255)
     n1=ddec(9)
     n2=ddec(10)
     iu(3)=ishft(iand(n1,255),24) + ishft(iand(n2,255),16)
  else if(cmode.eq.'JT64' .or. cmode.eq.'ISCAT') then
     nerase=0
     kk=5
     if(nbit.eq.48) kk=8
     if(nbit.eq.78) kk=13
     call krsdecode(gsym,kk,era,nerase,dat4,ncount)
     dbits=0
     call unpackbits(dat4,13,6,dbits)
     call packbits(dbits,3,32,iu)
  endif

  return
end subroutine chdec
