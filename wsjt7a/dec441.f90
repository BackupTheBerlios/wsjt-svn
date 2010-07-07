subroutine dec441(msg2,msg)

  character*28 msg2,msg
  character*4 tok(12)
  character*12 ctok
  integer ntok(12)
  integer len(12)

! Define tokens and token lengths
  data tok(1) /'CQ ' /,ntok(1)/3/                !001 1
  data tok(2) /'QRZ '/,ntok(2)/4/                !012 6
  data tok(3) /'TNX '/,ntok(3)/4/                !020 8
  data tok(4) /' 26' /,ntok(4)/3/                !033 _
  data tok(5) /' 27' /,ntok(5)/3/                !100 $
  data tok(6) /' 38' /,ntok(6)/3/                !113 G
  data tok(7) /' R26'/,ntok(7)/4/                !121 I
  data tok(8) /' R27'/,ntok(8)/4/                !132 N
  data tok(9) /' R38'/,ntok(9)/4/                !202 R
  data tok(10)/' RRR'/,ntok(10)/4/               !210 T
  data tok(11)/' 73' /,ntok(11)/3/               !223 O
  data tok(12)/'   ' /,ntok(12)/3/               !231 Z

! Permissible message lengths
  data len/4,7,9,11,13,14,15,17,19,21,23,28/

! Token indexes, as characters
  data ctok/'168_$GINRTOZ'/

  do i=28,1,-1
     if(msg2(i:i).ne.' ') go to 10
  enddo
10 len2=i                                !Encoded message length

  if(msg2(1:2).eq.'$!') then
     do i=1,12
        if(ctok(i:i).eq.msg2(3:3)) go to 12
     enddo
12   if(i.le.12) len2=len(i)
  endif

  len0=0
  do j=1,12
     if(msg2(3:3).eq.ctok(j:j)) len0=len(j)
  enddo

  iz=0
  do i=1,12
     if(msg2(4:4).eq.ctok(i:i)) iz=i
  enddo

  if(iz.eq.0) then
     msg=msg2
  else if(iz.le.3) then
     msg=tok(iz)(1:ntok(iz))//msg2(5:)
  else
     if(len2.eq.4) then
        msg=tok(iz)(2:)
     else
        do i=28,1,-1
           if(msg2(i:i).ne.' ') go to 20
        enddo
20      i2=i
        msg=msg2(5:i2)//tok(iz)
     endif
  endif

  return
end subroutine dec441

