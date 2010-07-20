subroutine dec441(msg1,n2,msg2)

! Decode the "++" enhancements from an FSK441 message

  character*28 msg1,msg2
  character*4 tok(12)
  character*12 ctok
  integer ntok(12)
  integer len(12)

  data tok(1) /'CQ ' /,ntok(1)/3/              !001 1  Tokens and token lengths
  data tok(2) /'QRZ '/,ntok(2)/4/              !012 6
  data tok(3) /'TNX '/,ntok(3)/4/              !020 8
  data tok(4) /' 26' /,ntok(4)/3/              !033 _
  data tok(5) /' 27' /,ntok(5)/3/              !100 $
  data tok(6) /' 38' /,ntok(6)/3/              !113 G
  data tok(7) /' R26'/,ntok(7)/4/              !121 I
  data tok(8) /' R27'/,ntok(8)/4/              !132 N
  data tok(9) /' R38'/,ntok(9)/4/              !202 R
  data tok(10)/' RRR'/,ntok(10)/4/             !210 T
  data tok(11)/' 73' /,ntok(11)/3/             !223 O
  data tok(12)/'   ' /,ntok(12)/3/             !231 Z
  data len/4,7,9,11,13,14,15,17,19,21,23,28/   !Permissible message lengths
  data ctok/'168_$GINRTOZ'/                    !Token indexes, as characters

  len0=0
  do j=1,12
     if(msg1(3:3).eq.ctok(j:j)) len0=len(j)
  enddo

  iz=0
  do i=1,12
     if(msg1(4:4).eq.ctok(i:i)) iz=i
  enddo
  if(msg1(4:4).eq.' ') iz=4

  if(iz.eq.0) then
     msg2=msg1
     if(msg1(1:2).eq.'$!') msg2=msg1(5:)
  else if(iz.le.3) then
     msg2=tok(iz)(1:ntok(iz))//msg1(5:)
  else
     if(n2.eq.4) then
        msg2=tok(iz)(2:)
     else
        msg2=msg1(5:n2)//tok(iz)
     endif
  endif

  return
end subroutine dec441

