subroutine token(c1,i1,tok1,n1)

  character c1*1,tok1*4
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

  i1=0
  do i=1,12
     if(ctok(i:i).eq.c1) go to 10
  enddo
  tok1='    '
  n1=0
  go to 900

10 i1=i
  tok1=tok(i1)
  n1=len(i1)

900 return
end subroutine token

