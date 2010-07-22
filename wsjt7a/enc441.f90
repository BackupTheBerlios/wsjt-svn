subroutine enc441(msg0,msg2,len2)

! Encode an FSK441++ message

  character*28 msg,msg0,msg2
  character*4 tok(12)
  character*12 ctok
  integer ntok(12)
  integer len(12)

  data tok(1) /'CQ ' /,ntok(1)/3/              !001 1   Tokens and token lengths
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
  data ctok/'168_$GINRTOZ'/                    !Token indexes as characters

  msg=msg0
  do i=28,1,-1
     if(msg(i:i).ne.' ') go to 10
  enddo
10 len0=min(i,24)                        !User's message length
  i1=0
  iz=0

  i2=max(len0+1,6)
  if(msg(1:3).eq.'26 ')  msg=msg0(4:i2)//'26'
  if(msg(1:3).eq.'27 ')  msg=msg0(4:i2)//'27'
  if(msg(1:3).eq.'38 ')  msg=msg0(4:i2)//'38'
  if(msg(1:4).eq.'R26 ') msg=msg0(5:i2)//'R26'
  if(msg(1:4).eq.'R27 ') msg=msg0(5:i2)//'R27'
  if(msg(1:4).eq.'R38 ') msg=msg0(5:i2)//'R38'
  if(msg(1:4).eq.'RRR ') msg=msg0(5:i2)//'RRR'
  if(msg(1:3).eq.'73 ')  msg=msg0(4:i2)//'73'

  do i=1,11
     i1=index(msg,tok(i)(1:ntok(i)))   !i1 marks start of token in msg
     if(i1.gt.0) go to 20
  enddo

20 iz=i                                !iz is token number (0 if none)

  if(i1.le.0) then                     !No token found
     jz=12
     do j=12,1,-1
        if(len(j).ge.len0+4) jz=j        !jz is index of msg length
     enddo
     len2=len(jz)
     msg2='$!'//ctok(jz:jz)//ctok(12:12)//msg
     go to 900
  endif

  jz=12
  do j=12,1,-1
     if(len(j).ge.len0+4-ntok(iz)) jz=j        !jz is index of msg length
  enddo
  len2=len(jz)

  if(iz.le.3 .or. len0.le.3) then
     msg2='$!'//ctok(jz:jz)//ctok(iz:iz)//msg(ntok(iz)+1:)
  else
     msg2='$!'//ctok(jz:jz)//ctok(iz:iz)//msg(1:i1-1)
  endif

900 return
end subroutine enc441

