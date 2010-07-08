program t1

  character*28 msg
  integer itone(645)

  call getarg(1,msg)
  call genscat6(msg,itone)

  write(*,3001) itone
3001 format(21i3)



end program t1
