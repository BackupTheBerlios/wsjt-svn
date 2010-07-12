program t74

! Tests experimental JT41 decoder

  parameter (NMAX=512*1024)
  real dat(NMAX)                          !Raw signal, 30 s at 11025 sps
  character arg*12                        !Command-line argument
  character cfile6*6                      !File time

  nargs=iargc()
  if(nargs.ne.1) then
     print*,'Usage: t74 nfile'
     go to 999
  endif
  call getarg(1,arg)
  read(arg,*) nfile
  open(74,file='dat.74',form='unformatted',status='old')
  
  do ifile=1,nfile
     read(74,end=999) jz,cfile6,(dat(j),j=1,jz)
     if(ifile.ne.nfile .and. nfile.ne.999) go to 900

     call jt41(dat,jz,cfile6)

900  continue
  enddo

999 end program t74
