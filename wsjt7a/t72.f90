program t72

! Tests experimental FSK441 decoder

  parameter (NMAX=512*1024)
  parameter (MAXFFT=8192)
  real dat0(NMAX)                          !Raw signal, 30 s at 11025 sps
  real dat(NMAX)                          !Raw signal, 30 s at 11025 sps
  character arg*12                        !Command-line argument
  character cfile6*6                      !File time
  character frag*28                       !Message fragment to be matched
  integer dftolerance
  logical pick
  real pingdat(3,100)                     !Detected pings
  character c*48
  common/scratch/work(NMAX)
  data c/' 123456789.,?/# $ABCD FGHIJKLMNOPQRSTUVWXY 0EZ  '/

  nargs=iargc()
  if(nargs.ne.2) then
     print*,'Usage: t72 nfile frag'
     go to 999
  endif
  call getarg(1,arg)
  read(arg,*) nfile
  call getarg(2,frag)
  open(72,file='dat.72',form='unformatted',status='old')

! Initialize variables
  minsigdb=2
  minwidth=40
  dftolerance=400
  pick=.false.
  nsps=25
  nsam=nsps*ndits
  xn=log(float(nsam))/log(2.0)
  n=xn
  if(xn-n .gt.0.001) n=n+1

  do ifile=1,nfile
     read(72,end=999) jz,nz,cfile6,(dat0(j),j=1,jz)
     if(ifile.ne.nfile .and. nfile.ne.999) go to 900

! If necessary, correct sample-rate errors
     if(ifile.eq.3 .or. ifile.eq.4 .or. ifile.eq.5) then
        j=0
        do i=1,jz
           j=j+1
           if(mod(i,147).eq.0) then
              j=j-1
           endif
           dat(i)=dat0(j)
        enddo
        dat(j:jz)=0.
     else
        dat(1:jz)=dat0(1:jz)
     endif

     call ping441(dat,jz,nz,MinSigdB,MinWidth,pick,pingdat,nping)   !Find pings

     do iping=1,nping                        !Process each ping
        tstart=pingdat(1,iping)
        width=pingdat(2,iping)
        call pp441(dat,jz,tstart,width,dftolerance,frag)
     enddo

900  continue
  enddo

999 end program t72
