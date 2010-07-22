subroutine match(s1,s2,nstart,nmatch)

  character*(*) s1,s2

  nstart=-1
  nmatch=0
  n1=len_trim(s1)+1
  n2=len(s2)
  if(n2.ge.n1) then
     do j=1,n2
        n=0
        do i=1,n1
           k=j+i-1
           if(k.gt.n2) k=k-n2
           if(s2(k:k).eq.s1(i:i)) n=n+1
        enddo
        if(n.gt.nmatch) then
           nmatch=n
           nstart=j
        endif
     enddo
  endif

  return
end subroutine match
        
