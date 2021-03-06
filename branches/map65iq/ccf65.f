      subroutine ccf65(ss,nhsym,sync1,dt1,flipk,mode65,syncshort,
     +     snr2,dt2)

C  Check one frequency bin for JT65 sync pattern or shorthand pattern.

      parameter (NFFT=512,NH=NFFT/2)
      parameter (LAGMAX=37)
      real ss(322)
                   !Input: half-symbol powers
      real s(NFFT)                     !CCF = ss*pr
      complex cs(0:NH)                 !Complex FT of s
      real s2(NFFT)                    !CCF = ss*pr2
      complex cs2(0:NH)                !Complex FT of s2
      real pr(NFFT)                    !JT65 pseudo-random sync pattern
      complex cpr(0:NH)                !Complex FT of pr
      real pr2(NFFT)                   !JT65 shorthand pattern
      complex cpr2(0:NH)               !Complex FT of pr2
      real tmp1(322)
      real tmp2(322)
      real ccf(-LAGMAX:LAGMAX)
      logical first
      integer npr(126)
      data first/.true./
      equivalence (s,cs),(pr,cpr),(s2,cs2),(pr2,cpr2)
      save

C  The JT65 pseudo-random sync pattern:
      data npr/
     + 1,0,0,1,1,0,0,0,1,1,1,1,1,1,0,1,0,1,0,0,
     + 0,1,0,1,1,0,0,1,0,0,0,1,1,1,0,0,1,1,1,1,
     + 0,1,1,0,1,1,1,1,0,0,0,1,1,0,1,0,1,0,1,1,
     + 0,0,1,1,0,1,0,1,0,1,0,0,1,0,0,0,0,0,0,1,
     + 1,0,0,0,0,0,0,0,1,1,0,1,0,0,1,0,1,1,0,1,
     + 0,1,0,1,0,0,1,1,0,0,1,0,0,1,0,0,0,0,1,1,
     + 1,1,1,1,1,1/

      if(first) then
C  Initialize pr, pr2; compute cpr, cpr2.
         fac=1.0/NFFT
         do i=1,NFFT
            pr(i)=0.
            k=2*mod((i-1)/8,2)-1
            pr2(i)=fac*k
         enddo
         do i=1,126
            j=2*i
            pr(j)=fac*(2*npr(i)-1)
         enddo
         call four2a(pr,NFFT,1,-1,0)
         call four2a(pr2,NFFT,1,-1,0)
         first=.false.
      endif

C  Look for JT65 sync pattern and shorthand square-wave pattern.
      ccfbest=0.
      ccfbest2=0.

      do i=1,nhsym              ! ?? nhsym-1 ??
!         s(i)=min(4.0,ss(i)+ss(i+1))
         s(i)=ss(i)+ss(i+1)
      enddo
      do i=nhsym+1,NFFT         ! ?? nhsym ??
         s(i)=0.
      enddo
      call four2a(s,NFFT,1,-1,0) !Real-to-complex FFT
      do i=0,NH
         cs2(i)=cs(i)*conjg(cpr2(i)) !Mult by complex FFT of pr2
         cs(i)=cs(i)*conjg(cpr(i)) !Mult by complex FFT of pr
      enddo
      call four2a(cs,NFFT,1,1,-1) !Complex-to-real inv-FFT
      call four2a(cs2,NFFT,1,1,-1) !Complex-to-real inv-FFT

      do lag=-LAGMAX,LAGMAX             !Check for best JT65 sync
         i=lag+28
         if(i.lt.1) i=i+NFFT
         ccf(lag)=s(i)
         if(abs(ccf(lag)).gt.ccfbest) then
            ccfbest=abs(ccf(lag))
            lagpk=lag
            flipk=1.0
            if(ccf(lag).lt.0.0) flipk=-1.0
         endif
      enddo

      do lag=-8,7               !Check for best shorthand
         ccf2=s2(lag+28)
         if(ccf2.gt.ccfbest2) then
            ccfbest2=ccf2
            lagpk2=lag
         endif
      enddo


C  Find rms level on baseline of "ccfblue", for normalization.
      sum=0.
      do lag=-LAGMAX+1,LAGMAX-1
         if(abs(lag-lagpk).gt.1) sum=sum + ccf(lag)
      enddo
      base=sum/50.0
      sq=0.
      do lag=-LAGMAX+1,LAGMAX-1
         if(abs(lag-lagpk).gt.1) sq=sq + (ccf(lag)-base)**2
      enddo
      rms=sqrt(sq/49.0)
      sync1=ccfbest/rms - 4.0
      dt1=3.5 + lagpk*(2048.0/11025.0)

C  Find base level for normalizing snr2.
      do i=1,nhsym
         tmp1(i)=ss(i)
      enddo
      call pctile(tmp1,tmp2,nhsym,40,base)
      snr2=0.398107*ccfbest2/base                !### empirical
      syncshort=0.5*ccfbest2/rms - 4.0           !### better normalizer than rms?
      dt2=2.5 + lagpk2*(2048.0/11025.0)

      return
      end
