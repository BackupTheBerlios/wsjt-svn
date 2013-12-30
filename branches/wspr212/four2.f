      SUBROUTINE FOUR2a (DATA,N,NDIM,ISIGN,IFORM)

C     Cooley-Tukey fast Fourier transform in USASI basic Fortran.
C     multi-dimensional transform, each dimension a power of two,
C     complex or real data.

C     TRANSFORM(K1,K2,...) = SUM(DATA(J1,J2,...)*EXP(ISIGN*2*PI*SQRT(-1)
C     *((J1-1)*(K1-1)/N(1)+(J2-1)*(K2-1)/N(2)+...))), summed for all
C     J1 and K1 from 1 to N(1), J2 and K2 from 1 TO N(2),
C     etc, for all NDIM subscripts.  NDIM must be positive and
C     each N(IDIM) must be a power of two.  ISIGN is +1 or -1.
C     Let NTOT = N(1)*N(2)*...*N(NDIM).  Then a -1 transform
C     followed by a +1 one (or vice versa) returns NTOT
C     times the original data.  

C     IFORM = 1, 0 or -1, as data is
C     complex, real, or the first half of a complex array.  Transform
C     values are returned in array DATA.  They are complex, real, or
C     the first half of a complex array, as IFORM = 1, -1 or 0.

C     The transform of a real array (IFORM = 0) dimensioned N(1) by N(2)
C     by ... will be returned in the same array, now considered to
C     be complex of dimensions N(1)/2+1 by N(2) by ....  Note that if
C     IFORM = 0 or -1, N(1) must be even, and enough room must be
C     reserved.  The missing values may be obtained by complex conjuga-
C     tion.  

C     The reverse transformation of a half complex array dimensioned
C     N(1)/2+1 by N(2) by ..., is accomplished by setting IFORM
C     to -1.  In the N array, N(1) must be the true N(1), not N(1)/2+1.
C     The transform will be real and returned to the input array.

C     Running time is proportional to NTOT*LOG2(NTOT), rather than
C     the naive NTOT**2.  Furthermore, less error is built up.

C     Written by Norman Brenner of MIT Lincoln Laboratory, January 1969.
C     See IEEE Audio Transactions (June 1967), Special issue on FFT.

      parameter(NMAX=2048*1024)
      DIMENSION DATA(NMAX), N(1)
      NTOT=1
      DO 10 IDIM=1,NDIM
 10   NTOT=NTOT*N(IDIM)
      IF (IFORM) 70,20,20
 20   NREM=NTOT
      DO 60 IDIM=1,NDIM
      NREM=NREM/N(IDIM)
      NPREV=NTOT/(N(IDIM)*NREM)
      NCURR=N(IDIM)
      IF (IDIM-1+IFORM) 30,30,40
 30   NCURR=NCURR/2
 40   CALL BITRV (DATA,NPREV,NCURR,NREM)
      CALL COOL2 (DATA,NPREV,NCURR,NREM,ISIGN)
      IF (IDIM-1+IFORM) 50,50,60
 50   CALL FIXRL (DATA,N(1),NREM,ISIGN,IFORM)
      NTOT=(NTOT/N(1))*(N(1)/2+1)
 60   CONTINUE
      RETURN
 70   NTOT=(NTOT/N(1))*(N(1)/2+1)
      NREM=1
      DO 100 JDIM=1,NDIM
      IDIM=NDIM+1-JDIM
      NCURR=N(IDIM)
      IF (IDIM-1) 80,80,90
 80   NCURR=NCURR/2
      CALL FIXRL (DATA,N(1),NREM,ISIGN,IFORM)
      NTOT=NTOT/(N(1)/2+1)*N(1)
 90   NPREV=NTOT/(N(IDIM)*NREM)
      CALL BITRV (DATA,NPREV,NCURR,NREM)
      CALL COOL2 (DATA,NPREV,NCURR,NREM,ISIGN)
 100  NREM=NREM*N(IDIM)
      RETURN
      END
      SUBROUTINE BITRV (DATA,NPREV,N,NREM)
C     SHUFFLE THE DATA BY BIT REVERSAL.
C     DIMENSION DATA(NPREV,N,NREM)
C     COMPLEX DATA
C     EXCHANGE DATA(J1,J4REV,J5) WITH DATA(J1,J4,J5) FOR ALL J1 FROM 1
C     TO NPREV, ALL J4 FROM 1 TO N (WHICH MUST BE A POWER OF TWO), AND
C     ALL J5 FROM 1 TO NREM.  J4REV-1 IS THE BIT REVERSAL OF J4-1.  E.G.
C     SUPPOSE N = 32.  THEN FOR J4-1 = 10011, J4REV-1 = 11001, ETC.
      parameter(NMAX=2048*1024)
      DIMENSION DATA(NMAX)
      IP0=2
      IP1=IP0*NPREV
      IP4=IP1*N
      IP5=IP4*NREM
      I4REV=1
C     I4REV = 1+(J4REV-1)*IP1
      DO 60 I4=1,IP4,IP1
C     I4 = 1+(J4-1)*IP1
      IF (I4-I4REV) 10,30,30
 10   I1MAX=I4+IP1-IP0
      DO 20 I1=I4,I1MAX,IP0
C     I1 = 1+(J1-1)*IP0+(J4-1)*IP1
      DO 20 I5=I1,IP5,IP4
C     I5 = 1+(J1-1)*IP0+(J4-1)*IP1+(J5-1)*IP4
      I5REV=I4REV+I5-I4
C     I5REV = 1+(J1-1)*IP0+(J4REV-1)*IP1+(J5-1)*IP4
      TEMPR=DATA(I5)
      TEMPI=DATA(I5+1)
      DATA(I5)=DATA(I5REV)
      DATA(I5+1)=DATA(I5REV+1)
      DATA(I5REV)=TEMPR
 20   DATA(I5REV+1)=TEMPI
C     ADD ONE WITH DOWNWARD CARRY TO THE HIGH ORDER BIT OF J4REV-1.
 30   IP2=IP4/2
 40   IF (I4REV-IP2) 60,60,50
 50   I4REV=I4REV-IP2
      IP2=IP2/2
      IF (IP2-IP1) 60,40,40
 60   I4REV=I4REV+IP2
      RETURN
      END
      SUBROUTINE COOL2 (DATA,NPREV,N,NREM,ISIGN)
C     DISCRETE FOURIER TRANSFORM OF LENGTH N.  IN-PLACE COOLEY-TUKEY
C     ALGORITHM, BIT-REVERSED TO NORMAL ORDER, SANDE-TUKEY PHASE SHIFTS.
C     DIMENSION DATA(NPREV,N,NREM)
C     COMPLEX DATA
C     DATA(J1,K4,J5) = SUM(DATA(J1,J4,J5)*EXP(ISIGN*2*PI*I*(J4-1)*
C     (K4-1)/N)), SUMMED OVER J4 = 1 TO N FOR ALL J1 FROM 1 TO NPREV,
C     K4 FROM 1 TO N AND J5 FROM 1 TO NREM.  N MUST BE A POWER OF TWO.
C     METHOD--LET IPREV TAKE THE VALUES 1, 2 OR 4, 4 OR 8, ..., N/16,
C     N/4, N.  THE CHOICE BETWEEN 2 OR 4, ETC., DEPENDS ON WHETHER N IS
C     A POWER OF FOUR.  DEFINE IFACT = 2 OR 4, THE NEXT FACTOR THAT
C     IPREV MUST TAKE, AND IREM = N/(IFACT*IPREV).  THEN--
C     DIMENSION DATA(NPREV,IPREV,IFACT,IREM,NREM)
C     COMPLEX DATA
C     DATA(J1,J2,K3,J4,J5) = SUM(DATA(J1,J2,J3,J4,J5)*EXP(ISIGN*2*PI*I*
C     (K3-1)*((J3-1)/IFACT+(J2-1)/(IFACT*IPREV)))), SUMMED OVER J3 = 1
C     TO IFACT FOR ALL J1 FROM 1 TO NPREV, J2 FROM 1 TO IPREV, K3 FROM
C     1 TO IFACT, J4 FROM 1 TO IREM AND J5 FROM 1 TO NREM.  THIS IS
C     A PHASE-SHIFTED DISCRETE FOURIER TRANSFORM OF LENGTH IFACT.
C     FACTORING N BY FOURS SAVES ABOUT TWENTY FIVE PERCENT OVER FACTOR-
C     ING BY TWOS.  DATA MUST BE BIT-REVERSED INITIALLY.
C     IT IS NOT NECESSARY TO REWRITE THIS SUBROUTINE INTO COMPLEX
C     NOTATION SO LONG AS THE FORTRAN COMPILER USED STORES REAL AND
C     IMAGINARY PARTS IN ADJACENT STORAGE LOCATIONS.  IT MUST ALSO
C     STORE ARRAYS WITH THE FIRST SUBSCRIPT INCREASING FASTEST.
      parameter(NMAX=2048*1024)
      DIMENSION DATA(NMAX)

      real*8 twopi,wstpr,wstpi,wr,wi,w2r,w2i,w3r,w3i,wtempr

      TWOPI=6.2831853072*FLOAT(ISIGN)
      IP0=2
      IP1=IP0*NPREV
      IP4=IP1*N
      IP5=IP4*NREM
      IP2=IP1
C     IP2=IP1*IPROD
      NPART=N
 10   IF (NPART-2) 60,30,20
 20   NPART=NPART/4
      GO TO 10
C     DO A FOURIER TRANSFORM OF LENGTH TWO
 30   IF (IP2-IP4) 40,160,160
 40   IP3=IP2*2
C     IP3=IP2*IFACT
      DO 50 I1=1,IP1,IP0
C     I1 = 1+(J1-1)*IP0
      DO 50 I5=I1,IP5,IP3
C     I5 = 1+(J1-1)*IP0+(J4-1)*IP3+(J5-1)*IP4
      I3A=I5
      I3B=I3A+IP2
C     I3 = 1+(J1-1)*IP0+(J2-1)*IP1+(J3-1)*IP2+(J4-1)*IP3+(J5-1)*IP4
      TEMPR=DATA(I3B)
      TEMPI=DATA(I3B+1)
      DATA(I3B)=DATA(I3A)-TEMPR
      DATA(I3B+1)=DATA(I3A+1)-TEMPI
      DATA(I3A)=DATA(I3A)+TEMPR
 50   DATA(I3A+1)=DATA(I3A+1)+TEMPI
      IP2=IP3
C     DO A FOURIER TRANSFORM OF LENGTH FOUR (FROM BIT REVERSED ORDER)
 60   IF (IP2-IP4) 70,160,160
 70   IP3=IP2*4
C     IP3=IP2*IFACT
C     COMPUTE TWOPI THRU WR AND WI IN DOUBLE PRECISION, IF AVAILABLE.
      THETA=TWOPI/FLOAT(IP3/IP1)
      SINTH=SIN(THETA/2)
      WSTPR=-2*SINTH*SINTH
      WSTPI=SIN(THETA)
      WR=1.
      WI=0.
      DO 150 I2=1,IP2,IP1
C     I2 = 1+(J2-1)*IP1
      IF (I2-1) 90,90,80
 80   W2R=WR*WR-WI*WI
      W2I=2*WR*WI
      W3R=W2R*WR-W2I*WI
      W3I=W2R*WI+W2I*WR
 90   I1MAX=I2+IP1-IP0
      DO 140 I1=I2,I1MAX,IP0
C     I1 = 1+(J1-1)*IP0+(J2-1)*IP1
      DO 140 I5=I1,IP5,IP3
C     I5 = 1+(J1-1)*IP0+(J2-1)*IP1+(J4-1)*IP3+(J5-1)*IP4
      I3A=I5
      I3B=I3A+IP2
      I3C=I3B+IP2
      I3D=I3C+IP2
C     I3 = 1+(J1-1)*IP0+(J2-1)*IP1+(J3-1)*IP2+(J4-1)*IP3+(J5-1)*IP4
      IF (I2-1) 110,110,100
C     APPLY THE PHASE SHIFT FACTORS
 100  TEMPR=DATA(I3B)
      DATA(I3B)=W2R*DATA(I3B)-W2I*DATA(I3B+1)
      DATA(I3B+1)=W2R*DATA(I3B+1)+W2I*TEMPR
      TEMPR=DATA(I3C)
      DATA(I3C)=WR*DATA(I3C)-WI*DATA(I3C+1)
      DATA(I3C+1)=WR*DATA(I3C+1)+WI*TEMPR
      TEMPR=DATA(I3D)
      DATA(I3D)=W3R*DATA(I3D)-W3I*DATA(I3D+1)
      DATA(I3D+1)=W3R*DATA(I3D+1)+W3I*TEMPR
 110  T0R=DATA(I3A)+DATA(I3B)
      T0I=DATA(I3A+1)+DATA(I3B+1)
      T1R=DATA(I3A)-DATA(I3B)
      T1I=DATA(I3A+1)-DATA(I3B+1)
      T2R=DATA(I3C)+DATA(I3D)
      T2I=DATA(I3C+1)+DATA(I3D+1)
      T3R=DATA(I3C)-DATA(I3D)
      T3I=DATA(I3C+1)-DATA(I3D+1)
      DATA(I3A)=T0R+T2R
      DATA(I3A+1)=T0I+T2I
      DATA(I3C)=T0R-T2R
      DATA(I3C+1)=T0I-T2I
      IF (ISIGN) 120,120,130
 120  T3R=-T3R
      T3I=-T3I
 130  DATA(I3B)=T1R-T3I
      DATA(I3B+1)=T1I+T3R
      DATA(I3D)=T1R+T3I
 140  DATA(I3D+1)=T1I-T3R
      WTEMPR=WR
      WR=WSTPR*WTEMPR-WSTPI*WI+WTEMPR
 150  WI=WSTPR*WI+WSTPI*WTEMPR+WI
      IP2=IP3
      GO TO 60
 160  RETURN
      END
      SUBROUTINE FIXRL (DATA,N,NREM,ISIGN,IFORM)
C     FOR IFORM = 0, CONVERT THE TRANSFORM OF A DOUBLED-UP REAL ARRAY,
C     CONSIDERED COMPLEX, INTO ITS TRUE TRANSFORM.  SUPPLY ONLY THE
C     FIRST HALF OF THE COMPLEX TRANSFORM, AS THE SECOND HALF HAS
C     CONJUGATE SYMMETRY.  FOR IFORM = -1, CONVERT THE FIRST HALF
C     OF THE TRUE TRANSFORM INTO THE TRANSFORM OF A DOUBLED-UP REAL
C     ARRAY.  N MUST BE EVEN.
C     USING COMPLEX NOTATION AND SUBSCRIPTS STARTING AT ZERO, THE
C     TRANSFORMATION IS--
C     DIMENSION DATA(N,NREM)
C     ZSTP = EXP(ISIGN*2*PI*I/N)
C     DO 10 I2=0,NREM-1
C     DATA(0,I2) = CONJ(DATA(0,I2))*(1+I)
C     DO 10 I1=1,N/4
C     Z = (1+(2*IFORM+1)*I*ZSTP**I1)/2
C     I1CNJ = N/2-I1
C     DIF = DATA(I1,I2)-CONJ(DATA(I1CNJ,I2))
C     TEMP = Z*DIF
C     DATA(I1,I2) = (DATA(I1,I2)-TEMP)*(1-IFORM)
C 10  DATA(I1CNJ,I2) = (DATA(I1CNJ,I2)+CONJ(TEMP))*(1-IFORM)
C     IF I1=I1CNJ, THE CALCULATION FOR THAT VALUE COLLAPSES INTO
C     A SIMPLE CONJUGATION OF DATA(I1,I2).
      parameter(NMAX=4048*1024)
      DIMENSION DATA(NMAX)
      TWOPI=6.283185307*FLOAT(ISIGN)
      IP0=2
      IP1=IP0*(N/2)
      IP2=IP1*NREM
      IF (IFORM) 10,70,70
C     PACK THE REAL INPUT VALUES (TWO PER COLUMN)
 10   J1=IP1+1
      DATA(2)=DATA(J1)
      IF (NREM-1) 70,70,20
 20   J1=J1+IP0
      I2MIN=IP1+1
      DO 60 I2=I2MIN,IP2,IP1
      DATA(I2)=DATA(J1)
      J1=J1+IP0
      IF (N-2) 50,50,30
 30   I1MIN=I2+IP0
      I1MAX=I2+IP1-IP0
      DO 40 I1=I1MIN,I1MAX,IP0
      DATA(I1)=DATA(J1)
      DATA(I1+1)=DATA(J1+1)
 40   J1=J1+IP0
 50   DATA(I2+1)=DATA(J1)
 60   J1=J1+IP0
 70   DO 80 I2=1,IP2,IP1
      TEMPR=DATA(I2)
      DATA(I2)=DATA(I2)+DATA(I2+1)
 80   DATA(I2+1)=TEMPR-DATA(I2+1)
      IF (N-2) 200,200,90
 90   THETA=TWOPI/FLOAT(N)
      SINTH=SIN(THETA/2.)
      ZSTPR=-2.*SINTH*SINTH
      ZSTPI=SIN(THETA)
      ZR=(1.-ZSTPI)/2.
      ZI=(1.+ZSTPR)/2.
      IF (IFORM) 100,110,110
 100  ZR=1.-ZR
      ZI=-ZI
 110  I1MIN=IP0+1
      I1MAX=IP0*(N/4)+1
      DO 190 I1=I1MIN,I1MAX,IP0
      DO 180 I2=I1,IP2,IP1
      I2CNJ=IP0*(N/2+1)-2*I1+I2
      IF (I2-I2CNJ) 150,120,120
 120  IF (ISIGN*(2*IFORM+1)) 130,140,140
 130  DATA(I2+1)=-DATA(I2+1)
 140  IF (IFORM) 170,180,180
 150  DIFR=DATA(I2)-DATA(I2CNJ)
      DIFI=DATA(I2+1)+DATA(I2CNJ+1)
      TEMPR=DIFR*ZR-DIFI*ZI
      TEMPI=DIFR*ZI+DIFI*ZR
      DATA(I2)=DATA(I2)-TEMPR
      DATA(I2+1)=DATA(I2+1)-TEMPI
      DATA(I2CNJ)=DATA(I2CNJ)+TEMPR
      DATA(I2CNJ+1)=DATA(I2CNJ+1)-TEMPI
      IF (IFORM) 160,180,180
 160  DATA(I2CNJ)=DATA(I2CNJ)+DATA(I2CNJ)
      DATA(I2CNJ+1)=DATA(I2CNJ+1)+DATA(I2CNJ+1)
 170  DATA(I2)=DATA(I2)+DATA(I2)
      DATA(I2+1)=DATA(I2+1)+DATA(I2+1)
 180  CONTINUE
      TEMPR=ZR-.5
      ZR=ZSTPR*TEMPR-ZSTPI*ZI+ZR
 190  ZI=ZSTPR*ZI+ZSTPI*TEMPR+ZI
C     RECURSION SAVES TIME, AT A SLIGHT LOSS IN ACCURACY.  IF AVAILABLE,
C     USE DOUBLE PRECISION TO COMPUTE ZR AND ZI.
 200  IF (IFORM) 270,210,210
C     UNPACK THE REAL TRANSFORM VALUES (TWO PER COLUMN)
 210  I2=IP2+1
      I1=I2
      J1=IP0*(N/2+1)*NREM+1
      GO TO 250
 220  DATA(J1)=DATA(I1)
      DATA(J1+1)=DATA(I1+1)
      I1=I1-IP0
      J1=J1-IP0
 230  IF (I2-I1) 220,240,240
 240  DATA(J1)=DATA(I1)
      DATA(J1+1)=0.
 250  I2=I2-IP1
      J1=J1-IP0
      DATA(J1)=DATA(I2+1)
      DATA(J1+1)=0.
      I1=I1-IP0
      J1=J1-IP0
      IF (I2-1) 260,260,230
 260  DATA(2)=0.
 270  RETURN
      END
