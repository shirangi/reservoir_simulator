	subroutine random1(N)

      INTEGER N
	DIMENSION RNSEQ1(N)
	INTEGER IDUM1
	OPEN (110, FILE='randomdev.dat', STATUS='UNKNOWN')
	OPEN (105, FILE='seed.dat', STATUS='UNKNOWN')
	READ(105,*) IDUM1
	OPEN (106, FILE='usedseeds.dat',ACCESS='APPEND')
	WRITE (106,*) IDUM1
	DO 5 I=1,N
	 RNSEQ1(I) = 0.0
    5 CONTINUE
	DO 7 I=1,N
       RNSEQ1(I) = GASDEV1(IDUM1)
       WRITE(110,*) RNSEQ1(I)  
    7 CONTINUE
	IF(IDUM1.GT.0) IDUM1 = - IDUM1
	REWIND(105)
	WRITE (105,*) IDUM1
	CLOSE(105)
      CLOSE(106)
	CLOSE(110)
	end subroutine random1
C
C  FUNCTION FOR GENERATING A NORMALLY (GAUSSIAN) DISTRIBUTED DEVIATE
C  WITH ZERO MEAN AND UNIT STANDARD DEVIATION,  USES FUNCTION RAN2
C  AS THE SOURCE OF THE UNIFORM DEVIATES
C
	FUNCTION gasdev1(idum)
      INTEGER idum
      REAL gasdev1
      INTEGER iset
	REAL fac,gset,rsq,v1,v2,ran2
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*ran3(idum)-1.
	  v2=2.*ran3(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev1=v2*fac
        iset=1
      else
        gasdev1=gset
        iset=0
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software
      FUNCTION ran3(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran3,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran3=min(AM*iy,RNMX)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software