      SUBROUTINE NLPMMA (IALG,
     &                   VAR,RESL,RESU,FF,G,DF,DG,
     &                   U,SCX,SCG,WA,ACTIVE,
     &                   NUMVAR,M,ME,
     &                   NMAX,MMAX,LWA,
     &                   MAXFUN,MAXIT,ISVAN,MIXU,MIXL,IAPPR,ITSUB,
     &                   IPRINT,IFAIL,FNAME,IFNSIZE,
     &                   ACC,SCBOU,ALM,SAU,SBU,SAL,SBL,
     &                   ASSCL,ASSCLU,ASSCLL,FU,FL,SA,SB,SC,DSTEP)
C     **********************************************************
C     *                                                        *
C     *    MAIN PROGRAM FOR METHODS OF MOVING ASYMPTOTES       *
C     *                                                        *
C     **********************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER FNAME*(*),CDUM*1
C
C     ******************************************************************
C     -                                                                -
C     -    INTEGER VARAIBLES                                           -
C     -    MAXFUN  --->  MAX. NUMBER OF FUNC CALLS IN LINESEARCH       -
C     -    MAXIT   --->  NUMB. OF ITERATIONS FOR MMA                   -
C     -    ISVAN   --->  TYPE OF ADAPTION RULE                         -
C     -                  0 : RULE OF FLEURY/KURITZ                     -
C     -                  1 : RULE OF SVANBERG                          -
C     -                  2 : LINEAR ADAPTION U=X/RS, L=X*RS            -
C     -                  3 : FIXED ASYMPTOTES                          -
C     -                  4 : MIX OF ALL RULES                          -
C     -    MIXU    --->  MIXRULES FOR ASYMPTOTES                       -
C     -    IAPR    --->  TYPE OF POLYNOMICAL APPR. IN LINE SEARCH      -
C     -                  2 : QUADRATICAL APPROXIMATION                 -
C     -                  3 : CUBIC APPROXIMATION                       -
C     -    ITSUB   --->  MAX. NUMBER OF CYCLES IN SUBPROBLEM           -
C     -    REAL VARIABLES  :                                           -
C     -    ACC     --->  ACCURACY IN PROGRAMM MMA                      -
C     -    SCBOU   --->  UPPER BOUND FOR AUTOM. SCALING OF OBJEKTIVE   -
C     -    ALM     --->  FACTOR FOR LINE SEARCH ADJUSTMENT             -
C     -                                                                -
C     -    VARIABLES FOR ASYMPTOTE UPDATE :                            -
C     -    SA      --->  ASYMPTOTES ADAPTION FACTOR IF ITER < 2        -
C     -    SB      --->  ASYMPTOTES ADAPTION FACTOR IF ITER => 2       -
C     -                  CONTRACTION                                   -
C     -    SC      --->  ASYMPTOTES ADAPTION FACTOR IF ITER => 2       -
C     -                  EXPANSION                                     -
C     -    DSTEP   --->  STEP SIZE IN LOCAL SUBPROBLEM                 -
C     -    SAU     --->  UPPER ASYMPTOTES ADAPTION FACTOR IF ITER < 2  -
C     -    SBU     --->  UPPER ASYMPTOTES ADAPTION FACTOR IF ITER => 2 -
C     -    SAL     --->  LOWER ASYMPTOTES ADAPTION FACTOR IF ITER < 2  -
C     -    SBL     --->  LOWER ASYMPTOTES ADAPTION FACTOR IF ITER => 2 -
C     -    ASSCL   --->  LINEAR ADAPTION FACTOR FOR ALL ASYMPTOTES     -
C     -    ASSCLU  --->  LINEAR ADAPTION FACTOR FOR UPPER ASYMPTOTES   -
C     -    ASSCLL  --->  LINEAR ADAPTION FACTOR FOR LOWER ASYMPTOTES   -
C     -    FIXUP   --->  FIX VALUE FOR UPPER ASYMPTOTE                 -
C     -    FIXLOW  --->  FIX VALUE FOR LOWER ASYMPTOTE                 -
C     -                                                                -
C     -    CHARACTER VARIABLE  :                                       -
C     -    PRINT   ---> PRINT MODE                                     -
C     -                  1 : LAST STATION                              -
C     -                  2 : ALL ITERATIONS                            -
C     -                  3 : ONLY OBJECTIVE FUNCTION                   -
C     ******************************************************************
C
C
C------------------------------- OPEN OUTPUT FILE AND FORWARD TO THE END
C
      IOUT= 45
      OPEN (IOUT,FILE=FNAME(1:IFNSIZE))
    1 READ (IOUT,1000,END=2) CDUM
      GOTO 1
    2 CONTINUE
C      
 1000 FORMAT(A1) 
C
      WRITE(IOUT,*) '%2 START OF OPTIMIZATION STRATEGY MMA '
      WRITE(IOUT,*) ' ==================================='
C
C------------------------------------------------------------CALL OF MMA
C
      NROPT=0
      ISCVAR=0
C
      CALL MMA1 (IALG,M,ME,MMAX,NUMVAR,NMAX,VAR,FF,G,DF,DG,RESL,RESU,
     *          U,SCX,SCG,WA,ACTIVE,
     *          ACC,SCBOU,SA,SB,SC,DSTEP,ISVAN,MAXIT,ITSUB,
     *          IPRINT,IFAIL,LWA,NROPT,ALM,MAXFUN,
     *          IAPPR,ISCVAR,IOUT,FU,FL,ASSCL,ASSCLU,
     *          ASSCLL,MIXU,MIXL,SAL,SBL,SAU,SBU)
C
C----------------------------------------------------END OF OPTIMIZATION
C
900   RETURN
C
      END
C=======================================================================
      SUBROUTINE MMA1 (IALG,M,ME,MMAX,N,NMAX,X,FF,G,DF,DG,XL,XU,
     1                 U,SCX,SCG,WA,ACTIVE,
     2                 ACC,SCBOU,SA,SB,SC,DSTEP,ISVAN,MAXIT,ITSUB,
     3                 IPRINT,IFAIL,LWA,NROPT,ALM,MAXFUN,
     4                 IAPPR,ISCVAR,IOUT,FU,FL,ASSCL,ASSCLU,
     5                 ASSCLL,MIXU,MIXL,SAL,SBL,SAU,SBU)
C-----------------------------------------------------------------------
C     METHOD OF MOVING ASYMPTOTES (MMA)
C-----------------------------------------------------------------------
C     DUAL MATHEMATICAL PROGRAMMING METHOD USING AN MODIFIED
C     VERSION OF FLETCHER/REEVES' CONJUGATE GRADIENT METHOD AS
C     OPTIMIZER.
C-----------------------------------------------------------------------
C     THIS PROGRAM SOLVES NONLINEAR PROGRAMMING PROBLEMS.
C     IT IS AN IMPLEMENTATION OF SVANBERG'S METHOD OF MOVING ASYMPTOTES
C     WHICH IS AN DUAL MATHEMATICAL PROGRAMMING ALGORITHM.
C-----------------------------------------------------------------------
C
C      MINIMIZE  ANY  NONLINEAR FUNCTION   F = F(X)   -->  MIN
C
C      SUBJECT TO NONLINEAR CONSTRAINTS    G(X) >= 0
C
C      AND THE RESTRICTIONS                XL <= X <= XU
C
C      (THE ALGORITHM IS NOT YET TESTED FOR THE USE OF EQUALITY CONSTR.)
C-----------------------------------------------------------------------
C     LITERATURE:
C     KRISTER SVANBERG:  THE METHOD OF MOVING ASYMPTOTES - A NEW METHOD
C             FOR STRUCTURAL OPTIMIZATION;
C             INTERNATION JOURNAL FOR NUMERICAL METHODS IN ENGINEERING,
C             VOL. 24, P. 359-373, 1987
C-----------------------------------------------------------------------
C     M      ... TOTAL NUMBER OF CONSTRAINTS
C     ME     ... NUMBER OF EQUALITY CONSTRAINTS
C     MMAX   ... MAXIMUM NUMBER OF CONSTRAINTS (MMAX >= M)
C     N      ... NUMBER OF VARIABLES
C     NMAX   ... MAXIMUM NUMBER OF VARIABLES (NMAX = N+1)
C     X      ... VARIABLES
C     F      ... OBJECTIVE FUNCTION
C     G      ... CONSTRAINTS
C     DF     ... OBJECTIVE DERIVATIVES WITH RESPECT TO X
C     DG     ... CONSTRAINT DERIVATIVES WITH RESPECT TO X
C     U      ... LAGRANGE MULTIPLIERS
C     XL     ... LOWER BOUNDS OF VARIABLES
C     XU     ... UPPER BOUNDS OF VARIABLES
C     SCX    ... VARIABLE SCALING FACTORS
C     SCG    ... CONSTRAINT SCALING FACTORS
C     SCF    ... OBJECTIVE SCALING FACTOR
C     ACTIVE ... INDICATOR OF ACTIVE CONSTRAINTS
C     ACC    ... REQUIRED ACCURACY OF CALCULATION
C     SCBOU  ... SCALE OBJECTIVE IF F > SCBOU
C     SA     ... ASYMPTOTES ADAPTATION FACTOR
C     SB     ... ASYMPTOTES ADAPTATION FACTOR
C     SC     ... ASYMPTOTES ADAPTATION FACTOR
C     DSTEP  ... STEP SIZE IN LOCAL SUBPROBLEM		  -
C     ISVAN  ... TYPE OF ADAPTATION RULES
C     MAXIT  ... MAXIMUM NUMBER OF ITERATIONS
C     ITSUB  ... MAXIMUM NUMBER OF CYCLES IN SUBPROBLEM
C     IPRINT ... PRINT CONTROLL
C     IFAIL  ... RETURN CODE
C                IFAIL = 0 :  REGULAR END OF ALGORITHM
C                IFAIL = 1 :  MAXIMUM NUMBER OF ITERATIONS REACHED
C                IFAIL = 2 :  SEARCH DIRECTION NOT PROFITABLE
C                IFAIL = 3 :  PREVENTION OF INCONSISTENCY FAILED
C     WA     ... REAL WORKING ARRAY
C     LWA    ... LENGTH OF WORKING ARRAY >WA<; LWA = 8*N+6*MMAX+N*MMAX
C     NROPT  ... TOTAL NUMBER OF OPTIMIZATION STEPS
C     ALM    ... INITIAL STEP LENGTH FACTOR IN LINE SEARCH
C     MAXFUN ... MAXIMUM FUNCTION CALLS IN LINE SEARCH
C     IAPPR  ... TYPE OF POLYNOMINAL APPROXIMATION IN LINE SEARCH
C                IAPPR = 2 :  QUADRATIC
C                IAPPR = 3 :  CUBIC
C     ISCVAR ... TYPE OF VARIABLE SCALING
C                ISCVAR = 0 :  NO SCALING
C                ISCVAR = 1 :  SCALING WITH RECIPROCAL INITIAL VALUES
C                ISCVAR = 2 :  SCALING WITH HEESIAN MATRIX (NOT IN MMA)
C     IOUT   ... OUTPUT DEVICE
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL*4 ACTIVE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,BARL=-1.0D20,BARU=1.0D20)
      DIMENSION X(N),G(MMAX),DF(N),DG(MMAX,N),XL(N),XU(N)
      DIMENSION U(MMAX),SCX(N),SCG(MMAX),WA(LWA),ACTIVE(MMAX)
C
      scf = 1.0d0
c
      do i=1,n
        SCX(i)=1.0d0
      enddo
c
      do i=1,mmax
        SCG(i)=1.0d0
      enddo
c
      do i=1,lwa
        wa(i)=0.0d0
      enddo
C
C-------------------------------------------------------------- ADRESSES
C
C
C     NX1  ... STATE OF VARAIABLES ONE ITERATION  BEFORE; X1(N)
C     NX2  ... STATE OF VARAIABLES TWO ITERATIONS BEFORE; X2(N)
C     NDA  ... UPDATE VALUE OF ASYMPTOTES; DA(N)
C     NAL  ... LOWER ASYMPTOTES; ASL(N)
C     NAU  ... UPPER ASYMPTOTES; ASU(N)
C     NALF ... LOWER MOVE LIMITS; ALF(N)
C     NBET ... UPPER MOVE LIMITS; BET(N)
C     NP0  ... FACTORS OF APPROXIMATED OBJECTIVE; P0(N)
C     NP   ... FACTORS OF APPROXIMATED CONSTRAINTS; P(M,N)
C     NR   ... FACTORS OF APPROXIMATED CONSTRAINTS; R(M)
C     NDU  ... INCREMENT OF LAGRANGE MULTIPLIERS; DU(M)
C     NDW  ... DERIVATIVES OF DUAL OBJECTIVE FUNCTION; W(M)
C     NS   ... SEARCH DIRECTION; S(M)
C     NZ   ... ARTIFICIAL VARIABLES TO PREVENT INCONSISTENCY; Z(M)
C     NPEN ... PENALTY FACTORS; PEN(M)
C
      NX1   = 1
      NX2   = NX1  + N
      NDA   = NX2  + N
      NAL   = NDA  + N
      NAU   = NAL  + N
      NALF  = NAU  + N
      NBET  = NALF + N
      NP0   = NBET + N
      NP    = NP0  + N
      NR    = NP   + MMAX*N
      NDU   = NR   + MMAX
      NDW   = NDU  + MMAX
      NS    = NDW  + MMAX
      NZ    = NS   + MMAX
      NPEN  = NZ   + MMAX
      NTEST = NPEN + MMAX - 1
C
C--------------------------------------------------------- PRINT HEADING
C
      IF (IPRINT.GT.0) WRITE (IOUT,2000) N,ME,M-ME,MAXIT,ITSUB,
     *                             ACC,SCBOU,ISVAN,MIXU,MIXL,FU,FL
      IF (IPRINT.GT.0) WRITE (IOUT,2001) ASSCL,ASSCLU,ASSCLL,SA,SB,SC,
     *                         DSTEP,SAU,SBU,SAL,SBL,ALM,MAXFUN,IAPPR
C
C--------------------------------------------------- CHECK WORKING ARRAY
C
      IF (NTEST.GT.LWA) THEN
        WRITE (IOUT,2010) LWA, NTEST
        GOTO 900
      ENDIF
C
C-------------------------------- CHECK FEASIBILITY OF ASYMPTOTE FACTORS
C
      IF (SA.GT.ONE.OR.SA.LT.ZERO) THEN
          IFAIL = 4
           WRITE (IOUT,2020)
           GOTO 900
      ELSE IF (SB.GT.ONE.OR.SB.LT.ZERO) THEN
          IFAIL = 4
           WRITE (IOUT,2020)
           GOTO 900
      ELSE IF (SBU.GT.ONE.OR.SBU.LT.ZERO) THEN
          IFAIL = 4
           WRITE (IOUT,2020)
           GOTO 900
      ELSE IF (SBL.GT.ONE.OR.SBL.LT.ZERO) THEN
          IFAIL = 4
           WRITE (IOUT,2020)
           GOTO 900
      ELSE IF (SAL.GT.ONE.OR.SAL.LT.ZERO) THEN
          IFAIL = 4
           WRITE (IOUT,2020)
           GOTO 900
      ELSE IF (SAU.GT.ONE.OR.SAU.LT.ZERO) THEN
          IFAIL = 4
           WRITE (IOUT,2020)
           GOTO 900
      ELSE IF (ASSCL.GT.ONE.OR.ASSCL.LT.ZERO) THEN
          IFAIL = 4
           WRITE (IOUT,2020)
           GOTO 900
      ELSE IF (ASSCLU.GT.ONE.OR.ASSCLU.LT.ZERO) THEN
          IFAIL = 4
           WRITE (IOUT,2020)
           GOTO 900
      ELSE IF (ASSCLL.GT.ONE.OR.ASSCLL.LT.ZERO) THEN
          IFAIL = 4
           WRITE (IOUT,2020)
           GOTO 900
      ENDIF
C
C------------------ CHECK FEASIBILITY OF UPPER AND LOWER FIXE ASYMPTOTES
C
      IF (ISVAN.EQ.3) THEN
        ICHK = 0
        FMIN = BARU
        DO 10 I=1,N
          IF (FL/SCX(I).GT.XL(I)) ICHK = ICHK+1
          FMIN = MIN(FMIN,XL(I)*SCX(I))
   10   CONTINUE
C
        IF (ICHK.GT.0) THEN
          WRITE (IOUT,2030) ICHK,FMIN
          IFAIL = 4
          GOTO 900
        ENDIF
      ENDIF
C
C-------------------------------------------------- PERFORM OPTIMIZATION
C
      CALL MMA2 (IALG,M,ME,MMAX,N,NMAX,X,FF,G,DF,WA(NALF),WA(NBET),
     1           DG,U,XL,XU,SCX,SCG,SCF,WA(NX1),WA(NX2),WA(NAL),WA(NAU),
     2           WA(NDA),WA(NP0),WA(NP),WA(NR),WA(NDU),
     3           WA(NDW),WA(NS),WA(NZ),WA(NPEN),ACTIVE,
     4           ACC,SCBOU,SA,SB,SC,DSTEP,ISVAN,MAXIT,ITSUB,IPRINT,
     5           IFAIL,NROPT,ALM,MAXFUN,IAPPR,ISCVAR,IOUT,FU,FL,ASSCL,
     6           ASSCLU,ASSCLL,MIXU,MIXL,SAL,SBL,SAU,SBU)
C
C-----------------------------------------------------------------------
  900 CONTINUE
C
      RETURN
C-----------------------------------------------------------------------
 2000 FORMAT(4X,' START OF METHOD OF MOVING ASYMPTOTES (MMA)'/
     1       4X,' =========================================='/
     2       4X,' NUMBER OF VARIABLES                :      N= ',I12/
     3       4X,' NUMBER OF EQUALITY CONSTRAINTS     :     ME= ',I12/
     4       4X,' NUMBER OF NON-EQUALITY CONSTRAINTS :      M= ',I12/
     5       4X,' MAXIMUM NUMBER OF ITERATIONS       :  MAXIT= ',I12/
     5       4X,' MAXIMUM NUMBER OF SUB CYCLES       :  ITSUB= ',I12/
     6       4X,' ACCURACY                           :    ACC= ',E12.5/
     7       4X,' SCALING CONTROLL                   :  SCBOU= ',E12.5/
     8       4X,' TYPE OF ASYMPTOTE ADAPTION         :  ISVAN= ',I12/
     9       4X,' 0=KURRITZ 1=SVANBERG 2=LINEAR 3= FIX 4=MIX   ',/
     1       4X,' TYPE OF UPPER ASYMTOTE ADAPTION    :   MIXU= ',I12/
     2       4X,' TYPE OF LOWER ASYMTOTE ADAPTION    :   MIXL= ',I12/
     3       4X,' 0=KURRITZ 1=SVANBERG 2=LINEAR 3= FIX         ',/
     4       4X,' FIXED VALUE FOR UPPER ASYMTOTE     :     FU= ',E12.5/
     5       4X,' FIXED VALUE FOR LOWER ASYMTOTE     :     FL= ',E12.5)
2001  FORMAT(4X,' GLOBAL ASYMTOTE SCALING FACTOR     :  ASSCL= ',E12.5/
     1       4X,' UPPER  ASYMTOTE SCALING FACTOR     : ASSCLU= ',E12.5/
     2       4X,' LOWER  ASYMTOTE SCALING FACTOR     : ASSCLL= ',E12.5/
     3       4X,' ASYMPTOTES ADAPTATION FACTOR       :     SA= ',E12.5/
     4       4X,' ASYMPTOTES ADAPTATION FACTOR       :     SB= ',E12.5/
     4       4X,' ASYMPTOTES ADAPTATION FACTOR       :     SC= ',E12.5/
     4       4X,' STEP SIZE IN LOCAL PROBLEM         :  DSTEP= ',E12.5/
     5       4X,' UPPER ASYMTOTE ADAPTION FACTOR     :    SAU= ',E12.5/
     6       4X,' UPPER ASYMTOTE ADAPTION FACTOR     :    SBU= ',E12.5/
     7       4X,' LOWER ASYMTOTE ADAPTION FACTOR     :    SAL= ',E12.5/
     8       4X,' LOWER ASYMTOTE ADAPTION FACTOR     :    SBL= ',E12.5/
     9       4X,' INITIAL STEP LENGTH                :    ALM= ',E12.5/
     1       4X,' MAX FUNCTION CALLS IN LINE SEARCH  : MAXFUN= ',I12/
     2       4X,' TYPE OF POLYNOMINAL APPROXIMATION  :  IAPPR= ',I12/)
C
 2010 FORMAT (' OPTMMA: WORKING ARRAY TOO SMALL: ',I9/
     1        '         SHOULD BE AT LEAST     : ',I9)
 2020 FORMAT (4X,' OPTMMA: ',
     1           'ASYMPTOTE SCALING FACTORS LESS ZERO OR GREATER ONE')
 2030 FORMAT (4X,' OPTMMA: LOWER BOUNDS OF ',I4,' VARIABLES'/
     1        4X,'         ARE IN CONFLICT WITH FIXED ASYMPTOTES'/
     2        4X,'         LOWER ASYMPTOTE < ',E12.5)
      END
C=======================================================================
      SUBROUTINE MMA2 (IALG,M,ME,MMAX,N,NMAX,X,FF,G,DF,ALF,BET,DG,U,
     &                 XL,XU,
     1                 SCX,SCG,SCF,X1,X2,ASL,ASU,DELA,P0,P,R,DELU,
     2                 DW,S,Z,PEN,ACTIVE,
     3                 ACC,SCBOU,SA,SB,SC,DSTEP,ISVAN,MAXIT,ITSUB,
     4                 IPRINT,IFAIL,NROPT,ALM,MAXFUN,IAPPR,ISCVAR,IOUT,
     5                 FU,FL,ASSCL,ASSCLU,ASSCLL,MIXU,MIXL,SAL,SBL,
     6                 SAU,SBU)
C-----------------------------------------------------------------------
C     METHOD OF MOVING ASYMPTOTES, MAIN ITERATION ROUTINE
C-----------------------------------------------------------------------
C     M      ... TOTAL NUMBER OF CONSTRAINTS
C     ME     ... NUMBER OF EQUALITY CONSTRAINTS
C     MMAX   ... MAXIMUM NUMBER OF CONSTRAINTS (MMAX >= M)
C     N      ... NUMBER OF VARIABLES
C     NMAX   ... MAXIMUM NUMBER OF VARIABLES (NMAX = N+1)
C     X      ... VARIABLES
C     F      ... OBJECTIVE FUNCTION
C     G      ... CONSTRAINTS
C     DF     ... OBJECTIVE DERIVATIVES WITH RESPECT TO X
C     DG     ... CONSTRAINT DERIVATIVES WITH RESPECT TO X
C     U      ... LAGRANGE MULTIPLIERS
C     XL     ... LOWER BOUNDS OF VARIABLES
C     XU     ... UPPER BOUNDS OF VARIABLES
C     SCX    ... VARIABLE SCALING FACTORS
C     SCG    ... CONSTRAINT SCALING FACTORS
C     SCF    ... OBJECTIVE SCALING FACTOR
C     X1     ... STATE OF VARIABLES ONE ITERATION  BEFORE
C     X2     ... STATE OF VARIABLES TWO ITERATIONS BEFORE
C     ASL    ... LOWER ASYMPTOTES
C     ASU    ... UPPER ASYMPTOTES
C     DELA   ... ASYMPTOTES IMPROVEMENT
C     P0     ... P- AND Q-FACTORS OF OBJECTIVE FUNCTION
C     P      ... P- AND Q-FACTORS OF CONSTRAINTS
C     R      ... R-FACTORS OF CONSTRAINTS
C     DELU   ... INCREMENT OF LAGRANGE MULTIPLIERS
C     DW     ... DERIVATIVES OF DUAL OBJECTIVE WITH RESPECT TO U
C     S      ... SEARCH DIRECTION OF FLETCHER/REEVES CONJUGATE GRADIENT
C     Z      ... ARTIFICIAL VARIABLES TO PREVENT INCONSISTENCY
C     PEN    ... PENALTY FACTORS
C     ACTIVE ... INDICATOR OF ACTIVE CONSTRAINTS
C     ACC    ... REQUIRED ACCURACY OF CALCULATION
C     SCBOU  ... SCALE OBJECTIVE IF F > SCBOU
C     SA     ... ASYMPTOTES ADAPTATION FACTOR
C     SB     ... ASYMPTOTES ADAPTATION FACTOR
C     SC     ... ASYMPTOTES ADAPTATION FACTOR
C     DSTEP  ... STEP SIZE IN LOCAL SUBPROBLEM		  -
C     ISVAN  ... TYPE OF ADAPTATION RULE
C     MAXIT  ... MAXIMUM NUMBER OF ITERATIONS
C     ITSUB  ... MAXIMUM NUMBER OF CYCLES IN SUBPROBLEM
C     IPRINT ... PRINT CONTROLL
C     IFAIL  ... RETURN CODE
C                IFAIL = 0 :  REGULAR END OF ALGORITHM
C                IFAIL = 1 :  MAXIMUM NUMBER OF ITERATIONS REACHED
C                IFAIL = 2 :  SEARCH DIRECTION NOT PROFITABLE
C                IFAIL = 3 :  LAGRANGE-FUNCTION LESSER THAN ZERO
C                IFAIL = 4 :  ASYMPTOTE SCALING FACTORS WRONG
C     IGRAD  ... TYPE OF GRADIENT CALCULATION
C                IGRAD = 0 :  NUMERICAL
C                IGRAD = 1 :  ANALYTICAL
C     NROPT  ... TOTAL NUMBER OF OPTIMIZATION STEPS
C     ALM    ... INITIAL STEP LENGTH FACTOR IN LINE SEARCH
C     MAXFUN ... MAXIMUM FUNCTION CALLS IN LINE SEARCH
C     IAPPR  ... TYPE OF POLYNOMINAL APPROXIMATION IN LINE SEARCH
C                IAPPR = 2 :  QUADRATIC
C                IAPPR = 3 :  CUBIC
C     ISCVAR ... TYPE OF VARIABLE SCALING
C                ISCVAR = 0 :  NO SCALING
C                ISCVAR = 1 :  SCALING WITH RECIPROCAL INITIAL VALUES
C     IOUT   ... OUTPUT DEVICE
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER TEXT*80
      LOGICAL*4 ACTIVE
      PARAMETER (ZERO=0.0D0,ONE=1.D0,TWO=2.0D0,FACMIN=1.0D-3,EPS=1.0D-1)
      PARAMETER (SEVEN=7.0D0,TEN=10.0D0,TWENTY=20.0D0)
C
      DIMENSION X(N),G(MMAX),DF(N),DG(MMAX,N),U(MMAX),XL(N),XU(N),SCX(N)
      DIMENSION SCG(MMAX),X1(N),X2(N),P(MMAX,N),P0(N),ALF(N),BET(N)
      DIMENSION R(MMAX),S(MMAX),ASL(N),ASU(N),ACTIVE(MMAX),DW(MMAX)
      DIMENSION DELA(N),DELU(MMAX),PEN(MMAX),Z(MMAX)
C
C---------------------------------------------------- SET INITIAL VALUES
C
      NRSOPT=NROPT
C
      FOLD   = ZERO
      ITER   = 0
      NFUNC  = 0
      ICON   = 0
      IPRSOL = 0
C
      IF (IPRINT.EQ.2) IPRSOL = 1
      IF (IPRINT.EQ.4) IPRSOL = 2
C
      W      = ZERO
      SUM0   = ZERO
      SRES0  = ZERO
      ACCSOL = ACC*ACC
C
      DO 10 J=1,M
        U(J) = ZERO
   10 CONTINUE
C
C---------------------------------------- CHECK FEASIBILITY OF VARIABLES
C
      DO 20 I=1,N
        X(I)  = MAX(XL(I),X(I))
        X(I)  = MIN(X(I),XU(I))
        X1(I) = X(I)
        X2(I) = X(I)
   20 CONTINUE
C
C----------------- DETERMINE OBJECTIVE AND CONSTRAINTS AT STARTING POINT
C
      NFUNC = NFUNC + 1
C
      DO 30 J=1,M
        ACTIVE(J) = .TRUE.
   30 CONTINUE
C
      CALL FUNC(IALG,ITER)
C
      F=FF
C
      CALL GRAD(IALG,ACTIVE)
C
C-------------------------- NORMALIZE ACTIVE CONSTRAINTS AND DERIVATIVES
C
      SCF = ONE
      IF (SCBOU.GT.ZERO) THEN
        IF (ABS(F).GT.SCBOU) SCF = ONE/SQRT(ABS(F))
        F = F*SCF
        DO 40 I=1,N
          DF(I) = DF(I)*SCF
   40   CONTINUE
      ENDIF
      DO 60 J=1,M
        SCG(J) = ONE
        IF (SCBOU.GT.ZERO) THEN
          IF (ABS(G(J)).GT.SCBOU) SCG(J) = ONE/SQRT(ABS(G(J)))
          G(J) = G(J)*SCG(J)
          IF (ACTIVE(J)) THEN
            DO 50 I=1,N
              DG(J,I) = DG(J,I)*SCG(J)
   50       CONTINUE
          ENDIF
        ENDIF
   60 CONTINUE
C
C***********************************************************************
C                                         MAIN ITERATION LOOP, LABEL 100
C***********************************************************************
C
  100 CONTINUE
C
       NROPT=NRSOPT+ITER
C
C-------------------------------------------- PRINT INTERMEDIATE ITERATE
C
      WRITE(IOUT,*)
      WRITE(IOUT,*)
      WRITE(IOUT,'(A,I5)') 'ITERATION ',ITER
      WRITE(IOUT,*)
      WRITE(IOUT,'(A)') 'VARIABLES:'
      DO IX=1,N
       WRITE (IOUT,'(A,I5,3E25.6)') 'VARIABLE ',IX,X(IX),ASL(IX),ASU(IX)
      ENDDO      
      WRITE(IOUT,*)      
      WRITE(IOUT,'(A,E20.10)') 'OBJECTIVE:',F
        IF (IPRINT.GT.1) THEN
          DO IX=1,N
            WRITE (IOUT,'(A,I5,E20.10)') 'GRADIENT VARIABLE ',IX,DF(IX)
          ENDDO
 	ENDIF
      WRITE(IOUT,*)      
      DO JX=1,M
        WRITE(IOUT,'(A,I5,2E20.10)') 'CONSTRAINT:',JX,G(JX),U(JX)
        IF (IPRINT.GT.1) THEN
          DO IX=1,N
          WRITE (IOUT,'(A,I5,E20.10)') 'GRADIENT VARIABLE ',IX,DG(JX,IX)
          ENDDO
	ENDIF
      ENDDO      
      WRITE(IOUT,*)      
C
C---------------------------------------- CHECK MAXIMUM ITERATION NUMBER
C
      IF (ITER.GE.MAXIT) THEN
        WRITE (IOUT,2050)
        IFAIL = 1
        GOTO 900
      ENDIF
C
      ITER = ITER + 1
      ISET = ISET + 1
C
C------------------------------------------------------- FORM SUBPROBLEM
C
      CALL MMASET (X,XL,XU,X1,X2,ASL,ASU,DELA,F,G,DF,DG,Z,FARS,U,P0,R0,
     *          P,R,PEN,ACTIVE,SCF,SCG,SA,SB,SC,DSTEP,ISVAN,ACC,N,ME,M,
     *          MMAX,ITER,FU,FL,ASSCL,ALF,BET,SCX,ASSCLU,ASSCLL,
     *          MIXU,MIXL,SAL,SBL,SAU,SBU,IOUT,IPRINT)
C
      DO 110 I=1,N
        X2(I) = X1(I)
        X1(I) = X(I)
  110 CONTINUE
C
C------------------------------------------------------ SOLVE SUBPROBLEM
C
      CALL MMASUB (X,X1,XL,XU,ASL,ASU,U,DELU,W,DW,S,P0,R0,FARS,
     *             P,R,Z,PEN,ACTIVE,ALM,ACCSOL,MAXFUN,IAPPR,ITSUB,
     *             N,ME,M,MMAX,IPRSOL,IFSOL,ITSOL,IOUT,ALF,BET)
C
      IF (IPRINT.EQ.2 .OR. IPRINT.EQ.4) THEN
        WRITE (IOUT,2060) ITSOL,ITSOL/MMAX
        IF (IFSOL.EQ.1) WRITE (IOUT,2070)
      ENDIF
C
C----------------------------------------------- PROJECTION OF VARIABLES
      DO 115 I=1,N
        X(I) = MAX(XL(I),X(I))
        X(I) = MIN(X(I),XU(I))
  115 CONTINUE
C
C--------------------------------------- DUAL WEIGHT WITH ACTUAL X AND U
      CALL MMADWE (X,ASL,ASU,U,P0,R0,P,R,W,ACTIVE,N,ME,M,MMAX)
C
C-------------------------------------- CHECK CONSISTENCY OF CONSTRAINTS
C
      CONSI = ZERO
      DO 120 J=1,M
        IF (ACTIVE(J)) CONSI = CONSI + Z(J)
  120 CONTINUE
      IF (CONSI.EQ.ZERO) THEN
        ICON = 0
      ELSE
        ICON = ICON+1
      ENDIF
C
      IF (ICON.EQ.1) THEN
        IF (IPRINT.GT.0) WRITE (IOUT,2100)
      ELSE IF (ICON.GT.1) THEN
        IFAIL = 3
        IF (IPRINT.GT.0) THEN
          WRITE (IOUT,2110)
          DO 130 J=1,M
            IF (ACTIVE(J) .AND. Z(J).NE.ZERO) WRITE (IOUT,2120) J,Z(J)
  130     CONTINUE
        ENDIF
        GOTO 900
      ENDIF
C
C------------------------------------------------- CHECK FOR CONVERGENCE
C
      DFDEL = ZERO
      DO 150 I=1,N
        DFDEL = DFDEL + DF(I)*(X(I)-X1(I))
  150 CONTINUE
C
      SRES1 = ZERO
      SUM1 = ABS(DFDEL)
C
      IF(M.GT.0) THEN
        DO 160 J=1,M
        UAD = ABS(G(J))
        IF(J.LE.ME.OR. (ACTIVE(J).AND.G(J).LT.ZERO)) SRES1 = SRES1+UAD
  160   SUM1 = SUM1 + ABS(U(J)*G(J))
      ENDIF
      SUM1 = SUM1/SCF
C
      DSUM = ABS(SUM0-SUM1)
      DRES = ABS(SRES0-SRES1)
      IF ((SUM1.LT.ACC .AND. SRES1.LE.SQRT(ACC)) .OR.
     *    (DSUM.LT.ACC .AND.  DRES.LE.SQRT(ACC))) THEN
        IFAIL = 0
        F = F/SCF
        W = W/SCF
        DO 170 J=1,M
          U(J) = U(J)*SCG(J)/SCF
          G(J) = G(J)/SCG(J)
  170   CONTINUE
        GOTO 900
      ENDIF
      SUM0 = SUM1
      SRES0 = SRES1
C
C------------- DETERMINE OBJECTIVE FUNCTION, CONSTRAINTS AND DERIVATIVES
C
      NFUNC = NFUNC+1
C
      CALL FUNC(IALG,ITER)
C
      F=FF
c
C------------------------------------------------------------ ACTIVE SET
c
      DO 180 J=ME+1,M
        ACTIVE(J) = .TRUE.
        IF (G(J).GT.EPS .AND. U(J).EQ.ZERO) ACTIVE(J)=.FALSE.
  180 CONTINUE
C
      CALL GRAD (IALG,ACTIVE)
C
C--------------------------------------- SCALE OBJECTIVE AND CONSTRAINTS
C
      IF (SCBOU.GT.ZERO) THEN
        F = F*SCF
        DO 190 I=1,N
          DF(I) = DF(I)*SCF
  190   CONTINUE
        DO 210 J=1,M
          G(J) = G(J)*SCG(J)
          IF (ACTIVE(J)) THEN
            DO 200 I=1,N
              DG(J,I) = DG(J,I)*SCG(J)
  200       CONTINUE
          ENDIF
  210   CONTINUE
      ENDIF
C
C-------------------------------------------- CONTINUE ITERATION PROCESS
C
      GO TO 100
C
C---------------------------------------------- FINAL CONVERGENCE REPORT
C
  900 CONTINUE
      IF (IPRINT.GT.0) THEN
        WRITE (IOUT,2080)
        WRITE (IOUT,2010) F,W
        WRITE (IOUT,2020) (X(I)*SCX(I),I=1,N)
        IF (M.GT.0) THEN
          WRITE (IOUT,2030) (U(J),J=1,M)
          WRITE (IOUT,2040) (G(J),J=1,M)
          WRITE (IOUT,2045) (Z(I),I=1,M)
        ENDIF
        WRITE (IOUT,2090) NFUNC
      ENDIF
C
      RETURN
C--------------------------------------------------------------- FORMATS
 2000 FORMAT (/4X,'ITERATION NO.',I5)
 2005 FORMAT (/4X,'STARTING POINT:',I5)
 2010 FORMAT (7X,'F(X) = ',E16.8,'  W(U) = ',E16.8)
 2015 FORMAT (4X,'ITERATION NO.',I4,' F(X) = ',E16.8,'  W(U) = ',E16.8)
 2020 FORMAT (7X,'PRIMAL VARIABLES :',/,(8X,4E16.8))
 2030 FORMAT (7X,'MULTIPLIERS:',/,(8X,4E16.8))
 2040 FORMAT (7X,'CONSTRAINTS:',/,(8X,4E16.8))
 2045 FORMAT (7X,'INDICATOR OF INCONSISTENCY :',/,(8X,4E16.8))
 2050 FORMAT (4X,'*** WARNING ***  MAXIMUM NUMBER OF ITERATIONS ',
     *            'REACHED')
 2060 FORMAT (7X,'SUBPROBLEM: ',I5,' ITERATIONS IN ',I5,' CYCLES')
 2070 FORMAT (4X,'*** WARNING ***  MAXIMUM NUMBER OF ITERATIONS ',
     *            'REACHED IN SUBPROBLEM')
 2080 FORMAT (/4X,'FINAL CONVERGENCE REPORT:')
 2090 FORMAT (4X,'NUMBER OF FUNCTION CALLS: ',I5)
 2100 FORMAT (4X,'*** WARNING ***  PREVENTION OF INCONSISTENCY')
 2110 FORMAT (4X,'*** WARNING ***  PREVENTION OF INCONSISTENCY FAILED'/
     *        7X,'CHECK THE FOLLOWING CONSTRAINTS FOR CONSISTENCY:')
 2120 FORMAT (7X,'CONSTRAINT NO.',I5,'   ARTIFICICAL VARIABLE : ',E12.5)
C
      END
C=======================================================================
      SUBROUTINE MMASET (X,XL,XU,X1,X2,ASL,ASU,DELA,F,G,DF,DG,Z,FARS,
     *                   U,P0,R0,P,R,PEN,ACTIVE,SCF,SCG,SA,SB,SC,DSTEP,
     *                   ISVAN,ACC,N,ME,M,MMAX,ITER,FU,FL,ASSCL,ALF,
     *                   BET,SCX,ASSCLU,ASSCLL,MIXU,MIXL,SAL,SBL,SAU,
     *                   SBU,IOUT,IPRINT)
C-----------------------------------------------------------------------
C     DETERMINE PARAMETERS TO FORM THE DUAL SUBPROBLEM
C-----------------------------------------------------------------------
C     INSTITUT FUER BAUSTATIK  *   STUTTGART   *  MAR. 19TH 1987  *  KUB
C-----------------------------------------------------------------------
C     X     ...  VARIABLES
C     XL    ...  LOWER BOUNDS OF VARIABLES
C     XU    ...  UPPER BOUNDS OF VARIABLES
C     X1    ...  STATE OF VARIABLES ONE ITERATION  BEFORE
C     X2    ...  STATE OF VARIABLES TWO ITERATIONS BEFORE
C     ASL   ...  LOWER ASYMPTOTES
C     ASU   ...  UPPER ASYMPTOTES
C     DELA  ...  ASYMPTOTE IMPROVEMENTS
C     F     ...  OBJECTIVE FUNCTION
C     G     ...  CONSTRAINTS
C     DF    ...  OBJECTIVE DERIVATIVES WITH RESPECT TO X
C     DG    ...  CONSTRAINT DERIVATIVES WITH RESPECT TO X
C     U     ...  LAGRANGE MULTIPLIERS
C     P0    ...  P- AND Q-FACTOR OF OBJECTIVE FUNCTION
C     R0    ...  R-FACTOR OF OBJECTIVE FUNCTION
C     P     ...  P- AND Q-FACTORS OF CONSTRAINTS
C     R     ...  R-FACTORS OF CONSTRAINTS
C     PEN   ...  PENALTY FACTOR OF ADDITIONAL VARIABLE
C     ACTIVE...  INDICATOR OF ACTIVE CONSTRAINTS
C     SCF   ...  OBJECTIVE SCALING FACTOR
C     SCG   ...  CONSTRAINT SCALING FACTORS
C
C     VARAIBLES FOR ASYMPTOTE ADAPTION RULES :
C
C     ISVAN ...  TYPE OF ASYMPTOTES UPDATE
C                ISVAN = 0 :  RULES BY KURITZ/FLEURY
C                ISVAN = 1 :  RULES BY SVANBERG
C                ISVAN = 2 :  LINEAR ADAPTION RULE
C                ISVAN = 3 :  FIX VALUES FOR LOWER AND UPPER ASYMPTOTES
C                ISVAN = 4 :  MIX BETWEEN ALL RULES
C     SA    ...  ASYMPTOTES ADAPTION FACTOR
C     SB    ...  ASYMPTOTES ADAPTION FACTOR
C     SC     ... ASYMPTOTES ADAPTATION FACTOR
C     DSTEP  ... STEP SIZE IN LOCAL SUBPROBLEM		  -
C     FU    ...  FIX VALUE FOR UPPER ASYMPTOTE
C     FL    ...  FIX VALUE FOR LOWER ASYMPTOTE
C     ASSCL ...  LINEAR ADAPTION FACTOR FOR ASYMPTOTES
C     SAU   ...  UPPER ASYMPTOTES ADAPTION FACTOR IF ITER < 2  -
C     SBU   ...  UPPER ASYMPTOTES ADAPTION FACTOR IF ITER => 2 -
C     SAL   ...  LOWER ASYMPTOTES ADAPTION FACTOR IF ITER < 2  -
C     SBL   ...  LOWER ASYMPTOTES ADAPTION FACTOR IF ITER => 2 -
C     ASSCLU...  LINEAR ADAPTION FACTOR FOR UPPER ASYMPTOTES   -
C     ASSCLL...  LINEAR ADAPTION FACTOR FOR LOWER ASYMPTOTES   -
C
C     ACC   ...  REQUIRED ACCURACY
C     N     ...  NUMBER OF VARIABLES
C     ME    ...  NUMBER OF EQUALITY CONSTRAINTS
C     M     ...  TOTAL NUMBER OF CONSTRAINTS
C     MMAX  ...  MAXIMUM NUMBER OF CONSTRAINTS (MMAX >= M)
C     ITER  ...  ACTUAL ITERATION STEP
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL*4 ACTIVE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,FAC=0.9D0,PENINI=1.0D3,PFAC=1.D03)
      PARAMETER (BAR=1.0D-6)
C
      DIMENSION X(N),XU(N),XL(N),G(MMAX),DF(N),DG(MMAX,N),U(MMAX),ASU(N)
      DIMENSION P(MMAX,N),R(MMAX),X1(N),X2(N),ASL(N),P0(N),DELA(N)
      DIMENSION PEN(MMAX),ACTIVE(MMAX),SCG(MMAX),ALF(N),BET(N),SCX(N)
      DIMENSION Z(MMAX)
C
C---------------------------------------------------- SET INITIAL VALUES
C
      R0 = F
      DO 10 J=1,M
        IF (ACTIVE(J)) THEN
          R(J)  = -G(J)
          IF (ABS(U(J)).GT.ZERO) THEN
            PEN(J) = ABS(U(J))*PFAC
          ELSE
            PMAX = PENINI/SCG(J)*SCF
            DFLEN = ZERO
            DGLEN = ZERO
            DO 5 I=1,N
              DFLEN = DFLEN + DF(I)*DF(I)
              DGLEN = DGLEN + DG(J,I)*DG(J,I)
    5       CONTINUE
            IF (DGLEN.GT.BAR) PMAX = SQRT(DFLEN/DGLEN)*PFAC
            PEN(J) = PMAX
          ENDIF
        ENDIF
   10 CONTINUE
C
C---------------------------------------------------------- MAIN LOOP 40
C
      DO 30 I=1,N
C
C---------------------------------- DETERMINE LOWER AND UPPER ASYMPTOTES
C
C------------------------------   INITIALISATION OF CONVEX LINEARISATION
C
          IF(ABS(FU).LT.ACC) THEN
             FU1 = 10.0D0*XU(I)
          ELSE
             FU1 = FU*XU(I)
             IF (FU1.EQ.ZERO) FU1 = FU/SCX(I)
           ENDIF
C
C------------------------------------------ UPDATE RULE OF KURITZ/FLEURY
C
        IF (ISVAN.EQ.0) THEN
          IF (ITER.LE.2) THEN
            DELA(I) = (XU(I)-XL(I))*SA
          ELSE
            CHSIGN = (X(I)-X1(I))*(X1(I)-X2(I))
            IF (CHSIGN.LT.SQRT(ACC)) THEN
              DELA(I) = DELA(I)*SB
            ELSE
              DELA(I) = DELA(I)*SC
            ENDIF
          ENDIF
          ASL(I) = X(I)-DELA(I)
          ASU(I) = X(I)+DELA(I)
        ENDIF
C
C--------------------------------------------------- SVANBERG'S FORMULAE
C
        IF (ISVAN.EQ.1) THEN
          IF (ITER.LE.2) THEN
            ASL(I) = X(I) - (XU(I)-XL(I))*SA
            ASU(I) = X(I) + (XU(I)-XL(I))*SA
          ELSE
            CHSIGN = (X(I)-X1(I))*(X1(I)-X2(I))
            IF (CHSIGN.LT.SQRT(ACC)) THEN
              ASL(I) = X(I) - (X1(I)-ASL(I))*SB
              ASU(I) = X(I) + (ASU(I)-X1(I))*SB
            ELSE
              ASL(I) = X(I) - (X1(I)-ASL(I))*SC
              ASU(I) = X(I) + (ASU(I)-X1(I))*SC
            ENDIF
          ENDIF
        ENDIF
C
C---------------------------------------------------- LINEAR UPDATE RULE
C
       IF (ISVAN .EQ.2 ) THEN
C
          ASL(I) = X(I)*ASSCL
          ASU(I) = X(I)/ASSCL
C
        ENDIF
C
C--------------------------------------------- FIX VALUES FOR ASYMPTOTES
C
        IF (ISVAN .EQ.3 ) THEN
C
          ASL(I) = FL/SCX(I)
          ASU(I) = FU1
C
        ENDIF
C
C------------------------------------- MIXED UPDATE RULES FOR ASYMPTOTES
        IF (ISVAN.EQ.4) THEN
C----------------------------------------- KURRITZ UPDATE RULE FOR UPPER
          IF (MIXU.EQ.0) THEN
            IF (ITER.LE.2) THEN
              DELA(I) = (XU(I)-XL(I))*SAU
            ELSE
              CHSIGN = (X(I)-X1(I))*(X1(I)-X2(I))
                IF (CHSIGN.LT.SQRT(ACC)) THEN
                  DELA(I) = DELA(I)*SBU
                ELSE
                  DELA(I) = DELA(I)/SBU
                ENDIF
            ENDIF
          ENDIF
          ASU(I) = X(I)+DELA(I)
C---------------------------------------- SVANBERG UPDATE RULE FOR UPPER
          IF (MIXU .EQ.1) THEN
            IF (ITER.LE.2) THEN
              ASU(I) = X(I) + (XU(I)-XL(I))*SAU
            ELSE
              CHSIGN = (X(I)-X1(I))*(X1(I)-X2(I))
                IF (CHSIGN.LT.SQRT(ACC)) THEN
                  ASU(I) = X(I) + (ASU(I)-X1(I))*SBU
                ELSE
                  ASU(I) = X(I) + (ASU(I)-X1(I))/SBU
                ENDIF
            ENDIF
          ENDIF
C------------------------------------------ LINEAR UPDATE RULE FOR UPPER
          IF (MIXU .EQ. 2) THEN
             ASU(I) = X(I)/ASSCLU
          ENDIF
C------------------------------------------- FIXED UPDATE RULE FOR UPPER
          IF (MIXU .EQ. 3) THEN
             ASU(I) = FU1
          ENDIF
C----------------------------------------- KURRITZ UPDATE RULE FOR LOWER
          IF (MIXL.EQ.0) THEN
            IF (ITER.LE.2) THEN
              DELA(I) = (XU(I)-XL(I))*SAL
            ELSE
              CHSIGN = (X(I)-X1(I))*(X1(I)-X2(I))
                IF (CHSIGN.LT.SQRT(ACC)) THEN
                  DELA(I) = DELA(I)*SBL
                ELSE
                  DELA(I) = DELA(I)/SBL
                ENDIF
            ENDIF
          ENDIF
          ASL(I) = X(I)+DELA(I)
C---------------------------------------- SVANBERG UPDATE RULE FOR LOWER
          IF (MIXL .EQ.1) THEN
            IF (ITER.LE.2) THEN
              ASL(I) = X(I) + (XU(I)-XL(I))*SAL
            ELSE
              CHSIGN = (X(I)-X1(I))*(X1(I)-X2(I))
                IF (CHSIGN.LT.SQRT(ACC)) THEN
                  ASL(I) = X(I) + (ASL(I)-X1(I))*SBL
                ELSE
                  ASL(I) = X(I) + (ASL(I)-X1(I))/SBL
                ENDIF
            ENDIF
          ENDIF
C------------------------------------------ LINEAR UPDATE RULE FOR LOWER
          IF (MIXL .EQ. 2) THEN
             ASL(I) = X(I)/ASSCLL
          ENDIF
C------------------------------------------- FIXED UPDATE RULE FOR LOWER
          IF (MIXL .EQ. 3) THEN
             ASL(I) = FL
          ENDIF
        ENDIF
C---------------------------------- UPPER AND LOWER BOUNDS FOR ASYMPTOTE
C	
        IF (ASL(I) .GT. X(I) ) THEN
	  WRITE(*,*) 'MMA - WARNING: ASL(I) > X(I)'
	  CALL FLUSH(6)
	  ASL(I)=XL(I)-1.0*SCX(I)
	ENDIF
        IF (ASU(I) .LT. X(I) ) THEN
	  WRITE(*,*) 'MMA - WARNING: ASU(I) < X(I)'
	  CALL FLUSH(6)
	  ASU(I)=XU(I)+1.0*SCX(I)
	ENDIF
C
        IF (ISVAN.NE.3) THEN
           ASL(I) = MAX(ASL(I), X(I) - 10.00D0 * (XU(I)-XL(I)) )
           ASL(I) = MIN(ASL(I), X(I) -  0.01D0 * (XU(I)-XL(I)) )
           ASU(I) = MAX(ASU(I), X(I) +  0.01D0 * (XU(I)-XL(I)) )
           ASU(I) = MIN(ASU(I), X(I) + 10.00D0 * (XU(I)-XL(I)) )
        ENDIF
C
        DIFL = X(I) - ASL(I)
        DIFU = ASU(I) - X(I)
C--------------------------------------------------- DETERMINE ALFA,BETA
C
        ALF(I) = MAX(0.9D0*ASL(I) + 0.1D0*X(I),XL(I))
        BET(I) = MIN(0.9D0*ASU(I) + 0.1D0*X(I),XU(I))
C
        ALF(I) = MAX(ALF(I),X(I) - DSTEP*(XU(I)-XL(I)))
	BET(I) = MIN(BET(I),X(I) + DSTEP*(XU(I)-XL(I)))
C	
C------------------------------------------- DETERMINE P,Q AND R FACTORS
C
        IF (DF(I).GT.ZERO) THEN
          P0(I) = DIFU*DIFU*DF(I)
          R0 = R0 - P0(I)/DIFU
        ELSE
          P0(I) = DIFL*DIFL*DF(I)
          R0 = R0 + P0(I)/DIFL
        ENDIF
C
        DO 20 J=1,M
          IF (ACTIVE(J)) THEN
            IF (DG(J,I).GT.ZERO) THEN
              P(J,I) = -DIFL*DIFL*DG(J,I)
              R(J) = R(J) + P(J,I)/DIFL
            ELSE
              P(J,I) = -DIFU*DIFU*DG(J,I)
              R(J) = R(J) - P(J,I)/DIFU
            ENDIF
          ENDIF
   20   CONTINUE
   30 CONTINUE
C
C----------------------------------- CALCULATION OF ARTIFICIAL VARIABLES
C
      Z0   = ZERO
      FARS = ZERO
      SUMU = ZERO
      FARS = R0
C      
71    CONTINUE
C
C------------------------------------------ WRITE ASYMTOTES ON FILE IOUT
      IF (IPRINT.EQ.2 .OR. IPRINT.EQ.4) THEN
        WRITE (IOUT,1000) (ASL(I)*SCX(I),I=1,N)
        WRITE (IOUT,1010) (ASU(I)*SCX(I),I=1,N)
        WRITE (IOUT,1020) (Z(J)/SCG(J),J=1,M)
      ENDIF
C
C-----------------------------------------------------------------------
      RETURN
C-----------------------------------------------------------------------
 1000 FORMAT (7X,'LOWER ASYMPTOTES :'/(8X,4E16.8))
 1010 FORMAT (7X,'UPPER ASYMPTOTES :'/(8X,4E16.8))
 1020 FORMAT (7X,'ARTIFICIAL VARIABLES :'/(8X,4E16.8))
      END

C=======================================================================

      SUBROUTINE MMADUA (IDERI,X,X1,XL,XU,ASL,ASU,U,ALF,BET,IOUT,
     *              P0,R0,P,R,W,DW,Z,PEN,ACTIVE,N,ME,M,MMAX,IRET,
     *              FARS)
C-----------------------------------------------------------------------
C     DUAL OBJECTIVE FUNCTION AND ITS DERIVATIVES WITH RESPECT
C     TO THE LAGRANGE MULTIPLIERS AS DUAL VARIABLES
C     PRIMAL VARIABLES AS A FUNCTION OF THE DUAL ONES
C-----------------------------------------------------------------------
C     IDERI ...  DETERMINE GRADIENT DW (0=NO,1=YES)
C     X     ...  VARIABLES
C     X1    ...  VARIABLES ONE ITERATION STEP BEFORE
C     XL    ...  LOWER BOUNDS OF VARIABLES
C     XU    ...  UPPER BOUNDS OF VARIABLES
C     ASL   ...  LOWER ASYMPTOTES
C     ASU   ...  UPPER ASYMPTOTES
C     U     ...  LAGRANGE MULTIPLIERS
C     P0    ...  P-AND Q-FACTOR OF OBJECTIVE FUNCTION
C     R0    ...  R-FACTOR OF OBJECTIVE FUNCTION
C     P     ...  P- AND Q-FACTORS OF CONSTRAINTS
C     R     ...  R-FACTORS OF CONSTRAINTS
C     W     ...  DUAL OBJECTIVE FUNCTION
C     DW    ...  DERIVATIVES OF DUAL OBJECTIVE WITH RESPECT TO U
C     Z     ...  ARTIFICIAL VARIABLES TO PREVENT INCONSISTENCY
C     PEN   ...  PENALTY FACTORS OF ARTIFICIAL VARIABLES
C     ACTIVE...  INDICATOR OF ACTIVE CONSTRAINTS
C     N     ...  NUMBER OF VARIABLES
C     ME    ...  NUMBER OF EQUALITY CONSTRAINTS
C     M     ...  TOTAL NUMBER OF CONSTRAINTS
C     MMAX  ...  MAXIMUM NUMBER OF CONSTRAINTS (MMAX >= M)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL*4 ACTIVE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,FAC=0.9D0,PFAC=0.1D0)
C
      DIMENSION X(N),X1(N),XU(N),XL(N),U(MMAX),P0(N),PEN(MMAX),Z(MMAX)
      DIMENSION P(MMAX,N),R(MMAX),DW(MMAX),ASL(N),ASU(N),ACTIVE(MMAX)
      DIMENSION ALF(N),BET(N)
C
C---------------------------------------------------- SET INITIAL VALUES
C
C------------------------------------------------------- SET RETURN CODE
      IRET=0
C
      W  = ZERO
C      
      IF (IDERI.GT.0) THEN
        DO 10 J=1,M
          DW(J) = ZERO
   10   CONTINUE
      ENDIF
C
C---------------------------------------------------------- MAIN LOOP 40
C
      DO 40 I=1,N
C
        AL = ASL(I)
        AU = ASU(I)
C
        IF (P0(I).GE.ZERO) THEN
          AP =  P0(I)
          AQ =  ZERO
        ELSE
          AP =  ZERO
          AQ = -P0(I)
        ENDIF
        DO 20 J=1,M
          IF (ACTIVE(J)) THEN
            IF (P(J,I).GE.ZERO) THEN
              AP = AP + U(J)*P(J,I)
            ELSE
              AQ = AQ - U(J)*P(J,I)
            ENDIF
          ENDIF
   20   CONTINUE
C
        DLA = AP/((AU-ALF(I))*(AU-ALF(I)))-AQ/((ALF(I)-AL)*(ALF(I)-AL))
        DLB = AP/((AU-BET(I))*(AU-BET(I)))-AQ/((BET(I)-AL)*(BET(I)-AL))
C
C------------------------------------------ DETERMINE X AS FUNCTION OF U
C
        IF (DLA.GE.ZERO) THEN
          X(I) = ALF(I)
C         IRET = 1
        ELSE IF (DLB.LE.ZERO) THEN
          X(I) = BET(I)
C         IRET = 1
        ELSE IF (DLA.LT.ZERO .AND. DLB.GT.ZERO) THEN
          X(I) = (SQRT(AP)*AL+SQRT(AQ)*AU)/(SQRT(AP)+SQRT(AQ))
        ELSE
          PRINT*,' CHECK YOUR PROGRAM'
          STOP
        ENDIF
C
        DIFU = AU - X(I)
        DIFL = X(I) - AL
C--------------------------------- DETERMINE DUAL WEIGHT AND DERIVATIVES
C
        W = W + AP/DIFU + AQ/DIFL
        IF (IDERI.GT.0) THEN
          DO 30 J=1,M
            IF (ACTIVE(J)) THEN
              IF (P(J,I).GE.ZERO) THEN
                DW(J) = DW(J) + P(J,I)/DIFU
              ELSE
                DW(J) = DW(J) - P(J,I)/DIFL
              ENDIF
            ENDIF
   30     CONTINUE
        ENDIF
C
   40 CONTINUE
C
C------------------------ LAST CORRECTION OF DUAL WEIGHT AND DERIVATIVES
C
      W = W + R0
      DO 60 J=1,M
        IF (ACTIVE(J)) THEN
          W     = W + U(J)*R(J)
          IF (IDERI.GT.0) DW(J) = DW(J) + R(J)
        ENDIF
   60 CONTINUE
C
C---------------------------- CORRECT DUAL WEIGHT AND DERIVATIVES WITH Z
      DO 80 J=1,M
        Z(J) = ZERO
        IF (ACTIVE(J)) THEN
          Z(J) = MAX(ZERO,(U(J)-PEN(J))/(TWO*PEN(J)))
C---------------------------- DETERMINE DUAL WEIGHT WITH ARTIFICIAL VAR.
          W     = W - U(J)*Z(J) + Z(J)*PEN(J)*(ONE+Z(J))
          IF (IDERI.GT.0) DW(J) = DW(J) - Z(J)
C
        ENDIF
   80 CONTINUE
C
      RETURN
      END
C=======================================================================
      SUBROUTINE MMADWE (X,ASL,ASU,U,P0,R0,P,R,W,ACTIVE,N,ME,M,MMAX)
C-----------------------------------------------------------------------
C     DUAL OBJECTIVE FUNCTION WITH GIVEN PRIMAL AND DUAL VARIABLES
C-----------------------------------------------------------------------
C     X     ...  VARIABLES
C     ASL   ...  LOWER ASYMPTOTES
C     ASU   ...  UPPER ASYMPTOTES
C     U     ...  LAGRANGE MULTIPLIERS
C     P0    ...  P-AND Q-FACTOR OF OBJECTIVE FUNCTION
C     R0    ...  R-FACTOR OF OBJECTIVE FUNCTION
C     P     ...  P- AND Q-FACTORS OF CONSTRAINTS
C     R     ...  R-FACTORS OF CONSTRAINTS
C     W     ...  DUAL OBJECTIVE FUNCTION
C     ACTIVE...  INDICATOR OF ACTIVE CONSTRAINTS
C     N     ...  NUMBER OF VARIABLES
C     ME    ...  NUMBER OF EQUALITY CONSTRAINTS
C     M     ...  TOTAL NUMBER OF CONSTRAINTS
C     MMAX  ...  MAXIMUM NUMBER OF CONSTRAINTS (MMAX >= M)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL*4 ACTIVE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,FAC=0.9D0)
C
      DIMENSION X(N),U(MMAX),P0(N)
      DIMENSION P(MMAX,N),R(MMAX),ASL(N),ASU(N),ACTIVE(MMAX)
C---------------------------------------------------- SET INITIAL VALUES
C
      W  = R0
      DO 10 J=1,M
        IF (ACTIVE(J)) W = W + U(J)*R(J)
C       W = W + U(J)*R(J)
   10 CONTINUE
C
      DO 40 I=1,N
C
        IF (P0(I).GE.ZERO) THEN
          AP =  P0(I)
          AQ =  ZERO
        ELSE
          AP =  ZERO
          AQ = -P0(I)
        ENDIF
        DO 20 J=1,M
          IF (ACTIVE(J)) THEN
            IF (P(J,I).GE.ZERO) THEN
              AP = AP + U(J)*P(J,I)
            ELSE
              AQ = AQ - U(J)*P(J,I)
            ENDIF
          ENDIF
   20   CONTINUE
C
        DIFU = ASU(I) - X(I)
        DIFL = X(I) - ASL(I)
C------------------------------------------------- DETERMINE DUAL WEIGHT
C
        W = W + AP/DIFU + AQ/DIFL
C
   40 CONTINUE
C
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
      SUBROUTINE MMASUB (X,X1,XL,XU,ASL,ASU,U,UOLD,W,DW,S,P0,R0,FARS,
     *                   P,R,Z,PEN,ACTIVE,ALM,ACC,MAXFUN,IAPPR,ITSUB,
     *                   N,ME,M,MMAX,IPRINT,IFAIL,ITER,IOUT,ALF,BET)
C-----------------------------------------------------------------------
C     SOLUTION OF THE DUAL SUBPROBLEM VIA THE CONJUGATE GRADIENT METHOD
C     BY FLETCHER AND REEVES SLIGHTLY MODIFIED WITH RESPECT TO THE FACT
C     THAT THE MULTIPLIERS MUST NOT BE NEGATIV AT ANY STATE OF
C     CALCULATION
C-----------------------------------------------------------------------
C     INSTITUT FUER BAUSTATIK  *   STUTTGART   *  MAR. 19TH 1987  *  KUB
C-----------------------------------------------------------------------
C     X     ...  VARIABLES
C     X1    ...  VARIABLES OF ONE ITERATION STEP BEFORE
C     XL    ...  LOWER BOUNDS OF VARIABLES
C     XU    ...  UPPER BOUNDS OF VARIABLES
C     ASL   ...  LOWER ASYMPTOTES
C     ASU   ...  UPPER ASYMPTOTES
C     U     ...  LAGRANGE MULTIPLIERS
C     UOLD  ...  LAGRANGE MULTIPLIERS OF PREVIOUS ITERATION
C     W     ...  DUAL WEIGHT
C     DW    ...  GRADIENT OF DUAL OBJECTIVE FUNCTION WITH RESPECT TO U
C     S     ...  SEARCH DIRECTION
C     P0    ...  P- AND Q-FACTORS OF OBJECTIVE FUNCTION
C     R0    ...  R-FACTOR OF OBJECTIVE FUNCTION
C     P     ...  P- AND Q-FACTORS OF CONSTRAINTS
C     R     ...  R-FACTORS OF CONSTRAINTS
C     Z     ...  ARTIFICIAL VARIABLES TO PREVENT INCONSISTENCY
C     PEN   ...  PENALTY FACTOR OF ARTIFICIAL VARIABLES
C     ACTIVE...  INDICATOR OF ACTIVE CONSTRAINTS
C     ALM   ...  INITIAL STEP IN LINE SEARCH
C     ACC   ...  REQUIRED ACCURACY
C     MAXFUN...  MAX NUMBER OF FUNCTION CALLS IN LINE SEARCH
C     IAPPR ...  TYPE OF POLYNOMINAL APPROXIMATION IN LINE SEARCH
C     ITSUB ...  MAX NUMBER OF ITERATION CYCLES IN SUBPROBLEM
C     N     ...  NUMBER OF VARIABLES
C     ME    ...  NUMBER OF EQUALITY CONSTRAINTS
C     M     ...  TOTAL NUMBER OF CONSTRAINTS
C     MMAX  ...  MAXIMUM NUMBER OF CONSTRAINTS (MMAX >= M)
C     IPRINT...  OUTPUT CONTROLL
C     IFAIL ...  RETURN CODE
C     ITER  ...  NUMBER OF ITERATIONS IN SUBPROBLEM
C     IOUT  ...  OUTPUT DEVICE
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL*4 ACTIVE
      PARAMETER (ZERO=0.0D0)
C
      DIMENSION X(N),X1(N),XU(N),XL(N),U(MMAX),ASU(N),ASL(N),DW(MMAX)
      DIMENSION P(MMAX,N),R(MMAX),S(MMAX),P0(N),UOLD(MMAX)
      DIMENSION PEN(MMAX),Z(MMAX),ACTIVE(MMAX),ALF(N),BET(N)
C
C---------------------------------------------------- SET INITIAL VALUES
C
      ITER = 0
      ISET = 0
      IFAIL = 0
      MAXIT = M*ITSUB
C
C--------------- DETERMINE INITIAL VALUES OF DUAL WEIGHT AND DERIVATIVES
C
      CALL MMADUA (1,X,X1,XL,XU,ASL,ASU,U,ALF,BET,IOUT,P0,R0,P,R,W,DW,
     *             Z,PEN,ACTIVE,N,ME,M,MMAX,IRET,FARS)
      W = -W
C
C------------------------------------ DETERMINE INITIAL SEARCH DIRECTION
C
      MDIM = M
      GRN0 = ZERO
      DO 20 J=1,M
        IF (ACTIVE(J)) THEN
          IF (J.GT.ME .AND. U(J).EQ.ZERO) THEN
            DW(J) = MAX(ZERO,DW(J))
            IF (DW(J).EQ.ZERO) MDIM = MDIM-1
          ENDIF
          S(J) = DW(J)
          GRN0 = GRN0 + S(J)*S(J)
          UOLD(J) = ZERO
        ELSE
          MDIM = MDIM-1
        ENDIF
   20 CONTINUE
C
C------------------------------- RETURN IF THERE IS NO ACTIVE CONSTRAINT
C
      IF (MDIM.EQ.0) GOTO 900
C
C---------------------------------------- MAIN ITERATION LOOP, LABEL 100
C
  100 CONTINUE
C
      ITER = ITER + 1
      ISET = ISET + 1
C
C------------------------------ DETERMINE GRADIENT IN THE DIRECTION OF S
C
      GRADS = ZERO
      DO 110 J=1,M
        IF (ACTIVE(J) .AND. S(J).NE.ZERO) GRADS = GRADS - S(J)*DW(J)
  110 CONTINUE
C
      IF (GRADS.GT.ZERO) THEN
        IF (IPRINT.GT.0) WRITE (IOUT,2000) GRADS
        MDIM = M
        GRN0 = ZERO
        DO 120 J=1,M
          IF (ACTIVE(J)) THEN
            IF (J.GT.ME .AND. U(J).EQ.ZERO) THEN
              DW(J) = MAX(ZERO,DW(J))
              IF (DW(J).EQ.ZERO) MDIM = MDIM-1
            ENDIF
            S(J) = DW(J)
            GRN0 = GRN0 + S(J)*S(J)
          ELSE
            MDIM = MDIM-1
          ENDIF
  120   CONTINUE
        GRADS = -GRN0
      ENDIF
C
C----------------------------------------------------------- LINE SEARCH
C
      IF (MDIM.GT.0) THEN
        CALL MMASEA (X,X1,XL,XU,ASL,ASU,U,P0,R0,P,R,S,W,DW,Z,PEN,ACTIVE,
     *              GRADS,ALM,ACC,MAXFUN,IAPPR,N,ME,M,MMAX,IPRINT,IOUT,
     *              ALF,BET,FARS)
      ENDIF
C
C------------------------------------------ PROJECTION OF DUAL VARIABLES
C
      DO 130 J=1,M
        IF (J.GT.ME .AND. ACTIVE(J)) U(J) = MAX(ZERO,U(J))
  130 CONTINUE
C
C----- DETERMINE VALUES OF DUAL WEIGHT AND DERIVATIVES AFTER LINE SEARCH
C
      CALL MMADUA (1,X,X1,XL,XU,ASL,ASU,U,ALF,BET,IOUT,P0,R0,P,R,W,DW,
     *             Z,PEN,ACTIVE,N,ME,M,MMAX,IRET,FARS)
      W = -W
C
C------------------------------------------------- CHECK FOR CONVERGENCE
C
      CHK = ZERO
      DO 140 J=1,M
        IF (ACTIVE(J)) CHK = CHK+ABS(DW(J)*(U(J)-UOLD(J)))
  140 CONTINUE
C
      IF (CHK.LE.ACC) GOTO 900
C
C------------------------------------ CHECK MAXIMUM NUMBER OF ITERATIONS
C
      IF (ITER.GE.MAXIT) THEN
        IFAIL = 1
        GOTO 900
      ENDIF
C
C---------------------------------------- DETERMINE NEW SEARCH DIRECTION
C
      MCHK = M
      GRN1 = ZERO
      DO 155 J=1,M
        IF (ACTIVE(J)) THEN
          IF (J.GT.ME .AND. U(J).EQ.ZERO) THEN
            IF (DW(J).LE.ZERO) THEN
              DW(J) = ZERO
              MCHK = MCHK-1
            ENDIF
          ENDIF
          IF (DW(J).NE.ZERO) GRN1 = GRN1 + DW(J)*DW(J)
          UOLD(J) = U(J)
        ELSE
          MCHK = MCHK-1
        ENDIF
  155 CONTINUE
C
C------------------------------------------------------ STEEPEST DESCENT
C
      IF (ISET.GE.MDIM .OR. MCHK.NE.MDIM) THEN
        ISET = 0
        MDIM = MCHK
        GRN0 = GRN1
        DO 160 J=1,M
          IF (ACTIVE(J)) S(J) = DW(J)
  160   CONTINUE
C
C--------------------------------------------------- CONJUGATE GRADIENTS
C
      ELSE
        QUOT = GRN1/GRN0
        GRN0 = GRN1
C
        DO 180 J=1,M
          IF (ACTIVE(J)) S(J) = DW(J) + QUOT*S(J)
  180   CONTINUE
C
      ENDIF
C
C------------------------------------------------- END OF ITERATION LOOP
C
      GO TO 100
C
  900 CONTINUE
C
      RETURN
C-----------------------------------------------------------------------
 2000 FORMAT (4X,'*** WARNING ***  SEARCH DIRECTION NOT PROFITABLE'/
     *        20X,'PROD = ',E16.8/
     *        20X,'NEGATIV GRADIENT INSTEAD')
      END
C=======================================================================
      SUBROUTINE MMASEA (X,X1,XL,XU,ASL,ASU,U,P0,R0,P,R,S,W,DW,Z,PEN,
     *        ACTIVE,GRADS,ALM,ACC,MAXFUN,IAPPR,N,ME,M,MMAX,IPRINT,IOUT,
     *        ALF,BET,FARS)
C-----------------------------------------------------------------------
C     LINE SEARCH ; BILINEAR INTERPOLATION  ................ IAPPR = 1
C                   QUADRATIC INTERPOLATION ................ IAPPR = 2
C                   CUBIC INTERPOLATION (HERMITE) .......... IAPPR = 3
C-----------------------------------------------------------------------
C     LIST OF PARAMETERS IS THE SAME AS EVERYWHERE
C-----------------------------------------------------------------------
C     NOTICE: THE SUCCESS OF THE WHOLE METHOD IS HEAVILEY DEPENDENT
C             ON AN EXACT SOLUTION OF THE LINE SEARCH PROCEDURE.
C             THEREFORE THE CONVERGENCE CHECK IN THIS ROUTINE IS DONE
C             WITH A FIXED ACCURACY OF VALUE 1.0E-14.
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL*4 ACTIVE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,HA=.5D0)
      PARAMETER (FAC=10.D0,BAR=1.0D-15,UBAR=1.0D15,UF=1.0D-14)
      PARAMETER (ONFO=0.25D0,ONEI=0.125D0,GOLD1=.618D0,GOLD2=.382D0)
      DIMENSION X(N),X1(N),XU(N),XL(N),U(MMAX),ASU(N),ASL(N),DW(MMAX)
      DIMENSION P(MMAX,N),R(MMAX),S(MMAX),P0(N),PEN(MMAX),Z(MMAX)
      DIMENSION ACTIVE(MMAX),BET(N),ALF(N)
C
C
C-------------------------------------------- SET LINE SEARCH PARAMETERS
C
      ILINE = 0
      F0   = W
      DF0  = GRADS
      DF1  = DF0
      ALF0 = ZERO
      ALFL = ALF0
      ALF1 = ABS(ALM)
      ACL  = UF
      SQACL= SQRT(ACL)
      AC2  = MAX(UF,ACL*SQACL)
C
      IF (ABS(GRADS).LT.ACL) GOTO 900
C***********************************************************************
C                                                    SEARCH FOR INTERVAL
C***********************************************************************
   10 CONTINUE
      ILINE = ILINE+1
   15 CONTINUE
C
C---------------------------------- CHECK IF MORE THAN MAXFUN FUNC-CALLS
C
      IF (ILINE.GT.MAXFUN) THEN
        IF (IPRINT.GT.0) WRITE (IOUT,2000) MAXFUN
        GOTO 900
      ENDIF
C
      ALFMAX = UBAR
      DALF  = ALF1 - ALF0
      DALF1 = ALF1 - ALFL
      DO 20 J=1,M
        IF (ACTIVE(J) .AND. S(J).NE.ZERO) U(J) = U(J) + DALF1*S(J)
        IF (U(J).LT.ZERO) THEN
          ALFMX = ALF1-U(J)/S(J)
          IF (ALFMX.LT.ALFMAX) THEN
            ALFMAX = ALFMX
            JMAX   = J
          ENDIF
        ENDIF
   20 CONTINUE
      ALFL = ALF1
C
C-----------------------------------------------------------------------
C                   SEARCH DIRECTION REACHES BOUNDARY OF FEASIBLE REGION
C-----------------------------------------------------------------------
      IF (ALFMAX.LT.UBAR) THEN
        ALF1  = ALFMAX
        DALF  = ALF1 - ALF0
        DALF1 = ALF1 - ALFL
        DO 25 J=1,M
          IF (ACTIVE(J) .AND. S(J).NE.ZERO) U(J) = U(J) + DALF1*S(J)
          IF (ABS(U(J)).LT.ACL) U(J) = ZERO
   25   CONTINUE
        U(JMAX) = ZERO
        ALFL = ALF1
C
C------------------------------- CHANGE SEARCH DIRECTION AND START AGAIN
C                                    IF BEGIN OF INTERVAL IS AT BOUNDARY
        IF (ABS(DALF).LT.SQACL) THEN
          DF0  = ZERO
          DO 27 J=1,M
            IF (ACTIVE(J)) THEN
              IF (U(J).EQ.ZERO.AND.S(J).LT.ZERO) S(J)=ZERO
              DF0 = DF0 - S(J)*DW(J)
            ENDIF
   27     CONTINUE
          ALF1 = ALF1 + ABS(ALM)
          GOTO 15
        ENDIF
      ENDIF
C-----------------------------------------------------------------------
C
      CALL MMADUA (1,X,X1,XL,XU,ASL,ASU,U,ALF,BET,IOUT,P0,R0,P,R,W,DW,
     *             Z,PEN,ACTIVE,N,ME,M,MMAX,IRET,FARS)
      W = -W
      F1 = W
C
      GRADS = ZERO
      DO 30 J=1,M
        IF (ACTIVE(J) .AND. S(J).NE.ZERO) GRADS = GRADS - S(J)*DW(J)
   30 CONTINUE
      DDF1 = GRADS-DF1
      DF1  = GRADS
C
      DDF = DF1-DF0
      IF (ABS(DF1).LT.ACL .OR.ABS(DALF1).LT.AC2) GOTO 900
C
C------------------------------------------------------ IMPROVE INTERVAL
      DDF1 = DDF1/DALF1
      IF (DF1.LT.ZERO .OR. (DF1.GT.ZERO .AND.
     *      DF1.LT.SQACL .AND. ABS(DDF1).GT.SQACL)) THEN
        IF (ABS(DDF).LT.ACL) THEN
          ALFM = ALF1*FAC
        ELSE
          ALFM = ALF0 - DF0*(ALF1-ALF0)/DDF
        ENDIF
        IF (ALFM.LE.ALF0) THEN
          ALFM = HA*(ALF0+ALF1)
        ELSE IF (ALFM.LT.ALF1) THEN
          ALFM = ALF1+DALF
        ENDIF
        IF (ALFM.GT.ALF1*FAC) ALFM = ALF1*FAC
        IF (DF1.LT.ZERO) THEN
          ALF0 = ALF1
          F0 = F1
          DF0 = DF1
        ENDIF
        ALF1 = ALFM
        GOTO 10
      ENDIF
C
C***********************************************************************
C                                DETERMINE MINIMUM WITHIN GIVEN INTERVAL
C***********************************************************************
C
   40 CONTINUE
      ILINE = ILINE+1
C
C---------------------------------- CHECK IF MORE THAN MAXFUN FUNC-CALLS
C
      IF (ILINE.GT.MAXFUN) THEN
        IF (IPRINT.GT.0) WRITE (IOUT,2000) MAXFUN
        GOTO 900
      ENDIF
C
      DALF = ALF1-ALF0
      GDIF = DF1-DF0
      IF (DALF.LT.AC2) GOTO 900
C***********************************************************************
C                                           GOLDEN SECTION INTERPOLATION
C***********************************************************************
      IF (ABS(DF0).LT.ABS(DF1)) THEN
        ALFM = ALF0 + GOLD2*DALF
      ELSE
        ALFM = ALF0 + GOLD1*DALF
      ENDIF
C
C***********************************************************************
C                                                     DIRECTION GRADIENT
C***********************************************************************
      DALF1 = ALFM - ALFL
      DO 45 J=1,M
        IF (ACTIVE(J) .AND. S(J).NE.ZERO) U(J) = U(J) + DALF1*S(J)
        IF (U(J) .LT. ZERO ) U(J)=ZERO
   45 CONTINUE
      ALFL = ALFM
C
      CALL MMADUA (1,X,X1,XL,XU,ASL,ASU,U,ALF,BET,IOUT,P0,R0,P,R,W,DW,
     *             Z,PEN,ACTIVE,N,ME,M,MMAX,IRET,FARS)
      W = -W
      FM = W
C
      GRADS = ZERO
      DO 47 J=1,M
        IF (ACTIVE(J) .AND. S(J).NE.ZERO) GRADS = GRADS - S(J)*DW(J)
   47 CONTINUE
      DFM = GRADS
C
      IF (ABS(DFM).LT.ACL .OR. ABS(DALF1).LT.AC2) GOTO 900
C
      IF (DFM.LT.ZERO) THEN
        ALF0 = ALFM
        F0 = FM
        DF0 = DFM
      ELSE IF (DFM.GT.ZERO) THEN
        ALF1 = ALFM
        F1 = FM
        DF1 = DFM
      ENDIF
C
      DALF = ALF1-ALF0
      GDIF = DF1-DF0
      IF (DALF.LT.AC2) GOTO 900
C***********************************************************************
C                                                 BILINEAR APPROXIMATION
C***********************************************************************
      IF (IAPPR.EQ.1) THEN
        ALFM = (F0-F1-ALF0*DF0+ALF1*DF1)/GDIF
C       IF ((ALFM-ALF0).LT.(ONFO*DALF)) ALFM = ALF0+ONFO*DALF
C       IF ((ALF1-ALFM).LT.(ONFO*DALF)) ALFM = ALF1-ONFO*DALF
C***********************************************************************
C                                                QUADRATIC INTERPOLATION
C***********************************************************************
      ELSE IF (IAPPR.EQ.2) THEN
        A1 = (-DF0*DALF-F0+F1)/(DALF*DALF)
        IF (ABS(A1).GT.AC2) THEN
          ALFM = ALF0 - HA*DF0/A1
        ELSE
          ALFM = HA*(ALF0+ALF1)
        ENDIF
        IF (ALFM.GE.ALF1) THEN
          IF (ABS(GDIF).LE.BAR) THEN
            ALFM = HA*(ALF0+ALF1)
          ELSE
            ALFM = ALF0 - DF0/(DF1-DF0)*DALF
          ENDIF
        ENDIF
        IF ((ALFM-ALF0).LT.(ONFO*DALF)) ALFM = ALF0+ONFO*DALF
        IF ((ALF1-ALFM).LT.(ONFO*DALF)) ALFM = ALF1-ONFO*DALF
C***********************************************************************
C                                                    CUBIC APPROXIMATION
C***********************************************************************
      ELSE IF (IAPPR.EQ.3) THEN
        DAL2 = DALF*DALF
        A1 = DF0/DALF + DF1/DALF +TWO*(F0-F1)/DAL2
        B1 = -HA*(DF0 - DF1 - THREE*A1*(ALF0*ALF0-ALF1*ALF1)/DALF )
        C1 = -(DF1*ALF0 - DF0*ALF1 - THREE*A1*ALF0*ALF1)
        ARG = B1*B1-THREE*A1*C1
        IF (ARG.GT.UBAR .OR. A1.GT.UBAR) THEN
          ALFM = HA*(ALF0+ALF1)
        ELSE IF (ARG.GE.ZERO .AND. ABS(A1).GT.ACL) THEN
          ALFM = (SQRT(ARG)-B1)/(THREE*A1)
        ELSE
          A1 = (-DF0*DALF-F0+F1)/(DALF*DALF)
          IF (ABS(A1).GT.AC2) THEN
            ALFM = ALF0 - HA*DF0/A1
          ELSE
            ALFM = HA*(ALF0+ALF1)
          ENDIF
        ENDIF
      ENDIF
C***********************************************************************
C                                                     DIRECTION GRADIENT
C***********************************************************************
C
      IF (ALFM.LT.ALF0 .OR. ALFM.GT.ALF1) ALFM = HA*(ALF0+ALF1)
C
      DALF1 = ALFM - ALFL
      DO 50 J=1,M
        IF (ACTIVE(J) .AND. S(J).NE.ZERO) U(J) = U(J) + DALF1*S(J)
        IF (U(J) .LT. ZERO ) U(J)=ZERO
   50 CONTINUE
      ALFL = ALFM
C
      CALL MMADUA (1,X,X1,XL,XU,ASL,ASU,U,ALF,BET,IOUT,P0,R0,P,R,W,DW,
     *             Z,PEN,ACTIVE,N,ME,M,MMAX,IRET,FARS)
      W = -W
      FM = W
C
      GRADS = ZERO
      DO 60 J=1,M
        IF (ACTIVE(J) .AND. S(J).NE.ZERO) GRADS = GRADS - S(J)*DW(J)
   60 CONTINUE
      DFM = GRADS
C
      IF (ABS(DFM).LT.ACL .OR. ABS(DALF1).LT.AC2) GOTO 900
C
      IF (DFM.LT.ZERO) THEN
        ALF0 = ALFM
        F0 = FM
        DF0 = DFM
      ELSE IF (DFM.GT.ZERO) THEN
        ALF1 = ALFM
        F1 = FM
        DF1 = DFM
      ENDIF
C
      GOTO 40
C
C***********************************************************************
C                                             REGULAR END OF LINE SEARCH
C***********************************************************************
C
  900 CONTINUE
      IF (IPRINT.GT.1) WRITE (IOUT,2010) ILINE
C
      RETURN
C-----------------------------------------------------------------------
 2000 FORMAT (4X,'*** WARNING ***  MORE THAN MAXFUN =',I3,
     *                                     ' FUNC-CALLS IN LINE SEARCH')
 2010 FORMAT (7X,I4,' FUNC-CALLS IN LINE SEARCH')
      END
