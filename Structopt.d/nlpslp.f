      SUBROUTINE NLPSLP (IALG,
     &                   VAR,RESL,RESU,FF,G,DF,DG,
     &                   DMIN,DMAX,SCX,SCG,Y,IWA,ACTIVE,
     &                   NUMVAR,M,ME,
     &                   IMAX,MODE,INCON,ICONV,ICYCLE,
     &                   IADAPT,ILISE,IPRINT,FNAME,IFNSIZE,
     &                   QMOVE,FMOVE,FBACK,TOL,OBJTOL,AMIJO,
     &                   IBL)
 
C-----------------------------------------------------------------------
C     SEQUENTIAL LINEAR PROGRAMMING (SLP): MAIN ROUTINE
C     ******************************************************************
C     *     INTEGER PARAMETER                                          *
C     *     IMAX   --->  MAXIMUM NUMBER OF ITERATIONS FOR SLP          *
C     *     INCON  --->  = 1: DOUBLING OF MOVE-LIMITS IN CASE OF       *
C     *                       INCONSISTANCE                            *
C     *                  = 2: ONE ARTIFICAL VARIABLE WILL BE ADDED     *
C     *                       IN CASE OF INCONSISTANCE                 *
C     *     ICONV  --->  = 0: OBJTOL WILL BE CHECKED                   *
C     *                  = 1:    TOL WILL BE CHECKED (KUHN-TUCKER)     *
C     *     ICYCLE --->  = 0: NEW MOVE LIMITS IN CASE OF CYCLING       *
C     *                       DMIN(I) = FMOVE*DMIN(I)                  *
C     *                    1: ONE STEP BACK IN CASE OF CYCLING WITH    *
C     *                       DX(I) = -FBACK*DX(I)                     *
C     *                    2: A LINESEARCH WITH DIRECTION DX WILL BE   *
C     *                       STARTED IN CASE OF CYCLING               *
C     *     IADAPT --->  = 0: NO AUTOMATICAL ADAPTION OF MOVE-LIMITS   *
C     *                  = 1: AUTOMATICAL ADAPTION OF MOVE-LIMITS      *
C     *                       ALF(K+1) = ALF(K)/(ALF(K) + 1); KNEPPE   *
C     *                  = 2: AUTOMATICAL ADAPTION OF MOVE-LIMITS      *
C     *                       ALF(K+1) = ALF(K)/(ALF(K)**2 + 1); BAIER *
C     *                  = 3: AUTOMATICAL ADAPTION OF MOVE-LIMITS      *
C     *                       ALF(K+1) = ALF(K)/2 ; REINSCHMIDT        *
C     *     ILISE  --->  = 2: QUADRATIC LINE SEARCH                    *
C     *                  = 3: CUBIC LINE SEARCH                        *
C     *                                                                *
C     *     REAL PARAMETER                                             *
C     *     OBJTOL --->  ACCURACY OF CHANGE FOR OBJ.FUNCTION           *
C     *     TOL    --->  ACCURACY OF SCHITTKOWSKI                      *
C     *     QMOVE  --->  INITIAL MOVE-LIMIT = (XU -XL)*QMOVE           *
C     *     FMOVE  --->  ADAPTION FACTOR FOR MOVE-LIMITS               *
C     *     FBACK  --->  ICYCLE = 1,FACTOR FOR A STEP BACKWARDS        *
C     *                  IN CASE OF CYCLING                            *
C     *     AMIJO  --->  FACTOR OF AMIJO CHECK IN LINE SEARCH          *
C     *                                                                *
C     *     CHARACTER PARAMETER                                        *
C     *     PRINT  --->  PRINT MODE                                    *
C     *                  1  --->  LAST ITERATION                       *
C     *                  2  --->  ALL ITERATIONS                       *
C     *                  3  --->  ONLY OBJECTIVE FUNCTION              *
C     ******************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)

      LOGICAL*4 ACTIVE

      CHARACTER FNAME*(*),CDUM*1
C
C--------------------------------------------------- STRATEGY PARAMETERS
C
      MMAX   = MAX(1,M)
      MNN    = M +2*NUMVAR
      LWA    = ((MNN+3)*(NUMVAR+2) 
     &          + 5*NUMVAR + MNN + 1)*2 + (MNN + NUMVAR + 5)*1
C
C-------------------------------- OPEN OUTPUT FILE AND FORWARD TO THE END
C
      IOUT= 45
      OPEN (IOUT,FILE=FNAME(1:IFNSIZE))
    1 READ (IOUT,1000,END=2) CDUM
      GOTO 1
    2 CONTINUE
C      
 1000 FORMAT(A1) 
C
      WRITE(IOUT,*) ' SEQUENTIAL LINEAR PROGRAMMING (SLP) '
C
C-------------------------------------------------- PERFORM OPTIMIZATION
C
      CALL SLP1 (IALG,M,MMAX,ME,NUMVAR,LWA,VAR,RESL,RESU,FF,G,DF,DG,
     &           DMIN,DMAX,SCX,SCG,Y,IWA,ACTIVE,
     &           IMAX,MODE,INCON,ICONV,ICYCLE,IADAPT,ILISE, IPRINT,IOUT,
     &           QMOVE,FMOVE,FBACK,TOL,OBJTOL,AMIJO,
     &           IBL,ITRACT,NROPT)
C
C--------------------------------------------------------- PRINT RESULTS
c
c      CALL SLPPRT (IBL,NROPT,TSLP,N,F,L(NVAR),L(NVSCL),IPRINT)
C
      RETURN
      END
C=======================================================================

      SUBROUTINE SLP1 (IALG, M, MMAX, ME, N, LWA, X,XL,XU,FF,G,DF,DG,
     &           DMIN,DMAX,SCX,SCG,Y,IWA,ACTIVE,
     &           IMAX,MODE,INCON,ICONV,ICYCLE,IADAPT,ILISE, IPRINT,IOUT,
     &           QMOVE,FMOVE,FBACK,TOL,OBJTOL,AMIJO,
     &           IBL,ITRACT,NROPT)

C-----------------------------------------------------------------------
C     SEQUENTIAL LINEAR PROGRAMMING (SLP): DEFINITION OF WORKING ARRAY
C-----------------------------------------------------------------------
C     M      ... NUMBER OF CONSTRAINTS
C     MMAX   ... MAX (M,1)
C     MEQ    ... NUMBER OF EQUALITY CONSTRAINTS
C     N      ... NUMBER OF VARIABLES
C     IMAX   ... MAX. NUMBER OF ITERATIONS
C     TOL    ... TOLERANCE OF CONVERGENCE
C     OBJTOL ... TOLERANCE OF OBJECTIVE CONVERGENCE
C     LWA    ... LENGTH OF WORKING ARRAY
C     X      ... OPTIMIZATION VARIABLES
C     DMIN   ... LOWER MOVE LIMITS
C     DMAX   ... UPPER MOVE LIMITS
C     G      ... CONSTRAINTS
C     XL     ... LOWER BOUNDS OF VARIABLES
C     XU     ... UPPER BOUNDS
C     DF     ... OBJECTIVE FUNCTION GRADIENT
C     DG     ... GRADIENTS OF CONSTRAINTS
C     ACTIVE ... INDICATOR OF ACTIVE CONSTRAINTS
C     SCX    ... VARIABLE SCALING FACTORS
C     SCG    ... CONSTRAINTS SCALING FACTORS
C     Y      ... LAGRANGE MULTIPLICATORS
C     IWA    ... WORKING ARRAY
C     IOUT   ... OUTPUT DEVICE
C     IBL    ... ERROR CODE
C     NROPT  ... NUMBER OF ITERAIONS WITHIN >OPTSLP<
C     ITRACT ... TOTAL NUMBER OF ITERATIONS
C     F      ... VALUE OF OBJECTIVE FUNCTION
C     INCON  ... INCONSISTENCY FLAG
C     ICYCLE ... CYCLING FLAG
C     ICONV  ... CONVERGENCE FLAG
C     MODE   ... RESTART MODE
C     IOPDOC ... DOCUMENTAION FLAG
C     QMOVE  ... INITIAL MOVE LIMIT FACTOR
C     FMOVE  ... MOVE LIMIT UPDATE FACTOR
C     IPRINT ... PRINT FLAG
C     FBACK  ... BACK SETTING FACTOR (IN CASE OF ICYCLE = 1)
C     IADAPT ... TYPE OF MOVE LIMIT ADAPTATION
C     ILISE  ... TYPE OF LINE SEARCH (QUADRATIC OR CUBIC)
C     AMIJO  ... FACTOR FOR AMIJO LINE SEARCH CHECK
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)

      LOGICAL*4 ACTIVE

      DIMENSION X(N),DMIN(N),DMAX(N),G(MMAX),XL(N),XU(N),IWA(LWA),DF(N)
      DIMENSION SCX(N),SCG(MMAX), ACTIVE(MMAX),Y(M+2*N+2), DG(MMAX,N)
C
C ---------------------------------- CALCULATION OF ADRESSES FOR WA(LWA)
C
      N1  = N + 1
      MNN = M + 2*N
C
      NS     = 1
      NWORK  = NS + ((MNN+3)*(N+2))*2
      NNR    = NWORK + N*2
      NNC    = NNR+ (MNN+3)*1
      NDX    = NNC + (N+2)*1
      NDMIN1 = NDX + (N+1)*2
      NDMAX1 = NDMIN1 + N*2
      NDXOLD = NDMAX1 + N*2
      NRPEN  = NDXOLD + N*2
      NSUM   = NRPEN  + MNN*2 - 1
C
      IF(NSUM.GT.LWA) THEN
        IF(IPRINT.GT.0) WRITE(IOUT,1000) LWA,NSUM
        GOTO 900
      ENDIF
C
C     DUMMY SCALING OF VARIABLES AND CONSTRAINTS
C     SINCE THIS IS HANDLED BEFORE
C
      DO I=1,N
       SCX(I)=1.0D0
      ENDDO
C      
      DO I=1,MMAX
       SCG(I)=1.0D0
      ENDDO
C
      DO I=1,m+2*n+2
        Y(I)=0.0D0
      ENDDO
C
      CALL SLP2 (IALG, M, MMAX, ME, N, MNN,N1,FF,
     &           X(1)       , G(1)       , XL(1)      , XU(1)     , 
     &           DF(1)      , DG(1,1)    , ACTIVE(1)  , SCX(1)    , 
     &           SCG(1)     , IWA(NS)    , Y(1)       , IWA(NNR)  ,
     &           IWA(NNC)   , IWA(NDX)   , DMIN       , DMAX      , 
     &           IWA(NDMIN1), IWA(NDMAX1), IWA(NDXOLD), IWA(NWORK), 
     &           IWA(NRPEN),
     &           IMAX,MODE,INCON,ICONV,ICYCLE,IADAPT,ILISE, IPRINT,IOUT,
     &           QMOVE,FMOVE,FBACK,TOL,OBJTOL,AMIJO,
     &           IBL,ITRACT,NROPT)

C
  900 CONTINUE
      RETURN
C
 1000 FORMAT(' *** WORKING ARRAY TOO SMALL  :',I9/
     *       '     SHOULD BE AT LEAST       :',I9)
      END

C=======================================================================
      SUBROUTINE SLP2 (IALG,M,MMAX,ME,N,MNN,N1,FF,
     &           X,G,XL,XU,DF,DG,ACTIVE,SCX,SCG,S,
     &           Y,NR,NC,DX,DMIN,DMAX,DMIN1,DMAX1,DXOLD,
     &           WORK,RPEN,
     &           IMAX,MODE,INCON,ICONV,ICYCLE,IADAPT,ILISE, IPRINT,IOUT,
     &           QMOVE,FMOVE,FBACK,TOL,OBJTOL,AMIJO,
     &           IBL,ITRACT,NROPT)

C-----------------------------------------------------------------------
C     SEQUENTIAL LINEAR PROGRAMMING (SLP): MAIN ITERATION ROUTINE
C-----------------------------------------------------------------------
C     M      ... NUMBER OF CONSTRAINTS
C     MMAX   ... MAX (M,1)
C     ME     ... NUMBER OF EQUALITY CONSTRAINTS
C     N      ... NUMBER OF VARIABLES
C     MNN    ... MNN = M + N + N
C     N1     ... N1  = N + 1
C     IMAX   ... MAX. NU,BER OF ITERATIONS
C     TOL    ... TOLERANCE OF CONVERGENCE
C     OBJTOL ... TOLERANCE OF OBJECTIVE CONVERGENCE
C     IOUT   ... OUTPUT DEVICE
C     MODE   ... RESTART MODE
C     X      ... OPTIMIZATION VARIABLES
C     G      ... CONSTRAINTS
C     XL     ... LOWER BOUNDS OF VARIABLES
C     XU     ... UPPER BOUNDS
C     DF     ... OBJECTIVE FUNCTION GRADIENT
C     DG     ... GRADIENTS OF CONSTRAINTS
C     ACTIVE ... INDICATOR OF ACTIVE CONSTRAINTS
C     SCX    ... VARIABLE SCALING FACTORS
C     SCG    ... CONSTRAINTS SCALING FACTORS
C     S      ... SIMPLEX TABLEAU
C     Y      ... LAGRANGE MULTIPLICATORS
C     NR     ... SIMPLEX ROW INDEX
C     NC     ... SIMPLEX COLUMN INDEX
C     DX     ... VARIABLE INCREMENT
C     DMIN   ... LOWER MOVE LIMITS
C     DMAX   ... UPPER MOVE LIMITS
C     DMIN1  ... LOWER MOVE LIMITS REFERENCE VALUE
C     DMAX1  ... UPPER MOVE LIMITS REFERENCE VALUE
C     DXOLD  ... VARIABLE INCREMENT OF PRIOR ITERATION STEP
C     WORK   ... WORKING ARRAY
C     RPEN   ... PENALTY FACTORS
C     IBL    ... ERROR CODE
C     ITER   ... NUMBER OF ITERATIONS IN >SLP2<
C     ITRACT ... TOTAL NUMBER OF ITERATIONS
C     FF     ... VALUE OF OBJECTIVE FUNCTION
C     INCON  ... INCONSISTENCY FLAG
C     ICYCLE ... CYCLING FLAG
C     ICONV  ... CONVERGENCE FLAG
C     IOPDOC ... DOCUMENTAION FLAG
C     QMOVE  ... INITIAL MOVE LIMIT FACTOR
C     FMOVE  ... MOVE LIMIT UPDATE FACTOR
C     IPRINT ... PRINT FLAG
C     IGRAD  ... TYPE OF GRADIENT EVALUATION (IMPLICIT,ANALYTICAL,SEMI)
C     FBACK  ... BACK SETTING FACTOR (IN CASE OF ICYCLE = 1)
C     IADAPT ... TYPE OF MOVE LIMIT ADAPTATION
C     ILISE  ... TYPE OF LINE SEARCH (QUADRATIC OR CUBIC)
C     AMIJO  ... FACTOR FOR AMIJO LINE SEARCH CHECK
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(ZERO=0.D0,ONE=1.0D0,TWO=2.0D0,EPS=1.D-6)
      LOGICAL *4 ACTIVE
      CHARACTER TEXT*80
      INTEGER simptyp
      DIMENSION G(MMAX),DF(N),DG(MMAX,N),SCX(N),ACTIVE(MMAX),DX(N+1)
      DIMENSION X(N), S(MNN+1,N+1), Y(MNN+2), NR(MNN+3), NC(N+2)
      DIMENSION DMIN(N), DMAX(N), DMIN1(N),DMAX1(N),DXOLD(N)
      DIMENSION XL(N), XU(N), WORK(N), RPEN(MNN)
c
      simptyp = 2
c
C
C-----------------------------------------------------------------------
C
      MN = M+N
      MS = MNN+1
      NS = N+1
      NF = N
      NFUNC = 0
      NGRAD = 0
      NSEN  = 0
      MAXCYC= 0
      NLIS  = 0
C
C-----------------------------------------------------------------------
      ITER   = 0
      ITR    = 0
      EXFAC  = ONE
      IFINAL = 0
      NCYCLE = 0
C-----------------------------------------------------------------------
C
      IF(IPRINT.GT.0) THEN
        WRITE(IOUT,9000) MODE,IMAX,IPRINT,ICONV,INCON,ICYCLE,ILISE,
     *                   IGRAD,IADAPT,QMOVE,FMOVE,FBACK,AMIJO,TOL,OBJTOL
      ENDIF
C                                              =========================
C--------------------------------------------- CHECK OF INITIAL VARIABLE
C                                              =========================
      DO 10 I=1,N
         X(I) = MAX(XL(I),X(I))
         X(I) = MIN(XU(I),X(I))
C
C                =======================================================
C--------------- RELATIVE MOVELIMITS  DMIN <=   DX = (X - X0)    <= DMAX
C                =======================================================
C
        DMIN1(I) = -(XU(I)-XL(I))*QMOVE
        DMAX1(I) =  (XU(I)-XL(I))*QMOVE
	DX(I)    =  ZERO
   10 CONTINUE
C
C                                                             ==========
C------------------------------------------------------------ ACTIVE SET
C                                                             ==========
      DO 20 I=1,M
       ACTIVE(I) = .TRUE.
       RPEN(I)   = ZERO
   20 CONTINUE
C                                                      =================
C----------------------------------------------------- SCALE MULTIPLIERS
C                                                      =================
      DO 30 I=1,N
        RPEN(M+I)  = ZERO
        RPEN(MN+I) = ZERO
        Y(M+I)     = Y(M+I) *SCX(I)
        Y(MN+I)    = Y(MN+I)*SCX(I)
   30 CONTINUE
C
C                                                      =================
C----------------------------------------------------- ASK FOR GRADIENTS
C                                                      =================
C
      NSEN = M

      CALL FUNC (IALG,1)
      FO=FF
      
      write(*,*) 'SLP:  FF:',FF,'  FO:',FO
      
      CALL GRAD (IALG,ACTIVE)

      NGRAD = NGRAD + 1
      FOPREV = FO
C
C***********************************************************************
C                                         MAIN ITERATION LOOP,LABEL 100
C***********************************************************************
  100 CONTINUE
C
C
      CALL SLPPSI (DF,DG,X,Y,G,N,M,MMAX,ME,TOL,FO,ISCHIT,IOUT,FLAG,SSUM,
     *             SRES,DLAN)
C
C                                                  =====================
C------------------------------------------------- CALUCLATE MOVE-LIMITS
C                                                  =====================
      CALL SLPMOV (N, DMIN, DMAX,X,XL,XU,DX,DXOLD,ITER,
     *                DMIN1,DMAX1,FMOVE,IADAPT,MAXCYC)
C
C                                                  =====================
C ------------------------------------------------ WRITE ACTUAL SOLUTION
C                                                  =====================
      IF (IPRINT.GT.1) THEN
        DO 110 I=1,N
          Y(M+I)  = Y(M+I) /SCX(I)
          Y(MN+I) = Y(MN+I)/SCX(I)
  110   CONTINUE
C
        CALL SLPERR (X,Y,G,SCX,XL,XU,N,M,MMAX,MNN,FO,ITER,FLAG,
     *              IOUT,SSUM,SRES,DLAN,IFINAL,NFUNC,NGRAD,IPRINT,
     *              NSEN,MAXCYC,NLIS)
C
        DO 120 I=1,N
          Y(M+I)  = Y(M+I) *SCX(I)
          Y(MN+I) = Y(MN+I)*SCX(I)
  120   CONTINUE
      ENDIF
C
C                                                  =====================
C------------------------------------------------- CHECK FOR CONVERGENCE
C                                                  =====================
C
      Q = (FOPREV - FO)
      IF (ICONV.EQ.1 .OR. ITER.EQ.0) THEN
        IF (ISCHIT.EQ.0) GOTO 900
      ELSE
        IF (ABS(Q).LE.OBJTOL) GOTO 900
      ENDIF
C
      ITER   = ITER + 1
      ITRACT = ITRACT + 1
      ITR    = 0
      EXFAC  = ONE
C                                     ==================================
C------------------------------------ CHECK FOR MAXIMUM ITERATION NUMBER
C                                     ==================================
C     IF(ITER.GT.10) ICYCLE = 0
      IF(ITER.GT.IMAX) THEN
         IBL = 6
         ITER = ITER - 1
         IF(IPRINT.EQ.0) GOTO 900
         WRITE(IOUT,9200)
         GOTO 900
      ELSE
        WRITE (IOUT,'(A,I3)') '>>> OPTSLP: ITERATION NO. ',ITER
      ENDIF
C
      FOPREV = FO
      CALL DCOPY(N,DX,1,DXOLD,1)
C                                    ===================================
C-------------------------------------- FORMING OF THE LINEAR SUBPROBLEM
C                                   ====================================
C------------------------------------& SOLUTION OF THE LINEAR SUBPROBLEM
C                                   ====================================
c
      if (simptyp.eq.1) then
c
        CALL SLPST1 (ITER,M,MMAX,ME,N,MNN,FO,G,DF,DG,X,SCX,ACTIVE,
     1               DMIN,DMAX,S,IOUT)
C
        CALL SLPSIM (S,DX,Y,F,NR,NC,MNN,ME,N,NF,IBL,EPS)
c
      else
c
        CALL NEWSIMP(M,ME,N,G,DF,DG,DX,DMIN,DMAX,IBL)
c
      endif
C
C                                                   ====================
C-------------------------------------------------- CHECK FOR ERROR FLAG
C                                                   ====================
  130 CONTINUE
      IF (IBL.EQ.0) GOTO 170
C
      IF(IPRINT.GT.0) THEN
        IF(IBL.EQ.1) WRITE(IOUT,5100)
        IF(IBL.EQ.2) WRITE(IOUT,5200)
        IF(IBL.EQ.3) WRITE(IOUT,5300)
        IF(IBL.EQ.4) WRITE(IOUT,5400)
        IF(IBL.EQ.5) WRITE(IOUT,5500)
      ENDIF
      IF (IBL.NE.3) GOTO 900
C
C***********************************************************************
C                                       CASE OF INCONSISTENT CONSTRAINTS
C***********************************************************************
C
      IF (INCON.EQ.1) THEN
C------------------------------------------------ EXTEND FEASIBLE REGION
        ITR = ITR + 1
        IF (ITR.GT.5) THEN
          WRITE(IOUT,5600)
          GOTO 900
        ENDIF
        EXFAC = EXFAC/TWO
        DO 140 I=1,N
          DMIN(I) = TWO*DMIN(I)
          DMAX(I) = TWO*DMAX(I)
  140   CONTINUE
C
C                                    ===================================
C-------------------------------------- FORMING OF THE LINEAR SUBPROBLEM
C                                   ====================================
C------------------------------------& SOLUTION OF THE LINEAR SUBPROBLEM
C                                   ====================================
c
        if (simptyp.eq.1) then
C
          CALL SLPST1 (ITER,M,MMAX,ME,N,MNN,FO,G,DF,DG,X,SCX,ACTIVE,
     1                 DMIN,DMAX,S,IOUT)
c
          CALL SLPSIM (S,DX,Y,F,NR,NC,MNN,ME,N,NF,IBL,EPS)
c
        else
c
c
          CALL NEWSIMP(M,ME,N,G,DF,DG,DX,DMIN,DMAX,IBL)
c
        endif
C
        IF (IBL.GT.0) GOTO 130
C
        DO 150 I=1,N
          DMIN(I) = EXFAC*DMIN(I)
          DMAX(I) = EXFAC*DMAX(I)
  150   CONTINUE
C
      ELSE
C
C                                    ===================================
C-------------------------------------- FORMING OF THE LINEAR SUBPROBLEM
C                                   ====================================
C------------------------------------& SOLUTION OF THE LINEAR SUBPROBLEM
C                                   ====================================
c
        if (simptyp.eq.1) then
c
          CALL SLPST2 (ITER,M,MMAX,ME,N,MNN,FO,G,DF,DG,X,SCX,ACTIVE,
     1                 DMIN,DMAX,S)
c
          CALL SLPSIM (S,DX,Y,F,NR,NC,MNN+2,ME,N+1,NF+1,IBL,EPS)
c
        else
c
          stop 'new simplex method not implemented for inconsistency'

c          CALL NEWSIMP(ITER,M,MMAX,ME,N,MNN,FO,G,DF,DG,X,SCX,ACTIVE,
c     1                 DMIN,DMAX,S,IOUT)
c
        endif
C
C--------------------------------------------------- REORDER MULTIPLIERS
        DO 160 I=1,N
          Y(MN+I) = Y(MN+I+1)
  160   CONTINUE
C
        IF (IBL.EQ.3) THEN
          WRITE(IOUT,5600)
          GOTO 900
        ENDIF
      ENDIF
C
C***********************************************************************
C                                                   CONTINUE CALCULATION
C***********************************************************************
  170 CONTINUE
C                                                  =====================
C------------------------------------------------- CHECK FOR MOVE LIMITS
C                                                  =====================
      DO 180 I=1,N
        DX(I) = MAX(DMIN(I),DX(I))
        DX(I) = MIN(DMAX(I),DX(I))
        IF (X(I)+DX(I).GT.XL(I)) Y(M+I)   = ZERO
        IF (X(I)+DX(I).LT.XU(I)) Y(M+N+I) = ZERO
  180 CONTINUE
C
C                                                      =================
C----------------------------------------------------- CHECK FOR CYCLING
C                                                      =================
C
      ILIS = 1
      IF (ICYCLE.GT.0) THEN
        NCYCLE = 0
C
        DO 190 I=1,N
          CHANGE = DX(I)*DXOLD(I)
          IF (CHANGE.LT.ZERO) THEN
            NCYCLE = NCYCLE+1
            IF (ICYCLE.EQ.1) DX(I) = (ONE-FBACK)*DX(I)
          ENDIF
  190   CONTINUE
C
        IF (ICYCLE.EQ.2 .AND. NCYCLE.GT.0) THEN
C                                             ==========================
C-------------------------------------------- LINESEARCH IN DIRECTION DX
C                                             ==========================
          NLIS = NLIS + 1
          CALL SLPLIS (DX,DXOLD,X,DMIN,DMAX, DF,N,FO,FOPREV,ILIS,IOUT,
     *                  XL,XU,Y,RPEN,DG,G,ACTIVE,SCX,M,MMAX,ME,MNN,
     *                  ITER,NFUNC,NGRAD,NSEN,IGRAD,ILISE,AMIJO,
     *                  IALG,FF)
C
        ENDIF
      ENDIF
C
C***********************************************************************
C                                                  END OF ITERATION STEP
C***********************************************************************
C
C
      IF (ILIS.EQ.1) THEN
C
        DO 200 I=1,N
          X(I) = X(I)+DX(I)
  200   CONTINUE
C
        CALL FUNC (IALG,ITER)
        FO=FF
        NFUNC = NFUNC + 1
      ENDIF
C                                                             ==========
C------------------------------------------------------------ ACTIVE SET
C                                                             ==========
      IF (ILIS.GE.0) THEN
        NSEN = NSEN + M
        DO 210 I=ME+1,M
          ACTIVE(I) = .TRUE.
          IF(Y(I).EQ.ZERO.AND.G(I).GT.EPS) THEN
            ACTIVE(I)= .FALSE.
            NSEN = NSEN - 1
          ENDIF
  210   CONTINUE
C
        CALL GRAD(IALG,ACTIVE)
        NGRAD = NGRAD + 1
      ENDIF
C
      GOTO 100
C
C***********************************************************************
C                                                           EXIT OF SLP2
C***********************************************************************
C
  900 CONTINUE
C                                                      =================
C----------------------------------------------------- SCALE MULTIPLIERS
C                                                      =================
      DO 910 I=1,N
        Y(M+I)  = Y(M+I) /SCX(I)
        Y(MN+I) = Y(MN+I)/SCX(I)
  910 CONTINUE
C
      IFINAL = 1
      IF (IPRINT.GT.0)
     *       CALL SLPERR (X,Y,G,SCX,XL,XU,N,M,MMAX,MNN,FO,ITER,FLAG,
     *                    IOUT,SSUM,SRES,DLAN,IFINAL,NFUNC,NGRAD,IPRINT,
     *                    NSEN,MAXCYC,NLIS)
C
C-----------------------------------------------------------------------
      RETURN
C
 9000 FORMAT(4X,
     * '-------------------------------------------------------------',
     * /4X,'START OF SLP',/4X,
     * '-------------------------------------------------------------',
     *      //5X,'PARAMETERS',
     *      //8X,'MODE    = ',I3,
     *       /8X,'IMAX    = ',I3,
     *       /8X,'IPRINT  = ',I3,
     *       /8X,'ICONV   = ',I3,
     *       /8X,'INCON   = ',I3,
     *       /8X,'ICYCLE  = ',I3,
     *       /8X,'ILISE   = ',I3,
     *       /8X,'GRAD    = ',I3,
     *       /8X,'IADAPT  = ',I3,
     *       /8X,'QMOVE   = ',E13.7,
     *       /8X,'FMOVE   = ',E13.7,
     *       /8X,'FBACK   = ',E13.7,
     *       /8X,'AMIJO   = ',E13.7,
     *       /8X,'TOL     = ',E13.7,
     *       /8X,'OBJTOL  = ',E13.7)
 9200 FORMAT(/5X,'*** MAXIMUM NUMBER OF ITERATIONS REACHED')
C
 5100 FORMAT(5X,'*** WRONG DIMENSIONS')
 5200 FORMAT(5X,'*** LINEAR DEPENDENT CONSTRAINTS')
 5300 FORMAT(5X,'*** INCONSISTENT CONSTRAINTS, SLP CONTINUED')
 5400 FORMAT(5X,'*** UNBOUNDED OBJECTIVE FUNCTION')
 5500 FORMAT(5X,'*** UNDERFLOW IN PIVOT ELEMENT')
 5600 FORMAT(5X,'*** STILL INCONSISTENT CONSTRAINTS, SLP STOPPED',
     *      /5X,'    NO FEASIBLE SOLUTION IN SUBPROBLEM !')
      END
C=======================================================================
      SUBROUTINE SLPPSI (DF,DG,X,Y,G,N,M,MMAX,ME,TOL,FO,ISCHIT,IOUT,
     *                   FLAG,SSUM,SRES,DLAN)
C-----------------------------------------------------------------------
C     SEQUENTIAL LINEAR PROGRAMMING (SLP): CONVERGENCE CRITERION
C-----------------------------------------------------------------------
C     DF     ... OBJECTIVE FUNCTION GRADIENT
C     DG     ... GRADIENTS OF CONTRAINTS
C     X      ... VARIABLES
C     Y      ... LAGRANGE MULTIPLICATORS
C     G      ... CONSTRAINTS
C     N      ... NUMBER OF VARIABLES
C     M      ... NUMBER OF CONSTRAINTS
C     MMAX   ... MAX (M,1)
C     ME     ... NUMBER OF EQUALITY CONSTRAINTS
C     TOL    ... TOLERANCE OF CONVERGENCE
C     FO     ... OBJECTIVE FUNCTION
C     ISCHIT ... TYPE OF CONVERGENCE CHECK
C     IOUT   ... OUTPUT DEVICE
C     FLAG   ... VALUE OF LAGRANGE FUNCTION
C     SSUM   ... SUM OF WEIGHTED RESIDUUM
C     SRES   ... SUM OF VIOLATED CONSTRAINTS
C     DLAN   ... NORM OF LAGRANGE FUNCTION GRADIENT
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ZERO=0.D0,ONE=1.D0)
      DIMENSION DF(N),DG(MMAX,N),X(N),Y(MMAX),G(MMAX)
C
C-------------------------------------------- CONVERGENCE-CRITERION NO 1
      SSUM = ZERO
      IF (M.GT.0) THEN
        DO 100 I=1,M
          SSUM = SSUM + ABS(Y(I)*G(I))
  100   CONTINUE
      ENDIF
C-------------------------------------------- CONVERGENCE-CRITERION NO 2
      SRES = ZERO
      IF (ME.GT.0) THEN
        DO 200 I=1,ME
          SRES = SRES + ABS(G(I))
  200   CONTINUE
      ENDIF
      IF (M.GT.0) THEN
        DO 210 I=ME+1,M
          GMIN = MIN(ZERO,G(I))
          SRES = SRES + ABS(GMIN)
  210   CONTINUE
      ENDIF
C
C-------------------------------------------- CONVERGENCE-CRITERION NO 3
      DLAN = ZERO
      DO 300 I=1,N
        SUM1 = ZERO
        IF (M.GT.0) THEN
          DO 310 J=1,M
            SUM1 = SUM1 + Y(J)*DG(J,I)
  310     CONTINUE
        ENDIF
        SUM1   = SUM1 + Y(M+I)
        SUM1   = SUM1 - Y(M+N+I)
        XNABLA = DF(I) - SUM1
        DLAN   = DLAN + XNABLA*XNABLA
  300 CONTINUE
      DLAN = SQRT(DLAN)
C
C----------------------------------------- VALUE OF THE LAGRANGEFUNCTION
      SUM = ZERO
      IF (M.GT.0) THEN
        DO 400 I=1,M
          SUM = SUM + Y(I) * G(I)
  400   CONTINUE
      ENDIF
      FLAG = FO - SUM
C
      TOL1   = SQRT(TOL)
      TOL2   = SQRT(TOL1)
      ISCHIT = 1
      IF (SSUM.LT.TOL.AND.SRES.LT.TOL1.AND.DLAN.LT.TOL2) ISCHIT = 0
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
      SUBROUTINE SLPERR (X,Y,G,SCX,XL,XU,N,M,MMAX,MNN,FO,ITER,FLAG,
     *                   IOUT,SSUM,SRES,DLAN,IFINAL,NFUNC,NGRAD,IPRINT,
     *                   NSEN,MAXCYC,NLIS)
C-----------------------------------------------------------------------
C     SEQUENTIAL LINEAR PROGRAMMING (SLP): ITERATION HISTORY
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(N),XL(N),XU(N),Y(MMAX),G(MMAX),SCX(N)
      PARAMETER (ZERO=0.D0)
C
C--------------------------------------------------WRITE ACTUAL SOLUTION
      IF(IFINAL.EQ.1) GOTO 555
      IF (IPRINT.EQ.1) THEN
        IF (ITER.GT.0) GOTO 999
      ELSE IF (IPRINT.EQ.2) THEN
        WRITE(IOUT,2000) ITER,FO,FLAG,SSUM,SRES,DLAN
      ELSE IF (IPRINT.EQ.3) THEN
        WRITE(IOUT,2001) ITER,FO,FLAG
        GOTO 999
      ENDIF
      WRITE(IOUT,2100) (X(I)*SCX(I),I=1,N)
      IF (M.GT.0) THEN
        WRITE(IOUT,2200) (Y(I),I=1,MNN)
        WRITE(IOUT,2300) (G(I),I=1,M)
      ENDIF
      GOTO 999
C
C-------------------------------------------------- WRITE FINAL SOLUTION
  555 CONTINUE
      WRITE(IOUT,3000)
      WRITE(IOUT,3100) FO
      WRITE(IOUT,3200) (X(I)*SCX(I),I=1,N)
      IF (M.GT.0) THEN
        WRITE(IOUT,3300) (Y(I),I=1,MNN)
        WRITE(IOUT,3400) (G(I),I=1,M)
      ENDIF
      WRITE(IOUT,3500) ((X(I)-XL(I))*SCX(I),I=1,N)
      WRITE(IOUT,3600) ((XU(I)-X(I))*SCX(I),I=1,N)
      WRITE(IOUT,3700) ITER,NFUNC,NGRAD,NSEN,MAXCYC,NLIS
C-----------------------------------------------------------------------
 999  CONTINUE
      RETURN
C
 2000 FORMAT(//4X,'ITERATION',I3,
     *        /7X,'FUNCTION VALUE:  F(X) =',E16.8,' L(X,U) =',E16.8,
     *        /7X,'                  SUM =',E16.8,
     *        /7X,'                 SRES =',E16.8,
     *        /7X,'                 DLAN =',E16.8)
 2001 FORMAT(4X,'ITERATION',I3,'  F(X) =',E16.8,' L(X,U) =',E16.8)
 2100 FORMAT(7X,'VARIABLE X =',/,(8X,4E16.8) )
 2200 FORMAT(7X,'MULTIPLIERS U =',/,(8X,4E16.8) )
 2300 FORMAT(7X,'CONSTRAINTS G(X) =',/,(8X,4E16.8) )
 3000 FORMAT(//4X,'FINAL CONVERGENCE ANALYSIS')
 3100 FORMAT(/7X,'OBJECTIVE FUNCTION VALUE:    F(X) =',E16.8)
 3200 FORMAT(7X,'APPROXIMATION OF SOLUTION:      X =',/,(8X,4E16.8) )
 3300 FORMAT(7X,'APPROXIMATION OF MULTIPLIERS:   U =',/,(8X,4E16.8) )
 3400 FORMAT(7X,'CONSTRAINT VALUES:           G(X) =',/,(8X,4E16.8) )
 3500 FORMAT(7X,'DISTANCE FROM LOWER BOUND: XL - X =',/,(8X,4E16.8) )
 3600 FORMAT(7X,'DISTANCE FROM UPPER BOUND: XU - X =',/,(8X,4E16.8) )
 3700 FORMAT(7X,'NUMBER OF ITERATIONS   :    ITER  = ',I4,
     *      /7X,'NUMBER OF FUNC CALLS   :    NFUNC = ',I4,
     *      /7X,'NUMBER OF GRAD CALLS   :    NGRAD = ',I4,
     *      /7X,'NUMBER OF SENS.ANALYSIS:    NSEN  = ',I4,
     *      /7X,'NUMBER OF CYCLINGS     :    MAXCYC= ',I4,
     *      /7X,'NUMBER OF LINESEARCHES :    NLIS  = ',I4)
C
      END
C=======================================================================
      SUBROUTINE SLPST1 (ITER,M,MMAX,ME,N,MNN,FO,G,DF,DG,X,SCX,ACTIVE,
     1                    DMIN,DMAX,S,IOUT)
C-----------------------------------------------------------------------
C     SEQUENTIAL LINEAR PROGRAMMING (SLP): SET SUBPROBLEM
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (ZERO=0.D0,ONE=1.D0)
      LOGICAL   ACTIVE
      DIMENSION G(MMAX),DF(N),DG(MMAX,N),SCX(N),ACTIVE(MMAX)
      DIMENSION X(N), S(MNN+1,N+1)
      DIMENSION DMIN(N), DMAX(N)
C
**********************************************************************
C     SLPST1 FORMS THE LINEAR SUBPROBLEM
**********************************************************************
      MN = M+N
      MS = MNN + 1
      NS = N + 1
C                                                         ==============
C-------------------------------------------------------- CLEAR MATRIX S
C                                                         ==============
      DO 55 J=1,N+1
       DO 66 I=1,MNN+1
        S(I,J) = ZERO
   66 CONTINUE
   55 CONTINUE
C
C                                            ===========================
C ------------------------------------------ INPUT OF OBJECTIVE FUNCTION
C                                            ===========================
      DO 100 I=1,N
         S(MNN+1,I) =  DF(I)
  100 CONTINUE
         S(MNN+1,N+1) = ZERO
C                                                   ====================
C ------------------------------------------------- INPUT OF CONSTRAINTS
C                                                   ====================
      IF(M.GT.0) THEN
        DO 200 I=1,N
          DO 210 J=1,M
             S(J,I)= -DG(J,I)
  210 CONTINUE
  200 CONTINUE
        DO 300 I=1,M
           S(I,N+1) = G(I)
  300 CONTINUE
      ENDIF
C
C                                            ===========================
C------------------------------------------- INPUT OF LOWER RESTRICTIONS
C                                            ===========================
      I = 0
      DO 500 J=1,N
         I = I + 1
         S(M+J,I) = -ONE
  500 CONTINUE
      DO 510 J=1,N
         SUM=ABS(DF(J))
         DO 520 K=1,M
  520    SUM=SUM+ABS(DG(K,J))
         IF(SUM.NE.ZERO) THEN
           S(M+J,N+1) = -DMIN(J)
         ELSE
           S(M+J,N+1) = ZERO
         ENDIF
  510 CONTINUE
C                                            ===========================
C------------------------------------------- INPUT OF UPPER RESTRICTIONS
C                                            ===========================
      I = 0
      DO 600 J=1,N
         I = I + 1
         S(MN+J,I) = ONE
  600 CONTINUE
      DO 610 J=1,N
         SUM=ABS(DF(J))
         DO 620 K=1,M
  620    SUM=SUM+ABS(DG(K,J))
         IF(SUM.NE.ZERO) THEN
           S(MN+J,N+1) = DMAX(J)
         ELSE
           S(MN+J,N+1) = ZERO
         ENDIF
  610 CONTINUE
C
C-----------------------------------------------------------------------
      RETURN
C
 9008 FORMAT(/5X,'INITIAL MATRIX')
 9095 FORMAT(/5X,'UPPER MOVE-LIMITS')
 9096 FORMAT(/5X,'LOWER MOVE-LIMITS')
      END
C=======================================================================
      SUBROUTINE SLPST2 (ITER,M,MMAX,ME,N,MNN,FO,G,DF,DG,X,SCX,ACTIVE,
     1                    DMIN,DMAX,S)
C-----------------------------------------------------------------------
C     SEQUENTIAL LINEAR PROGRAMMING (SLP): SET MODIFIED SUBPROBLEM
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(ZERO=0.D0,ONE=1.D0)
      LOGICAL ACTIVE
      DIMENSION G(MMAX),DF(N),DG(MMAX,N),SCX(N),ACTIVE(MMAX)
      DIMENSION X(N), S(MNN+3,N+2)
      DIMENSION DMIN(N), DMAX(N)
C
C
**********************************************************************
*     SLPST2 FORMS THE LINEAR SUBPROBLEM, ONE ARTIFICAL VARIABLE N+1
*     AND 2 RESTRICTIONS,MNN+1,MNN+2, ARE ADDED THO THE MATRIX S, TO
*     MAKE THE PROBLEM FEASIBLE
**********************************************************************
C
      MN = M + N + 1
      MS = MNN + 3
      NS = N + 2
C                                                     ==================
C---------------------------------------------------- CALCULATION OF RHO
C                                                     ==================
      SUM = ZERO
      DO 11 I=1,N
   11   SUM = SUM + DF(I)*DF(I)
      RHO = SQRT(SUM)/N
C
C                                                         ==============
C-------------------------------------------------------- CLEAR MATRIX S
C                                                         ==============
      DO 55 J=1,N+2
       DO 66 I=1,MNN+3
        S(I,J) = ZERO
   66 CONTINUE
   55 CONTINUE
C
C                                            ===========================
C ------------------------------------------ INPUT OF OBJECTIVE FUNCTION
C                                            ===========================
      DO 100 I=1,N
         S(MNN+3,I) =  DF(I)
  100 CONTINUE
         S(MNN+3,N+2) = ZERO
         S(MNN+3,N+1) = RHO
C                                                   ====================
C ------------------------------------------------- INPUT OF CONSTRAINTS
C                                                   ====================
      IF(M.GT.0) THEN
        DO 200 I=1,N
          DO 210 J=1,M
             S(J,I)= -DG(J,I)
  210 CONTINUE
  200 CONTINUE
C                                            ===========================
C ------------------------------------------ INPUT OF ARTIFICAL VARIABLE
C                                            ===========================
        DO 300 J=N+1,N+2
          DO 310 I=1,M
            S(I,J) = G(I)
  310 CONTINUE
  300 CONTINUE
      ENDIF
C
C                                            ===========================
C------------------------------------------- INPUT OF LOWER RESTRICTIONS
C                                            ===========================
      I = 0
      DO 500 J=1,N+1
         I = I + 1
         S(M+J,I) = -ONE
  500 CONTINUE
      DO 510 J=1,N
         S(M+J,N+2) = -DMIN(J)
  510 CONTINUE
         S(M+N+1,N+2) = ZERO
C
C                                            ===========================
C------------------------------------------- INPUT OF UPPER RESTRICTIONS
C                                            ===========================
      I = 0
      DO 600 J=1,N+1
         I = I + 1
         S(MN+J,I) = ONE
  600 CONTINUE
      DO 610 J=1,N
         S(MN+J,N+2) = DMAX(J)
  610 CONTINUE
         S(MN+N+1,N+2) = ONE
C
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
      SUBROUTINE SLPMOV (N, DMIN, DMAX,X,XL,XU,DXNEW,DXOLD,ITE,
     1                   DMIN1,DMAX1,FMOVE,IADAPT,MAXCYC)
C-----------------------------------------------------------------------
C     SEQUENTIAL LINEAR PROGRAMMING (SLP): SET MOVE LIMITS
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H, O-Z)
      PARAMETER(ZERO=0.D0,HALF=.5D0,ONE=1.D0)
      DIMENSION DMIN(N), DMAX(N),X(N),XL(N),XU(N),DXNEW(N),DXOLD(N)
      DIMENSION DMIN1(N),DMAX1(N)
      SAVE ADAPT
C
C                                                 ======================
C------------------------------------------------ ASSIGN NEW MOVE-LIMITS
C                                                 ======================
C
      IF(ITE.EQ.0) THEN
        ADAPT = ONE
        GOTO 111
      ENDIF
C
      DO 10 I=1,N
        CHANGE = DXNEW(I)*DXOLD(I)
        IF (CHANGE.LT.ZERO) MAXCYC = MAXCYC + 1
  10  CONTINUE
C
      IF (IADAPT.EQ.0) ADAPT = FMOVE
      IF (IADAPT.EQ.1) ADAPT = ONE/(ONE + ADAPT)
      IF (IADAPT.EQ.2) ADAPT = ONE/(ONE + ADAPT*ADAPT)
      IF (IADAPT.EQ.3) ADAPT = FMOVE*ADAPT
      DO 20 I=1,N
        DMIN1(I) = ADAPT*DMIN1(I)
        DMAX1(I) = ADAPT*DMAX1(I)
  20  CONTINUE
C
  111 CONTINUE
      DO 30 I=1,N
        DMIN(I) = MAX(DMIN1(I),XL(I)-X(I))
        DMAX(I) = MIN(DMAX1(I),XU(I)-X(I))
   30 CONTINUE
C
C-----------------------------------------------------------------------
      RETURN
      END
C=======================================================================
      SUBROUTINE SLPLIS (DX,DXOLD,X,DMIN,DMAX, DF,N,FO,FOPREV,ILIS,
     *                   IOUT,XL,XU,Y,RPEN,DG,G,ACTIVE,SCX,M,MMAX,ME,
     *                   MNN,ITER,NFUNC,NGRAD,NSEN,IGRAD,ILISE,AMIJO,
     *                   IALG,FF)
C-----------------------------------------------------------------------
C     SEQUENTIAL LINEAR PROGRAMMING (SLP): LINESEARCH
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H, O-Z)
      PARAMETER (ZERO=0.D0,ONE=1.D0,TWO=2.D0,EPS=1.D-6,R01=.1D0)
      PARAMETER (THREE=3.0D0,HALF=0.5D0,UF=1.0D-15)
      LOGICAL ACTIVE
      DIMENSION DX(N),DXOLD(N),DF(N),DMIN(N),DMAX(N),X(N),XL(N),XU(N)
      DIMENSION Y(MNN),DG(MMAX,N),G(MMAX),ACTIVE(MMAX),SCX(N),RPEN(MNN)
C
C-----------------------------------------------------------------------
      ILIS = 0
      RHO  = ONE
C                                        ===============================
C--------------------------------------- CALCULATIONS OF PENALTY FACTORS
C                                        ===============================
      DO 10 I=1,MNN
        RPEN(I) = MAX(ABS(Y(I)),HALF*(RPEN(I)+ABS(Y(I))))
   10 CONTINUE
C                               ========================================
C----------------------------- CALCULATION OF GRADIENT OF MERIT FUNCTION
C                               ========================================
      DFZERO = ZERO
C
      DO 30 I=1,N
        DFZERO = DFZERO + DF(I)*DX(I)
        IF(M.GT.0) THEN
          SUM = ZERO
          DO 20 J=1,M
            IF(G(J).GT.ZERO.AND.J.LE.ME) SUM = SUM + RPEN(J)*DG(J,I)
            IF(G(J).LE.ZERO) SUM = SUM - RPEN(J)*DG(J,I)
   20     CONTINUE
          DFZERO = DFZERO + SUM*DX(I)
        ENDIF
        IF (X(I).LE.XL(I)) DFZERO = DFZERO - RPEN(M+I)  *DX(I)
        IF (X(I).GE.XU(I)) DFZERO = DFZERO + RPEN(M+N+I)*DX(I)
   30 CONTINUE
C                                              =========================
C--------------------------------------------- CHECK GRADIENT AT RHO = 0
C                                              =========================
      IF (DFZERO.GT.ZERO) THEN
        WRITE (IOUT,1000)
        GOTO 950
      ENDIF
C
C                               ========================================
C------------------------------ CALCULATION OF MERIT FUNCTION AT RHO = 0
C                               ========================================
C
      FZERO = FOPREV
      IF(M.GT.0) THEN
        DO 40 J=1,M
        IF(J.LE.ME .OR. G(J).LT.ZERO) FZERO = FZERO + RPEN(J)*ABS(G(J))
   40   CONTINUE
      ENDIF
C                               ========================================
C------------------------------ CALCULATION OF MERIT FUNCTION AT RHO = 1
C                               ========================================
C
      DO 50 I=1,N
        X(I) = X(I) + DX(I)
   50 CONTINUE
C
      CALL FUNC (IALG,ITER)
      FO=FF
      NFUNC = NFUNC + 1
C
      FONE = FO
      IF(M.GT.0) THEN
        DO 60 J=1,M
        IF(J.LE.ME .OR. G(J).LT.ZERO) FONE = FONE + RPEN(J)*ABS(G(J))
   60   CONTINUE
      ENDIF
C                                     ==================================
C------------------------------------ CALCULATION OF GRADIENT AT RHO = 1
C                                     ==================================
      IF (ILISE.EQ.3) THEN
        ILIS = -1
        IF (M.GT.0) THEN
          NSEN = NSEN + M
          DO 70 I=ME+1,M
            ACTIVE(I) = .TRUE.
            IF(Y(I).EQ.ZERO.AND.G(I).GT.EPS) THEN
              ACTIVE(I)= .FALSE.
              NSEN = NSEN - 1
            ENDIF
   70     CONTINUE
        ENDIF
C
        CALL GRAD(IALG,ACTIVE)
        NGRAD = NGRAD + 1
C
        DFONE = ZERO
        DO 90 I=1,N
          DFONE = DFONE + DF(I)*DX(I)
          IF(M.GT.0) THEN
            SUM = ZERO
            DO 80 J=1,M
              IF(G(J).GT.ZERO.AND.J.LE.ME) SUM = SUM + RPEN(J)*DG(J,I)
              IF(G(J).LE.ZERO) SUM = SUM - RPEN(J)*DG(J,I)
   80       CONTINUE
            DFONE = DFONE + SUM*DX(I)
          ENDIF
          IF (X(I).LE.XL(I)) DFONE = DFONE - RPEN(M+I)  *DX(I)
          IF (X(I).GE.XU(I)) DFONE = DFONE + RPEN(M+N+I)*DX(I)
   90   CONTINUE
      ENDIF
C
C----------------------------------------------------------- AMIJO CHECK
C
      DIFF = FZERO+AMIJO*DFZERO-FONE
      IF (DIFF.GT.ZERO) GOTO 900
C                                                  =====================
C------------------------------------------------- QUADRATIC LINE SEARCH
C                                                  =====================
      IF (ILISE.EQ.2) THEN
        A = FONE-DFZERO-FZERO
        B = DFZERO
        C = FZERO
        IF (A.LT.UF) GOTO 900
        RHO = -B/(TWO*A)
C                                                      =================
C----------------------------------------------------- CUBIC LINE SEARCH
C                                                      =================
      ELSE IF (ILISE.EQ.3) THEN
        A = FZERO
        B = DFZERO
        C = THREE*(FONE-FZERO)-TWO*DFZERO-DFONE
        D = TWO*(FZERO-FONE)+DFZERO+DFONE
C
        IF (D.LT.EPS)  GOTO 900
        TT1 = -C/(3*D)
        TT2 = B/(3*D)
        DIFF  = (TT1*TT1 - TT2)
        IF (DIFF.LT.ZERO) GOTO 900
        RHO1  =  TT1 + SQRT(DIFF)
        RHO2  =  TT1 - SQRT(DIFF)
        IF (RHO1.GT.TT1) THEN
          RHO = RHO1
        ELSE IF (RHO2.GT.TT1) THEN
          RHO = RHO2
        ENDIF
      ENDIF
C
      RHO = MIN(RHO,ONE)
      RHO = MAX(RHO,R01)
C                                                  =====================
C------------------------------------------------- CALCULATION OF NEW DX
C                                                  =====================
      IF (RHO.EQ.ONE) GOTO 900
      ILIS = 1
      RHOX = RHO-ONE
      DO 100 I=1,N
         DX(I) = RHOX*DX(I)
  100 CONTINUE
C-----------------------------------------------------------------------
  900 CONTINUE
      IF (ILIS.LE.0) THEN
        IF (RHO.EQ.ONE) WRITE (IOUT,1100)
        IF (RHO.LT.ONE) WRITE (IOUT,1200) RHO
      ELSE
        IF (RHO.EQ.ONE) WRITE (IOUT,1300)
        IF (RHO.LT.ONE) WRITE (IOUT,1400) RHO
      ENDIF
  950 CONTINUE
      RETURN
C
 1000 FORMAT(7X,'*** LINE SEARCH NOT POSSIBLE, GRADIENT GT. ZERO')
 1100 FORMAT(7X,'*** LINE SEARCH,  ALPHA =  1.0')
 1200 FORMAT(7X,'*** LINE SEARCH,  ALPHA = ',E16.8)
 1300 FORMAT(7X,'*** LINE SEARCH (ADD. FUNCT. CALL), ALPHA =  1.0')
 1400 FORMAT(7X,'*** LINE SEARCH (ADD. FUNCT. CALL), ALPHA = ',E16.8)
      END
C=======================================================================
      SUBROUTINE SLPSIM (S,X,Y,F,NR,NC,M,ME,N,NF,IBL,EPS)
C-----------------------------------------------------------------------
C   @ COPYRIGHT ADAM BORKOWSKI 1985
C
C
C
C     SOLVES LP-PROBLEM:
C
C        MAX( C'*X )  SUBJECT TO   A*X#B ,  WHERE '#' IS '=' OR '<='
C
C     AND ITS DUAL:
C
C        MIN( B'*Y )  SUBJECT TO  A'*Y#C ,  WHERE '#' IS '=' OR '>='
C
C
C     PARAMETERS:
C
C       S(MNN+1,N+1) - SIMPLEX MATRIX
C       X(N)         - PRIMAL SOLUTION
C       Y(MNN)       - DUAL SOLUTION
C       F            - VALUE OF COST FUNCTION
C       NR(MNN)      - INDICES OF ROWS
C       NC(N)        - INDICES OF COLUMNS
C
C       M  - NUMBER OF CONSTRAINTS
C       ME - NUMBER OF EQUALITIES
C       MI - NUMBER OF INEQUALITIES
C       MNN- TOTAL NUMBER OF CONSTRAINTS + RESTRICTIONS
C       N  - NUMBER OF VARIABLES
C       NF - NUMBER OF FREE VARIABLES
C
C       IBL - ERROR FLAG ; IBL=0 - NO ERROR
C                          IBL=1 - WRONG DIMENSIONS
C                          IBL=2 - LINEAR DEPENDENT CONSTRAINTS
C                          IBL=3 - INCONSISTENT CONSTRAINTS
C                          IBL=4 - UNBOUNDED OBJECTIVE FUNCTION
C                          IBL=5 - UNDERFLOW IN PIVOT ELEMENT
C
C       EPS - NUMERICAL TRESHOLD
C
C       IPRIN - OUTPUT LEVEL
C               IPRIN=0 - FINAL RESULTS ONLY
C               IPRIN=1 - SIMPLEX MATRIX AFTER EACH ITERATION
C
C       IPROT - PROTOKOL FLAG
C               IPROT=0 - RESULTS DISPLAYED ON SCREEN
C               IPROT=1 - RESULTS SAVED ON DISK
C
C     INPUT MODES:
C          K - DATA COMING FROM KEYBOARD
C          D - DATA COMING FROM DISK
C
C     BACKUP FILES:
C          DATAS.DAT - INPUT DATA (NAME, M, ME, N, NF,
C                      EPS, IPRIN, IPROT)
C          RESUS.DAT - RESULTS (F, X, Y)
C
C     INITIAL SIMPLEX MATRIX:
C
C          !  A   B  !
C      S = !         !
C          ! -C'  0  !
C
C     REMARK:
C        EQUALITIES AND FREE VARIABLES (IF PRESENT)
C        MUST PRECEDE INEQUALITIES AND NON-NEGATIVE
C        VARIABLES
C
C
C***************************************************************
C
C
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION S(1),X(1),Y(1),NR(1),NC(1)
C
      IPRIN=0
      IPROT=1
C
C     TEST OF PARAMETERS
      IF(N.GT.0.AND.M.GE.N.AND.ME.GE.0
     *   .AND.NF.GE.0.AND.NF.LE.N) GO TO 1
         IBL= 1
         WRITE(6,9000)
         GO TO 999
C
C     INITIAL VALUES
 1    IBL= 0
      DO 100 I= 1,M
 100     NR(I)= -I
      DO 200 I= 1,N
 200     NC(I)= I
C
C     CALCULATION
      IF(ME.EQ.0) GO TO 2
         CALL SLPEQS (S,NR,NC,M,ME,N,NF,IBL,EPS)
         IF(IBL.GT.0) GO TO 999
 2          IF(NF.EQ.0) GO TO 3
               CALL SLPFRE (S,NR,NC,M,ME,N,NF,IBL,EPS)
               IF(IBL.GT.0) GO TO 999
 3                CALL SLPFEA (S,NR,NC,M,ME,N,NF,IBL,EPS)
                  IF(IBL.GT.0) GO TO 999
                     CALL SLPOPT (S,NR,NC,M,ME,N,NF,IBL,EPS)
                     IF(IBL.GT.0) GO TO 999
                        CALL SLPSOL (S,X,Y,NR,NC,M,N,F)
 999  CONTINUE
      RETURN
 9000 FORMAT(/1X,'*** SLPSIM: BAD PARAMETERS ***')
      END
C=======================================================================
      SUBROUTINE SLPEQS (S,NR,NC,M,ME,N,NF,IBL,EPS)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (ONE=1.0D0)
      DIMENSION S(1),NR(1),NC(1)
      IPRIN=0
      IPROT=1
C
C     ELIMINATES EQUATIONS
C
      MS= M+1
      NS= N+1
      DO 100 I= 1,ME
C
C     SEARCH FOR PIVOT COLUMN
         P= -ONE
         L= 0
         DO 101 J= 1,N
            IJ= (J-1)*MS+I
            R= ABS(S(IJ))
            IF(R.LE.P.OR.NC(J).LT.0) GO TO 101
               P= R
               L= J
 101     CONTINUE
C
C     TEST OF LINEAR DEPENDENCE
         IF(L.NE.0) GO TO 1
            IBL= 2
            WRITE(16,9000)
            GO TO 999
C
 1       CALL SLPJOR (S,NR,NC,MS,NS,I,L,IBL,EPS)
         IF(IBL.GT.0) GO TO 999
 100  CONTINUE
C
C     PERMUTATION OF COLUMNS
      DO 200 I= 1,N
         IF(NC(I).LT.0) GO TO 200
            DO 201 J= I,N
               IF(NC(J).GT.0) GO TO 201
                  CALL SLPWEK (S,NC,MS,NS,I,J)
 201        CONTINUE
 200  CONTINUE
C
 999  RETURN
 9000 FORMAT(/1X,'*** SLPEQS: LINEARLY DEPENDENT ROWS ***')
      END
C=======================================================================
      SUBROUTINE SLPFRE (S,NR,NC,M,ME,N,NF,IBL,EPS)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (ONE=1.0D0)
      DIMENSION S(1),NR(1),NC(1)
      IPRIN=0
      IPROT=1
C
C     ELIMINATES FREE VARIABLES
C
      MS= M+1
      NS= N+1
      IC= ME+1
      DO 100 I= IC,NF
C
C     SEARCH FOR PIVOT ROW
         P= -ONE
         K= 0
         DO 101 J= IC,M
            JI= (I-1)*MS+J
            R= ABS(S(JI))
            IF(R.LE.P.OR.NR(J).GT.0) GO TO 101
               P= R
               K= J
 101     CONTINUE
C
C     TEST OF LINEAR DEPENDENCE
         IF(K.GT.0) GO TO 1
            IBL= 2
            WRITE(16,9000)
            GO TO 999
C
 1       CALL SLPJOR (S,NR,NC,MS,NS,K,I,IBL,EPS)
         IF(IBL.GT.0) GO TO 999
 100  CONTINUE
C
C     PERMUTATION OF ROWS
      DO 200 I= IC,M
         IF(NR(I).GT.0) GO TO 200
            DO 201 J= I,M
            IF(NR(J).LT.0) GO TO 201
               CALL SLPSER (S,NR,MS,NS,I,J)
               GO TO 200
 201        CONTINUE
 200  CONTINUE
C
 999  RETURN
 9000 FORMAT(/1X,'*** SLPFRE: LINEARLY DEPENDENT ROWS ***')
      END
C=======================================================================
      SUBROUTINE SLPFEA (S,NR,NC,M,ME,N,NF,IBL,EPS)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (UBAR=1.0D10)
      DIMENSION S(1),NR(1),NC(1)
      IPRIN=0
      IPROT=1
C
C     FINDS FEASIBLE SOLUTION
C
      ITER= 0
      MS= M+1
      NS= N+1
      IC= ME+1
      JC= NF+1
      EPS1= -EPS
C
C     START OF MAIN LOOP
 1    CONTINUE
      IF(ITER.GT.( M+2*N) ) GOTO 999
C
C     TEST OF FEASIBILITY
      DO 100 I= JC,M
         INS= MS*N+I
         IF(S(INS).LT.EPS1) GO TO 2
 100  CONTINUE
C        WRITE(16,9000) ITER
         GO TO 999
C
C     SEARCH OF NEGATIVE FREE TERM
 2    K= I
      R= UBAR
      DO 200 I= IC,N
         KI= (I-1)*MS+K
         IF(S(KI).GT.EPS1) GO TO 200
            L= I
            GO TO 3
 200  CONTINUE
C
C     CONTRADICTORY CONSTRAINTS
      IBL= 3
      WRITE(16,9001)
      GO TO 999
C
C     SEARCH OF PIVOT ROW
 3    DO 300 I=JC,M
         IL= (L-1)*MS+I
         IF(ABS(S(IL)).LE.EPS) GO TO 300
            INS= N*MS+I
            P= S(INS)/S(IL)
            IF(P.LE.EPS1.OR.P.GE.R) GO TO 300
               IF(P.LE.EPS) GO TO 4
                  R= P
                  K= I
                  GO TO 300
 4             IF(S(IL).LE.EPS) GO TO 300
                  K= I
                  GO TO 5
 300  CONTINUE
C
 5    CALL SLPJOR (S,NR,NC,MS,NS,K,L,IBL,EPS)
      ITER= ITER+1
      GO TO 1
C     END OF MAIN LOOP
C
 999  RETURN
 9000 FORMAT(/1X,'FEASIBLE SOLUTION FOUND,        ITER= ',I4)
 9001 FORMAT(/1X,'*** SLPFEA: CONTRADICTORY CONSTRAINTS ***')
      END
C=======================================================================
      SUBROUTINE SLPOPT (S,NR,NC,M,ME,N,NF,IBL,EPS)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (UBAR=1.0D10)
      DIMENSION S(1),NR(1),NC(1)
      IPRIN=0
      IPROT=1
C
C     FINDS OPTIMUM SOLUTION
C
      ITER= 0
      MS= M+1
      NS= N+1
      IC= ME+1
      JC= NF+1
      EPS1= -EPS
C
C     START OF MAIN LOOP
 1    CONTINUE
      IF(ITER.GT.(M+2*N) ) GOTO 999
C
C     TEST OF OPTIMALITY
      DO 100 I= IC,N
         MSI= MS*I
         IF(S(MSI).LT.EPS1) GO TO 2
 100  CONTINUE
C        WRITE(16,9000) ITER
         GO TO 999
C
C     SEARCH FOR PIVOT COLUMN
 2    L= I
      R= UBAR
      K= 0
      DO 200 I= JC,M
         IL= (L-1)*MS+I
         INS= N*MS+I
         IF(ABS(S(IL)).LE.EPS) GO TO 200
            P= S(INS)/S(IL)
            IF(P.LE.EPS1.OR.P.GE.R) GO TO 200
               IF(P.GE.EPS) GO TO 3
                  IF(S(IL).GT.EPS) GO TO 4
                     GO TO 200
 3             R= P
               K= I
               GO TO 200
 4                K= I
                  GO TO 5
 200  CONTINUE
C
C     TEST OF BOUNDED COST FUNCTION
      IF(K.GT.0) GO TO 5
         IBL= 4
         WRITE(16,9001)
         GO TO 999
C
 5    CALL SLPJOR (S,NR,NC,MS,NS,K,L,IBL,EPS)
      ITER= ITER+1
      GO TO 1
C     END OF MAIN LOOP
C
 999  RETURN
 9000 FORMAT(/1X,'OPTIMUM SOLUTION FOUND,       ITER=',I4)
 9001 FORMAT(/1X,'*** SLPOPT: UNBOUNDED COST FUNCTION ***')
      END
C=======================================================================
      SUBROUTINE SLPSOL (S,X,Y,NR,NC,M,N,F)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (ZERO=0.0D0)
      DIMENSION S(1),X(1),Y(1),NR(1),NC(1)
C
C     GENERATES SOLUTIONS X, Y
C
      MS= M+1
      NS= N+1
C
C     CLEARING X, Y
      DO 100 I= 1,N
 100     X(I)= ZERO
      DO 200 I= 1,M
 200     Y(I)= ZERO
C
C     GENERATING X
      DO 300 I= 1,M
         IF(NR(I).LT.0) GO TO 300
            J= NR(I)
            INS= N*MS+I
            X(J)= S(INS)
 300  CONTINUE
C
C     GENERATING Y
      DO 400 I= 1,N
         IF(NC(I).GT.0) GO TO 400
            J= -NC(I)
            MSI= MS*I
            Y(J)= S(MSI)
 400  CONTINUE
C
C     OPTIMUM VALUE OF COST FUNCTION
      I= MS*NS
      F= S(I)
C
      RETURN
      END
C=======================================================================
      SUBROUTINE SLPJOR (A,NR,NC,M,N,K,L,IBL,EPS)
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER (ONE=1.0D0)
      DIMENSION A(1),NR(1),NC(1)
      IPRIN=0
      IPROT=1
C
C   PERFORMS JORDAN TRANSFORM
C   WITH PIVOT ELEMENT A(K,L)
C
      KL=(L-1)*M+K
C
C   TEST OF PIVOT ELEMENT
      IF(ABS(A(KL)).GT.EPS) GO TO 1
         IBL=5
C        WRITE(16,9000)
         GO TO 999
 1    P= ONE/A(KL)
C
C   TRANSF. OF ROWS 1 - (K-1)
      IF(K.EQ.1) GO TO 2
         I1= 1
         I2= K-1
         CALL SLPROW (A,M,N,K,L,I1,I2,P)
 2    CONTINUE
C
C   TRANSF. OF ROWS (K+1) - M
      IF(K.EQ.M) GO TO 3
         I1= K+1
         I2= M
         CALL SLPROW (A,M,N,K,L,I1,I2,P)
 3    CONTINUE
C
C   TRANSF. OF K-TH ROW
      DO 100 J=1,N
         KJ=(J-1)*M+K
         A(KJ)=A(KJ)*P
 100  CONTINUE
C
      A(KL)= P
      J= NR(K)
      NR(K)= NC(L)
      NC(L)= J
C
 999  RETURN
 9000 FORMAT(//5X,'*** SLPJOR: TOO SMALL PIVOT ELEMENT ***')
 9001 FORMAT(5X,'K=',I4,8X,'L=',I4)
 9003 FORMAT(/5X,'INDICES OF ROWS')
 9004 FORMAT(/5X,'INDICES OF COLUMNS')
      END
C=======================================================================
      SUBROUTINE SLPROW (A,M,N,K,L,I1,I2,P)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(1)
C
C   TRANSFORMS ROWS I1-I2 AFTER JORDAN
C
      DO 100 I= I1,I2
         IL= (L-1)*M+I
         AIL= A(IL)
         DO 200 J= 1,N
            J1= (J-1)*M
            IJ= J1+I
            KJ= J1+K
            A(IJ)= A(IJ)-A(KJ)*AIL*P
  200    CONTINUE
         A(IL)=-AIL*P
  100 CONTINUE
      RETURN
      END
C=======================================================================
      SUBROUTINE SLPWEK (S,NC,MS,NS,K,L)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION S(1),NC(1)
C
C     PERMUTES COLUMNS K, L
C
      DO 100 I= 1,MS
         IK= (K-1)*MS+I
         IL= (L-1)*MS+I
          Y= S(IK)
         S(IK)= S(IL)
         S(IL)= Y
 100  CONTINUE
      NC(NS)= NC(K)
      NC(K) = NC(L)
      NC(L) = NC(NS)
      RETURN
      END
C=======================================================================
      SUBROUTINE SLPSER (S,NR,MS,NS,K,L)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION S(1),NR(1)
C
C     PERMUTES ROWS K, L
C
      DO 100 I= 1,NS
         KI= (I-1)*MS+K
         LI= (I-1)*MS+L
          X= S(KI)
         S(KI)= S(LI)
         S(LI)= X
 100  CONTINUE
      NR(MS)= NR(K)
      NR(K) = NR(L)
      NR(L) = NR(MS)
      RETURN
      END
