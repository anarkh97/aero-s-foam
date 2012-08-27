      SUBROUTINE NLPOCM (IALG,
     &                   VAR,RESL,RESU,FF,G,DF,DG,
     &                   NUMVAR,M,ME,
     &                   ACC,BETA,DELTA,ETHA,XDGO,
     &                   MAXIT,LINES,FNAME,IFNSIZE,IFAIL)
C     ***************************************************************
C     *                                                             *
C     *           OC - Algorithm  -- main control routine           *
C     *                                                             *
C     ***************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D0,UF=1.0D-10)
C
      CHARACTER FNAME*(*),CDUM*1
      LOGICAL ACTIVE(M)
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
C------------------------------------------------------- CHECK PARAMETER
C
      IF (M.GT.ME .OR. ME.GT.1) THEN
       WRITE(*,*)    "ERROR IN NLPOCM"
       WRITE(*,*)    "OC-METHOD WORKS ONLY WITH ONE EQUALITY CONSTRAINT"
       WRITE(IOUT,*) "ERROR IN NLPOCM"
       WRITE(IOUT,*) "OC-METHOD WORKS ONLY WITH ONE EQUALITY CONSTRAINT"
       STOP "EXIT"
      ENDIF
C
C----------------------------------------- INITIALIZE ACTIVE CONSTRAINTS
C
      DO I=1,M
       ACTIVE(I) = .TRUE.
      ENDDO
C
C-------------------------------------------- SOLVE OPTIMIZATION PROBLEM
C
      CALL OCMSOLVE(IALG,
     &              VAR,RESL,RESU,FF,G,DF,DG,ACTIVE,
     &              NUMVAR,M,ME,
     &              ACC,BETA,DELTA,ETHA,XDGO,
     &              MAXIT,LINES,FNAME,IFNSIZE,IFAIL,IOUT)
C
      RETURN
      END

C=======================================================================
      SUBROUTINE OCMSOLVE (IALG,
     &                     VAR,RESL,RESU,FF,G,DF,DG,ACTIVE,
     &                     NUMVAR,M,ME,
     &                     ACC,BETA,DELTA,ETHA,XDGO,
     &                     MAXIT,LINES,FNAME,IFNSIZE,IFAIL,IOUT)
C     ***************************************************************
C     *                                                             *
C     *           OC - Algorithm  -- SOLUTION ITERATION             *
C     *                                                             *
C     ***************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D0,UF=1.0D-10)
C
      CHARACTER FNAME*(*)
C
      DIMENSION G(1)
      DIMENSION VAR(NUMVAR),VAROLD(NUMVAR)
C
      LOGICAL ACTIVE(M)
C
      NUMCON=ME
C
C-------------------------------------------------------- ITERATION LOOP
C
      DO 200 IT=1,MAXIT
C
C----------------- DETERMINE OBJECTIVE,CONSTRAINTS AND THEIR DERIVATIVES        
C
        CALL FUNC(IALG,IT)
C
        F=FF
C
        CALL GRAD(IALG,ACTIVE)
C
C----------------------------------------------------------- LINE SEARCH
C
        IF(LINES.GT.0) THEN
C
          IF(IT.GT.1) THEN
C
             CALL NLPLINE (IALG,IT,
     &                     VAR,VAROLD,FOLD,FF,F,G,ACTIVE,
     &                     NUMVAR,LINES,IOUT)
C
           ELSE
C
             FOLD=F 
C
           ENDIF
C
        ENDIF
c
C--------------------------------------------------------- PRINT ITERATE
C
        CALL NLPOCMP (IT,VAR,F,G,DF,DG,
     &                NUMVAR,NUMCON,IOUT,
     &                ACC,MAXIT,LINES,BETA,DELTA,ETHA,XDGO)
C
C--------------------------------------------- AND CHECK FOR CONVERGENCE
C
        IF(IT.GT.1) THEN 
          ISTOP = IOCMCONV(F,FOLD,ACC)
          IF (ISTOP .GT.0) GOTO 300
        ENDIF
C          
C-------------------------------------------------------- CALL OC-METHOD 
C
        DO I=1,NUMVAR
         VAROLD(I)=VAR(I)
        ENDDO
C
        CALL NLPOCM1 (VAR,RESU,RESL,DF,DG,
     &                NUMVAR,IOUT,ACC,BETA,
     & 	              DELTA,G(1),ETHA,XDGO)
C
C------------------------------------------------- END OF ITERATION LOOP
C
  200 CONTINUE
C
      IFAIL=1
c
      RETURN   
C
C-----------------------------------------------------------------------
C
  300 CONTINUE
c
      IFAIL=0
C
      RETURN
      END

C=======================================================================
C
      INTEGER FUNCTION IOCMCONV(F,FOLD,ACC)
C
C-----------------------------------------------------------------------
C
C     CHECK FOR CONVERGENCE
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER (ZERO=0.0D0,UF=1.0D-10)
C
      IOCMCONV = 0
C
      DIFF = ABS(F-FOLD)/MAX(UF,ABS(FOLD))
C
      IF (DIFF.LT.ACC) IOCMCONV=1
C
      FOLD = F
C   
      RETURN
      END
C
C=======================================================================
      SUBROUTINE  NLPOCMP (ITER,VAR,F,G,DF,DG,
     &                     NUMVAR,NUMCON,IOUT,
     &                     ACC,MAXIT,LINES,BETA,DELTA,ETHA,XDGO)
C-----------------------------------------------------------------------
C
C     print results for each iteration step
C
C-----------------------------------------------------------------------
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION VAR(NUMVAR),G(NUMCON)      
      DIMENSION DF(NUMVAR),DG(NUMCON,NUMVAR)
C
      IF (ITER.GT.1) GOTO 100

      WRITE(IOUT,1000)
      WRITE(IOUT,1100) MAXIT,LINES,ACC,DELTA,BETA,INT(ETHA),XDGO
      CALL FLUSH(IOUT)
C
 1000 FORMAT(////,
     & '------------------------------------------------------------',
     & /,/,
     & '              START OF OPTIMALITY CRITERIA METHOD           ',
     & /,/,
     & '------------------------------------------------------------',
     & ////)
C     
 1100 FORMAT(/,
     &     5X,'PARAMETERS:',
     &     //,
     &     8X,'MAXIT =',I13,/,
     &     8X,'LINES =',I13,/,
     &     8X,'ACC   =',D13.4,/,
     &     8X,'DELTA =',D13.4,/,
     &     8X,'BETA  =',D13.4,/,
     &     8X,'SHIFT =',I13,/,
     &     8X,'MCOR  =',D13.4,
     &     //)
C     
  100 CONTINUE   
C
      WRITE(IOUT,*)
      WRITE(IOUT,*)
      WRITE(IOUT,'(A,I5)') 'ITERATION ',ITER
      WRITE(IOUT,*)
      WRITE(IOUT,'(A)') 'VARIABLES:'
      DO IX=1,NUMVAR
        WRITE (IOUT,'(A,I5,E20.10)') 'VARIABLE ',IX,VAR(IX)
      ENDDO      
      WRITE(IOUT,*)      
      WRITE(IOUT,'(A,E20.10)') 'OBJECTIVE:',F
      DO IX=1,NUMVAR
        WRITE (IOUT,'(A,I5,E20.10)') 'GRADIENT VARIABLE ',IX,DF(IX)
      ENDDO
      WRITE(IOUT,*)      
      DO JX=1,NUMCON
        WRITE(IOUT,'(A,I5,E20.10)') 'CONSTRAINT:',JX,G(JX)
        DO IX=1,NUMVAR
          WRITE (IOUT,'(A,I5,E20.10)') 'GRADIENT VARIABLE ',IX,DG(JX,IX)
        ENDDO
      ENDDO      
      WRITE(IOUT,*)      
C
      return
      end

C=======================================================================
      SUBROUTINE  NLPOCM1 (VAR,RESU,RESL,DF,DG,
     &                     NUMVAR,NERR,ACCIT,BETA,
     & 	                   DELTA,CONACT,ETHA,XDGO)
C-----------------------------------------------------------------------
C
C     OC algorithm for optimization problems 
c     with one equality constraint
c
c     based on Ph.D. thesis of Kurt Maute 1998
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
c
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=0.5D0)
      PARAMETER (PI=3.1415927D0,PIHF=1.5707963D0) 
      PARAMETER (UP=1.0D+08,UL=1.0D-08)
c
      DIMENSION VAR(NUMVAR),RESU(NUMVAR),RESL(NUMVAR)
c      
      DIMENSION DF(NUMVAR),DG(NUMVAR),ETAI(NUMVAR)
      DIMENSION VARUP(NUMVAR),VARLO(NUMVAR)
C
C------------------------------------------------------ IN LINE FUNCTION
C
      POW(A,B) = EXP(B*LOG(A))
C
C------------------------------------------------------- INITIALIZE DATA
C
      ITRMAX=NUMVAR*10
C
C------------------------------------ GET CURRENT UPPER AND LOWER BOUNDS
C
      DO NVAR=1,NUMVAR            
C
        V1 = VAR(NVAR)*(ONE-DELTA)
        V2 = VAR(NVAR)*(ONE+DELTA)
        VARUP(NVAR) = MIN(RESU(NVAR),MAX(V1,V2))
        VARLO(NVAR) = MAX(RESL(NVAR),MIN(V1,V2))
C
      ENDDO
C                                                                       
C------------------------------------------------- SHIFT UPDATE GRADIENT 
C
C     ETHA = 0 .... NO SHIFT OF GRADIENTS
C     ETHA > 0 .... SHIFT OF GRADIENTS WITH CORRECTION OF DG
C     ETHA < 0 .... SHIFT OF GRADIENTS ONLY  
C
      SCLA=ZERO
C      
      IF(ETHA.NE.ZERO) THEN
C      
        XUDF = ZERO
        XLDF = ZERO
C
        IF(ETHA.GT.ZERO .AND. XDGO.LT.UL) 
     &    STOP 'NLPOCM: ERROR: ORTHOTROPIC MASS CORRECTION = 0 !'
C   
        DO NVAR=1,NUMVAR
C
          DGX  = DG(NVAR)
C
          IF (ABS(DGX).GT.UL) THEN
C
            IF(ETHA.GT.ZERO) THEN
              DGF  = XDGO * (ETHA + (ONE-ETHA)*(DGX/XDGO))
            ELSE
              DGF  = DGX
            ENDIF
C
            FAC  = DF(NVAR)/DGF
            XUDF = MAX(FAC,XUDF)
            XLDF = MIN(FAC,XLDF)
C
          ENDIF
C
        ENDDO
C      
        IF(XUDF.GT.ZERO) THEN
C        
          SCLA= XUDF + UL 
C
          WRITE(NERR,*) "" 
          WRITE(NERR,*) "NLPOCM: MODIFIKATION OF GRADIENTS IN NLPOCM"
          WRITE(NERR,*) "        XUDF:",XUDF
          WRITE(NERR,*) "        XLDF:",XLDF
C
        ENDIF
C
      ENDIF             
C                                                                       
C------------------------------------------- DETERMINE GRADIENT QUOTIENT 
C
      IF (ETHA.NE.ZERO) THEN
C
        DO NVAR=1,NUMVAR
C
          DGX  = DG(NVAR)
C
          IF(ABS(DGX).GT.ZERO) THEN
C
            IF(ETHA.GT.ZERO) THEN
              DGF  = XDGO * (ETHA + (ONE-ETHA)*(DGX/XDGO))
            ELSE
              DGF  = DGX
            ENDIF
C
            FAC  = SCLA - DF(NVAR)/DGF
C
            IF(FAC.LT.ZERO) THEN
              WRITE(NERR,*) "NLPOCM: ETHA > 0 , UPDATE RULE < 0"
              FAC = ZERO
            ENDIF  
C
          ELSE
C
            IF(DF(NVAR).GT.ZERO) THEN
              FAC = VARLO(NVAR)/VAR(NVAR)
            ELSEIF(DF(NVAR).LT.ZERO) THEN
              FAC = VARUP(NVAR)/VAR(NVAR)
            ELSE
              FAC = ONE
            ENDIF
C
          ENDIF
C
          FAC        =  POW(FAC,BETA) 
          ETAI(NVAR) =  FAC

        ENDDO
C
      ELSE 
C
        DO NVAR=1,NUMVAR
C
          IF(ABS(DG(NVAR)).GT.UL) THEN
C
            FAC = - DF(NVAR) / DG(NVAR) 
            IF(FAC.LT.ZERO) THEN
              WRITE(NERR,*) "NLPOCM: ETHA = 0 , UPDATE RULE < 0"
              FAC = ZERO
            ENDIF
C
          ELSE
C
            IF(DF(NVAR).GT.ZERO) THEN
              FAC = VARLO(NVAR)/VAR(NVAR)
            ELSEIF(DF(NVAR).LT.ZERO) THEN
              FAC = VARUP(NVAR)/VAR(NVAR)
            ELSE
              FAC = ONE
            ENDIF  
C
          ENDIF
C
          FAC        =  POW(FAC,BETA) 
          ETAI(NVAR) =  FAC
C
        ENDDO
C
      ENDIF
C
C--------------------------------------------- DETERMINE MIN/MAX OF RLAM      
C
      RLAMIN =  ONE/UL
      RLAMAX = -ONE/UL
C      
      DO NVAR=1,NUMVAR
        FAC    = VAR(NVAR)*ETAI(NVAR)
        IF(ABS(FAC).GT.UL) THEN
          FAC1   = VARUP(NVAR)/FAC
          FAC2   = VARLO(NVAR)/FAC
        ELSE
          FAC1   = UP
          FAC2   = UP
        ENDIF
        RLAMIN = MIN(RLAMIN,FAC1,FAC2)
        RLAMAX = MAX(RLAMAX,FAC1,FAC2)
      ENDDO
C                     
C--------------- CHECK IF CONSTRAINT CAN BE SATISFIED BY RLAMIN / RLAMAX
C
      SRMIN = ZERO
      SRMAX = ZERO
C      
      DO NVAR=1,NUMVAR
        IF(ABS(DG(NVAR)).GT.UL) THEN
          X = RLAMIN*VAR(NVAR)*ETAI(NVAR)
          X = MIN(VARUP(NVAR),MAX(VARLO(NVAR),X))
          SRMIN = SRMIN + DG(NVAR)*(X-VAR(NVAR))
          X = RLAMAX*VAR(NVAR)*ETAI(NVAR)
          X = MIN(VARUP(NVAR),MAX(VARLO(NVAR),X))
          SRMAX = SRMAX + DG(NVAR)*(X-VAR(NVAR))
        ENDIF
      ENDDO
C      
      SRMIN = SRMIN + CONACT 
      SRMAX = SRMAX + CONACT
      SRQ   = SRMIN/SRMAX
C      
      IF(SRQ.GT.ZERO) THEN 
        WRITE(NERR,*) 'CONSTRAINTS CAN NOT BE SATISFIED'
        IF(SRQ.LT.ONE) THEN
          WRITE(NERR,*) '- ALL VARIABLES SET ON CURRENT LOWER BOUNDS'
          RLAM = RLAMIN
          GOTO 900
        ELSE
          WRITE(NERR,*) '- ALL VARIABLES SET ON CURRENT UPPER BOUNDS'
          RLAM = RLAMAX
          GOTO 900
        ENDIF
      ENDIF  
C
C-------------------------------------- INITIALIZE LOOP TO DETERMIN RLAM 
C
      ITR   = 0
      SRDIF = TWO*ACCIT
C
      DO WHILE (ABS(SRDIF).GT.ACCIT .AND. ITR.LT.ITRMAX)
        ITR   = ITR + 1
        SR    = ZERO
        RLAM  = HALF*(RLAMIN+RLAMAX)
        DO NVAR=1,NUMVAR
          IF(ABS(DG(NVAR)).GT.UL) THEN 
            X  = RLAM*VAR(NVAR)*ETAI(NVAR)
            X  = MIN(VARUP(NVAR),MAX(VARLO(NVAR),X))
            SR = SR + DG(NVAR)*(X-VAR(NVAR))
          ENDIF 
        ENDDO
        SRDIF = SR + CONACT
        IF((SRDIF*SRMIN).GE.0 .AND. (SRDIF*SRMAX).LE.0 ) THEN
          RLAMIN = RLAM
          SRMIN  = SRDIF
        ELSEIF((SRDIF*SRMAX).GE.0 .AND. (SRDIF*SRMIN).LE.0 ) THEN
          RLAMAX = RLAM
          SRMAX  = SRDIF
        ELSE
          STOP 'ERROR IN NLPOCM: RLAM CAN NOT BE DETERMINED'
        ENDIF
      ENDDO
C
C----------------------------------------------------- CHECK CONVERGENCE
C
      IF(ITR.GE.ITRMAX) THEN
        WRITE(NERR,*) 'CONSTRAINTS CAN NOT BE SATISFIED - MAXIMUM NUMBE'
     *               ,'R OF ITERATION REACHED'
        WRITE(NERR,*) 'ERROR IN CONSTRAINT :',SRDIF
      ENDIF
C
C---------------------------------------------------------- FINAL UPDATE  
C
  900 CONTINUE
C
      DO NVAR=1,NUMVAR
        IF(ABS(DG(NVAR)).GT.UL) THEN 
          X         = RLAM*VAR(NVAR)*ETAI(NVAR)
          VAR(NVAR) = MIN(VARUP(NVAR),MAX(VARLO(NVAR),X))
        ELSE
          X         = VAR(NVAR)*ETAI(NVAR)
          VAR(NVAR) = MIN(VARUP(NVAR),MAX(VARLO(NVAR),X))
        ENDIF 
      ENDDO
C
C------------------------------------------------------------------- END
C      
      RETURN
      END
C=======================================================================
      SUBROUTINE  NLPLINE (IALG,IT,
     &                     VAR,VAROLD,FOLD,FF,F,G,ACTIVE,
     &                     NUMVAR,LINES,IOUT)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
c
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,HALF=0.5D0)
      PARAMETER (PI=3.1415927D0,PIHF=1.5707963D0) 
      PARAMETER (UP=1.0D+08,UL=1.0D-08)
C
      LOGICAL ACTIVE
c
      DIMENSION G(1)
C
      DIMENSION VAR(NUMVAR),VAROLD(NUMVAR),VARINC(NUMVAR)
c      
      IF (F.GT.FOLD) THEN
C
        MAXLINE = LINES
        ALM	= 1.0
        ALINE	= 0.5
C
        rkl	= SQRT(Abs(FOLD/G(1)))
C
        DO I=1,NUMVAR
          VARINC(I)=VAR(I)-VAROLD(I)
        ENDDO
C
        DO ILINE=1,MAXLINE
C
          DO I=1,NUMVAR
            VAR(I)=VAROLD(I)+ALINE*VARINC(I)
          ENDDO
C
          CAll OCMFUNC(IALG,IT,FF,F,G)
C
          XMERIT=F+rkl*ABS(G(1))
C          
          IF(XMERIT.LT.ALM*FOLD)  GOTO 111
C
          ALINE=ALINE*0.5
C         
        ENDDO
C
        WRITE(IOUT,*)'LINE SEARCH STOPPED WITH MAX. NUMBER OF STEPS'
C
c        DO I=1,NUMVAR
c          VAR(I)=VAROLD(I)+VARINC(I)
c        ENDDO
C
  111   WRITE(IOUT,*)'STEPS IN LINE SEARCH :',ILINE
        WRITE(IOUT,*)'STEP SIZE IN LINE SEARCH:',ALINE
C
      ENDIF                       
C
      RETURN
      END
C=======================================================================
      SUBROUTINE OCMFUNC(IALG,IT,FF,F,G)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C-----------------------------------------------------------------------
      CALL FUNC(IALG,IT)
      F=FF
      RETURN
      END
C=======================================================================
      SUBROUTINE OCMGRAD(IALG,DF,DG,ACTIVE)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C-----------------------------------------------------------------------
      LOGICAL ACTIVE
C
      CALL GRAD(IALG,ACTIVE)
C
      RETURN
      END
