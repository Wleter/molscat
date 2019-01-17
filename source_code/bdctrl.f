      SUBROUTINE BDCTRL(N, NSQ, MXLAM, NPOTL, Y1, Y2,
     1                  U, PSIMID, SUMPSI, VL, IV,
     2                  EINT, CENT, P, NODEC, ERED,
     3                  RMLMDA, EIGMIN, WAVE, IPRINT)
C  Copyright (C) 2018 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  ---------------------------------------------------------------
C  ROUTINE TO PERFORM A BOUND-STATE PROPAGATION. ON EXIT,
C  EIGMIN CONTAINS THE SMALLEST EIGENVALUE OF THE MATCHING MATRIX.
C  Y1 CONTAINS THE MATCHING MATRIX FOR USE BY ZERVEC
C  ---------------------------------------------------------------
C  THIS VERSION ORIGINALLY BY JM Hutson, 1988.
C  MODIFIED BY AE Thornley JULY 94 TO ALLOW WAVEFUNCTION PROPAGATION
C  FOR LDMD PROPAGATOR (ONLY)
C  EXTENDED BY JM Hutson 2006 TO ALLOW HYBRID LDMD/AIRY (ENERGIES ONLY)
C  MODIFIED SEPT 2012 TO ELIMINATE 3RD PROPAGATION SEGMENT AND
C  USE THE NUMBER OF NEGATIVE EIGENVALUES OF THE MATCHING MATRIX
C
C  03-12-15 CR Le Sueur:
C  3RD PROPAGATION SEGMENT BROUGHT BACK IN ORDER TO ALLOW
C  DECOUPLING OF MATCHING POINT FROM CHANGEOVER POINT BETWEEN
C  PROPAGATORS
C  CODE FOR PROPAGATORS EXTENDED TO CODE LONG AND SHORT RANGE
C  SEPARATELY -  VALUES PASSED IN PRPDTA COMMON BLOCK IN IFLG
C
C  DIMENSION STATEMENTS FOR ARGUMENT LIST
      DIMENSION U(NSQ),Y1(NSQ),Y2(NSQ),PSIMID(N),SUMPSI(N)
      DIMENSION P(MXLAM),VL(2),IV(2),EINT(N),CENT(N)
      LOGICAL WAVE
C
C  DYNAMIC MEMORY COMMON BLOCK
      COMMON /MEMORY/ MX,IXNEXT,NIPR,IDUMMY,X(1)
C
      DIMENSION STPSEG(3),CAYSEG(3),TOLSEG(3),DRSEG(3),IPRSEG(3),
     1          RBSEG(3),RESEG(3),NSTEPS(3),POWSEG(3)
C
C  COMMON BLOCK FOR CONTROL OF USE OF PROPAGATION SCRATCH FILE
      LOGICAL IREAD,IWRITE
      COMMON /PRPSCR/ ESHIFT,ISCRU,IREAD,IWRITE
C
C  COMMON BLOCK FOR CONTROL OF PROPAGATION SEGMENTS
      COMMON /RADIAL/ RMNINT,RMXINT,RMID,RMATCH,DRS,DRL,STEPS,STEPL,
     1                POWRS,POWRL,TOLHIS,TOLHIL,CAYS,CAYL,UNSET,
     2                IPROPS,IPROPL,NSEG
C
C  COMMON BLOCK FOR INPUT/OUTPUT CHANNEL NUMBERS
      LOGICAL PSIFMT
      COMMON /IOCHAN/ IPSISC,IWAVSC,IPSI,NWVCOL,PSIFMT

      CHARACTER(2) NCOL
      CHARACTER(40) F2000
C
      IF (WAVE) THEN
        CALL WVPROP(N)
        WRITE(NCOL,'(I2)') MIN(N,NWVCOL)
        F2000='(E15.7,'//NCOL//'E15.7)'

        IF (IPRINT.GE.30) WRITE(6,50)
  50    FORMAT('  WRITING OUT WAVEFUNCTION AT EVERY STEP')

      ENDIF

      NODEC=0
      NODES=-1
C
      IF (IREAD .OR. IWRITE) REWIND ISCRU

      IF (WAVE) SUMPSI=0D0
      IWREC=1
      NTSTPS=0

      CALL BDPSET(IPRSEG,RBSEG,RESEG,DRSEG,
     1            STPSEG,TOLSEG,CAYSEG,POWSEG,ICHNGE)
C
C  NOW LOOP OVER THE SEGMENTS.
C
      DO 100 ISEG=1,NSEG
        NUSED=0
        IPROP=IPRSEG(ISEG)
        CAY=CAYSEG(ISEG)
        DRT=DRSEG(ISEG)
        STEP=STPSEG(ISEG)
        TOLHIT=TOLSEG(ISEG)
        RSTART=RBSEG(ISEG)
        RSTOP=RESEG(ISEG)
        POW=POWSEG(ISEG)
        IF (IREAD) THEN
          READ(ISCRU) RSTART,RSTOP,EFIRST,ISTART,NSTEP,DR
          ESHIFT=ERED-EFIRST

        ELSE
          ISTART=1
          IF (ISEG.EQ.1 .OR. ISEG.EQ.ICHNGE+1) ISTART=0

          CALL DRSET(RSTART,RSTOP,-10.D0,TOLHIT,0.D0,DRT,NSTEP,DR,
     1               UNSET,IPROP,POW)
          IF (IWRITE) WRITE(ISCRU) RSTART,RSTOP,ERED,ISTART,NSTEP,DR
        ENDIF
C
C  INITIALISE Y MATRIX (HERE CALLED Y2)
C
        IT1=IXNEXT  ! DIAG
        IT2=IT1+N   ! EVAL
        IT3=IT2+N   ! EVECS
        IXNEXT=IT3+NSQ
        CALL CHKSTR(NUSED)
        IF (ISTART.EQ.0)
     1    CALL YINIT(Y2,U,VL,IV,P,CENT,EINT,X(IT1),X(IT2),
     2               X(IT3),N,MXLAM,NPOTL,
     3               ERED,RSTART,RMLMDA,RSTART.LT.RSTOP,
     4               IPRINT)
        IXNEXT=IT1

        IF (IPROP.EQ.5) THEN
C  JOHNSON'S LOG-DERIVATIVE PROPAGATOR
          IT1=IXNEXT   ! DIAG
          IXNEXT=IT1+N
          CALL CHKSTR(NUSED)
          CALL LDPROP(N,MXLAM,NPOTL,
     1                Y2,U,VL,IV,EINT,CENT,P,X(IT1),
     3                RSTART,RSTOP,NSTEP,DR,NODES,
     4                ERED,RMLMDA,IPRINT)
          IXNEXT=IT1
C
        ELSEIF (IPROP.EQ.6) THEN
C  MANOLOPOULOS' DIABATIC LOG-DERIVATIVE PROPAGATOR
          IT1=IXNEXT   ! Y14
          IT2=IT1+N    ! Y23
          IT3=IT2+N    ! DIAG
C  *** EXTRA STORAGE REQUIRED FOR WAVEFUNCTIONS
          IT4=IT3+N    ! W
          IT5=IT4+NSQ  ! W2
          IT6=IT5+NSQ  ! W3
          IF (WAVE) THEN
            IXNEXT=IT6+NSQ
          ELSE
            IXNEXT=IT4
          ENDIF
          CALL CHKSTR(NUSED)
          CALL MDPROP(N,MXLAM,NPOTL,
     1                Y2,U,VL,IV,EINT,CENT,P,
     2                X(IT1),X(IT2),X(IT3),X(IT4),X(IT5),X(IT6),
     3                RSTART,RSTOP,NSTEP,DR,NODES,IWREC,WAVE,
     4                ERED,RMLMDA,IPRINT)

          IXNEXT=IT1
C
        ELSEIF (IPROP.EQ.7) THEN
C  MANOLOPOULOS' QUASI-ADIABATIC LOG-DERIVATIVE PROPAGATOR
          IT1=IXNEXT  ! T
          IT2=IT1+NSQ ! Q
          IT3=IT2+NSQ ! W
          IT4=IT3+NSQ ! EIVAL
          IT5=IT4+N   ! Y1
          IT6=IT5+N   ! Y2
          IT7=IT6+N   ! Y3
          IT8=IT7+N   ! Y4
          IT9=IT8+N   ! DIAG
          IXNEXT=IT9+N
          CALL CHKSTR(NUSED)
          CALL MAPROP(N,NSQ,MXLAM,NPOTL,
     1                Y2,X(IT1),U,VL,IV,EINT,CENT,P,
     2                X(IT2),X(IT3),X(IT4),X(IT5),X(IT6),X(IT7),X(IT8),
     3                X(IT9),
     4                RSTART,RSTOP,NSTEP,DR,NODES,
     5                ERED,RMLMDA,IPRINT)
          IXNEXT=IT1
C
        ELSEIF (IPROP.EQ.8) THEN
C  MANOLOPOULOS+GRAY SYMPLECTIC PROPAGATOR
          IT1=IXNEXT
          IXNEXT=IT1+N
          CALL CHKSTR(NUSED)
          CALL MGPROP(N,MXLAM,NPOTL,
     1                Y2,U,VL,IV,EINT,CENT,P,X(IT1),
     3                RSTART,RSTOP,NSTEP,DR,NODES,
     4                ERED,RMLMDA,IPRINT)
          IXNEXT=IT1
C
        ELSEIF (IPROP.EQ.9) THEN
C  AIRY PROPAGATOR (LONG-RANGE ONLY)
            IT1=IXNEXT    ! W
            IT2=IT1+NSQ   ! Y1
            IT3=IT2+N     ! Y2
            IT4=IT3+N     ! Y3
            IT5=IT4+N     ! Y4
            IT6=IT5+N     ! VECNOW
            IT7=IT6+NSQ   ! VECNEW
            IT8=IT7+NSQ   ! EIGOLD
            IT9=IT8+N     ! EIGNOW
            IT10=IT9+N    ! HP
            IXNEXT=IT10+N
            CALL CHKSTR(NUSED)
            CALL AIPROP(N,MXLAM,NPOTL,
     1                  Y2,X(IT1),U,VL,IV,EINT,CENT,P,
     2                  X(IT2),X(IT3),X(IT4),X(IT5),X(IT6),X(IT7),
     3                  X(IT8),X(IT9),X(IT10),
     4                  RSTART,RSTOP,NSTEP,DR,POW,TOLHIT,NODES,
     5                  ERED,RMLMDA,IPRINT)
          IXNEXT=IT1
C
        ELSE
          WRITE(6,110) IPROP
  110     FORMAT(/'  *** ERROR IN BDCTRL. NO IMPLEMENTATION FOR ',
     1           'PROPAGATOR CODE =',I3)
          STOP
        ENDIF
C
C  ---------------------------------------------------------------
C  END OF PROPAGATION
C
        IF (IPRINT.GE.8) THEN
          IF (ISEG.EQ.1) WRITE(6,*)
          WRITE(6,120) RSTART,RSTOP,NSTEP,NODES
        ENDIF
  120   FORMAT('  LOG-DERIVATIVE MATRIX PROPAGATED FROM ',
     1         F12.4,'  TO ',F12.4,'  IN ',I6,'  STEPS.',
     2         I6,' NODES FOUND.')
C
        NODEC=NODEC+NODES
        NSTEPS(ISEG)=NSTEP
        IF (ISEG.EQ.ICHNGE .AND. .NOT.WAVE) THEN
          DO I=1,NSQ
            Y1(I)=Y2(I) ! Y1 IS NOW THE END OF THE OUTWARDS PROPAGATION
          ENDDO
        ELSEIF (ISEG.EQ.1 .AND. WAVE) THEN
C  Y1 IS NO LONGER NEEDED, SO IS USED TO STORE PSIMID FOR BOTH PROPAGATIONS
C  AWAY FROM RMATCH
          DO I=1,N
            Y1(I)=PSIMID(I)
          ENDDO
        ENDIF
C
        IF (WAVE .AND. (ISEG.EQ.NSEG .OR. ISEG.EQ.ICHNGE)) THEN
C  REVERSE PROPAGATION FROM RMATCH TO GET WAVEFUNCTION
          IF (ISEG.EQ.ICHNGE) THEN
C  KREC IS USED TO SAVE THE VALUE OF IWREC FOR SECOND PART OF WAVEFUNCTION
C  PROPAGATION.
C  IWREC IS WHERE TO READ RAB MATRIX FROM (SO GOES FROM LARGE NUMBER BACK
C  TOWARDS 1 IN BOTH CALLS TO EFPROP).
            KREC=IWREC
            ISLST=1
          ELSEIF (ISEG.EQ.NSEG) THEN
            ISLST=ICHNGE+1
          ENDIF
C  JPREC IS WHERE TO WRITE PSI TO (SO IN FIRST USE OF EFPROP GOES FROM
C  LARGE NUMBER BACK TOWARDS 1, AND IN SECOND USE OF EFPROP GOES FROM
C  LARGE NUMBER+1 FORWARDS TO EVEN BIGGER NUMBER).
          JPREC=KREC

          DO JSEG=ISEG,ISLST,-1
            IF (IPRINT.GE.8) WRITE(6,130) RESEG(JSEG),RBSEG(JSEG),
     1                                    NSTEPS(JSEG)
  130       FORMAT('  PROPAGATING WAVEFUNCTION FROM ',F10.4,' TO ',
     1             F10.4,' IN ',I6, ' STEPS'/)
C  Y2 AND U ARE NOW NO LONGER NEEDED SO THEY ARE REUSED FOR WORKSPACE HERE
            PSIMID=Y1(1:N)
            RSTART=RBSEG(ISEG)
            RSTOP=RESEG(ISEG)
            CALL EFPROP(N,RSTART,RSTOP,NSTEPS(JSEG),Y2,
     1                  U,PSIMID,IWREC,SUMPSI,IPRINT,JPREC)

          ENDDO

C  COLLECT THE CONTRIBUTIONS TO THE NORMALISATION CONSTANT
          IF (ISEG.EQ.NSEG) THEN
            IF (IPRINT.GE.6) WRITE(6,140)
 140        FORMAT(/'  NORMALISATION FACTOR CALCULATED BY ',
     1             'SIMPSON''S RULE FOR EACH CHANNEL')
            ANORM=0.D0
            DO I=1,N
              ANORM=ANORM+SUMPSI(I)
              IF (IPRINT.GE.8) WRITE(6,150) I,SUMPSI(I)
 150          FORMAT(I4,5X,G25.15)
            ENDDO
            IF (IPRINT.GE.6) WRITE(6,160) ANORM
 160        FORMAT(/'  TOTAL NORMALISATION FACTOR',5X,G18.10)
          ENDIF

        ENDIF

        IF (IPRINT.GE.15) CALL MATPRN(6,Y2,N,N,N,2,Y2,' LOGD MATRIX',1)
        NTSTPS=NTSTPS+NSTEP
  100 CONTINUE
C  END OF LOOP OVER SEGMENTS
C
C  WRITE OUT WAVEFUNCTION TO IPSI IN ORDER RMIN TO RMAX
C
      IF (WAVE) THEN
        CALL WVSTPS(NTSTPS+1)
        DO I=1,NTSTPS+1
          READ(IPSISC,REC=I) R,PSIMID
          IF (PSIFMT) THEN
            WRITE(IPSI,FMT=F2000) R,PSIMID/SQRT(ANORM)
          ELSE
            WRITE(IPSI) R,PSIMID/SQRT(ANORM)
          ENDIF 
        ENDDO
        CLOSE (IPSISC)
        CLOSE (IWAVSC)
        RETURN
      ENDIF

      DO I=1,NSQ
        Y2(I)=-Y2(I)+Y1(I) ! MATCHING MATRIX IS DEFINED AS INWARD-OUTWARD
        Y1(I)=Y2(I)        ! EVMTCH NEEDS A COPY OF THE MATCHING MATRIX
      ENDDO
C
      IF (IPRINT.GE.12) CALL MATPRN(6,Y2,N,N,N,2,Y2,
     1                              ' MATCHING MATRIX',1)
C
C  DIAGONALISE THE MATCHING MATRIX.
C
      IFAIL=0
      CALL F02AAF(Y2, N, N, U, X(IT1), IFAIL)
      IF (IFAIL.NE.0) GOTO 900
C
      NEGCNT=0
      EIGMIN=1.D30
      DO I=1,N
        IF (U(I).LE.0.D0) NEGCNT=NEGCNT+1
        IF (ABS(U(I)).GT.ABS(EIGMIN)) CYCLE
        EIGMIN=U(I)
      ENDDO
C  INCLUDE NEGCNT INSTEAD OF DOING 3RD PROPAGATION (SEPT 2012)
      NODEC=NODEC+NEGCNT
      IF (IPRINT.GE.8) WRITE(6,70) NEGCNT
  70  FORMAT('  COUNT OF NEGATIVE EIGENVALUES OF MATCHING MATRIX',
     1       ' GIVES',I6,' NODES.')
C
      IF (IPRINT.GE.9) WRITE(6,80) (I,U(I),I=1,N)
  80  FORMAT(/'  EIGENVALUES OF MATCHING MATRIX ARE:'/
     1       (2X,5(0P,I4,1P,1X,G15.7)))
C
C     IF (IPRINT.GE.8) WRITE(6,90) EIGMIN
  90  FORMAT(/'  SMALLEST EIGENVALUE OF MATCHING MATRIX IS ',1P,G11.4)
      RETURN
C
  900 WRITE(6,910) IFAIL
  910 FORMAT(/'  *** ERROR IN BDCTRL: EIGENVALUE ROUTINE FAILS WITH ',
     1       'IFAIL =',I3)
      STOP
      END
