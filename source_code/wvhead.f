      SUBROUTINE WVHEAD(IPROGM,LABEL,ITYPE,URED,IPRINT)
C  Copyright (C) 2018 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
      USE efvs
      IMPLICIT NONE
      SAVE
C  THIS SUBROUTINE AND ITS ENTRY POINTS WRITE OUT INFORMATION REQUIRED
C  FOR WAVEFUNCTIONS
C
C  WRITTEN BY CRLS 11-11-2016

C  INPUTS FOR WVHEAD
      INTEGER, INTENT(IN)::IPROGM,ITYPE,IPRINT
      CHARACTER(80), INTENT(IN):: LABEL
      DOUBLE PRECISION, INTENT(IN):: URED

C  INPUTS FOR WVINFO
      INTEGER, INTENT(IN):: JTOT,M,ICHND,N,NQN,NSTATE,JSTATE,JSINDX,L
      DOUBLE PRECISION, INTENT(IN):: ENERGY,EREF,EFACT
      CHARACTER(8), INTENT(IN):: UNAME
      DIMENSION JSTATE(NSTATE,NQN),JSINDX(N),L(N)

C  INPUTS FOR WVSTPS
      INTEGER, INTENT(IN):: NTSTPS

C  COMMON BLOCK FOR PROPAGATOR CONTROL
      INTEGER MXSEG,IDIR,IFLG,NSEGS,NSTEPS
      DOUBLE PRECISION RST,REND,DRPR
      PARAMETER (MXSEG=3)
      COMMON/PRPDTA/RST(MXSEG),REND(MXSEG),DRPR(MXSEG),
     1              IDIR(MXSEG),IFLG(MXSEG),NSEGS,NSTEPS(MXSEG)

C  COMMON BLOCK FOR INPUT/OUTPUT CHANNEL NUMBERS
      COMMON /IOCHAN/ IPSISC,IWAVSC,IPSI,NWVCOL,PSIFMT
      LOGICAL PSIFMT
      INTEGER IPSISC,IWAVSC,IPSI,NWVCOL

C  COMMON BLOCK TO DESCRIBE WHICH DRIVER IS USED
      COMMON /CNTROL/ CDRIVE
      CHARACTER(1) CDRIVE

C  INTERNAL VARIABLES
      CHARACTER(2) NCOL
      CHARACTER(40) F100,F101
      CHARACTER(7) NPROGM

      INTEGER JQN,ILEV,ISEG,IEFV,LENEFV

C  SCATTERING WAVEFUNCTIONS ARE COMPLEX, SO HALVE THE NUMBER OF COLUMNS
C  IN THE WAVEFUNCTION FILE
      IF (CDRIVE.EQ.'M') THEN
        NPROGM='MOLSCAT'
        NWVCOL=10
      ELSEIF (CDRIVE.EQ.'F') THEN
        NPROGM='FIELD'
        NWVCOL=20
      ELSEIF (CDRIVE.EQ.'B') THEN
        NPROGM='BOUND'
        NWVCOL=20
      ENDIF

C  OPEN UNIT FOR STORING WAVEFUNCTIONS
      IF (PSIFMT) THEN
        CLOSE(IPSI)
        OPEN(UNIT=IPSI,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     1       FORM='FORMATTED')
      ELSE
        CLOSE(IPSI)
        OPEN(UNIT=IPSI,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     1       FORM='UNFORMATTED')
      ENDIF
      WRITE(6,1000) IPSI
 1000 FORMAT(/'  WAVEFUNCTION REQUESTED, WRITING TO FILE ON UNIT ',I3)

C  WRITE HEADER WHICH GOES ONCE AT THE TOP OF THE WAVEFUNCTION FILE
      IF (PSIFMT) THEN
        WRITE(IPSI,'(A,A,A,I3)') '# CALCULATION FROM ',NPROGM,
     1                           ' OUTPUT FORMAT ',IPROGM
        WRITE(IPSI,1010) ' # LABEL ON CALCULATION WAS: ',TRIM(LABEL)
 1010   FORMAT(A,A)
        WRITE(IPSI,*) '# ITYPE AND REDUCED MASS ARE ',ITYPE,URED
      ELSE
        WRITE(IPSI) IPROGM
        WRITE(IPSI) LABEL
        WRITE(IPSI) ITYPE,URED
      ENDIF

      RETURN
C =========================================================== END OF WVHEAD
      ENTRY WVINFO(JTOT,M,ICHND,N,NQN,NSTATE,JSTATE,
     1             JSINDX,L,ENERGY,EREF,EFACT,UNAME)

      WRITE(NCOL,'(I2)') MIN(N,NWVCOL)
      IF (CDRIVE.NE.'M') THEN
        F100='(" # QNAME ",I2,14X,'//NCOL//'(I3,12X))'
        F101='(" # L     ",16X,'//NCOL//'(I3,12X))'
      ELSE
        F100='(" # QNAME ",I2,14X,'//NCOL//'(I3,27X))'
        F101='(" # L     ",16X,'//NCOL//'(I3,27X))'
      ENDIF

      LENEFV=MAX(NEFVP-IEFVST+1,0)
      IF (IPRINT.GE.1 .AND. LENEFV.GT.0)
     1  WRITE(6,1020) 'HEADER ON WAVEFUNCTION FILE CONTAINS ',
     2                'THE FOLLOWING VARIABLES: ',
     3                (TRIM(EFVNAM(IEFV)),IEFV=IEFVST,NEFVP)
 1020   FORMAT(2X,A,A,(A:,', '))

C  WRITE A HEADER FOR THIS PARTICULAR WAVEFUNCTION (QUANTUM NUMBERS,
C  EXTERNAL FIELDS)
      IF (PSIFMT) THEN
        WRITE(IPSI,1030) JTOT,M,(TRIM(EFVNAM(IEFV)),EFV(IEFV),
     2                           IEFV=IEFVST,NEFVP)
 1030   FORMAT(' # JTOT AND M ARE ',2I4:/ ' # EXTERNAL VARIABLES ARE: ',
     3         10(A,' = ',G20.13))

        IF (N.GT.NWVCOL) THEN
          WRITE(6,*) ' **** WAVEFUNCTION COLUMNS ARE SPLIT OVER MORE ',
     1               'THAN ONE ROW'
          WRITE(IPSI,*) '# WAVEFUNCTIONS ARE SPLIT OVER MORE THAN ONE ',
     1                  'ROW - BEWARE!'
        ENDIF
        WRITE(IPSI,*) '# NUMBER OF QUANTUM LABELS, SIZE OF BASIS,',
     1                ' AND TOTAL SIZE OF PRIMITIVE BASIS'
        WRITE(IPSI,*) '# ',NQN-1,N,NSTATE

C  WRITE THE BASIS FUNCTION LABELS USED FOR THIS WAVEFUNCTION
        DO JQN=1,NQN-1
          WRITE(IPSI,FMT=F100) JQN,(JSTATE(JSINDX(ILEV),JQN),ILEV=1,N)
        ENDDO
        WRITE(IPSI,FMT=F101) (L(ILEV),ILEV=1,N)

C  INFORMATION ABOUT THE PARTICULAR WAVEFUNCTION (ENERGY AND EITHER
C  INCOMING CHANNEL, OR NODE NUMBER)
        IF (CDRIVE.EQ.'M') THEN
          WRITE(IPSI,1040) ' # SCATTERING CHANNEL IS ',ICHND,
     1                     ' AT ENERGY ',(ENERGY-EREF)/EFACT,UNAME,
     2                     ' REFERRED TO REFERENCE ENERGY = ',
     3                     EREF/EFACT,UNAME
        ELSE
          WRITE(IPSI,1040) ' # WAVEFUNCTION FOR NODE ',ICHND,
     1                     ' WITH BINDING ENERGY ',(ENERGY-EREF)/EFACT,
     2                     UNAME,' REFERRED TO REFERENCE ENERGY = ',
     3                     EREF/EFACT,UNAME
 1040     FORMAT(A,I6,A,G20.13,1X,A,A,G20.13,1X,A)
        ENDIF

      ELSE
         WRITE(IPSI) JTOT,M,ENERGY,EREF,EFACT,
     1               (EFV(IEFV),IEFV=IEFVST,NEFVP)

        WRITE(IPSI) NQN-1,NSTATE,N
        DO JQN=1,NQN-1
          WRITE(IPSI) (JSTATE(JSINDX(ILEV),JQN),ILEV=1,N)
        ENDDO
        WRITE(IPSI) (L(ILEV),ILEV=1,N)
        WRITE(IPSI) ICHND
      ENDIF

      RETURN
C =========================================================== END OF WVINFO
      ENTRY WVPROP(N)

C  OPEN SCRATCH FILES USED FOR STORING WAVEFUNCTION AND PROPAGATION
C  MATRICES
      OPEN(UNIT=IPSISC,FORM='UNFORMATTED',ACCESS='DIRECT',
     1     RECL=8*(N+1),STATUS='SCRATCH')
      OPEN(UNIT=IWAVSC,FORM='UNFORMATTED',ACCESS='DIRECT',
     1     RECL=8*(N*N+1),STATUS='SCRATCH')
C  INFORMATION ABOUT THE PROPAGATION

      IF (PSIFMT) THEN
        WRITE(IPSI,*) '# NUMBER OF PROPAGATIONS IS ',NSEGS
        WRITE(IPSI,1050)
 1050   FORMAT(' #',10X,'RST',17X,'REND',16X,'DRPR',9X,'INTFLG IDIR')
        WRITE(IPSI,1060) (RST(ISEG),REND(ISEG),DRPR(ISEG),
     1                    IFLG(ISEG),IDIR(ISEG),ISEG=1,NSEGS)
 1060   FORMAT((' #',3(G20.13,1X),I3,5X,I3))
      ELSE
        WRITE(IPSI) NSEGS
        WRITE(IPSI) (RST(ISEG),REND(ISEG),DRPR(ISEG),
     1               IFLG(ISEG),IDIR(ISEG),ISEG=1,NSEGS)
      ENDIF

      RETURN
C =========================================================== END OF WVPROP
      ENTRY WVSTPS(NTSTPS)
C  INFORMATION ABOUT THE TOTAL NUMBER OF STEPS

      IF (PSIFMT) THEN
        WRITE(IPSI,*) '# TOTAL NUMBER OF STEPS IS ',NTSTPS
      ELSE
        WRITE(IPSI) NTSTPS
      ENDIF

      END