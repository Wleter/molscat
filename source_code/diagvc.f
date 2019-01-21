      SUBROUTINE DIAGVC(A, LD, N, EVAL, EVEC)
C  Copyright (C) 2018 J. M. Hutson & C. R. Le Sueur
C  Distributed under the GNU General Public License, version 3
C
C  THIS SUBROUTINE DIAGONALISES THE N X N MATRIX A WITH LAPACK 
C  USING THE RELATIVELY ROBUST REPRESENTATION (RRR) ALGORITHM.
C  THE EIGENVALUES ARE RETURNED IN EVAL. 
C  THE EIGENVECTORS ARE RETURNED IN EVEC
C  J M HUTSON, JANUARY 2019
C
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: LD, N
      DOUBLE PRECISION, INTENT(INOUT) :: A(LD,N)
      DOUBLE PRECISION, INTENT(OUT)   :: EVEC(LD,N), EVAL(N)
      INTEGER IL, IU, LWORK, LIWORK, M, INFO
      DOUBLE PRECISION WORKSZ(1), VL, VU, ABSTOL, DLAMCH
      DOUBLE PRECISION, ALLOCATABLE :: E(:), WORK(:)
      INTEGER, ALLOCATABLE :: ISUPPZ(:), IWORK(:)
C
      ALLOCATE(ISUPPZ(2*N))
      ALLOCATE(E(N))
C
C  SOME CALLS IN MOLSCAT/BOUND MAY REQUIRE A TO BE PRESERVED.
C  DSYEVR DESTROYS THE LOWER TRIANGLE INCLUDING THE DIAGONAL.
C  SAVE DIAGONAL ELEMENTS IN E
C  LDW AND N ARE ALWAYS THE SAME IN MOLSCAT/BOUND 
C  (AND ACTUALLY THE SAME AS N)
C
      CALL DCOPY(N,A,LD+1,E,1)
C
C  ASK DSYEVR FOR OPTIMAL SIZES FOR WORK AND IWORK ARRAYS
C
      LWORK=-1
      ABSTOL=DLAMCH('S')
      CALL DSYEVR('V', 'A', 'L', N, A, LD, VL, VU, IL, IU, ABSTOL, M,
     1   EVAL, EVEC, LD, ISUPPZ, WORKSZ, LWORK, LIWORK, LWORK, INFO)
      LWORK=WORKSZ(1)
      ALLOCATE(IWORK(LIWORK))
      ALLOCATE(WORK(LWORK))
C
C  COMPUTE EIGENVALUES AND EIGENVECTORS
C
      CALL DSYEVR('V', 'A', 'L', N, A, LD, VL, VU, IL, IU, ABSTOL, M,
     1   EVAL, EVEC, LD, ISUPPZ, WORK, LWORK, IWORK, LIWORK, INFO)
C
      IF(INFO.NE.0) THEN
        WRITE(6,*)' *** ERROR: INFO =',INFO,' RETURNED FROM DSYEVR',
     1            ' IN DIAGONALISER'
        STOP
      ENDIF
C
C  RESTORE LOWER TRIANGLE OF A FROM UNCHANGED UPPER TRIANGLE
C  AND DIAGONAL FROM E
C
      CALL DSYFIL('L',N,A,LD)
      CALL DCOPY(N,E,1,A,LD+1)
C
      DEALLOCATE(E,ISUPPZ,WORK,IWORK)
C
      RETURN
      END