      module BASE9_SUITE

      IMPLICIT NONE

c  nqns and nlabvs are used to set variables nqn and nlabv that are
c  returned from set9
c  nqns is one greater than the number of quantum labels for a pair state
c  nlabvs is the number of indices for each term in the potential expansion
c  nqn is passed back into various routines in the basis-set suite
c  nlabvs is used internally in potin9
      integer, parameter        :: nqns=2, nlabvs=1

      end module BASE9_SUITE

C=========================================================================

      SUBROUTINE BAS9IN(PRTP, IBOUND, IPRINT)
      USE potential, ONLY: NCONST, NDGVL, NEXTMS, NEXTRA, NRSQ,
     1                     VCONST
      USE efvs  ! only needed if efvs are included in the calculation
      USE base9_suite
c     USE basis_data
      IMPLICIT NONE

      CHARACTER(32), INTENT(OUT)   :: PRTP

      INTEGER,       INTENT(INOUT) :: IBOUND

      INTEGER,       INTENT(IN)    :: IPRINT

      integer :: iconst,iextra,iefv


      RETURN
      END SUBROUTINE BAS9IN

C=========================================================================

      SUBROUTINE SET9(LEVIN, EIN, NSTATE, JSTATE, NQN, QNAME, NBLOCK,
     1                NLABV, IPRINT)
      USE basis_data ! needed only if H_intl is diagonal
      USE potential, ONLY: NCONST
      USE base9_suite
      IMPLICIT NONE

      INTEGER, INTENT(OUT)      :: NQN, NBLOCK, NLABV, NSTATE,
     1                             JSTATE(*)
c  JSTATE is conceptually dimensioned (NSTATE,NQN), but NSTATE is not
c  known on entry, so it cannot be given those dimensions explicitly.
c  It is instead indexed as a 1-dimensional array

      CHARACTER(8), INTENT(OUT) :: QNAME(10)

      LOGICAL, INTENT(IN)       :: LEVIN, EIN

      INTEGER, INTENT(IN)       :: IPRINT

      integer :: iqn,iqn1,iqn1mn,iqn1mx,iqn1st,ilevel,iloop,istate,
     1           nqlev


      RETURN
      END SUBROUTINE SET9

C=========================================================================

      SUBROUTINE BASE9(LCOUNT, N, JTOT, IBLOCK, JSTATE, NSTATE, NQN,
     1                 JSINDX, L, IPRINT)
      USE base9_suite
      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: N

      INTEGER, INTENT(OUT)   :: JSINDX(N), L(N)

      LOGICAL, INTENT(IN)    :: LCOUNT

      INTEGER, INTENT(IN)    :: JTOT, IBLOCK, NSTATE, NQN,
     1                          JSTATE(NSTATE,NQN), IPRINT

      integer :: ifunc,istate,iqn1,lmin,lmax,ll


      RETURN
      END SUBROUTINE BASE9

C=========================================================================

      SUBROUTINE POTIN9(ITYPP, LAM, MXLAM, NPTS, NDIM, XPT, XWT, MXPT,
     1                  IVMIN, IVMAX, L1MAX, L2MAX, MXLMB, XFN, MX,
     2                  IXFAC)
      USE base9_suite
      IMPLICIT NONE

      INTEGER, INTENT(INOUT)          :: ITYPP, MXLAM

      INTEGER, INTENT(OUT)            :: LAM(*)

      LOGICAL :: LVRTP
      DATA LVRTP/.FALSE./

      INTEGER, INTENT(IN)             :: IVMIN, IVMAX, L1MAX, L2MAX

      INTEGER :: L1

c  the remaining quantities are used only if quadrature is to be used
c  to project out potential expansion coefficients

      INTEGER, INTENT(INOUT)          :: NDIM, NPTS(NDIM), IXFAC, MX

      DOUBLE PRECISION, INTENT(INOUT) :: XFN(*)

      DOUBLE PRECISION, INTENT(OUT)   :: XPT(MXPT,NDIM), XWT(MXPT,NDIM)

      INTEGER, INTENT(IN)             :: MXLMB, MXPT

      double precision, allocatable :: fn1(:,:) !,fn2(:,:) etc
      integer :: il,npttot,nfun,ipt,ipt1tm,ipt1,i,ix
      double precision pt1,wt1


      RETURN
      END SUBROUTINE POTIN9

C=========================================================================

      SUBROUTINE CPL9(N, IBLOCK, NHAM, LAM, MXLAM, NSTATE, JSTATE,
     1                JSINDX, L, JTOT, VL, IV, CENT, DGVL, IBOUND,
     2                IEXCH, IPRINT)
      USE base9_suite
      USE potential, ONLY: NCONST, NDGVL, NRSQ, NVLBLK
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT) :: VL(NVLBLK,N*(N+1)/2), CENT(N),
     1                                 DGVL(N,NDGVL)

      INTEGER, INTENT(OUT)          :: IV(NVLBLK,N*(N+1)/2)

      INTEGER, INTENT(IN)           :: N, IBLOCK, NHAM, MXLAM,
     1                                 NSTATE, JSTATE(NSTATE,NQNS),
     2                                 JSINDX(N), L(N), JTOT, IBOUND,
     3                                 IEXCH, IPRINT, LAM(*)

      integer :: irc,icol,irow,i1col,i1row,lcol,lrow,iham,idgvl
      logical livuse
      data livuse/.false./


      RETURN
      END SUBROUTINE CPL9
C=========================================================================
