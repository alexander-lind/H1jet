*
* $Id: dgset.F,v 1.1.1.1 1996/04/01 15:02:15 mclareni Exp $
*
* $Log: dgset.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:15  mclareni
* Mathlib gen
*
*
C#include "gen/pilot.h"
*#if defined(CERNLIB_DOUBLE)
      SUBROUTINE DGSET(A,B,N,X,W)
C
C#include "gen/imp64.inc"
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      EXTERNAL DGQUAD
      DIMENSION X(*),W(*)

      CALL D107D1(2,DGQUAD,A,B,N,X,W)
      RETURN
      END
*#endif
