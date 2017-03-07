#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: fabm_pclake_utility
!
! !INTERFACE:
   module pclake_utility
!
!  Originally implemented by Fenjuan Hu (fenjuan@bios.au.dk)
!
! !USES:
   use fabm_types
   use fabm_expressions
   implicit none
!
!  default: all is private.
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public uFunTmAbio
   public uFunTmBio
   public uFunTmVeg
!
! !PRIVATE DATA MEMBERS:
   real(rk), parameter                 :: cTmRef=20.0_rk
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   pure real(rk) function uFunTmAbio(uTm,cTheta)

! !INPUT PARAMETERS:
   real(rk), intent(in)                :: uTm,cTheta
!
! !REVISION HISTORY:
!  Original author(s): Fenjuan Hu
!EOP
!-----------------------------------------------------------------------
!BOC
   uFunTmAbio = cTheta**(uTm-cTmRef)
   end function uFunTmAbio
!EOC

!-----------------------------------------------------------------------
!BOP
!
!
! !INTERFACE:
   pure real(rk) function uFunTmBio(uTm,cSigTm,cTmOpt)
!
!
! !INPUT PARAMETERS:
   real(rk), intent(in)                :: uTm,cSigTm,cTmOpt
!
! !REVISION HISTORY:
!  original author(s): Fenjuan Hu
!
!EOP
!-----------------------------------------------------------------------
!BOC
   uFunTmBio = exp(-0.5_rk/cSigTm**2 *((uTm-cTmOpt)**2- (cTmRef -cTmOpt)**2))
   end function uFunTmBio
!EOC

!-----------------------------------------------------------------------
!BOP
!
!
! !INTERFACE:
   pure real(rk) function uFunTmVeg(uTm,cQ10Veg)
!

! !INPUT PARAMETERS:
   real(rk), intent(in)                :: uTm,cQ10Veg
!
! !REVISION HISTORY:
!  Original author(s): Fenjuan Hu
!
!EOP
!-----------------------------------------------------------------------
!BOC
   uFunTmVeg = ((cQ10Veg )** (0.1_rk * (uTm - cTmRef)))
   end function uFunTmVeg
!EOC

!-----------------------------------------------------------------------

   end module pclake_utility

!------------------------------------------------------------------------------
! Copyright by the FABM_PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------
