#include "fabm_driver.h"
!-----------------------------------------------------------------------
!BOP
!
! !INTERFACE:
   module foraging_arena

! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type, extends(type_base_model),public :: type_foraging_arena
!  state dependencies identifiers
!  defined prey and predator (starting from both bottom variables first, will improve later)
   type (type_bottom_state_variable_id)            :: id_prey, id_predator
!  Diagnostic variables for local output: vulnerable prey concentrtion, predation rate
   type (type_horizontal_diagnostic_variable_id)       :: id_aVulPrey, id_tPred
!  foraging arena theory parameters
   real(rk)   :: cEffPred,cVulPrey

   contains

!  Model procedures
   procedure :: initialize
   procedure :: do_bottom

   end type type_foraging_arena

!  private data members(API0.92)
   real(rk),parameter :: secs_pr_day=86400.0_rk
   real(rk),parameter :: NearZero = 0.000000000000000000000000000000001_rk
   real(rk)           :: Pi=3.14159265358979_rk

!EOP
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:
!
! !INTERFACE:
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
   class (type_foraging_arena), intent(inout), target :: self
   integer,                     intent(in)            :: configunit
!EOP
!-----------------------------------------------------------------------
!BOC
!  Store parameter values in our own derived type
!  NB: all rates must be provided in values per day,
!  and are converted here to values per second.
!  Paramters regarding foraging arena theory(Ahrens, Walters, & Christensen, 2012). 
!  effective searching rate for predator: has units area or volume per time searched by the predator divided by the area or volume of the foraging arena
   call self%get_parameter(self%cEffPred,     'cEffPred',     'gDW-1 m2 day-1','effective searching rate for predator',                                    default=1.0_rk)
!  definition of vulnerable factor for prey: prey exchange between the vulnerable and invulnerable states at instantaneous rates v, so vulerable status prey gains 
!  biomass at rate of v*(N-V), where N is the total prey, and V is the vulnerable part of prey, here v = cVulPrey
   call self%get_parameter(self%cVulPrey,     'cVulPrey',     '/day',          'vulnerable factor of prey',                                                default=1.0_rk)
!  Register local state variable
!  register prey pointer and predator pointer
   call self%register_state_dependency(self%id_prey,     'forgaing_arena_prey',     'g m-2', 'forgaing arena prey')
   call self%register_state_dependency(self%id_predator, 'foraging_arena_predator', 'g m-2', 'foraging arena predator')
! register output diagnostic variables
   call self%register_diagnostic_variable(self%id_tPred,    'tPred',      'g m-2 day-1','foraging arena predation rate', output=output_instantaneous)
   call self%register_diagnostic_variable(self%id_aVulPrey, 'aVulPrey',   'g m-2',      'vulnerable prey',               output=output_instantaneous)
   
   return

   end subroutine initialize
!EOC

!-----------------------------------------------------------------------
!BOP

! !IROUTINE:
!
! !INTERFACE:
   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
!
!  INPUT PARAMETERS:
   class (type_foraging_arena), intent(in)    :: self
   _DECLARE_ARGUMENTS_DO_BOTTOM_
!  LOCAL VARIABLES:
!  carriers for prey and predator
   real(rk)     :: prey, predator
!  foraging arena theory vaiables
   real(rk)     :: aVulPrey,tPred

!EOP
!-----------------------------------------------------------------------
!BOC
!  Spatial loop
!   _LOOP_BEGIN_
   _FABM_HORIZONTAL_LOOP_BEGIN_
!  Retrieve current (local) state variable values.
   _GET_HORIZONTAL_(self%id_prey,Prey)
   _GET_HORIZONTAL_(self%id_predator,predator)
!---------------------------------------------------------------------------
! foraging arena predation process
!---------------------------------------------------------------------------
    aVulPrey = prey * self%cVulPrey/(2.0_rk * self%cVulPrey + self%cEffPred * predator)
    tPred = self%cEffPred * aVulPrey * predator/secs_pr_day
!    print *, prey, predator,tPred

!  update abiotic variables in water
   _SET_ODE_BEN_(self%id_prey, -tPred)
   _SET_ODE_BEN_(self%id_predator, tPred)
!-----------------------------------------------------------------------
!  output diagnostic variables for external links
!-----------------------------------------------------------------------
!!  Export diagnostic variables
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_aVulPrey,aVulPrey)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_tPred,tPred * secs_pr_day)

! Horizontal loop end
   _FABM_HORIZONTAL_LOOP_END_
!
!EOP
!-----------------------------------------------------------------------

   end subroutine do_bottom

!EOC
!-----------------------------------------------------------------------
   end module foraging_arena

!------------------------------------------------------------------------------
! Copyright by the FABM-PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------
