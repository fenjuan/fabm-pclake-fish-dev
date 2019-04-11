#include "fabm_driver.h"

module pclake_diel_vertical_movement

   use fabm_types
   use fabm_particle

   

   implicit none

!  default: all is private.
   private
   
   type,extends(type_particle_model),public :: type_pclake_diel_vertical_movement
!  register carrier for state variables that perform diel vertical movement
       type (type_state_variable_id) :: id_SwimVarD,id_SwimVarN,id_SwimVarP
!  register environment dependencies
       type (type_dependency_id)                :: id_par, id_I_0
       type (type_dependency_id)                :: id_temp, id_dens
       type (type_dependency_id)                :: id_depth   ! note it's the water depth, the current depth
!  register diagnostic variables
       type (type_diagnostic_variable_id)       :: id_swimspeed, id_rPredForg
!  parameters used
     !  diel vertical movement variables
      character(len=64) :: dvmvar
      integer    :: movement_methods
      real(rk)   :: vSwimVar
!  parameter for dvm case 3, STOKES
      real(rk)   :: cell_diameter, cell_density
!  parameters for dvm case 4, Ross et al, 1997
      real(rk)   ::  rCDSwimVar, cQSwimVarMax, cISwimVar
      real(rk)   :: f1, f2
!  parameters for dvm case 5, anti-predation window (Clark&Levy 1988 and Scheuerell&Schindle, 2003)
      real (rk)  :: h, prey_cell_diameter, prey_cell_density, pred_cell_density
      real(rk)   :: PredDeptMean, PredDeptVan, lowrisk, highrisk, fLife
      
       
   contains
      procedure :: initialize
      procedure :: get_vertical_movement
      
      
      
   end type type_pclake_diel_vertical_movement

!  private data members(API0.92)
   real(rk),parameter :: secs_pr_day=86400.0_rk
   real(rk),parameter :: NearZero = 0.000000000000000000000000000000001_rk

contains
      
      
      
      
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
      class (type_pclake_diel_vertical_movement),intent(inout),target :: self
      integer,                            intent(in)             :: configunit
!
! !REVISION HISTORY:
!

      
!EOP
!-----------------------------------------------------------------------
!BOC
!  pars for fish vertical movement
      call self%get_parameter(self%dvmvar,'dvmvar','[-]', 'variable that perform diel vertical movement',          default= 'FiJv' )
      call self%get_parameter(self%movement_methods,'movement_methods','[-]', '0 - no vertical movement, 1 - constant swimming speed, &
          &  2 - constant speed corrected by density, 3 - STOKKES, 4 - Ross&Sharples 2007',          default= 0 )
      call self%get_parameter(self%vSwimVar,'vSwimVar',  'm/s',   'constant swimming speed',     default=0.000012_rk)
!  pars for case 3
      call self%get_parameter(self%cell_density,'cell_density',  'kg/m3',   'cell density, max.1300, min.900',     default=1000.0_rk)
      call self%get_parameter(self%cell_diameter,'cell_diameter',  'm',   'cell size in diameter',     default = 0.00001_rk)
!  pars for case 4
      call self%get_parameter(self%f1,  'f1', '[-]',   'correction ratio of maximun nutrient quota for fish to swim up,0.6 < f1 < 0.8 ',   default= 0.675_rk)
      call self%get_parameter(self%f2,  'f2', '[-]',   'correction ratio of maximun nutrient quota for fish to swim down,0.6 < f2 < 0.9 ',   default= 0.75_rk)
      call self%get_parameter(self%rCDSwimVar, 'rCDSwimVar', '[-]',   'carbon/dry-weight ratio for swim variables, 0.4 - 0.5',   default=0.5_rk)
      call self%get_parameter(self%cQSwimVarMax,  'cQSwimVarMax', '[-]',   'maximun nutrient quota for swim variable',   default=0.28_rk)  ! need to adjust to fish
      call self%get_parameter(self%cISwimVar,  'cISwimVar', '',   'critical light intensity for swim variable',   default=0.1_rk)
!  extra pars for case 5
      call self%get_parameter(self%h,  'h', 's',   'average attack and handling time of fish',   default=1.8_rk)
      call self%get_parameter(self%prey_cell_density,'prey_cell_density',  'kg/m3',   'cell density of prey, zooplankton',     default=1000.0_rk)
      call self%get_parameter(self%prey_cell_diameter,'prey_cell_diameter',  'm',   'cell size in diameter, prey',     default = 0.00001_rk)
      call self%get_parameter(self%pred_cell_density,'pred_cell_density',  'kg/m3',   'cell density of predator,pisc.',     default=1000.0_rk)
      call self%get_parameter(self%PredDeptMean,'PredDeptMean',  'm',   'mean depth of predator',     default=19.7_rk)
      call self%get_parameter(self%PredDeptVan,'PredDeptVan',  'm**2',   'variance of predator depth',     default=7.17_rk)
      call self%get_parameter(self%lowrisk,'lowrisk',  '-',   'low risk ratio of u/f',     default=0.01_rk)
      call self%get_parameter(self%highrisk,'highrisk',  '-',   'high risk ratio of u/f',     default=0.5_rk)
!  this factor was used to low risk range, that how much the fish is trying to avoid its predator (between 0-1)
     call self%get_parameter(self%fLife,'fLife',  '-',   'life stage correction factor for avoiding predator,ranging 0 to 1',     default=0.5_rk)
      
      
      
!  register dependencies for swim variables
      call self%register_state_dependency(self%id_SwimVarD,'SwimVarD','g DW/m^3',   'swim variable in dry-weight')
      call self%register_state_dependency(self%id_SwimVarN,'SwimVarN','g N/m^3',    'swim variable in nitrogen')
      call self%register_state_dependency(self%id_SwimVarP,'SwimVarP','g P/m^3',    'swim variable in phosphorus')
!   the swim_var latter can be more general while PClake is modulized, that swim var can just use one
!   both for fish, zooplankton and phytoplankton. For now, I wrote it for current un-modulized fish
      call self%request_coupling_to_model(self%id_SwimVarD,'dvm_mod','sD'//trim(self%dvmvar))
      call self%request_coupling_to_model(self%id_SwimVarN,'dvm_mod','sN'//trim(self%dvmvar))
      call self%request_coupling_to_model(self%id_SwimVarP,'dvm_mod','sP'//trim(self%dvmvar))
!  environmental dependencies
      call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
      call self%register_dependency(self%id_I_0,standard_variables%downwelling_shortwave_flux)
      call self%register_dependency(self%id_dens,standard_variables%density)
      call self%register_dependency(self%id_temp,standard_variables%temperature)
      call self%register_dependency(self%id_depth,standard_variables%depth)
!  diagnostic variables
      call self%register_diagnostic_variable(self%id_swimspeed,  'swim_speed',      'm/s',     'vertical movement speed', missing_value = 0.0_rk, output=output_instantaneous)
      if (self%movement_methods == 5) then
          call self%register_diagnostic_variable(self%id_rPredForg,  'rPredForg',      '-',     'predation risk / foraging gain', missing_value = 0.0_rk, output=output_instantaneous)
     endif
   return
   
   
   end subroutine

   
!BOP
   subroutine get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)
   class (type_pclake_diel_vertical_movement),intent(in) :: self
   _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
   
   real(rk)    :: swim_speed
   real(rk)    :: par, I_0
   real(rk)    :: density, temperature
   real(rk)    :: depth
   real(rk)    :: SwimVarD, SwimVarN, SwimVarP
   real(rk)    :: QSwimVar, aNDSwimVar
   real(rk)    :: mu, mu20, dens20
!  for case 5, anti-predation window
   real(rk)    :: r, x, f
   real(rk)    :: PredDist, PredR
   real(rk)    :: u, rPredForg
   _LOOP_BEGIN_
!  retrieve environmental dependencies
   _GET_(self%id_par,par)
   _GET_(self%id_I_0,I_0)
   _GET_(self%id_dens,density)
   _GET_(self%id_temp,temperature)
   _GET_(self%id_depth, depth)
   
!  get the fish state
   _GET_(self%id_SwimVarD,SwimVarD)
   _GET_(self%id_SwimVarN,SwimVarN)
   _GET_(self%id_SwimVarP,SwimVarP)

!-------------------------------------------------------------------------------   
!  calculate water viscocity, if this is used in case 2, then will be moved to utility afterwards
!  adopted from aed2_util
!
! From Table A.1b, FLUID_MECHANICS With Engineering Applications
! by Robert L. Daugherty and Joseph B. Franzini,
! however, note these values are common in most fluid mechanics texts.
! NOTE: N s / m^2  = kg / m / s
!
!  Temp (C)     Viscosity (N s / m^2) x 10^3
!  --------     ---------
!      0          1.781
!      5          1.518
!     10          1.307
!     15          1.139
!     20          1.002
!     25          0.890
!     30          0.798
!     40          0.653
!     50          0.547
!     60          0.466
!     70          0.404
!     80          0.354
!     90          0.315
!    100          0.282
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!BEGIN
   !-- Check for non-sensical temperatures
   IF( temperature < 0.0_rk) temperature = 0.0
   IF( temperature > 100.0_rk ) temperature = 100.0

   IF( temperature <= 20.0 ) THEN
!  0C to 20C
!  y = 0.0008 * x^2 - 0.0556 * x + 1.7789
!  r^2 = 0.9999
     mu = 0.0008 * temperature**2. - 0.0556 * temperature + 1.7789
   ELSEIF(temperature <= 60) THEN
!  20C to 60C
!  y = 0.0002 * x^2 - 0.0323 * x + 1.5471
!  r^2 = 0.9997
     mu = 0.0002 * temperature**2. - 0.0323 * temperature + 1.5471
   ELSE
!  60C to 100C
!  y = 0.00006 * x^2 - 0.0141 * x + 1.1026
!  r^2 = 0.9995
     mu = 0.00006 *temperature**2. - 0.0141 * temperature + 1.1026
   ENDIF
!   Now convert to units of: N s / m^2
    mu = mu / 1e3

!  END OF WATER VISCOCITY CALCULATION
!-------------------------------------------------------------------------------


   SELECT CASE (self%movement_methods)
!  select the movement methods
!  case 0: no vertical movement
   case (0)
       swim_speed = 0.0_rk
!  case 1: constant movement
   case (1)
       swim_speed = self%vSwimVar
!  case 2: constant settling velocity at 20C corrected for density changes
   case (2)
        mu20 = 0.001002  ! N s/m2
        dens20 = 998.2000  ! kg/m3 (assuming freshwater)
        swim_speed = self%vSwimVar * mu20 * density / ( mu*dens20 )
!  case3: settling velocity based on Stokes Law calculation and cell density
   case (3)
!  Comments: the calculation of phytoplankton density will be implemented properly later
!   cell density should be a state variable
!           for now, use the simplified version
            swim_speed = -9.807*(self%cell_diameter**2.)*( self%cell_density - density ) / ( 18.* mu )
   case (4)
!!!! Keep in mind that the light units in this mehtods is uE m-2 s-1
! par and I_0 need to convert to uEm-2 -1
!  we will use W/m2 for fish condition
       ! strategy 1, adopted phytoplankton from the paper
!   the current N/dry-weight ratio of junenile fish 
     aNDSWimVar = SwimVarN / (SwimVarD +NearZero)
!   convert to N/C ratio, used in the equation adopted
     QSwimVar = aNDSwimVar * self%rCDSwimVar
!    calculate maximun nutrient qouta Qmax of juvenile fish
     if (QSwimVar < self%cQSwimVarMax * self%f1) then
         ! swimming down
        swim_speed = - self%vSwimVar
     elseif (QSwimVar > self%cQSwimVarMax * self%f2) then
         if (par > self%cISwimVar) then
             swim_speed = self%vSwimVar
         else    ! low light
             swim_speed = 0.0_rk  !  ?? if the light is critically low, shouldn't phy swim up?
         endif
     else
         swim_speed = 0.0_rk
     endif
     if (par > 0.95_rk * I_0)  swim_speed = 0.0_rk ! already at surface
   case (5)
! case 5: a simple light dependent vertical migration function for diurnal movememnt of zooplankton
       
!  case 6: anti-predation window theory, Brett et.al 1990,& ....2003
   case (6)
!  calculate light dependent reactive distance of fish to its prey(zooplankton)
     if (log10(par) < -3.6_rk) then
         r = 0.0_rk
     elseif (-3.6_rk < log10(par) < 0.377) then
         r = 0.0533*log10(par) + 0.199
     elseif (log10(par) > 0.377_rk) then
         r = 0.022_rk
     endif
!  calculate the foraging rate of fish on its prey(zooplankton)
!  here reuse pars cell diameter and cell density, it refers to prey (zooplankton, but not
!  the swim variable itself.
!  The first term of x is PI, I just chose this precision to be simple
     x = 3.1415926_rk *((r + self%prey_cell_diameter)**2) * self%vSwimVar * self%prey_cell_density
     f = x / (1 + self%h * x)  ! h--average attack and handling time (user defined)
!     print *, 'foraging gain', f
!  calculate the predator encounter probability P, i.e. the probability density function(pdf) of the preator depths
!  distribution times the predator temporal density (cell density) 
     PredDist = self%pred_cell_density * exp(-((depth - self%PredDeptMean)**2)/(2.0 * self%PredDeptVan))
!  calculate predator reactive distance, PredR
     if (par < 0.0757_rk) then
         PredR = 1.906*par**0.4747_rk
     else
         PredR = 1.012_rk
     endif
!  calculate the relative predation risk for an individual fish
     u = PredDist * PredR
!     print *, 'predation risk', u
!  calculate the ratio of predation risk to foraging gain( u :f)
     rPredForg= u / (f + NearZero)
!  important variable, make it output as diagnostic variable
     _SET_DIAGNOSTIC_(self%id_rPredForg, rPredForg)
     if (rPredForg < self%lowrisk) then
         swim_speed = self%vSwimVar ! safe from predator, approching to surface and get food
     elseif (self%lowrisk < rPredForg < self%highrisk) then
        swim_speed  = - self%fLife * self%vSwimVar   ! f is a life stage correction factor
     elseif ( rPredForg > self%highrisk) then
         swim_speed = -self%vSwimVar   ! too risky, going downward
     endif

!  case 6, the Ideal Free Distribution Model (IFD)

     

   END SELECT
   
!   print*, 'swim_speed=', swim_speed
   
!  Set vertical movement for juvenile fish
   _SET_VERTICAL_MOVEMENT_(self%id_SwimVarD,swim_speed)
   _SET_VERTICAL_MOVEMENT_(self%id_SwimVarN,swim_speed)
   _SET_VERTICAL_MOVEMENT_(self%id_SwimVarP,swim_speed)
!  updated diagnostic variables
    _SET_DIAGNOSTIC_(self%id_swimspeed, swim_speed)
    
    
   _LOOP_END_
   end subroutine get_vertical_movement
!EOC
!-----------------------------------------------------------------------

   end module pclake_diel_vertical_movement
!------------------------------------------------------------------------------
! Copyright by the FABM-PCLake-team under the GNU Public License - www.gnu.org
!------------------------------------------------------------------------------

