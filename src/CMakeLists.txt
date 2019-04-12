# Generate placeholder "pclake" module (versions all set to "unknown" - updated at build time if FABM_EMBED_VERSION is set)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/version.F90.in ${CMAKE_CURRENT_BINARY_DIR}/version.F90)

add_library(fabm_models_pclake OBJECT
            pclake_model_library.F90
            ${CMAKE_CURRENT_BINARY_DIR}/version.F90
            abiotic_sediment.F90
            abiotic_water.F90
            auxiliary.F90
            fish.F90
            macrophytes.F90
            phytoplankton_sediment.F90
            phytoplankton_water.F90
            utility.F90
            zoobenthos.F90
            zooplankton.F90
			pclake_dvm.F90
#			foraging_arena.F90
           )

add_dependencies(fabm_models_pclake fabm_base)

if(FABM_EMBED_VERSION)
  # Dynamically generate "pclake_version" module with git commit/branch information at every build.
  add_custom_target(fabm_models_pclake_version
     ${CMAKE_COMMAND} -DINFILE=${CMAKE_CURRENT_SOURCE_DIR}/version.F90.in -DOUTFILE=${CMAKE_CURRENT_BINARY_DIR}/version.F90 -P ${GET_GIT_INFO_SCRIPT}
     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
     COMMENT "Retrieving description of last ERSEM commit..." 
     SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/version.F90.in ${GET_GIT_INFO_SCRIPT}
     VERBATIM
  )
  add_dependencies(fabm_models_pclake fabm_models_pclake_version)
endif()


# added by F.Hu, for foraging arena test
# turn on foraging arena connection between piscvorous fish and zooplanktivorous and zoobenthivorous fish
add_definitions(-D_FORAGING_ARENA_Pisc_)
# turn on foraging arena connection between zooplankton and zooplanktivorous fish
add_definitions(-D_FORAGING_ARENA_FiJv_)
# turn on foraging arena connection between zoobenthos and benthivorous fish
add_definitions(-D_FORAGING_ARENA_FiAd_)
# turn on foraging arena connection between phytoplankton and zooplankton
add_definitions(-D_FORAGING_ARENA_Zoo_)
#-----------------------------------------------------------------------
# Copyright (C) 2017 - AU
#-----------------------------------------------------------------------
