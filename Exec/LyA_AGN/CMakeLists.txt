set(_sources     Prob.H)
set(_input_files inputs inputs.rt )
list(APPEND _input_files 32.nyx 64sssss_20mpc.nyx TREECOOL_middle)

# Required
if (NOT Nyx_AGN)
   message(STATUS "Setting Nyx_AGN=ON for ${CMAKE_CURRENT_LIST_DIR}")
   set(Nyx_AGN ON)
endif ()

if (NOT Nyx_HEATCOOL)
   message(STATUS "Setting Nyx_HEATCOOL=ON for ${CMAKE_CURRENT_LIST_DIR}")
   set(Nyx_HEATCOOL ON)
endif ()

message("======================= HEATCOOL = ${Nyx_HEATCOOL} =================")
nyx_setup_executable(_sources _input_files)

unset(_sources)
unset(_input_files)
