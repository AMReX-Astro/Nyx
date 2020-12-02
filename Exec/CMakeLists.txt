include(NyxSetupExecutable)

add_subdirectory(AMR-density)
add_subdirectory(AMR-zoom)
add_subdirectory(LyA)

if (Nyx_HEATCOOL)
  # add_subdirectory(LyA_AGN)  # This is broken
endif ()

add_subdirectory(LyA_Neutrinos)

if (Nyx_HYDRO)
   add_subdirectory(HydroTests)
   add_subdirectory(DrivenTurbulence)
endif ()

add_subdirectory(MiniSB)
