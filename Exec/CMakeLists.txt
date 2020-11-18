include(NyxSetupExecutable)

add_subdirectory(AMR-density)
add_subdirectory(AMR-zoom)
add_subdirectory(LyA)
# add_subdirectory(LyA_AGN)  # This is broken
add_subdirectory(LyA_Neutrinos)

if (Nyx_HYDRO)
   add_subdirectory(HydroTests)
   add_subdirectory(DrivenTurbulence)
endif ()

add_subdirectory(MiniSB)
