include(NyxSetupExecutable)

add_subdirectory(AMR-density)

add_subdirectory(LyA)

if (Nyx_HEATCOOL)
  # add_subdirectory(LyA_AGN)  # This is broken
endif ()


if (Nyx_HYDRO)
   add_subdirectory(AMR-zoom)
   add_subdirectory(LyA_Neutrinos)
   add_subdirectory(HydroTests)
   add_subdirectory(MiniSB)
   add_subdirectory(DrivenTurbulence)
else()
   message(WARNING "\nDisabling AMR-zoom, LyA_Neutrinos, HydroTests, MiniSB and DrivenTurbulence: "
      "re-configure with -DNyx_HYDRO=YES to enable")
endif ()
