include(NyxSetupExecutable)

add_subdirectory(AMR-density)
add_subdirectory(AMR-zoom)
add_subdirectory(LyA)
add_subdirectory(LyA_AGN)
message("======================= HEATCOOL OUT = ${Nyx_HEATCOOL} =================")
add_subdirectory(LyA_Neutrinos)
