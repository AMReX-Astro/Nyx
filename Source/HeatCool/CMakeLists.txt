if (Nyx_HEATCOOL)   # Nyx_HEATCOOL=ON

   target_include_directories(${EXE_NAME}
      PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

   target_sources( ${EXE_NAME}  PRIVATE f_rhs.H )

   if (Nyx_SDC)
      target_sources( ${EXE_NAME}  PRIVATE f_rhs_struct.H )
   endif()

   if (Nyx_CVODE)
         target_sources( ${EXE_NAME}
            PRIVATE
            integrate_state_vec_3d.cpp
            integrate_state_with_source_3d.cpp
            )
   else ()
      target_sources( ${EXE_NAME}
         PRIVATE
         integrate_state_box_3d_stubs.cpp
         integrate_state_with_source_3d_stubs.cpp
         )
   endif ()

else ()  # Nyx_HEATCOOL=OFF

   target_sources( ${EXE_NAME}
      PRIVATE
      integrate_state_box_3d_stubs.cpp
      integrate_state_exact_3d_stubs.cpp
      )

endif ()
