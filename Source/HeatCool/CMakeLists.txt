if (Nyx_HEATCOOL)   # Nyx_HEATCOOL=ON

   target_include_directories(nyx
      PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

   target_sources( nyx  PRIVATE f_rhs.H )

   if (Nyx_SDC)
      target_sources( nyx  PRIVATE f_rhs_struct.H )
   endif()

   if (Nyx_CVODE)
         target_sources( nyx
            PRIVATE
            integrate_state_vec_3d.cpp
            integrate_state_with_source_3d.cpp
            )
   else ()
      target_sources( nyx
         PRIVATE
         integrate_state_box_3d_stubs.cpp
         integrate_state_with_source_3d_stubs.cpp
         )
   endif ()

else ()  # Nyx_HEATCOOL=OFF

   target_sources( nyx
      PRIVATE
      integrate_state_box_3d_stubs.cpp
      integrate_state_exact_3d_stubs.cpp
      )

endif ()
