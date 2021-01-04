target_sources( nyxcore
   PRIVATE
   strang_reactions.cpp
   )

if (Nyx_HEATCOOL)   # Nyx_HEATCOOL=ON

   target_include_directories(nyxcore
      PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )

   target_sources( nyxcore  PRIVATE f_rhs.H )

   if (Nyx_SDC) # Nyx_SDC=ON

      target_sources( nyxcore
         PRIVATE
         f_rhs_struct.H
         sdc_reactions.cpp )

      if (Nyx_SUNDIALS)
         target_sources( nyxcore
            PRIVATE
            integrate_state_with_source_3d.cpp
            )
      else ()
         target_sources( nyxcore
            PRIVATE
            integrate_state_with_source_3d_stubs.cpp
            )
      endif ()

   endif ()

   if (Nyx_SUNDIALS)
      target_sources( nyxcore
         PRIVATE
         integrate_state_vec_3d.cpp
         )
   else ()
      target_sources( nyxcore
         PRIVATE
         integrate_state_vec_3d_stubs.cpp
         )
   endif ()


   else () # Nyx_HEATCOOL=OFF

   target_sources( nyxcore
      PRIVATE
      integrate_state_vec_3d_stubs.cpp
      )

endif ()
