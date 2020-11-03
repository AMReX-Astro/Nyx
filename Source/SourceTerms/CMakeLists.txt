target_sources( nyx
   PRIVATE
   Nyx_time_center_sources.cpp
   )

if (NOT Nyx_AGN)
   target_sources( nyx
      PRIVATE
      Nyx_sources.cpp
      )
endif ()
