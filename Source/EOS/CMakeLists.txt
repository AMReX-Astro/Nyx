target_include_directories(nyxcore
   PUBLIC
   $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}> )


if (Nyx_HEATCOOL)
   target_sources( nyxcore
      PRIVATE
      eos_hc.H
      atomic_rates.H
      atomic_rates_data.H
      atomic_rates_data.cpp
      )
else()
   target_sources( nyxcore
      PRIVATE
      eos_adiabatic.H
      )
endif()

target_sources( nyxcore PRIVATE reset_internal_e.H )
