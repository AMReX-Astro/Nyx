#ifndef _INDEX_DEFINES_H_
#define _INDEX_DEFINES_H_

#include "Nyx.H"
#include <AMReX_REAL.H>
#include <AMReX_Arena.H>

#define NUM_SPECIES 2

#define UFA 6

#define QFA 8

#define GDRHO 0
#define GDU 1
#define GDV 2
#define GDW 3
#define GDPRES 4
#define GDGAME 5

#ifndef PELEC_USE_AUX
#define NUM_ADV 0
#define NUM_AUX 0
#endif

#define UFS UFA + NUM_ADV
#define UFX UFS + NUM_SPECIES

#define QFS QFA + NUM_ADV
#define QFX QFS + NUM_SPECIES
#define NGDNV 6

#define NPASSIVE NUM_SPECIES + NUM_ADV + NUM_AUX + 3 - AMREX_SPACEDIM

// This value was hardcoded to 4 in get_method_params in the old Fortran code
#define NUM_GROW 4

#endif
