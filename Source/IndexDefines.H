#ifndef _INDEX_DEFINES_H_
#define _INDEX_DEFINES_H_

#include "Nyx.H"
#include <AMReX_REAL.H>
#include <AMReX_Arena.H>

#define GDRHO 0
#define GDU 1
#define GDV 2
#define GDW 3
#define GDPRES 4
#define GDGAME 5

#define NGDNV 6

// This value was hardcoded to 4 in get_method_params in the old Fortran code
#define NUM_GROW 4

#endif
