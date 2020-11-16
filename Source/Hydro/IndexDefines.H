#ifndef _INDEX_DEFINES_H_
#define _INDEX_DEFINES_H_

#include <AMReX_REAL.H>
#include <AMReX_Arena.H>

enum State_Type_Index {
  Density   = 0,
  Xmom      = 1,
  Ymom      = 2,
  Zmom      = 3,
  Eden      = 4,
  Eint      = 5,
#ifndef CONST_SPEC
  FirstSpec = 6
#endif
};

enum State_Type_UIndex {
  URHO     = 0,
  UMX      = 1,
  UMY      = 2,
  UMZ      = 3,
  UEDEN    = 4,
  UEINT    = 5 
};

#define QRHO   0
#define QU     1
#define QV     2
#define QW     3
#define QPRES  4
#define QREINT 5
#define NGDNV  6 // Note this never includes the species 

#ifdef CONST_SPECIES
#define QVAR   6 // This should be equal to the number of components of U
#else
#define QVAR   8 // This should be equal to the number of components of U
#endif

enum DiagEOS_Type_Index {
  Temp_comp = 0,
  Ne_comp   = 1,
  Zhi_comp  = 2
};

#define NUM_GROW 4


#endif
