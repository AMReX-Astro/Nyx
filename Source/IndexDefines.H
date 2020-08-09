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
  FirstSpec = 6

};

enum State_Type_UIndex {
  URHO     = 0,
  UMX      = 1,
  UMY      = 2,
  UMZ      = 3,
  UEDEN    = 4,
  UEINT    = 5,

};

enum State_Type_QIndex {
  QRHO     = 0,
  QU       = 1,
  QV       = 2,
  QW       = 3,
  QPRES    = 4,
  QREINT   = 5

};

enum DiagEOS_Type_Index {
  Temp_comp = 0,
  Ne_comp   = 1,
  Zhi_comp  = 2

};

#define GDRHO 0
#define GDU 1
#define GDV 2
#define GDW 3
#define GDPRES 4
#define GDGAME 5
#define NGDNV 6

#define NUM_GROW 4

#endif
