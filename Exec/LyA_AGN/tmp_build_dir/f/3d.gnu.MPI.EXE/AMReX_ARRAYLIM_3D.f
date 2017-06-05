

c ::: -----------------------------------------------------------
c ::: This routine sets the values for the lo() and hi() arrays
c ::: from the ARG_L1, ARG_H1, ... macros.  This is done since
c ::: it is more convenient to use the lo() and hi() arrays.
c :::
c ::: INPUTS/OUTPUTS:
c :::
c ::: holder_l1, holder_l2, holder_l3, holder_h1, holder_h2, holder_h3=> index extent of place holder array
c ::: lo(3)    <= lower index limits
c ::: hi(3)    <= upper index limits
c ::: -----------------------------------------------------------
      subroutine SET_LOHI(holder_l1, holder_l2, holder_l3, holder_h1, ho
     &lder_h2, holder_h3, lo, hi)

      implicit none

c
c     :::: Passed Variables ::::
c
      integer holder_l1, holder_l2, holder_l3, holder_h1, holder_h2, hol
     &der_h3
      integer lo(3), hi(3)

c
c     --------------------------------------
c     :::: Set Values for lo() and hi() ::::
c     --------------------------------------
c
      lo(1) = holder_l1
      hi(1) = holder_h1
      lo(2) = holder_l2
      hi(2) = holder_h2
      lo(3) = holder_l3
      hi(3) = holder_h3

      return
      end

c ::: -----------------------------------------------------------
c ::: This routine sets the values for the ARG_L1, ARG_H1, ... macros
c ::: from the lo() and hi() arrays.  This is done since
c ::: it is more convenient to use the macros to dimension arrays.
c :::
c ::: INPUTS/OUTPUTS:
c :::
c ::: FF_DIMS(holder) <=  index extent of place holder array
c ::: lo(3)         => lower index limits
c ::: hi(3)         => upper index limits
c ::: -----------------------------------------------------------
      subroutine SET_ARGS(holder_l1, holder_l2, holder_l3, holder_h1, ho
     &lder_h2, holder_h3, lo, hi)

      implicit none

c
c     :::: Passed Variables ::::
c
      integer holder_l1, holder_l2, holder_l3, holder_h1, holder_h2, hol
     &der_h3
      integer lo(3), hi(3)

c
c     --------------------------------------
c     :::: Set Values for lo() and hi() ::::
c     --------------------------------------
c
      holder_l1 = lo(1)
      holder_h1 = hi(1)
      holder_l2 = lo(2)
      holder_h2 = hi(2)
      holder_l3 = lo(3)
      holder_h3 = hi(3)

c
c
      return
      end

