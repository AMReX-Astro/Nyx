#!/bin/bash
make -j
/bin/rm -rf plt* out*
./Nyx3d.gnu.PROF.ex* inputs max_step=$1 | tee out.1

xstart=1;xend=$1;xstep=1
for (( x = $xstart; x <= $xend; x += $xstep));
do
    mv plt0000$x plt1_0000$x
done
./Nyx3d.gnu.PROF.ex* inputs max_step=$1 | tee out.0

#xstart=1;xend=9;xstep=1
for (( x = $xstart; x <= $xend; x += $xstep));
do
    mv plt0000$x plt0_0000$x
    ~/amrex/Tools/Postprocessing/F_Src/fcompare.Linux.gfortran.exe -z rho_e plt0_0000$x plt1_0000$x
done

meld out.0 out.1 &
