
 .. role:: cpp(code)
    :language: c++

.. _InSitu:

In situ analysis
=============

Halo finding
------

To find halos in situ while Nyx is running we use the Reeber package.
To compile with Reeber, add in GNUmakefile

REEBER = TRUE

and set the location of Boost library

BOOST_INLUDE_DIR := ${BOOST_ROOT}/include/boost

Note that these codes are in separate repositories and are not included with Nyx.

