# Example diagnostic/analysis tool

*****

## About

This is a simple example of opening an AMReX `plt` file, reading and manipulating data in it.
It is meant to be a starting place for building analysis tools for user's specific needs.
In this particular example, given a N-dimensional data set, the program computes a one-dimensional 
profile in direction *dir* where each value in the profile is the average of the data
in the non-*dir* directions.   This can work with multiple components at one time.

## Usage

This is a stand-alone application built on top of AMReX.  Set the compiler value in the `GNUmakefile` (variable `COMP`)
and make an executable with:
```sh
$ make -j 8
```

To run, you can execute with runtime parameters in the `inputs` file
```sh
$ AmrDerive3d.Linux.gcc.gfortran.MPI.ex inputs
```
where the file `inputs` contains, for example:
```sh
infile  = plt00000
dir     = 0
nComp   = 1
sComp   = 0
```

or you can specify runtime parameters via command-line:
```sh
AmrDerive3d.Linux.gcc.gfortran.MPI.ex infile=plt00000 dir=0 nComp=2 sComp=0
```

You can also set verbose = 1 (via command line or inputs file) in order to get a quite verbose output.
