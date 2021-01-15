#!/bin/bash
set -e # Exit with nonzero exit code if anything fails

# Doxygen
echo "Build the Doxygen documentation"
cd Docs/sphinx_documentation
doxygen Doxyfile
cd ../..

# sphinx
cd Docs/sphinx_documentation
echo "Build the Sphinx documentation for Nyx."
#Skip latexpdf since it fails
#make PYTHON="python3" latexpdf
make PYTHON="python3" html &> make_source_html.out
cd ../../

mkdir build
cd build
mkdir docs_html

# add sphinx
cp -rp ../Docs/sphinx_documentation/build/html/* docs_html/
