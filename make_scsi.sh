#!/bin/sh

dir=`pwd`

cd scsi_src/python
# Chang name of cvar object for global variables to gv.
swig -c++ -python -globals gv pyscsi.i
cd $dir

# Remove autoconf cashe.
rm -rf autom4te.cache
rm -rf aclocal.m4

mkdir -p config

# Configure libtool (for shared libraries).
libtoolize

./bootstrap
./configure --prefix=$dir/scsi_src

make install
