#!/bin/sh

dir=`pwd`

cd scsi_src/python
swig -c++ -python python_api.i
cd $dir

rm -rf autom4te.cache
rm -rf aclocal.m4

mkdir -p config

./bootstrap
./configure --prefix=$dir/scsi_src

make install
