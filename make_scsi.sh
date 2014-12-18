#!/bin/sh

dir=`pwd`

rm -rf autom4te.cache
rm -rf aclocal.m4

mkdir config

./bootstrap
./configure --prefix=$dir/scsi_src

make install
