#!/bin/sh

swig -c++ -python -globals gv pyscsi.i

gcc -shared -fPIC -I../inc -I/usr -I/usr/include/python2.6 \
    -L../lib -lscsi -lglps -lgsl -lstdc++ -lgslcblas pyscsi_wrap.cxx \
    -o _pyscsi.so

cp _pyscsi.so ../lib/.
