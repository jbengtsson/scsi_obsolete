#!/bin/sh

g++ -shared -fPIC -I../inc -I/usr -I/usr/include/python2.6 -L../lib -lscsi -lstdc++ -lgsl -lgslcblas pyscsi_wrap.cxx -o _pyscsi.so
