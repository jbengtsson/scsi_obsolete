#!/bin/sh

gcc -DHAVE_CONFIG_H -I. -I../.. -I../inc -I/usr -g -O0 -Wall -Wno-error=all -fno-implicit-templates -fPIC -MT scsi_lib.o -MD -MP -MF .deps/scsi_lib.Tpo -c -o libscsi.so scsi_lib.cc -shared
