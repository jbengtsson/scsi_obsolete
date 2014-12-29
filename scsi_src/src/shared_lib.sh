#!/bin/sh

gcc -shared -fPIC -I../inc -I/usr scsi_lib.cc -o libscsi.so

cp libscsi.so ../lib/.
