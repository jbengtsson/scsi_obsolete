#!/bin/sh

time=`date --iso-8601`

dir=`pwd`

tracy="scsi"

# set home directory
cd "$dir/$tracy"

# clean up directory
make distclean

rm -rf autom4te.cache
rm -rf aclocal.m4
rm -rf bin/*
rm -rf tracy/bin/*

cd $dir

fname="back_up_$tracy.$time.tgz"
cmd_str="tar -czf $fname "back_up_$tracy.sh" "make_$tracy.sh" $tracy"
echo "$cmd_str"
`$cmd_str`
