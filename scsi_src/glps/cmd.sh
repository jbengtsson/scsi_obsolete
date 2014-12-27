#!/bin/sh

flex -o scanner.c scanner.l
bison -t -v -d -y -o parser.c parser.y

gcc glpslat.cc glpsutil.cc elemtree.c glps.c glpserror.c scanner.c parser.c \
    -lstdc++ -lfl -o glpslat