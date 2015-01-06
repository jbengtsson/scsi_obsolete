# Shared library required for pytracy.
lib_LTLIBRARIES   = libglps.a libglps.la

BUILT_SOURCES     = parser.h

libglps_a_SOURCES  = glpsutil.cc elemtree.c glps.c glpserror.c \
                     scanner.l parser.y
libglps_la_SOURCES = glpsutil.cc elemtree.c glps.c glpserror.c \
                     scanner.l parser.y

libglps_la_LDFLAGS = -module -avoid-version -shared

LFLAGS   = -i # Ignore case in patterns.
YFLAGS   = -t -v -d -y
CXXFLAGS = -g -O0 -ggdb3 -Wall -Wno-error=all -fPIC
#CXXFLAGS = -pg -O2 -Wall -fno-implicit-templates -fPIC
#CXXFLAGS = -g -O2 -Wall -frepo
