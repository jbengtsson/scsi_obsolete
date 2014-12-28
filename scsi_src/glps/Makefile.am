lib_LIBRARIES     = libglps.a

BUILT_SOURCES     = parser.h

libglps_a_SOURCES = glpsutil.cc elemtree.c glps.c glpserror.c \
                    scanner.l parser.y

LFLAGS   = -i
YFLAGS   = -t -v -d -y
CXXFLAGS = -g -O0 -ggdb3 -Wall -Wno-error=all -fPIC
#CXXFLAGS = -pg -O2 -Wall -fno-implicit-templates -fPIC
#CXXFLAGS = -g -O2 -Wall -frepo
