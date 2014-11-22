lib_LIBRARIES         = libscsi.a libscsi_ptc.a

libscsi_a_SOURCES     = scsi_lib.cc

libscsi_ptc_a_SOURCES = scsi_ptc_lib.cc

INCLUDES = -I../inc -I$(NUM_REC)/inc
#INCLUDES = -I../inc -I$(NUM_REC)/inc \
#	   -I/usr/include/c++/3.3.6

CXXFLAGS = -g -O0 -Wall -Wno-error=all -fno-implicit-templates -fPIC
#CXXFLAGS = -pg -O2 -Wall -fno-implicit-templates -fPIC
#CXXFLAGS = -g -O2 -Wall -frepo
