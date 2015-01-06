# Test program.
bin_PROGRAMS = tst
main_SOURCES = tst.cc

LIBS = $(HOME)/git_repos/scsi/scsi_src/lib/libscsi.a \
       $(HOME)/git_repos/scsi/scsi_src/lib/libglps.a \
       -lstdc++ \
       -lgsl -lgslcblas \
       -lfl

# Shared library required for pytracy.
lib_LTLIBRARIES        = libscsi.a libscsi_ptc.a libscsi.la libscsi_ptc.la 
libscsi_a_SOURCES      = scsi_lib.cc
libscsi_la_SOURCES     = scsi_lib.cc
libscsi_ptc_a_SOURCES  = scsi_ptc_lib.cc
libscsi_ptc_la_SOURCES = scsi_ptc_lib.cc

libscsi_la_LDFLAGS     = -module -avoid-version -shared
libscsi_ptc_la_LDFLAGS = -module -avoid-version -shared

INCLUDES = -I../inc -I/usr

CXXFLAGS = -g -O0 -Wall -Wno-error=all -fPIC
#CXXFLAGS = -g -O0 -Wall -Wno-error=all -fno-implicit-templates -fPIC
#CXXFLAGS = -pg -O2 -Wall -fno-implicit-templates -fPIC
#CXXFLAGS = -g -O2 -Wall -frepo
