# Needs to be made after library has been made.
# Test program.
#bin_PROGRAMS = tst
#main_SOURCES = tst.cc

#LIBS = $(HOME)/git_repos/scsi/scsi_src/src/libscsi.a \
#       $(HOME)/git_repos/scsi/scsi_src/glps/libglps.a \
#       -lstdc++ \
#       -lgsl -lgslcblas \
#       -lfl

# Shared library required for pytracy.
lib_LTLIBRARIES        = libscsi.la libscsi_ptc.la 
libscsi_la_SOURCES      = scsi_lib.cc \
			 ../glps/glpsutil.cc ../glps/elemtree.c ../glps/glps.c \
			 ../glps/glpserror.c \
			 ../glps/scanner.l ../glps/parser.y
libscsi_ptc_la_SOURCES = scsi_ptc_lib.cc \
			 ../glps/glpsutil.cc ../glps/elemtree.c ../glps/glps.c \
			 ../glps/glpserror.c \
			 ../glps/scanner.l ../glps/parser.y

libscsi_la_LDFLAGS     = -module -avoid-version -shared
libscsi_ptc_la_LDFLAGS = -module -avoid-version -shared

INCLUDES = -I../inc -I/usr -I../glps

LFLAGS   = -i # Ignore case in patterns.
YFLAGS   = -t -v -d -y
CXXFLAGS = -g -O0 -Wall -Wno-error=all -fPIC
#CXXFLAGS = -g -O0 -Wall -Wno-error=all -fno-implicit-templates -fPIC
#CXXFLAGS = -pg -O2 -Wall -fno-implicit-templates -fPIC
#CXXFLAGS = -g -O2 -Wall -frepo
