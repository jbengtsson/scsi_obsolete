# Shared library required for pytracy.
lib_LTLIBRARIES        = libscsi.la libscsi_ptc.la 
libscsi_la_SOURCES     = scsi_lib.cc \
			 ../glps/glpsutil.cc ../glps/elemtree.c ../glps/glps.c \
			 ../glps/glpserror.c \
			 ../glps/scanner.l ../glps/parser.y
libscsi_ptc_la_SOURCES = scsi_ptc_lib.cc \
			 ../glps/glpsutil.cc ../glps/elemtree.c ../glps/glps.c \
			 ../glps/glpserror.c \
			 ../glps/scanner.l ../glps/parser.y

libscsi_la_LDFLAGS     = -module -avoid-version -shared
libscsi_ptc_la_LDFLAGS = -module -avoid-version -shared

# Include GSL directly.
LIBS = -lgsl -lgslcblas

INCLUDES = -I../inc -I/usr -I../glps

# Ignore case in patterns.
LFLAGS   = -i
# Flags: debug, verbose, produce header file, and emulate POSIX Yacc.
YFLAGS   = -t -v -d -y
# Flags: generate Position Independenct Code (PIC); for shared library.
CXXFLAGS = -g -O0 -Wall -Wno-error=all -fPIC

#CXXFLAGS = -g -O0 -Wall -Wno-error=all -fno-implicit-templates -fPIC
#CXXFLAGS = -pg -O2 -Wall -fno-implicit-templates -fPIC
#CXXFLAGS = -g -O2 -Wall -frepo
