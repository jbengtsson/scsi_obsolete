# Shared library.
lib_LTLIBRARIES   = _pyscsi.la

_pyscsi_la_SOURCES = pyscsi_wrap.cxx

_pyscsi_la_LDFLAGS = -module -avoid-version -shared

LIBS     = -L../src -lscsi \
           -L../glps -lglps \
	   -lgsl -lstdc++ -lgslcblas pyscsi_wrap.cxx

INCLUDES = -I../inc -I/usr -I/usr/include/python2.6

CXXFLAGS = -g -O0 -Wall -Wno-error=all -fPIC
#LDFLAGS = -shared -fPIC
#CXXFLAGS = -g -O0 -Wall -Wno-error=all -fno-implicit-templates -fPIC
#CXXFLAGS = -pg -O2 -Wall -fno-implicit-templates -fPIC
#CXXFLAGS = -g -O2 -Wall -frepo
