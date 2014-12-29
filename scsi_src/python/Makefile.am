lib_LIBRARIES     = _pyscsi.a
#lib_LTLIBRARIES   = _pyscsi.la

_pyscsi_a_SOURCES = pyscsi_wrap.cxx

LIBS     = -L$(HOME)/git_repos/scsi/scsi_src/lib -lscsi -lglps

INCLUDES = -I../inc -I/usr -I/usr/include/python2.6

CXXFLAGS = -g -O0 -Wall -Wno-error=all -fPIC
#CXXFLAGS = -g -O0 -Wall -Wno-error=all -fPIC -shared
#LDFLAGS = -shared -fPIC
#CXXFLAGS = -g -O0 -Wall -Wno-error=all -fno-implicit-templates -fPIC
#CXXFLAGS = -pg -O2 -Wall -fno-implicit-templates -fPIC
#CXXFLAGS = -g -O2 -Wall -frepo
