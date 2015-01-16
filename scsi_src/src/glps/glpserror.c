#include <stdio.h>
#include <stdlib.h>

#include "glpserror.h"


void
glps_lat_error (const char * reason, const char * file, int line,
                const char * src, int glps_errno)
{
    abort();
}

const char *
glps_strerror (const int glps_errno)
{
  switch (glps_errno) {
  case GLPS_SUCCESS:
      return "success" ;
  case GLPS_FAILURE:
      return "failure" ;
  case GLPS_CONTINUE:
      return "the iteration has not converged yet" ;
  case GLPS_EEOF:
      return "end of file" ;
  case GLPS_ENOLAT:
      return "lattice file is not found" ;
  case GLPS_ESYNTAX:
      return "lattice syntax error" ;
  case GLPS_EBLNODEF:
      return "beamline not defined" ;
  case GLPS_ENOBL:
      return "no beamline defined in lattice file";
  case GLPS_EFLAT:
      return "can not write flat output file";
  case GLPS_EXML:
      return "can not write XML file as output";
  case GLPS_EERR:
      return "unknown general error" ;
  default:
      return "unknown error code" ;
    }
}
