%module python_api
%{

// C standard library
#include <stdio.h>
#include <stddef.h>
#include <setjmp.h>
#include <time.h>
#include <memory.h>
#include <malloc.h>
//#include <execinfo.h>


// C++ standard library
#include <cstdlib>
#include <cfloat>
#include <cctype>
#include <cmath>
#include <complex>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>
using namespace std;

// GSL
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_multimin.h>
#include "gslport.h"

#include "../inc/gslport.h"

#define NO 1

// Friend declarations are ignored by SWIG.
void TPSAEps(const double);

#include "../inc/field.h"
#include "../inc/mathlib.h"

#include "../inc/tpsa_lin.h"
#include "../inc/tpsa_lin_pm.h"

#include "../inc/scsi.h"
#include "../inc/scsi_global.h"
#include "../inc/ety.h"
#include "../inc/eigenv.h"

#include "../inc/radia2scsi.h"
#include "../inc/pascalio.h"

#include "../inc/t2elem.h"
#include "../inc/t2cell.h"
#include "../inc/t2lat.h"
#include "../inc/t2ring.h"

#include "../inc/fft.h"

#include "../inc/physlib.h"
#include "../inc/nsls-ii_lib.h"

#include "../inc/lsoc.h"

#include "../inc/modnaff.h"

#include "../inc/naffutils.h"
#include "../inc/complexeheader_naff.h"

#include "../inc/soleillib.h"

#include "../inc/rdmfile.h"
#include "../inc/prtmfile.h"


extern const int  nv_tps, nd_tps, iref_tps;
extern int        no_tps, ndpt_tps;
extern double     eps_tps;

extern ElemFamType ElemFam[];

extern CellType Cell[];

extern globvalrec globval;

%}


%include "../inc/gslport.h"

%include "../inc/field.h"
%include "../inc/mathlib.h"

%include "../inc/tpsa_lin.h"
%include "../inc/tpsa_lin_pm.h"

%include "../inc/scsi.h"
%include "../inc/scsi_global.h"
%include "../inc/ety.h"
%include "../inc/eigenv.h"

%include "../inc/radia2scsi.h"
%include "../inc/pascalio.h"

%include "../inc/t2elem.h"
%include "../inc/t2cell.h"
%include "../inc/t2lat.h"
%include "../inc/t2ring.h"

%include "../inc/fft.h"

%include "../inc/physlib.h"
%include "../inc/nsls-ii_lib.h"

%include "../inc/lsoc.h"

%include "../inc/modnaff.h"

%include "../inc/naffutils.h"
%include "../inc/complexeheader_naff.h"

%include "../inc/soleillib.h"

%include "../inc/rdmfile.h"
%include "../inc/prtmfile.h"


export const int  nv_tps, nd_tps, iref_tps;
export int        no_tps, ndpt_tps;
export double     eps_tps;

export ElemFamType ElemFam[];

export CellType Cell[];

export globvalrec globval;
