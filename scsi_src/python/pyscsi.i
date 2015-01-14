%module pyscsi

%include "std_string.i"

%{

#define PYAPI 1
#define NO    1

// Include the header files in the wrapper code.
#include "scsi_lib.h"

// Friend declarations are ignored by SWIG.
void TPSAEps(const double);

// Not defined for linear TPSA.
int ndpt_tps = 2;

%}

#define PYAPI 1
#define NO    1

// Declared in "scsi_lib.h".

extern const int   nv_tps, nd_tps, iref_tps;

// Not defined in the C code...
extern int         ndpt_tps;
extern double      eps_tps;

extern ElemFamType ElemFam[];
extern CellType    Cell[];
extern globvalrec  globval;


// Parse the header files to generate the wrapper code.

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
//%include "../inc/t2lat.h"
%include "../inc/t3lat.h"
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

%include "typemaps.i"

%apply double *OUTPUT { double &bn, double &an };
%inline %{
  extern void get_bn_design_elem(const int Fnum, const int Knum,
				 const int n, double &bn, double &an);
%}
