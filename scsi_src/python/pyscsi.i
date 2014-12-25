%module pyscsi
%{

#define NO 1

// Include the header files in the wrapper code.
#include "scsi_lib.h"

int no_tps   = NO;

// Friend declarations are ignored by SWIG.
void TPSAEps(const double);

// Not defined for linear TPSA.
int ndpt_tps = 2;

double *SWIG_vec;

%}

// Declared in "scsi_lib.h".

extern const int   nv_tps, nd_tps, iref_tps;

// Not defined in C code...
extern int         no_tps, ndpt_tps;
extern double      eps_tps;

extern ElemFamType ElemFam[];
extern CellType    Cell[];
extern globvalrec  globval;


%inline %{
  struct globv_vect {
    double *vec;
    // Python method for array access.
    double __getitem__(int k) { return vec[k]; };
  };
%}

%extend globvalrec {
  // Python method for array access.
  double __getitem__(int k) {
    return SWIG_vec[k];
  };

  // Python method for attribute access are:
  //  __getattr__(char *attr)
  //  __getattribute__(char *attr)

  globvalrec* gvec(char *attr) {
    if (strcmp(attr, "TotalTune") == 0)
      SWIG_vec = self->TotalTune;
    else if (strcmp(attr, "Chrom") == 0)
      SWIG_vec = self->Chrom;
    else if (strcmp(attr, "CODvect") == 0)
      SWIG_vec = &self->CODvect[0];
  }
}

%extend CellType {
  // Python method for array access.
  CellType* __getitem__(int k) { return self+k; };
#  CellType* __setitem__(int k) { return self+k; };
}


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

