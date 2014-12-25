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
  struct vec {
    double *a;
    // Python method for array access.
    double __getitem__(int k) { return *(a+k); };
  };
%}

%inline %{
  struct mat {
    Matrix *a;
    // Python method for array access.
    vec __getitem__(int k) { 
      vec v; v.a = &(*a)[k][0];
      return v;
    };
  };
%}

%extend globvalrec {
  // Python method for attribute access are:
  //  __getattr__(char *attr)
  //  __getattribute__(char *attr)
  vec gvec(const char *attr) {
    vec v;
    if (strcmp(attr, "TotalTune") == 0)
      v.a = &self->TotalTune[0];
    else if (strcmp(attr, "Chrom") == 0)
      v.a = &self->Chrom[0];
    else if (strcmp(attr, "CODvect") == 0)
      v.a = &self->CODvect[0];
    return v;
  }

  mat gmat(const char *attr) {
    mat m;
    if (strcmp(attr, "OneTurnMat") == 0)
      m.a = &self->OneTurnMat;
    else if (strcmp(attr, "Ascr") == 0)
      m.a = &self->Ascr;
    return m;
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

