%module pyscsi
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
    else if (strcmp(attr, "wr") == 0)
      v.a = &self->wr[0];
    else if (strcmp(attr, "wi") == 0)
      v.a = &self->wi[0];
    else if (strcmp(attr, "alpha_rad") == 0)
      v.a = &self->alpha_rad[0];
    else if (strcmp(attr, "D_rad") == 0)
      v.a = &self->D_rad[0];
    else if (strcmp(attr, "J") == 0)
      v.a = &self->J[0];
    else if (strcmp(attr, "tau") == 0)
      v.a = &self->tau[0];
    else if (strcmp(attr, "D_IBS") == 0)
      v.a = &self->D_IBS[0];
    else if (strcmp(attr, "eps") == 0)
      v.a = &self->eps[0];
    else if (strcmp(attr, "epsp") == 0)
      v.a = &self->epsp[0];
    return v;
  }

  mat gmat(const char *attr) {
    mat m;
    if (strcmp(attr, "OneTurnMat") == 0)
      m.a = &self->OneTurnMat;
    else if (strcmp(attr, "Ascr") == 0)
      m.a = &self->Ascr;
    else if (strcmp(attr, "Ascrinv") == 0)
      m.a = &self->Ascrinv;
    else if (strcmp(attr, "Vr") == 0)
      m.a = &self->Vr;
    else if (strcmp(attr, "Vi") == 0)
      m.a = &self->Vi;
    return m;
  }
}

%extend CellType {
  // Python method for array access.
  CellType* __getitem__(int k) { return self+k; };
#  CellType* __setitem__(int k) { return self+k; };

  // Python method for attribute access are:
  //  __getattr__(char *attr)
  //  __getattribute__(char *attr)
  vec gvec(const char *attr) {
    vec v;
    if (strcmp(attr, "dS") == 0)
      v.a = &self->dS[0];
    else if (strcmp(attr, "dT") == 0)
      v.a = &self->dT[0];
    else if (strcmp(attr, "Nu") == 0)
      v.a = &self->Nu[0];
    else if (strcmp(attr, "Alpha") == 0)
      v.a = &self->Alpha[0];
    else if (strcmp(attr, "Beta") == 0)
      v.a = &self->Beta[0];
    else if (strcmp(attr, "Eta") == 0)
      v.a = &self->Eta[0];
    else if (strcmp(attr, "Etap") == 0)
      v.a = &self->Etap[0];
    else if (strcmp(attr, "BeamPos") == 0)
      v.a = &self->BeamPos[0];
    return v;
  }

  mat gmat(const char *attr) {
    mat m;
    if (strcmp(attr, "A") == 0)
      m.a = &self->A;
    else if (strcmp(attr, "sigma") == 0)
      m.a = &self->sigma;
    return m;
  }
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
