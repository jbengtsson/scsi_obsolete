%module pyscsi
%{

#define NO 1

#include "scsi_lib.h"

int no_tps   = NO;

// Friend declarations are ignored by SWIG.
void TPSAEps(const double);

// Not defined for linear TPSA.
int ndpt_tps = 2;

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
