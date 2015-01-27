/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -

*/

#include "scsi_lib.h"

#include "gslport.cc"

#include "field.cc"

#if NO == 1
  // linear TPSA
  #include "tpsa_lin.cc"
  #include "tpsa_lin_pm.cc"
#else
  // interface to M. Berz' TPSA
  #include "tpsa_for_pm.cc"
#endif

#include "mathlib.cc"

#include "ety.cc"
#include "eigenv.cc"

#if false
  #include "t2lat.cc"
#else
  #include "t3lat.cc"
#endif
#include "t2elem.cc"
#include "t2cell.cc"
#include "t2ring.cc"

#include "pascalio.cc"

#include "lsoc.cc"

#include "prtmfile.cc"
#include "rdmfile.cc"

#include "fft.cc"

#include "physlib.cc"

#include "naffutils.cc"

#include "modnaff.cc"
#include "radia2scsi.cc"
#include "soleillib.cc"

#include "nsls-ii_lib.cc"
#include "usr_lib.cc"


// Truncated Power Series Algebra (TPSA)
const int no_tps = NO,       // Order.
          nv_tps   = ss_dim, // Number of variables.
          nd_tps   = 3,      // Number of degrees of freedom.
          iref_tps = 0;      // File with resonances to be excluded from
			     // the map normal form: fort.7.

double eps_tps  = 1e-25;     // Floating point truncation level.


// Instantiate templates.

template class ss_vect<double>;

template class ss_vect<tps>;


template void LinearInterpolation2(double &, double &, double &, double &, double &,
				   CellType &, bool &, int);

template void LinearInterpolation2(tps &, tps &, tps &, tps &, tps &,
				   CellType &, bool &, int);

template void SplineInterpolation2(double &, double &, double &, double &,
				   CellType &, bool &);

template void SplineInterpolation2(tps &, tps &, tps &, tps &,
				   CellType &, bool &);

template void spline(const double [], const double [], int const,
		     double const, const double, double []);

template void spline(const double [], const tps [], int const,
		     double const, const double, tps []);

template void splint(const double[], const double [], const double [],
		     const int, const double &, double &);

template void splint(const double[], const double [], const double [],
		     const int, const tps &, tps &);

template void splint(const double[], const tps [], const tps [],
		     const int, const tps &, tps &);

template void splin2(const double [], const double [],
		     double **, double **, const int, const int,
		     const double &, const double &, double &);

template void splin2(const double [], const double [],
		     double **, double **, const int, const int,
		     const tps &, const tps &, tps &);


template void Elem_Pass(const long, ss_vect<double> &);

template void Elem_Pass(const long, ss_vect<tps> &);


template void Cell_Pass(const long, const long, ss_vect<double> &, long &);

template void Cell_Pass(const long, const long, ss_vect<tps> &, long &);


/* Global variable used through the code */
globvalrec globval;

statusrec status;
bool trace, traceID;
bool cellconcat;

/* Random stuff */
long rseed0, rseed;
double normcut_;

double d_sign(double a, double b)
{
  double x;

  x = (a >= 0 ? a : - a);
  return( b >= 0 ? x : -x);
}

int P_eof(FILE *f)
{
  register int ch;

  if (feof(f)) return 1;
  if (f == stdin) return 0; /* not safe to look-ahead on the keyboard! */
  ch = getc(f);
  if (ch == EOF) return 1;
  ungetc(ch, f);

  return 0;
}


/* Check if at end of line (or end of entire file). */

int P_eoln(FILE *f)
{
  register int ch;

  ch = getc(f);
  if (ch == EOF) return 1;
  ungetc(ch, f);
  return (ch == '\n');
}
