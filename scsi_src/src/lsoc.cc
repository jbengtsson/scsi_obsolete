/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/


static bool               first_h[] = {true, true}, first_v[] = {true, true};
       int                n_bpm_[2], n_corr_[2];
       long unsigned int  *bpms_[2], *corrs_[2];
static double             *w_lsoc[2], **A_lsoc[2], **U_lsoc[2], **V_lsoc[2];
static double             *w_lstc[2], **A_lstc[2], **U_lstc[2], **V_lstc[2];

gsl_vector *vw_lsoc[2], *S;
gsl_matrix *mA_lsoc[2], *mU_lsoc[2], *mV_lsoc[2];
gsl_vector *vw_lstc[2];
gsl_matrix *mA_lstc[2], *mU_lstc[2], *mV_lstc[2];


void zero_trims(void)
{
  int       j, k;
  long int  loc;

  for (k = 0; k < 2; k++)
    for (j = 1; j <= n_corr_[k]; j++) {
      loc = corrs_[k][j];
      set_bn_design_elem(Cell[loc].Fnum, Cell[loc].Knum, Dip, 0.0, 0.0);
    }
}


void prt_gcmat(const int plane)
{
  int     i, j, k;
  FILE    *outf = NULL;

  k = plane - 1;

  printf("\n");
  printf("no of bpms = %d, no of corrs = %d, plane = %d\n",
	 n_bpm_[k], n_corr_[k], plane);

  if (plane == 1)
    outf = file_write("svdh.out");
  else if (plane == 2)
    outf = file_write("svdv.out");
  else {
    printf("prt_gcmat: undefined plane %d\n", plane);
    exit_(1);
  }

  fprintf(outf,"# total no of monitors:                %d\n", n_bpm_[k]);

  if (plane == 1)
    fprintf(outf,"# total no of horizontal correctors: %d\n", n_corr_[k]);
  else
    fprintf(outf,"# total no of vertical correctors:   %d\n", n_corr_[k]);

  fprintf(outf, "# A[%d][%d] = \n", n_bpm_[k], n_corr_[k]);
  for (i = 1; i <= n_bpm_[k]; i++) {
    for (j = 1; j <= n_corr_[k]; j++)
      fprintf(outf, "% .13e ", A_lsoc[k][i][j]);
    fprintf(outf, "\n");
  }

  fprintf(outf, "# U[%d][%d] = \n", n_bpm_[k], n_corr_[k]);
  for (i = 1; i <= n_bpm_[k]; i++) {
    for (j = 1; j <= n_corr_[k]; j++)
      fprintf(outf, "% .13e ", U_lsoc[k][i][j]);
    fprintf(outf, "\n");
  }

  fprintf(outf, "# w[%d]    = \n", n_corr_[k]);
  for (j = 1; j <= n_corr_[k]; j++)
    fprintf(outf, "% .13e ", w_lsoc[k][j]);
  fprintf(outf, "\n");
  fprintf(outf, "# V[%d][%d] = \n", n_bpm_[k], n_corr_[k]);

  for (i = 1; i <= n_corr_[k]; i++) {
    for (j = 1; j <= n_corr_[k]; j++)
      fprintf(outf, "% .13e ", V_lsoc[k][i][j]);
    fprintf(outf, "\n");
  }

  fclose(outf);
}


void gcmat(const int plane)
{
  /* Get orbit response matrix

                -----------
              \/beta  beta
                    i     j
        A   = ------------- cos(nu pi - 2 pi|nu  - nu |)
         ij   2 sin(pi nu)                     i     j

  */

  int       i, j, k;
  long int  loc;
  double    nu, betai, betaj, nui, nuj, spiq;

  const double  eps = 1e-10;

  k = plane - 1;

  nu = globval.TotalTune[k]; spiq = sin(M_PI*nu);

  for (i = 1; i <= n_bpm_[k]; i++) {
    loc = bpms_[k][i]; betai = Cell[loc].Beta[k]; nui = Cell[loc].Nu[k];
    for (j = 1; j <= n_corr_[k]; j++) {
      loc = corrs_[k][j]; betaj = Cell[loc].Beta[k]; nuj = Cell[loc].Nu[k];
      A_lsoc[k][i][j] =
	sqrt(betai*betaj)/(2.0*spiq)*cos(nu*M_PI-fabs(2.0*M_PI*(nui-nuj)));
    }
  }

  for (i = 1; i <= n_bpm_[k]; i++)
    for (j = 1; j <= n_corr_[k]; j++)
      U_lsoc[k][i][j] = A_lsoc[k][i][j];

  gsl_linalg_SV_decomp (mU_lsoc[k], mV_lsoc[k], S, vw_lsoc[k]);

  printf("\n");
  printf("gcmat singular values:\n");
  for (j = 1; j <= n_corr_[k]; j++) {
    printf("%11.3e", w_lsoc[k][j]);
    if (w_lsoc[k][j] < eps) {
      w_lsoc[k][j] = 0.0;
      printf(" (zeroed)");
    }
    if (j % 5 == 0) printf("\n");
  }
  if (n_corr_[k] % 5 != 0) printf("\n");

  if (trace) prt_gcmat(plane);
}


void gcmat(const int n_bpm, const long int bpms[],
	   const int n_corr, const long int corrs[], const int plane,
	   const bool svd)
{
  bool  first;
  int   i, k;

  k = plane - 1;

  first = (plane == 1)? first_h[0] : first_v[0];
  if (first) {
    if (plane == 1)
      first_h[0] = false;
    else
      first_v[0] = false;

    bpms_[k] = gslport_lvector(1, n_bpm);
    corrs_[k] = gslport_lvector(1, n_corr);

    mA_lsoc[k] = gsl_matrix_alloc(n_bpm, n_corr);
    GSL2NRDM2(dmA_lsoc, mA_lsoc[k], A_lsoc[k],0);
    mU_lsoc[k] = gsl_matrix_alloc(n_bpm, n_corr);	
    GSL2NRDM2(dmU_lsoc, mU_lsoc[k], U_lsoc[k],0);
    vw_lsoc[k] = gsl_vector_alloc(n_corr);
    GSL2NRDV2(vw_lsoc[k],w_lsoc[k]);
    mV_lsoc[k] = gsl_matrix_alloc(n_corr, n_corr);	
    GSL2NRDM2(dmAV_lsoc, mV_lsoc[k], V_lsoc[k],0);
	
    S = gsl_vector_alloc(n_corr);
  }

  for (i = 1; i <= n_bpm; i++)
    bpms_[k][i] = bpms[i-1];

  for (i = 1; i <= n_corr; i++)
    corrs_[k][i] = corrs[i-1];

  n_bpm_[k] = n_bpm; n_corr_[k] = n_corr;

  if (svd) gcmat(plane);
}


void gcmat(const int bpm, const int corr, const int plane)
{
  int  i, k;

  k = plane - 1; n_bpm_[k] = GetnKid(bpm); n_corr_[k] = GetnKid(corr);

  long int  bpms[n_bpm_[k]], corrs[n_corr_[k]];

  for (i = 1; i <= n_bpm_[k]; i++)
    bpms[i-1] = Elem_GetPos(bpm, i);

  for (i = 1; i <= n_corr_[k]; i++)
    corrs[i-1] = Elem_GetPos(corr, i);

  gcmat(n_bpm_[k], bpms, n_corr_[k], corrs, plane, true);
}


void lsoc(const int plane)
{
  int       j, k;
  long int  loc;
  double    *b, *x;

  gsl_vector *vb;
  gsl_vector *vx;

  k = plane - 1;

  vb = gsl_vector_alloc(n_bpm_[k]);
  vx = gsl_vector_alloc(n_corr_[k]);
  GSL2NRDV2(vb, b);
  GSL2NRDV2(vx, x);

  for (j = 1; j <= n_bpm_[k]; j++) {
    loc = bpms_[k][j];
    b[j] = -Cell[loc].BeamPos[2*k] + Cell[loc].dS[k];
  }
      
  gsl_vector *vs = gsl_vector_alloc(n_corr_[k]);
  gsl_vector *work = gsl_vector_alloc(n_corr_[k]);
  gsl_linalg_SV_decomp (mU_lsoc[k], mV_lsoc[k], vs, work);
  gsl_linalg_SV_solve(mU_lsoc[k], mV_lsoc[k],vs,vb,vx);
  gsl_vector_free(vs);
  gsl_vector_free(work);

  for (j = 1; j <= n_corr_[k]; j++) {
    loc = corrs_[k][j];
    if (plane == 1)
      set_dbnL_design_elem(Cell[loc].Fnum, Cell[loc].Knum, Dip, -x[j], 0.0);
    else
      set_dbnL_design_elem(Cell[loc].Fnum, Cell[loc].Knum, Dip, 0.0, x[j]);
  }

  gsl_vector_free(vb);
  gsl_vector_free(vx);
}


void gtcmat(const int plane)
{
  /* Get trajectory response matrix

                -----------
        A   = \/beta  beta  sin(2 pi(nu  - nu ))
         ij         i     j            i     j

  */

  int       i, j, k;
  long int  loc_bpm, loc_corr;
  double    betai, betaj, nui, nuj;

  const double  eps = 1e-10;

  k = plane - 1;

  for (i = 1; i <= n_bpm_[k]; i++) {
    loc_bpm = bpms_[k][i];
    betai = Cell[loc_bpm].Beta[k]; nui = Cell[loc_bpm].Nu[k];
    for (j = 1; j <= n_corr_[k]; j++) {
      loc_corr = corrs_[k][j];
      betaj = Cell[loc_corr].Beta[k]; nuj = Cell[loc_corr].Nu[k];
      if (loc_bpm > loc_corr)
	A_lstc[k][i][j] = sqrt(betai*betaj)*sin(2.0*M_PI*(nui-nuj));
      else
	A_lstc[k][i][j] = 0e0;
    }
  }

  for (i = 1; i <= n_bpm_[k]; i++)
    for (j = 1; j <= n_corr_[k]; j++)
      U_lstc[k][i][j] = A_lstc[k][i][j];

  gsl_linalg_SV_decomp (mU_lstc[k], mV_lstc[k], S, vw_lstc[k]);

  printf("\n");
  printf("gcmat singular values:\n");
  for (j = 1; j <= n_corr_[k]; j++) {
    printf("%11.3e", w_lstc[k][j]);
    if (w_lstc[k][j] < eps) {
      w_lstc[k][j] = 0.0;
      printf(" (zeroed)");
    }
    if (j % 5 == 0) printf("\n");
  }
  if (n_corr_[k] % 5 != 0) printf("\n");

  if (trace) prt_gcmat(plane);
}


void gtcmat(const int n_bpm, const long int bpms[],
	    const int n_corr, const long int corrs[], const int plane,
	    const bool svd)
{
  bool  first;
  int   i, k;

  k = plane - 1;

  first = (plane == 1)? first_h[1] : first_v[1];
  if (first) {
    if (plane == 1)
      first_h[1] = false;
    else
      first_v[1] = false;

    bpms_[k] = gslport_lvector(1, n_bpm);
    corrs_[k] = gslport_lvector(1, n_corr);

    mA_lstc[k] = gsl_matrix_alloc(n_bpm, n_corr);
    GSL2NRDM2(dmA_lstc, mA_lstc[k], A_lstc[k],0);
    mU_lstc[k] = gsl_matrix_alloc(n_bpm, n_corr);	
    GSL2NRDM2(dmU_lstc, mU_lstc[k], U_lstc[k],0);
    vw_lstc[k] = gsl_vector_alloc(n_corr);
    GSL2NRDV2(vw_lstc[k],w_lstc[k]);
    mV_lstc[k] = gsl_matrix_alloc(n_corr, n_corr);	
    GSL2NRDM2(dmAV_lstc, mV_lstc[k], V_lstc[k],0);
	
    S = gsl_vector_alloc(n_corr);
  }

  for (i = 1; i <= n_bpm; i++)
    bpms_[k][i] = bpms[i-1];

  for (i = 1; i <= n_corr; i++)
    corrs_[k][i] = corrs[i-1];

  n_bpm_[k] = n_bpm; n_corr_[k] = n_corr;

  if (svd) gtcmat(plane);
}


void lstc(const int plane, const long int lastpos)
{
  int       j, k;
  long int  loc;
  double    *b, *x;

  gsl_vector *vb;
  gsl_vector *vx;

  k = plane - 1;

  vb = gsl_vector_alloc(n_bpm_[k]);
  vx = gsl_vector_alloc(n_corr_[k]);
  GSL2NRDV2(vb, b);
  GSL2NRDV2(vx, x);

  for (j = 1; j <= n_bpm_[k]; j++) {
    loc = bpms_[k][j];
    if (loc < lastpos)
      b[j] = -Cell[loc].BeamPos[2*k] + Cell[loc].dS[k];
    else
      b[j] = 0e0;

    if (trace) cout << scientific << setprecision(5)
		    << "b[" << setw(3) << j << "] = "
		    << setw(12) << b[j] << endl;
  }

  gsl_vector *vs = gsl_vector_alloc(n_corr_[k]);
  gsl_vector *work = gsl_vector_alloc(n_corr_[k]);
  gsl_linalg_SV_decomp (mU_lstc[k], mV_lstc[k], vs, work);
  gsl_linalg_SV_solve(mU_lstc[k], mV_lstc[k],vs,vb,vx);
  gsl_vector_free(vs);
  gsl_vector_free(work);

  for (j = 1; j <= n_corr_[k]; j++) {
    loc = corrs_[k][j];
    if (plane == 1) {
      if (trace) cout << scientific << setprecision(5)
		      << "(b_1L)[" << setw(3) << j << "] = "
		      << setw(12)<< -x[j] << endl;

      set_dbnL_design_elem(Cell[loc].Fnum, Cell[loc].Knum, Dip,
			   -x[j], 0.0);
    } else {
      if (trace) cout << scientific << setprecision(5)
		      << "(a_1L)[" << setw(3) << j << "] = "
		      << setw(12)<< x[j] << endl;

      set_dbnL_design_elem(Cell[loc].Fnum, Cell[loc].Knum, Dip,
			   0.0, x[j]);
    }
  }

  gsl_vector_free(vb);
  gsl_vector_free(vx);
}
