/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -        

*/


bool               first_h = true, first_v = true;
int                n_bpm_[2], n_corr_[2];
long unsigned int  *bpms_[2], *corrs_[2];
double             *w_lsoc[2], **A_lsoc[2], **U_lsoc[2], **V_lsoc[2];

//GSL add
gsl_vector		   *vw_lsoc[2], *S;
gsl_matrix         *mA_lsoc[2], *mU_lsoc[2], *mV_lsoc[2];
//end GSL add


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
      fprintf(outf, "% .3e ", A_lsoc[k][i][j]);
    fprintf(outf, "\n");
  }

  fprintf(outf, "# U[%d][%d] = \n", n_bpm_[k], n_corr_[k]);
  for (i = 1; i <= n_bpm_[k]; i++) {
    for (j = 1; j <= n_corr_[k]; j++)
      fprintf(outf, "% .3e ", U_lsoc[k][i][j]);
    fprintf(outf, "\n");
  }

  fprintf(outf, "# w[%d]    = \n", n_corr_[k]);
  for (j = 1; j <= n_corr_[k]; j++)
    fprintf(outf, "% .3e ", w_lsoc[k][j]);
  fprintf(outf, "\n");
  fprintf(outf, "# V[%d][%d] = \n", n_bpm_[k], n_corr_[k]);

  for (i = 1; i <= n_corr_[k]; i++) {
    for (j = 1; j <= n_corr_[k]; j++)
      fprintf(outf, "% .3e ", V_lsoc[k][i][j]);
    fprintf(outf, "\n");
  }

  fclose(outf);
}


void gcmat(const int plane)
{
  /* Get correlation matrix

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

  /*Orginal
  dsvdcmp(U_lsoc[k], n_bpm_[k], n_corr_[k], w_lsoc[k], V_lsoc[k]);
  */
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

  first = (plane == 1)? first_h : first_v;
  if (first) {
    if (plane == 1)
      first_h = false;
    else
      first_v = false;

	/*Orginal
    bpms_[k] = lvector(1, n_bpm); corrs_[k] = lvector(1, n_corr);
	*/
	
	//GSL add
	bpms_[k] = gslport_lvector(1, n_bpm);
	corrs_[k] = gslport_lvector(1, n_corr);
	//end GSL add

	/*Orginal
    A_lsoc[k] = dmatrix(1, n_bpm, 1, n_corr);
    U_lsoc[k] = dmatrix(1, n_bpm, 1, n_corr);
    w_lsoc[k] = dvector(1, n_corr);
    V_lsoc[k] = dmatrix(1, n_corr, 1, n_corr);
	*/
	
	//GSL add
	mA_lsoc[k] = gsl_matrix_alloc(n_bpm, n_corr);
	GSL2NRDM2(dmA_lsoc, mA_lsoc[k], A_lsoc[k],0);
    mU_lsoc[k] = gsl_matrix_alloc(n_bpm, n_corr);	
	GSL2NRDM2(dmU_lsoc, mU_lsoc[k], U_lsoc[k],0);
	vw_lsoc[k] = gsl_vector_alloc(n_corr);
	GSL2NRDV2(vw_lsoc[k],w_lsoc[k]);
    mV_lsoc[k] = gsl_matrix_alloc(n_corr, n_corr);	
	GSL2NRDM2(dmAV_lsoc, mV_lsoc[k], V_lsoc[k],0);
	
	S = gsl_vector_alloc(n_corr);
	//end GSL add
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


void gcmat1(const int bpm, const int corr, const int plane)
{
  /* Get correlation matrix

                -----------
              \/beta  beta
                    i     j
        A   = ------------- cos(nu pi - 2 pi|nu  - nu |)
         ij   2 sin(pi nu)                     i     j

  */

  bool      first;
  int       i, j, k;
  long int  loc;
  double    nu, betai, betaj, nui, nuj, spiq;

  const double  eps = 1e-10;

  k = plane - 1; n_bpm_[k] = GetnKid(bpm); n_corr_[k] = GetnKid(corr);

  first = (plane == 1)? first_h : first_v;
  if (first) {
    if (plane == 1)
      first_h = false;
    else
      first_v = false;
	  
	/*Original
    A_lsoc[k] = dmatrix(1, n_bpm_[k], 1, n_corr_[k]);
    U_lsoc[k] = dmatrix(1, n_bpm_[k], 1, n_corr_[k]);
    w_lsoc[k] = dvector(1, n_corr_[k]);
    V_lsoc[k] = dmatrix(1, n_corr_[k], 1, n_corr_[k]);
	*/
	
	//GSL add
	mA_lsoc[k] = gsl_matrix_alloc(n_bpm_[k], n_corr_[k]);
	GSL2NRDM2(dmA_lsoc, mA_lsoc[k], A_lsoc[k],0);
    mU_lsoc[k] = gsl_matrix_alloc(n_bpm_[k], n_corr_[k]);	
	GSL2NRDM2(dmU_lsoc, mU_lsoc[k], U_lsoc[k],0);
	vw_lsoc[k] = gsl_vector_alloc(n_corr_[k]);
	GSL2NRDV2(vw_lsoc[k],w_lsoc[k]);
    mV_lsoc[k] = gsl_matrix_alloc(n_corr_[k], n_corr_[k]);	
	GSL2NRDM2(dmAV_lsoc, mV_lsoc[k], V_lsoc[k],0);
	
	S = gsl_vector_alloc(n_corr_[k]);
	//end GSL add
  }

  nu = globval.TotalTune[k]; spiq = sin(M_PI*nu);

  for (i = 1; i <= n_bpm_[k]; i++) {
    loc = Elem_GetPos(bpm, i);
    betai = Cell[loc].Beta[k]; nui = Cell[loc].Nu[k];
    for (j = 1; j <= n_corr_[k]; j++) {
      loc = Elem_GetPos(corr, j);
      betaj = Cell[loc].Beta[k]; nuj = Cell[loc].Nu[k];
      A_lsoc[k][i][j] =
	sqrt(betai*betaj)/(2.0*spiq)*cos(nu*M_PI-fabs(2.0*M_PI*(nui-nuj)));
    }
  }

  for (i = 1; i <= n_bpm_[k]; i++)
    for (j = 1; j <= n_corr_[k]; j++)
      U_lsoc[k][i][j] = A_lsoc[k][i][j];

  /*Orgianl
  dsvdcmp(U_lsoc[k], n_bpm_[k], n_corr_[k], w_lsoc[k], V_lsoc[k]);
  */
  
  //GSL add
  gsl_linalg_SV_decomp (mU_lsoc[k], mV_lsoc[k], S, vw_lsoc[k]);
  //end GSL add
  
  for (j = 1; j <= n_corr_[k]; j++)
    if (w_lsoc[k][j] < eps) {
      printf("gcmat: singular beam response matrix"
	     " %12.3e, plane = %d, j = %d\n", w_lsoc[k][j], plane, j);
      prt_gcmat(plane);
      exit_(1);
    }
}


void lsoc(const int niter, const int plane)
{
  int       i, j, k;
  long int  loc;
  double    *b, *x;
  
  //GSL add
  gsl_vector *vb;
  gsl_vector *vx;
  //end GSL add

  k = plane - 1;

  /*Orginal
  b = dvector(1, n_bpm_[k]); x = dvector(1, n_corr_[k]);
  */
  
  //GSL add
  vb = gsl_vector_alloc(n_bpm_[k]);
  vx = gsl_vector_alloc(n_corr_[k]);
  GSL2NRDV2(vb, b);
  GSL2NRDV2(vx, x);
  //end GSL add

  for (i = 1; i <= niter; i++) {
    for (j = 1; j <= n_bpm_[k]; j++) {
      loc = bpms_[k][j];
      b[j] = -Cell[loc].BeamPos[2*k] + Cell[loc].dS[k];
    }
      
	/*Orginal
    dsvbksb(U_lsoc[k], w_lsoc[k], V_lsoc[k], n_bpm_[k], n_corr_[k], b, x);
	*/
	
	//GSL add
	gsl_vector *vs = gsl_vector_alloc(n_corr_[k]);
	gsl_vector *work = gsl_vector_alloc(n_corr_[k]);
	gsl_linalg_SV_decomp (mU_lsoc[k], mV_lsoc[k], vs, work);
	gsl_linalg_SV_solve(mU_lsoc[k], mV_lsoc[k],vs,vb,vx);
	gsl_vector_free(vs);
	gsl_vector_free(work);
	//end GSL add

    for (j = 1; j <= n_corr_[k]; j++) {
      loc = corrs_[k][j];
      if (plane == 1)
	set_dbnL_design_elem(Cell[loc].Fnum, Cell[loc].Knum, Dip, -x[j], 0.0);
      else
        set_dbnL_design_elem(Cell[loc].Fnum, Cell[loc].Knum, Dip, 0.0, x[j]);
    }
  }

  /*Orginal
  free_dvector(b, 1, n_bpm_[k]); free_dvector(x, 1, n_corr_[k]);
  */
  
  //GSL add
  gsl_vector_free(vb);
  gsl_vector_free(vx);
  //end GSL add
  
}


void lsoc(const int niter, const int bpm, const int corr, const int plane)
{

  lsoc(niter, plane);
}


void lsoc1(const int niter, const int bpm, const int corr, const int plane)
{
  int       i, j, k;
  long int  loc;
  double    *b, *x;
  
   //GSL add
  gsl_vector *vb;
  gsl_vector *vx;
  //end GSL add

  k = plane - 1;

  /*Orginal
  b = dvector(1, n_bpm_[k]); x = dvector(1, n_corr_[k]);
  */
  
  //GSL add
  vb = gsl_vector_alloc(n_bpm_[k]);
  vx = gsl_vector_alloc(n_corr_[k]);
  GSL2NRDV2(vb, b);
  GSL2NRDV2(vx, x);
  //end GSL add

  for (i = 1; i <= niter; i++) {
    for (j = 1; j <= n_bpm_[k]; j++) {
      loc = Elem_GetPos(bpm, j);
      b[j] = -Cell[loc].BeamPos[2*k] + Cell[loc].dS[k];
    }
	
    /*Orginal
    dsvbksb(U_lsoc[k], w_lsoc[k], V_lsoc[k], n_bpm_[k], n_corr_[k], b, x);
	*/
	
	//GSL add
	gsl_vector *vs = gsl_vector_alloc(n_corr_[k]);
	gsl_vector *work = gsl_vector_alloc(n_corr_[k]);
	gsl_linalg_SV_decomp (mU_lsoc[k], mV_lsoc[k], vs, work);
	gsl_linalg_SV_solve(mU_lsoc[k], mV_lsoc[k],vs,vb,vx);
	gsl_vector_free(vs);
	gsl_vector_free(work);
	//end GSL add

    for (j = 1; j <= n_corr_[k]; j++)
      if (plane == 1)
        SetdKpar(corr, j, Dip, -x[j]);
      else
        SetdKpar(corr, j, -Dip, x[j]);
  }

  /*Orginal
  free_dvector(b, 1, n_bpm_[k]); free_dvector(x, 1, n_corr_[k]);
  */
  
  //GSL add
  gsl_vector_free(vb);
  gsl_vector_free(vx);
  //end GSL add
  
}
