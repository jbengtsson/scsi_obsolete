/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -

                                                                              */

bool        first_FM = true;
double      c_1, d_1, c_2, d_2, cl_rad, q_fluct;
double      I2, I4, I5, dcurly_H, dI4, s_FM;
ElemFamType ElemFam[Elem_nFamMax];
CellType    Cell[Cell_nLocMax+1];
ofstream    outf_;

// For FieldMap.
bool  sympl             = true;
int   FieldMap_filetype = 2;


MarkerType::MarkerType(void) {}


DriftType::DriftType(void) {}


MpoleType::MpoleType(void)
{
  int j;

  method = Meth_Fourth; n = 0;
  for (j = 0; j <= 1; j++) {
    dSsys[j] = 0e0; dSrms[j] = 0e0; dSrnd[j] = 0e0;
  }
  rollpar = 0e0; rollsys = 0e0; rollrms = 0e0; rollrnd = 0e0;
  for (j = -HOMmax; j <= HOMmax; j++) {
    bn[j+HOMmax]    = 0e0; bnpar[j+HOMmax] = 0e0;
    bnsys[j+HOMmax] = 0e0; bnrms[j+HOMmax] = 0e0;
    bnrnd[j+HOMmax] = 0e0;
  }
  order = 0; n_design = 0; irho = 0e0; tx1 = 0e0; tx2 = 0e0; gap = 0e0;

  c0 = 0e0; c1 = 0e0; s1 = 0e0;
}


CavityType::CavityType() { volt = 0e0; freq = 0e0; phi = 0e0; h = 0; }


WigglerType::WigglerType(void)
{
  int j;

  method = Meth_Linear; n = 0;
  for (j = 0; j <= 1; j++) {
    dSsys[j] = 0e0; dSrnd[j] = 0e0;
  }
  rollpar = 0e0; rollsys = 0e0; rollrnd = 0e0; n_harm = 0; lambda = 0e0;
  for (j = 0; j < n_harm_max; j++) {
    BoBrhoV[j] = 0e0; BoBrhoH[j] = 0e0; kxV[j] = 0e0; kxH[j] = 0e0;
    phi[j] = 0e0;
  }
  for (j = 0; j <= HOMmax; j++)
    bn[j+HOMmax] = 0e0;
  order = 0;
}


InsertionType::InsertionType(void)
{
  int i, j;

  method = Meth_Linear; n = 0; nx = 0; nz = 0;

  /* Initialisation thetax and thetaz to 0*/

  // first order kick map
  if (firstorder){
    for (i = 0; i < IDZMAX; i++){
      for (j = 0; j < IDXMAX; j++) {
	thetax1[i][j] = 0e0; thetaz1[i][j] = 0e0; B2[i][j] = 0e0;
      }
    }
  }

  // second order kick map
  if (secondorder) {
    for (i = 0; i < IDZMAX; i++) {
      for (j = 0; j < IDXMAX; j++) {
          thetax[i][j] = 0e0; thetaz[i][j] = 0e0; B2[i][j] = 0e0;
      }
    }
  }

  // stuffs for interpolation
  for (j = 0; j < IDXMAX; j++)
    tabx[j] = 0e0;

  for (j = 0; j < IDZMAX; j++)
    tabz[j] = 0e0;

  // filenames
  strcpy(fname1,""); strcpy(fname2,"");

//  kx = 0e0;
  for (j = 0; j <= 1; j++) {
    dSsys[j] = 0e0; dSrnd[j] = 0e0;
  }
  rollpar = 0e0; rollsys = 0e0; rollrnd = 0e0;
//  for (j = 0; j <= HOMmax; j++)
//    bnW[j+HOMmax] = 0e0;
  order = 0;
}


FieldMapType::FieldMapType(void)
{

  n_step = 0; n[X_] = 0; n[Y_] = 0; n[Z_] = 0; scl = 1e0;
  phi = 0e0; Ld = 0e0; L1 = 0e0; cut = 0; x0 = 0e0;
}


SpreaderType::SpreaderType(void)
{
  int     k;

  for (k = 0; k < Spreader_max; k++)
    Cell_ptrs[k] = NULL;
}


RecombinerType::RecombinerType(void) {}


SolenoidType::SolenoidType(void)
{
  int j;

  n = 0;
  for (j = 0; j <= 1; j++) {
    dSsys[j] = 0e0; dSrms[j] = 0e0; dSrnd[j] = 0e0;
  }
  rollpar = 0e0; rollsys = 0e0; rollrnd = 0e0;
}


void Marker_Init(int Fnum)
{
  int         i;
  ElemFamType *elamfam;
  ElemType    *elem;
  CellType    *cell;

  elamfam = &ElemFam[Fnum-1];
  for (i = 0; i < elamfam->nKid; i++) {
    cell = &Cell[elamfam->KidList[i]];
    elem = dynamic_cast<ElemType*>(new MarkerType::MarkerType());
    // Copying of struct/class is provided by compiler.
    cell->Elem = elamfam->Elem;
    // Initialize Galilean misalignment group.
    cell->droll[0] = 1e0; cell->droll[1] = 0e0;
    cell->dS[0] = 0e0; cell->dS[1] = 0e0;
  }
}


void Drift_Init(int Fnum)
{
  int         i;
  ElemFamType *elemfam;
  ElemType    *elem;
  CellType    *cell;

  elemfam = &ElemFam[Fnum-1];
  for (i = 1; i <= elemfam->nKid; i++) {
    cell = &Cell[elemfam->KidList[i-1]];
    elem = dynamic_cast<ElemType*>(new DriftType::DriftType());
    // Copying of struct/class is provided by compiler.
    *elem = *elemfam->Elem;
    // Initialize Galilean misalignment group.
    cell->droll[0] = 1e0; cell->droll[1] = 0e0;
    cell->dS[0] = 0e0; cell->dS[1] = 0e0;
  }
}


static int Updateorder(MpoleType *M)
{
  int i, order;

  if (M->irho != 0e0) // Non zero curvature => bend.
    order = 1;
  else // Multipole.
    order = 0;
  for (i = -HOMmax; i <= HOMmax; i++)
    if (M->bn[i+HOMmax] != 0e0 && abs(i) > order) order = abs(i);

  return order;
}


void Mpole_Init(int Fnum)
{
  int         i;
  ElemFamType *elemfam;
  ElemType    *elem;
  CellType    *cell;

  elemfam = &ElemFam[Fnum-1];
  for (i = 1; i <= elemfam->nKid; i++) {
    cell = &Cell[elemfam->KidList[i-1]];
    elem = dynamic_cast<ElemType*>(new MpoleType::MpoleType());
    // Copying of struct/class is provided by compiler.
    *elem = *elemfam->Elem;
    // Initialize Galilean misalignment group.
    cell->droll[0] = 1e0; cell->droll[1] = 0e0;
    cell->dS[0] = 0e0; cell->dS[1] = 0e0;
  }
}


void Cav_Init(int Fnum)
{
  int         i;
  ElemFamType *elemfam;
  ElemType    *elem;
  CellType    *cell;

  elemfam = &ElemFam[Fnum-1];
  for (i = 1; i <= elemfam->nKid; i++) {
    cell = &Cell[elemfam->KidList[i-1]];
    elem = dynamic_cast<ElemType*>(new CavityType::CavityType());
    // Copying of struct/class is provided by compiler.
    *elem = *elemfam->Elem;
    // Initialize Galilean misalignment group.
    cell->droll[0] = 1e0; cell->droll[1] = 0e0;
    cell->dS[0] = 0e0; cell->dS[1] = 0e0;
  }
}


void Wiggler_Init(int Fnum)
{
  int         i;
  ElemFamType *elemfam;
  ElemType    *elem;
  CellType    *cell;

  elemfam = &ElemFam[Fnum-1];
  for (i = 1; i <= elemfam->nKid; i++) {
    cell = &Cell[elemfam->KidList[i-1]];
    elem = dynamic_cast<ElemType*>(new WigglerType::WigglerType());
    // Copying of struct/class is provided by compiler.
    *elem = *elemfam->Elem;
    // Initialize Galilean misalignment group.
    cell->droll[0] = 1e0; cell->droll[1] = 0e0;
    cell->dS[0] = 0e0; cell->dS[1] = 0e0;
  }
}


void Insertion_Init(int Fnum)
{
  int         i;
  ElemFamType *elemfam;
  ElemType    *elem;
  CellType    *cell;

  elemfam = &ElemFam[Fnum-1];
  for (i = 1; i <= elemfam->nKid; i++) {
    cell = &Cell[elemfam->KidList[i-1]];
    elem = dynamic_cast<ElemType*>(new InsertionType::InsertionType());
    // Copying of struct/class is provided by compiler.
    *elem = *elemfam->Elem;
    // Initialize Galilean misalignment group.
    cell->droll[0] = 1e0; cell->droll[1] = 0e0;
    cell->dS[0] = 0e0; cell->dS[1] = 0e0;
  }
}


void FieldMap_Init(int Fnum)
{
  int         i;
  ElemFamType *elemfam;
  ElemType    *elem;
  CellType    *cell;

  elemfam = &ElemFam[Fnum-1];
  for (i = 1; i <= elemfam->nKid; i++) {
    cell = &Cell[elemfam->KidList[i-1]];
    elem = dynamic_cast<ElemType*>(new FieldMapType::FieldMapType());
    // Copying of struct/class is provided by compiler.
    *elem = *elemfam->Elem;
    // Initialize Galilean misalignment group.
    cell->droll[0] = 1e0; cell->droll[1] = 0e0;
    cell->dS[0] = 0e0; cell->dS[1] = 0e0;
  }
}


void Spreader_Init(int Fnum)
{
  int         i;
  ElemFamType *elemfam;
  ElemType    *elem;
  CellType    *cell;

  elemfam = &ElemFam[Fnum-1];
  for (i = 1; i <= elemfam->nKid; i++) {
    cell = &Cell[elemfam->KidList[i-1]];
    elem = dynamic_cast<ElemType*>(new SpreaderType::SpreaderType());
    // Copying of struct/class is provided by compiler.
    *elem = *elemfam->Elem;
    // Initialize Galilean misalignment group.
    cell->droll[0] = 1e0; cell->droll[1] = 0e0;
    cell->dS[0] = 0e0; cell->dS[1] = 0e0;
  }
}


void Recombiner_Init(int Fnum)
{
  int         i;
  ElemFamType *elemfam;
  ElemType    *elem;
  CellType    *cell;

  elemfam = &ElemFam[Fnum-1];
  for (i = 1; i <= elemfam->nKid; i++) {
    cell = &Cell[elemfam->KidList[i-1]];
    elem = dynamic_cast<ElemType*>(new RecombinerType::RecombinerType());
    // Copying of struct/class is provided by compiler.
    *elem = *elemfam->Elem;
    // Initialize Galilean misalignment group.
    cell->droll[0] = 1e0; cell->droll[1] = 0e0;
    cell->dS[0] = 0e0; cell->dS[1] = 0e0;
  }
}


void Solenoid_Init(int Fnum)
{
  int         i;
  ElemFamType *elemfam;
  ElemType    *elem;
  CellType    *cell;

  elemfam = &ElemFam[Fnum-1];
  for (i = 1; i <= elemfam->nKid; i++) {
    cell = &Cell[elemfam->KidList[i-1]];
    elem = dynamic_cast<ElemType*>(new SolenoidType::SolenoidType());
    // Copying of struct/class is provided by compiler.
    *elem = *elemfam->Elem;
    // Initialize Galilean misalignment group.
    cell->droll[0] = 1e0; cell->droll[1] = 0e0;
    cell->dS[0] = 0e0; cell->dS[1] = 0e0;
  }
}


template<typename T>
void spline_(const double x[], const T y[], const int n, const double yp1,
	     const double ypn, T y2[])
{
  int     i,k;
  double  sig;
  T       p, u[n], qn, un;

  if (yp1 > 0.99e30)
    y2[1] = u[1] = 0e0;
  else {
    y2[1] = -0.5;
    u[1] = (3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }
  for (i = 2; i <= n-1; i++) {
    sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p = sig*y2[i-1]+2.0;
    y2[i] = (sig-1e0)/p;
    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn = un = 0e0;
  else {
    qn = 0.5;
    un = (3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n] = (un-qn*u[n-1])/(qn*y2[n-1]+1e0);
  for (k = n-1; k >= 1; k--)
    y2[k] = y2[k]*y2[k+1]+u[k];
}


void splie2_(const double x1a[], const double x2a[], double **ya,
	     const int m, const int n, double **y2a)
{
  int  j;

  for (j = 1; j <= m; j++)
    spline_(x2a,ya[j],n,1.0e30,1.0e30,y2a[j]);
}


template<typename T, typename U>
void splint_(const double xa[], const U ya[], const U y2a[],
	     const int n, const T &x, T &y)
{
  int     klo,khi,k;
  double  h;
  T       a,b;

  klo = 1;
  khi = n;
  while (khi-klo > 1) {
    k = (khi+klo) >> 1;
    if (xa[k] > x) khi = k;
    else klo = k;
  }
  h = xa[khi]-xa[klo];
  if (h == 0e0) {
    printf("Bad xa input to routine splint_\n");
    exit(1);
  }

  a = (xa[khi]-x)/h;
  b = (x-xa[klo])/h;
  y = a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}


template<typename T>
void splin2_(const double x1a[], const double x2a[], double **ya, double **y2a,
	     const int m, const int n, const T &x1, const T &x2, T &y)
{
  int  j;
  T    ytmp[m+1], yytmp[m+1];

  if ((x1 < x1a[1]) || (x1 > x1a[m])) {
    cout << fixed << setprecision(8)
	 << "splin2_: x undefined ["
	 << is_double<T>::cst(x1) << ", " << is_double<T>::cst(x2) << "] (["
	 << x1a[1] << ", " << x1a[m] << "])" << endl;

    y = NAN;

    return;
  }

  if ((x2 < x2a[1]) || (x2 > x2a[n])) {
    cout << fixed << setprecision(8)
	 << "splin2_: y undefined ["
	 << is_double<T>::cst(x1) << ", " << is_double<T>::cst(x2) << "] (["
	 << x2a[1] << ", " << x2a[n] << "])" << endl;

    y = NAN;

    return;
  }

  for (j = 1; j<= m; j++)
    splint_(x2a,ya[j],y2a[j],n,x2,yytmp[j]);
  spline_(x1a,yytmp,m,1.0e30,1.0e30,ytmp);
  splint_(x1a,yytmp,ytmp,m,x1,y);
}


void get_B_DIAMOND(const char *filename, FieldMapType *FM)
{
  char         line[max_str];
  int          i, j, n, ny;
  double       x0, y0, z0;
  ifstream     inf;
  ofstream     outf;

  const int     skip = 8;
  const double  Brho = globval.Energy*1e9/c0;

  cout << endl;
  cout << "get_B_DIAMOND: loading field map: " << filename << endl;

  file_rd(inf, filename);

  for (i = 1; i <= skip; i++)
    inf.getline(line, max_str);

  inf.getline(line, max_str);
  sscanf(line, "Starting point[cm] : %lf %lf %lf", &x0, &y0, &z0);
  inf.getline(line, max_str);
  sscanf(line, "Step size[cm]      : %lf %lf %lf",
	 &FM->dx[X_], &FM->dx[Y_], &FM->dx[Z_]);
  inf.getline(line, max_str);
  sscanf(line, "Number of points   : %d %d %d", &FM->n[X_], &ny, &FM->n[Z_]);

  // convert from cm to m
  x0 *= 1e-2; y0 *= 1e-2; z0 *= 1e-2;
  FM->dx[X_] *= 1e-2; FM->dx[Y_] *= 1e-2; FM->dx[Z_] *= 1e-2;
  FM->Lr = FM->dx[Z_]*(FM->n[Z_]-1);

  FM->n[Y_] = 2*ny - 1;

  FM->vx[X_] = gsl_vector_alloc(FM->n[X_]);
  GSL2NRDV2(FM->vx[X_],FM->x[X_]);

  FM->vx[Y_] = gsl_vector_alloc(FM->n[Y_]);
  GSL2NRDV2(FM->vx[Y_],FM->x[Y_]);

  FM->vx[Z_] = gsl_vector_alloc(FM->n[Z_]);
  GSL2NRDV2(FM->vx[Z_],FM->x[Z_]);

  FM->BoBrho[X_]  = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho[Y_]  = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho[Z_]  = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[X_] = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[Y_] = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[Z_] = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  FM->AoBrho[X_]  = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho[Y_]  = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho2[X_] = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho2[Y_] = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  for (i = 1; i <= 2; i++)
    inf.getline(line, max_str);

  for (n = 1; n <= FM->n[Z_]; n++)
    for (j = 1; j <= ny; j++)
      for (i = 1; i <= FM->n[X_]; i++) {
	inf.getline(line, max_str);
	sscanf(line, "%lf %lf %lf %lf %lf %lf",
	       &FM->x[X_][i], &FM->x[Y_][ny-1+j], &FM->x[Z_][n],
	       &FM->BoBrho[X_][n][i][ny-1+j],
	       &FM->BoBrho[Y_][n][i][ny-1+j],
	       &FM->BoBrho[Z_][n][i][ny-1+j]);

	// convert from cm to m
	FM->x[X_][i] *= 1e-2; FM->x[Y_][ny-1+j] *= 1e-2; FM->x[Z_][n] *= 1e-2;
	// convert from Gauss to Tesla
	FM->BoBrho[X_][n][i][ny-1+j] /= 1e+4*Brho;
	FM->BoBrho[Y_][n][i][ny-1+j] /= 1e+4*Brho;
	FM->BoBrho[Z_][n][i][ny-1+j] /= 1e+4*Brho;

	// Compute vector potential (axial gauge) by extended trapezodial rule
	if (n == 1) {
	  FM->AoBrho[X_][n][i][ny-1+j] =
	    -FM->BoBrho[Y_][n][i][ny-1+j]*FM->dx[Z_]/2.0;
	  FM->AoBrho[Y_][n][i][ny-1+j] =
	    FM->BoBrho[X_][n][i][ny-1+j]*FM->dx[Z_]/2.0;
	} else if (n == FM->n[Z_]) {
	  FM->AoBrho[X_][n][i][ny-1+j] =
	    FM->AoBrho[X_][n-1][i][ny-1+j]
	    - FM->BoBrho[Y_][n][i][ny-1+j]*FM->dx[Z_]/2.0;
	  FM->AoBrho[Y_][n][i][ny-1+j] =
	    FM->AoBrho[Y_][n-1][i][ny-1+j]
	    + FM->BoBrho[X_][n][i][ny-1+j]*FM->dx[Z_]/2.0;
	} else {
	  FM->AoBrho[X_][n][i][ny-1+j] =
	    FM->AoBrho[X_][n-1][i][ny-1+j]
	    - FM->BoBrho[Y_][n][i][ny-1+j]*FM->dx[Z_];
	  FM->AoBrho[Y_][n][i][ny-1+j] =
	    FM->AoBrho[Y_][n-1][i][ny-1+j]
	    + FM->BoBrho[X_][n][i][ny-1+j]*FM->dx[Z_];
	}
      }

  inf.close();

  cout << fixed << setprecision(5)
       << setw(10) << x0 << setw(10) << y0 << setw(10) << z0 << endl;
  cout << fixed << setprecision(5)
       << setw(10) << FM->dx[X_] << setw(10) << FM->dx[Y_]
       << setw(10) << FM->dx[Z_] << endl;
  cout << setw(10) << FM->n[X_] << setw(10) << ny
       << setw(10) << FM->n[Z_] << endl;
  cout << fixed << setprecision(3)
       << setw(10) << FM->x[X_][1] << setw(10) << FM->x[X_][FM->n[X_]]
       << setw(10) << FM->x[Y_][1] << setw(10) << FM->x[Y_][FM->n[Y_]]
       << setw(10) << FM->x[Z_][1] << setw(10) << FM->x[Z_][FM->n[Z_]] << endl;
  cout << fixed << setprecision(5)
       << "Magnet length [m]:" << setw(10) << FM->Lr << endl;

  for (j = 1; j <= ny-1; j++) {
    FM->x[Y_][j] = -FM->x[Y_][2*ny-j];
    for (i = 1; i <= FM->n[X_]; i++) {
      for (n = 1; n <= FM->n[Z_]; n++) {
	// B[X_] is antisymmetric in y (rot(A) = 0)
	FM->BoBrho[X_][n][i][j] = -FM->BoBrho[X_][n][i][2*ny-j];
	FM->BoBrho[Y_][n][i][j] =  FM->BoBrho[Y_][n][i][2*ny-j];
	// Bz is antisymmetric in y
	FM->BoBrho[Z_][n][i][j] = -FM->BoBrho[Z_][n][i][2*ny-j];

	FM->AoBrho[X_][n][i][j] = FM->AoBrho[X_][n][i][2*ny-j];
	// Ay is antisymmetric in y
	FM->AoBrho[Y_][n][i][j] = -FM->AoBrho[Y_][n][i][2*ny-j];
      }
    }
  }

  if (true) {
    file_wr(outf, "field_map.dat");
    cout << scientific << setprecision(3)
	 << setw(11) << FM->x[X_][7] << setw(11) << FM->x[Y_][7] << endl;
    for (i = 1; i <= FM->n[X_]; i++)
      outf << scientific << setprecision(3)
	   << setw(11) << FM->x[X_][i] << setw(11) << FM->BoBrho[Y_][7][i][7]
	   << endl;
    outf.close();
  }

  for (n = 1; n <= FM->n[Z_]; n++) {
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[X_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[X_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Y_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[Y_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Z_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[Z_][n]);

    splie2_(FM->x[X_], FM->x[Y_], FM->AoBrho[X_][n],
	    FM->n[X_], FM->n[Y_], FM->AoBrho2[X_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][n],
	    FM->n[X_], FM->n[Y_], FM->AoBrho2[Y_][n]);
  }

  cout << "field map loaded: " << filename << endl;

/*  free_dvector(FM->x[X_], 1, FM->n[X_]); free_dvector(FM->x[Y_], 1, FM->n[Y_]);
  free_dvector(FM->x[Z_], 1, FM->n[Z_]);

  free_df3tensor(FM->BoBrho[X_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho[Y_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho[Z_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[X_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[Y_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[Z_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  free_df3tensor(FM->AoBrho[X_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho[Y_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho2[X_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho2[Y_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);*/
}


void get_B_NSLS_II(const char *filename, FieldMapType *FM)
{
  char         line[max_str];
  int          i, j, n;
  double       x_min[3], x_max[3];
  ifstream     inf;

  const double  Brho = globval.Energy*1e9/c0;

  cout << endl;
  cout << "get_B_NSLS_II: loading field map: " << filename << endl;

  file_rd(inf, filename);

  inf.getline(line, max_str);
  inf.getline(line, max_str);
  sscanf(line, "#%lf <= x <= %lf, dx = %lf, nx = %d",
	 &x_min[X_], &x_max[X_], &FM->dx[X_], &FM->n[X_]);
  x_min[X_] *= 1e-2; x_max[X_] *= 1e-2; FM->dx[X_] *= 1e-2;
  inf.getline(line, max_str);
  sscanf(line, "#%lf <= y <= %lf, dy = %lf, ny = %d",
	 &x_min[Y_], &x_max[Y_], &FM->dx[Y_], &FM->n[Y_]);
  x_min[Y_] *= 1e-2; x_max[Y_] *= 1e-2; FM->dx[Y_] *= 1e-2;
  inf.getline(line, max_str);
  sscanf(line, "#%lf <= z <= %lf, dz = %lf, nz = %d",
	 &x_min[Z_], &x_max[Z_], &FM->dx[Z_], &FM->n[Z_]);
  x_min[Z_] *= 1e-2; x_max[Z_] *= 1e-2; FM->dx[Z_] *= 1e-2;

  FM->vx[X_] = gsl_vector_alloc(FM->n[X_]);
  GSL2NRDV2(FM->vx[X_],FM->x[X_]);

  FM->vx[Y_] = gsl_vector_alloc(FM->n[Y_]);
  GSL2NRDV2(FM->vx[Y_],FM->x[Y_]);

  FM->vx[Z_] = gsl_vector_alloc(FM->n[Z_]);
  GSL2NRDV2(FM->vx[Z_],FM->x[Z_]);

  FM->BoBrho[X_]  = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho[Y_]  = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho[Z_]  = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[X_] = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[Y_] = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[Z_] = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  FM->AoBrho[X_]  = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho[Y_]  = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho2[X_] = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho2[Y_] = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  for (i = 1; i <= FM->n[X_]; i++)
    for (j = 1; j <= FM->n[Y_]; j++)
      for (n = 1; n <= FM->n[Z_]; n++) {
	inf.getline(line, max_str);
	sscanf(line, "%lf %lf %lf %lf %lf %lf",
	       &FM->x[X_][i], &FM->x[Y_][j], &FM->x[Z_][n],
	       &FM->BoBrho[X_][n][i][j],
	       &FM->BoBrho[Y_][n][i][j],
	       &FM->BoBrho[Z_][n][i][j]);

	// convert from cm to m
	FM->x[X_][i] *= 1e-2; FM->x[Y_][j] *= 1e-2; FM->x[Z_][n] *= 1e-2;

	FM->BoBrho[X_][n][i][j] = Brho;
	FM->BoBrho[Y_][n][i][j] = Brho;
	FM->BoBrho[Z_][n][i][j] = Brho;

	// Compute vector potential (axial gauge) by extended trapezodial rule
 	if (n == 1) {
	  FM->AoBrho[X_][n][i][j] = -FM->BoBrho[Y_][n][i][j]*FM->dx[Z_]/2.0;
	  FM->AoBrho[Y_][n][i][j] =  FM->BoBrho[X_][n][i][j]*FM->dx[Z_]/2.0;
	} else if (n == FM->n[Z_]) {
	  FM->AoBrho[X_][n][i][j] =
	    FM->AoBrho[X_][n-1][i][j] - FM->BoBrho[Y_][n][i][j]*FM->dx[Z_]/2.0;
	  FM->AoBrho[Y_][n][i][j] =
	    FM->AoBrho[Y_][n-1][i][j] + FM->BoBrho[X_][n][i][j]*FM->dx[Z_]/2.0;
	} else {
	  FM->AoBrho[X_][n][i][j] =
	    FM->AoBrho[X_][n-1][i][j] - FM->BoBrho[Y_][n][i][j]*FM->dx[Z_];
	  FM->AoBrho[Y_][n][i][j] =
	    FM->AoBrho[Y_][n-1][i][j] + FM->BoBrho[X_][n][i][j]*FM->dx[Z_];
	}
     }

  inf.close();

  FM->Lr = FM->dx[Z_]*(FM->n[Z_]-1);

  cout << fixed << setprecision(5)
       << setw(10) << 1e3*FM->dx[X_] << setw(10) << 1e3*FM->dx[Y_]
       << setw(10) << 1e3*FM->dx[Z_] << endl;
  cout << setw(10) << FM->n[X_] << setw(10) << FM->n[Y_]
       << setw(10) << FM->n[Z_] << endl;
  cout << fixed << setprecision(3)
       << setw(10) << FM->x[X_][1] << setw(10) << FM->x[X_][FM->n[X_]]
       << setw(10) << FM->x[Y_][1] << setw(10) << FM->x[Y_][FM->n[Y_]]
       << setw(10) << FM->x[Z_][1] << setw(10) << FM->x[Z_][FM->n[Z_]] << endl;
  cout << fixed << setprecision(5)
       << "Magnet length [m]:" << setw(10) << FM->Lr << endl;

  for (n = 1; n <= FM->n[Z_]; n++) {
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[X_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[X_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Y_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[Y_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Z_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[Z_][n]);

    splie2_(FM->x[X_], FM->x[Y_], FM->AoBrho[X_][n],
	    FM->n[X_], FM->n[Y_], FM->AoBrho2[X_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][n],
	    FM->n[X_], FM->n[Y_], FM->AoBrho2[Y_][n]);
  }

  cout << "field map loaded: " << filename << endl;

/*  free_dvector(FM->x[X_], 1, FM->n[X_]);
  free_dvector(FM->x[Y_], 1, FM->n[Y_]);
  free_dvector(FM->x[Z_], 1, FM->n[Z_]);

  free_df3tensor(FM->BoBrho[X_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho[Y_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho[Z_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[X_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[Y_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[Z_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  free_df3tensor(FM->AoBrho[X_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho[Y_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho2[X_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho2[Y_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);*/
}


void get_B_Oleg1(const char *filename, FieldMapType *FM)
{
  char         line[max_str];
  int          i, j, n;
  double       x_min[3], x_max[3];
  ifstream     inf;

  const double  Brho = globval.Energy*1e9/c0;

  cout << endl;
  cout << "get_B_Oleg1: loading field map: " << filename << endl;

  file_rd(inf, filename);

  inf.getline(line, max_str);
  inf.getline(line, max_str);
  sscanf(line, "#%lf <= x <= %lf, dx = %lf, nx = %d",
	 &x_min[X_], &x_max[X_], &FM->dx[X_], &FM->n[X_]);
  FM->dx[X_] *= 1e-3;
  inf.getline(line, max_str);
  sscanf(line, "#%lf <= y <= %lf, dy = %lf, ny = %d",
	 &x_min[Y_], &x_max[Y_], &FM->dx[Y_], &FM->n[Y_]);
  FM->dx[Y_] *= 1e-3;
  inf.getline(line, max_str);
  sscanf(line, "#%lf <= z <= %lf, dz = %lf, nz = %d",
	 &x_min[Z_], &x_max[Z_], &FM->dx[Z_], &FM->n[Z_]);
  FM->dx[Z_] *= 1e-3;

  FM->vx[X_] = gsl_vector_alloc(FM->n[X_]);
  GSL2NRDV2(FM->vx[X_],FM->x[X_]);

  FM->vx[Y_] = gsl_vector_alloc(FM->n[Y_]);
  GSL2NRDV2(FM->vx[Y_],FM->x[Y_]);

  FM->vx[Z_] = gsl_vector_alloc(FM->n[Z_]);
  GSL2NRDV2(FM->vx[Z_],FM->x[Z_]);

  FM->BoBrho[X_]  = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho[Y_]  = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho[Z_]  = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[X_] = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[Y_] = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[Z_] = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  FM->AoBrho[X_]  = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho[Y_]  = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho2[X_] = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho2[Y_] = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  for (n = 1; n <= FM->n[Z_]; n++)
    for (j = 1; j <= FM->n[Y_]; j++)
      for (i = 1; i <= FM->n[X_]; i++) {
	inf.getline(line, max_str);
	sscanf(line, "%lf %lf %lf %lf %lf %lf",
	       &FM->x[X_][i], &FM->x[Y_][j], &FM->x[Z_][n],
	       &FM->BoBrho[X_][n][i][j],
	       &FM->BoBrho[Y_][n][i][j],
	       &FM->BoBrho[Z_][n][i][j]);

	// convert from mm to m
	FM->x[X_][i] *= 1e-3; FM->x[Y_][j] *= 1e-3; FM->x[Z_][n] *= 1e-3;

	FM->BoBrho[X_][n][i][j] /= Brho;
	FM->BoBrho[Y_][n][i][j] /= Brho;
	FM->BoBrho[Z_][n][i][j] /= Brho;

	// Compute vector potential (axial gauge) by extended trapezodial rule
 	if (n == 1) {
	  FM->AoBrho[X_][n][i][j] = -FM->BoBrho[Y_][n][i][j]*FM->dx[Z_]/2.0;
	  FM->AoBrho[Y_][n][i][j] =  FM->BoBrho[X_][n][i][j]*FM->dx[Z_]/2.0;
	} else if (n == FM->n[Z_]) {
	  FM->AoBrho[X_][n][i][j] =
	    FM->AoBrho[X_][n-1][i][j] - FM->BoBrho[Y_][n][i][j]*FM->dx[Z_]/2.0;
	  FM->AoBrho[Y_][n][i][j] =
	    FM->AoBrho[Y_][n-1][i][j] + FM->BoBrho[X_][n][i][j]*FM->dx[Z_]/2.0;
	} else {
	  FM->AoBrho[X_][n][i][j] =
	    FM->AoBrho[X_][n-1][i][j] - FM->BoBrho[Y_][n][i][j]*FM->dx[Z_];
	  FM->AoBrho[Y_][n][i][j] =
	    FM->AoBrho[Y_][n-1][i][j] + FM->BoBrho[X_][n][i][j]*FM->dx[Z_];
	}
     }

  inf.close();

  FM->Lr = FM->dx[Z_]*(FM->n[Z_]-1);

  cout << fixed << setprecision(5)
       << setw(10) << 1e3*FM->dx[X_] << setw(10) << 1e3*FM->dx[Y_]
       << setw(10) << 1e3*FM->dx[Z_] << endl;
  cout << setw(10) << FM->n[X_] << setw(10) << FM->n[Y_]
       << setw(10) << FM->n[Z_] << endl;
  cout << fixed << setprecision(3)
       << setw(10) << FM->x[X_][1] << setw(10) << FM->x[X_][FM->n[X_]]
       << setw(10) << FM->x[Y_][1] << setw(10) << FM->x[Y_][FM->n[Y_]]
       << setw(10) << FM->x[Z_][1] << setw(10) << FM->x[Z_][FM->n[Z_]] << endl;
  cout << fixed << setprecision(5)
       << "Magnet length [m]:" << setw(10) << FM->Lr << endl;

  for (n = 1; n <= FM->n[Z_]; n++) {
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[X_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[X_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Y_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[Y_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Z_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[Z_][n]);

    splie2_(FM->x[X_], FM->x[Y_], FM->AoBrho[X_][n],
	    FM->n[X_], FM->n[Y_], FM->AoBrho2[X_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][n],
	    FM->n[X_], FM->n[Y_], FM->AoBrho2[Y_][n]);
  }

  cout << "field map loaded: " << filename << endl;

/*  free_dvector(FM->x[X_], 1, FM->n[X_]);
  free_dvector(FM->x[Y_], 1, FM->n[Y_]);
  free_dvector(FM->x[Z_], 1, FM->n[Z_]);

  free_df3tensor(FM->BoBrho[X_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho[Y_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho[Z_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[X_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[Y_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->BoBrho2[Z_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  free_df3tensor(FM->AoBrho[X_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho[Y_],  1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho2[X_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  free_df3tensor(FM->AoBrho2[Y_], 1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);*/
}


void get_B_Oleg2(const char *filename, FieldMapType *FM)
{
  char         line[max_str];
  int          i, j, n;
  double       x_min[3];
  ifstream     inf;

  const double  Brho = globval.Energy*1e9/c0;

  cout << endl;
  cout << "get_B_Oleg2: loading field map: " << filename << endl;

  file_rd(inf, filename);

  inf.getline(line, max_str);

  inf.getline(line, max_str); sscanf(line, "#%lf", &x_min[X_]);
  inf.getline(line, max_str); sscanf(line, "#%lf", &FM->dx[X_]);
  inf.getline(line, max_str); sscanf(line, "#%d", &FM->n[X_]);

  inf.getline(line, max_str); sscanf(line, "#%lf", &x_min[Y_]);
  inf.getline(line, max_str); sscanf(line, "#%lf", &FM->dx[Y_]);
  inf.getline(line, max_str); sscanf(line, "#%d", &FM->n[Y_]);

  inf.getline(line, max_str); sscanf(line, "#%lf", &x_min[Z_]);
  inf.getline(line, max_str); sscanf(line, "#%lf", &FM->dx[Z_]);
  inf.getline(line, max_str); sscanf(line, "#%d", &FM->n[Z_]);

  cout << fixed << setprecision(5)
       << setw(10) << 1e3*FM->dx[X_] << setw(10) << 1e3*FM->dx[Y_]
       << setw(10) << 1e3*FM->dx[Z_] << endl;
  cout << setw(10) << FM->n[X_] << setw(10) << FM->n[Y_]
       << setw(10) << FM->n[Z_] << endl;
  cout << fixed << setprecision(3)
       << setw(10) << x_min[X_] << setw(10) << x_min[Y_]
       << setw(10) << x_min[Z_] << endl;

  FM->vx[X_] = gsl_vector_alloc(FM->n[X_]);
  GSL2NRDV2(FM->vx[X_],FM->x[X_]);

  FM->vx[Y_] = gsl_vector_alloc(FM->n[Y_]);
  GSL2NRDV2(FM->vx[Y_],FM->x[Y_]);

  FM->vx[Z_] = gsl_vector_alloc(FM->n[Z_]);
  GSL2NRDV2(FM->vx[Z_],FM->x[Z_]);

  FM->BoBrho[X_]  = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho[Y_]  = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho[Z_]  = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[X_] = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[Y_] = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->BoBrho2[Z_] = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  FM->AoBrho[X_]  = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho[Y_]  = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho2[X_] = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);
  FM->AoBrho2[Y_] = gslport_tensor(1, FM->n[Z_], 1, FM->n[X_], 1, FM->n[Y_]);

  for (n = 1; n <= FM->n[Z_]; n++) {
    FM->x[Z_][n] = (n == 1)? x_min[Z_] : FM->x[Z_][n-1] + FM->dx[Z_];

    for (j = 1; j <= FM->n[Y_]; j++) {
      FM->x[Y_][j] = (j == 1)? x_min[Y_] : FM->x[Y_][j-1] + FM->dx[Y_];

      for (i = 1; i <= FM->n[X_]; i++) {
	FM->x[X_][i] = (i == 1)? x_min[X_] : FM->x[X_][i-1] + FM->dx[X_];

	inf.getline(line, max_str);
	sscanf(line, "%lf %lf %lf",
	       &FM->BoBrho[X_][n][i][j], &FM->BoBrho[Y_][n][i][j],
	       &FM->BoBrho[Z_][n][i][j]);

	FM->BoBrho[X_][n][i][j] /= Brho;
	FM->BoBrho[Y_][n][i][j] /= Brho;
	FM->BoBrho[Z_][n][i][j] /= Brho;

	// Compute vector potential (axial gauge) by extended trapezodial rule
 	if (n == 1) {
	  FM->AoBrho[X_][n][i][j] = -FM->BoBrho[Y_][n][i][j]*FM->dx[Z_]/2.0;
	  FM->AoBrho[Y_][n][i][j] =  FM->BoBrho[X_][n][i][j]*FM->dx[Z_]/2.0;
	} else if (n == FM->n[Z_]) {
	  FM->AoBrho[X_][n][i][j] =
	    FM->AoBrho[X_][n-1][i][j] - FM->BoBrho[Y_][n][i][j]*FM->dx[Z_]/2.0;
	  FM->AoBrho[Y_][n][i][j] =
	    FM->AoBrho[Y_][n-1][i][j] + FM->BoBrho[X_][n][i][j]*FM->dx[Z_]/2.0;
	} else {
	  FM->AoBrho[X_][n][i][j] =
	    FM->AoBrho[X_][n-1][i][j] - FM->BoBrho[Y_][n][i][j]*FM->dx[Z_];
	  FM->AoBrho[Y_][n][i][j] =
	    FM->AoBrho[Y_][n-1][i][j] + FM->BoBrho[X_][n][i][j]*FM->dx[Z_];
	}
      }
    }
  }

  inf.close();

  FM->Lr = FM->dx[Z_]*(FM->n[Z_]-1);

  cout << fixed << setprecision(5)
       << setw(10) << 1e3*FM->dx[X_] << setw(10) << 1e3*FM->dx[Y_]
       << setw(10) << 1e3*FM->dx[Z_] << endl;
  cout << setw(10) << FM->n[X_] << setw(10) << FM->n[Y_]
       << setw(10) << FM->n[Z_] << endl;
  cout << fixed << setprecision(3)
       << setw(10) << FM->x[X_][1] << setw(10) << FM->x[X_][FM->n[X_]]
       << setw(10) << FM->x[Y_][1] << setw(10) << FM->x[Y_][FM->n[Y_]]
       << setw(10) << FM->x[Z_][1] << setw(10) << FM->x[Z_][FM->n[Z_]] << endl;
  cout << fixed << setprecision(5)
       << "Magnet length [m]:" << setw(10) << FM->Lr << endl;

  for (n = 1; n <= FM->n[Z_]; n++) {
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[X_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[X_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Y_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[Y_][n]);
    splie2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Z_][n],
	    FM->n[X_], FM->n[Y_], FM->BoBrho2[Z_][n]);
  }

  cout << "field map loaded: " << filename << endl;

}


void get_B(const char *filename, FieldMapType *FM)
{
  // Do not scale fieldmaps only Hamiltonians, i.e., the kick.  Note that RADIA
  // (2nd order) kick maps are quadratic in the field, and 1st order linear.

  switch (FieldMap_filetype) {
  case 1:
    get_B_DIAMOND(filename, FM);
    break;
  case 2:
    get_B_NSLS_II(filename, FM);
    break;
  case 3:
    get_B_Oleg1(filename, FM);
    break;
  case 4:
    get_B_Oleg2(filename, FM);
    break;
  default:
    cout << "get_B: unknown FieldMap type " << FieldMap_filetype << endl;
    break;
  }
}


template<typename T>
void GtoL(ss_vect<T> &X, Vector2 &S, Vector2 &R,
	  const double c0, const double c1, const double s1)
{
  ss_vect<T>  x1;

  /* Simplified rotated p_rot */
  X[px_] += c1; X[py_] += s1;
  /* Translate */
  X[x_] -= S[X_]; X[y_] -= S[Y_];
  /* Rotate */
  x1 = X;
  X[x_]  =  R[X_]*x1[x_]  + R[Y_]*x1[y_];
  X[px_] =  R[X_]*x1[px_] + R[Y_]*x1[py_];
  X[y_]  = -R[Y_]*x1[x_]  + R[X_]*x1[y_];
  X[py_] = -R[Y_]*x1[px_] + R[X_]*x1[py_] ;
  /* Simplified p_rot */
  X[px_] -= c0;
}


template<typename T>
void LtoG(ss_vect<T> &X, Vector2 &S, Vector2 &R,
	  double c0, double c1, double s1)
{
  ss_vect<T>  x1;

  /* Simplified p_rot */
  X[px_] -= c0;
  /* Rotate */
  x1 = X;
  X[x_]  = R[X_]*x1[x_]  - R[Y_]*x1[y_];
  X[px_] = R[X_]*x1[px_] - R[Y_]*x1[py_];
  X[y_]  = R[Y_]*x1[x_]  + R[X_]*x1[y_];
  X[py_] = R[Y_]*x1[px_] + R[X_]*x1[py_];
  /* Translate */
  X[x_] += S[X_]; X[y_] += S[Y_];
  /* p_rot rotated */
  X[px_] += c1; X[py_] += s1;
}


template<typename T>
inline T get_p_s(const ss_vect<T> &x)
{
  T  p_s, p_s2;

  if (!globval.H_exact)
    p_s = 1e0+x[delta_];
  else {
    p_s2 = sqr(1e0+x[delta_]) - sqr(x[px_]) - sqr(x[py_]);
    if (p_s2 >= 0e0)
      p_s = sqrt(p_s2);
    else {
      printf("get_p_s: *** Speed of light exceeded!\n");
      p_s = NAN;
    }
  }
  return(p_s);
}


template<typename T>
void Drift(double L, ss_vect<T> &x)
{
  T  u;

  if (!globval.H_exact) {
    u = L/(1e0+x[delta_]);
    x[ct_] += u*(sqr(x[px_])+sqr(x[py_]))/(2.0*(1e0+x[delta_]));
  } else {
    u = L/get_p_s(x);
    x[ct_] += u*(1e0+x[delta_]) - L;
  }
  x[x_] += x[px_]*u; x[y_] += x[py_]*u;
  if (globval.pathlength) x[ct_] += L;
}


template<typename T>
void Drift_Pass(CellType &Cell, ss_vect<T> &x) { Drift(Cell.Elem->L, x); }


void zero_mat(const int n, double** A)
{
  int  i, j;

  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      A[i][j] = 0e0;
}


void identity_mat(const int n, double** A)
{
  int  i, j;

  for (i = 1; i <= n; i++)
    for (j = 1; j <= n; j++)
      A[i][j] = (i == j)? 1e0 : 0e0;
}


double det_mat_gsl(const int n, gsl_matrix *tmpA)
{
  double det;
  int signum;
  gsl_permutation *p = gsl_permutation_alloc(n);

  gsl_linalg_LU_decomp(tmpA , p , &signum);
  det = gsl_linalg_LU_det(tmpA , signum);
  gsl_permutation_free(p);
  gsl_matrix_free(tmpA);

  return det;
}


double det_mat(const int n, double **A)
{
  gsl_matrix *tmpA = gsl_matrix_alloc(n,n);
  int ii,jj;
  for (ii=1; ii <=n; ii++) {
    for (jj=1; jj <=n; jj++){
      gsl_matrix_set (tmpA, ii-1, jj-1, A[ii][jj]);
    }
  }

  double d = det_mat_gsl(n,tmpA);
  return d;
}


double trace_mat(const int n, double **A)
{
  int     i;
  double  d;

  d = 0e0;
  for (i = 1; i <= n; i++)
    d += A[i][i];

  return d;
}

double trace_mat_gsl(const int n, gsl_matrix *m)
{
  int     i;
  double  d;

  d = 0e0;
  for (i = 0; i < n; i++)
    d += gsl_matrix_get(m, i, i);

  return d;
}


// partial template-class specialization
// primary version
template<typename T>
class is_tps { };

// partial specialization
template<>
class is_tps<double> {
 public:
  static inline void get_ps(const ss_vect<double> &x, CellType &Cell)
  { Cell.BeamPos = x; }

  static inline double set_prm(const int k) { return 1e0; }

  static inline double get_curly_H(const ss_vect<tps> &x)
    {
      cout << "get_curly_H: operation not defined for double" << endl;
      exit_(1);
      return 0e0;
    }

  static inline double get_dI4(const double h, const double b2, const double L,
			     const ss_vect<tps> &x)
    {
      cout << "get_dI4: operation not defined for double" << endl;
      exit_(1);
      return 0e0;
    }

  static inline void emittance(const double B2, const double u,
			       const double ps0, const ss_vect<double> &xp) { }

  static inline void diff_mat(const double B2, const double u,
			      const double ps0, const ss_vect<double> &xp) { }

};


// partial specialization
template<>
class is_tps<tps> {
 public:
  static inline void get_ps(const ss_vect<tps> &x, CellType &Cell)
  {
    Cell.BeamPos = x.cst(); getlinmat(6, x, Cell.A);
  }

  static inline tps set_prm(const int k) { return tps(0e0, k); }

  static inline double get_curly_H(const ss_vect<tps> &A)
  {
    int              j;
    double           curly_H[2];
    ss_vect<double>  eta;

    eta.zero();
    for (j = 0; j < 4; j++)
      eta[j] = A[j][delta_];

    get_twoJ(2, eta, A, curly_H);

    return curly_H[X_];
  }

  static inline double get_dI4(const ss_vect<tps> &A) { return A[x_][delta_]; }

  static inline void emittance(const tps &B2_perp, const tps &ds, const tps &ps0,
			       const ss_vect<tps> &A) {
    // M. Sands "The Physics of Electron Storage Rings" SLAC-121, p. 118.
    // d<delta^2>/ds = 3*C_U*C_gamma*h_bar*c*E_0^5*(1+delta)^4*(B_perp/(Brho))^3
    //                 /(4*pi*m_e^3)
    // A contains the eigenvectors.
    int           j;
    double        B_66;
    ss_vect<tps>  A_inv;

    if (B2_perp > 0e0) {
      B_66 = (q_fluct*pow(B2_perp.cst(), 1.5)*pow(ps0, 4)*ds).cst();
      A_inv = Inv(A);
      // D_11 = D_22 = curly_H_x,y * B_66 / 2,
      // curly_H_x,y = eta_Fl^2 + etap_Fl^2
      for (j = 0; j < 3; j++)
	globval.D_rad[j] +=
	  (sqr(A_inv[j*2][delta_])+sqr(A_inv[j*2+1][delta_]))*B_66/2.0;
    }
  }

  static inline void diff_mat(const tps &B2_perp, const tps &ds, const tps &ps0,
			      ss_vect<tps> &x) { }

};


template<typename T>
void get_B2(const double h_ref, const T B[], const ss_vect<T> &xp,
	    T &B2_perp, T &B2_par)
{
  // compute B_perp^2 and B_par^2
  T  xn, e[3];

  xn = 1e0/sqrt(sqr(1e0+xp[x_]*h_ref)+sqr(xp[px_])+sqr(xp[py_]));
  e[X_] = xp[px_]*xn; e[Y_] = xp[py_]*xn; e[Z_] = (1e0+xp[x_]*h_ref)*xn;

  // left-handed coordinate system
  B2_perp =
    sqr(B[Y_]*e[Z_]-B[Z_]*e[Y_]) + sqr(B[X_]*e[Y_]-B[Y_]*e[X_])
    + sqr(B[Z_]*e[X_]-B[X_]*e[Z_]);

//  B2_par = sqr(B[X_]*e[X_]+B[Y_]*e[Y_]+B[Z_]*e[Z_]);
}


template<typename T>
void radiate(ss_vect<T> &x, const double L, const double h_ref, const T B[])
{
  // M. Sands "The Physics of Electron Storage Rings" SLAC-121, p. 98.
  // ddelta/d(ds) = -C_gamma*E_0^3*(1+delta)^2*(B_perp/(Brho))^2/(2*pi)
  T           ps0, ps1, ds, B2_perp = 0e0, B2_par = 0e0;
  ss_vect<T>  xp;

  // large ring: x' and y' unchanged
  xp = x; ps0 = get_p_s(x); xp[px_] /= ps0; xp[py_] /= ps0;

  // H = -p_s => ds = H*L
  ds = (1e0+xp[x_]*h_ref+(sqr(xp[px_])+sqr(xp[py_]))/2.0)*L;
  get_B2(h_ref, B, xp, B2_perp, B2_par);

  if (globval.radiation) {
    x[delta_] -= cl_rad*sqr(ps0)*B2_perp*ds;
    ps1 = get_p_s(x); x[px_] = xp[px_]*ps1; x[py_] = xp[py_]*ps1;
  }

  if (globval.emittance) is_tps<T>::emittance(B2_perp, ds, ps0, xp);
}


template<typename T>
void radiate_ID(ss_vect<T> &x, const double L, const T &B2_perp)
{
  T           ps0, ps1, ds;
  ss_vect<T>  xp;

  // large ring: x' and y' unchanged
  xp = x; ps0 = get_p_s(x); xp[px_] /= ps0; xp[py_] /= ps0;

  // H = -p_s => ds = H*L
  ds = (1e0+(sqr(xp[px_])+sqr(xp[py_]))/2.0)*L;

  if (globval.radiation) {
    x[delta_] -= cl_rad*sqr(ps0)*B2_perp*ds;
    ps1 = get_p_s(x); x[px_] = xp[px_]*ps1; x[py_] = xp[py_]*ps1;
  }

  if (globval.emittance) is_tps<T>::emittance(B2_perp, ds, ps0, xp);
}


static double get_psi(double irho, double phi, double gap)
{
  /* Correction for magnet gap (longitudinal fringe field)

       irho h = 1/rho [1/m]
       phi  edge angle
       gap  full gap between poles

                                    2
                   K1*gap*h*(1 + sin phi)
            psi = ----------------------- * (1 - K2*g*gap*tan phi)
                        cos phi

            K1 is usually 1/2
            K2 is zero here                                                  */

  double  psi;

  const double k1 = 0.5, k2 = 0e0;

  if (phi == 0e0)
    psi = 0e0;
  else
    psi = k1*gap*irho*(1e0+sqr(sin(dtor(phi))))/cos(dtor(phi))
          *(1e0 - k2*gap*irho*tan(dtor(phi)));

  return psi;
}


template<typename T>
void thin_kick(int Order, double MB[], double L, double h_bend, double h_ref,
	       ss_vect<T> &x)
{
  /* The kick is given by

              e L       L delta    L x              e L
     Dp_x = - --- B_y + ------- - ----- ,    Dp_y = --- B_x
              p_0         rho     rho^2             p_0

    where

                           ====
                           \
      (B_y + iB_x) = B rho  >   (ia_n  + b_n ) (x + iy)^n-1
                           /
                           ====

    where

       e      1
      --- = -----
      p_0   B rho                                                            */

  int         j;
  T           BxoBrho, ByoBrho, ByoBrho1, B[3];
  ss_vect<T>  x0, cod;

  if ((h_bend != 0e0) || ((1 <= Order) && (Order <= HOMmax))) {
    x0 = x;
    /* compute field with Horner's rule */
    ByoBrho = MB[Order+HOMmax]; BxoBrho = MB[HOMmax-Order];
    for (j = Order-1; j >= 1; j--) {
      ByoBrho1 = x0[x_]*ByoBrho - x0[y_]*BxoBrho + MB[j+HOMmax];
      BxoBrho  = x0[y_]*ByoBrho + x0[x_]*BxoBrho + MB[HOMmax-j];
      ByoBrho  = ByoBrho1;
    }

    if (globval.radiation || globval.emittance) {
      B[X_] = BxoBrho; B[Y_] = ByoBrho + h_bend; B[Z_] = 0e0;
      radiate(x, L, h_ref, B);
    }

    if (h_ref != 0e0) {
      x[px_] -= L*(ByoBrho+(h_bend-h_ref)/2.0+h_ref*h_bend*x0[x_]
		-h_ref*x0[delta_]);
      x[ct_] += L*h_ref*x0[x_];
    } else
      x[px_] -= L*(ByoBrho+h_bend);
    x[py_] += L*BxoBrho;
  }
}


template<typename T>
static void EdgeFocus(double irho, double phi, double gap, ss_vect<T> &x)
{
  x[px_] += irho*tan(dtor(phi))*x[x_];
  if (false)
    // warning: => diverging Taylor map (see SSC-141)
    x[py_] -= irho*tan(dtor(phi)-get_psi(irho, phi, gap))*x[y_]
              /(1e0+x[delta_]);
  else
    x[py_] -= irho*tan(dtor(phi)-get_psi(irho, phi, gap))*x[y_];
}


template<typename T>
void p_rot(double phi, ss_vect<T> &x)
{
  T           c, s, ps, p;
  ss_vect<T>  x1;

  c = cos(dtor(phi)); s = sin(dtor(phi)); ps = get_p_s(x);

  if (!globval.H_exact) {
     x[px_] = s*ps + c*x[px_];
  } else {
    x1 = x; p = c*ps - s*x1[px_];
    x[x_] = x1[x_]*ps/p; x[px_] = s*ps + c*x1[px_];
    x[y_] += x1[x_]*x1[py_]*s/p;
    x[ct_] += (1e0+x1[delta_])*x1[x_]*s/p;
  }
}


template<typename T>
void bend_fringe(double hb, ss_vect<T> &x)
{
  T           coeff, u, ps, ps2, ps3;
  ss_vect<T>  x1;

  coeff = -hb/2.0; x1 = x; ps = get_p_s(x); ps2 = sqr(ps); ps3 = ps*ps2;
  u = 1e0 + 4.0*coeff*x1[px_]*x1[y_]*x1[py_]/ps3;
  if (u >= 0e0) {
    x[y_] = 2.0*x1[y_]/(1e0+sqrt(u));
    x[x_] = x1[x_] - coeff*sqr(x[y_])*(ps2+sqr(x1[px_]))/ps3;
    x[py_] = x1[py_] + 2.0*coeff*x1[px_]*x[y_]/ps;
    x[ct_] = x1[ct_] - coeff*x1[px_]*sqr(x[y_])*(1e0+x1[delta_])/ps3;
  } else {
    printf("bend_fringe: *** Speed of light exceeded!\n");
    x[x_] = NAN; x[px_] = NAN; x[y_] = NAN; x[py_] = NAN;
    x[delta_] = NAN; x[ct_] = NAN;
  }
}


template<typename T>
void quad_fringe(double b2, ss_vect<T> &x)
{
  T  u, ps;

  u = b2/(12.0*(1e0+x[delta_])); ps = u/(1e0+x[delta_]);
  x[py_] /= 1e0 - 3.0*u*sqr(x[y_]); x[y_] -= u*cube(x[y_]);
  if (globval.Cavity_on) x[ct_] -= ps*cube(x[y_])*x[py_];
  x[px_] /= 1e0 + 3.0*u*sqr(x[x_]);
  if (globval.Cavity_on) x[ct_] += ps*cube(x[x_])*x[px_];
  x[x_] += u*cube(x[x_]); u = u*3.0; ps = ps*3.0;
  x[y_] = exp(-u*sqr(x[x_]))*x[y_]; x[py_] = exp(u*sqr(x[x_]))*x[py_];
  x[px_] += 2.0*u*x[x_]*x[y_]*x[py_];
  if (globval.Cavity_on) x[ct_] -= ps*sqr(x[x_])*x[y_]*x[py_];
  x[x_] = exp(u*sqr(x[y_]))*x[x_]; x[px_] = exp(-u*sqr(x[y_]))*x[px_];
  x[py_] -= 2.0*u*x[y_]*x[x_]*x[px_];
  if (globval.Cavity_on) x[ct_] += ps*sqr(x[y_])*x[x_]*x[px_];
}


template<typename T>
void Mpole_Pass(CellType &Cell, ss_vect<T> &x)
{
  int        seg = 0;
  double     dL = 0e0, dL1 = 0e0, dL2 = 0e0;
  double     dkL1 = 0e0, dkL2 = 0e0, h_ref = 0e0;
  ElemType   *elem;
  MpoleType  *M;

  M = Cell.Elem;

  /* Global -> Local */
  GtoL(x, Cell.dS, Cell.droll, M->c0, M->c1, M->s1);

  if ((M->method == Meth_Second) || (M->method == Meth_Fourth)) {
    /* fringe fields */
    if (globval.quad_fringe && (M->bn[Quad+HOMmax] != 0e0))
      quad_fringe(M->bn[Quad+HOMmax], x);
    if (!globval.H_exact) {
      if (M->irho != 0e0) EdgeFocus(M->irho, M->tx1, M->gap, x);
    } else {
//      p_rot(elem->L*M->irho/2.0*180e0/M_PI, x);
      p_rot(M->tx1, x); bend_fringe(M->irho, x);
    }
  }

  if (M->thick == thick) {
//    if (!globval.H_exact || ((M->PTx1 == 0e0) && (M->PTx2 == 0e0))) {
    if (!globval.H_exact) {
      // polar coordinates
      h_ref = M->irho; dL = elem->L/M->n;
    } else {
      // Cartesian coordinates
      h_ref = 0e0;
      if (M->irho == 0e0)
	dL = elem->L/M->n;
      else
//	dL = 1e0/M->irho*(sin(dtor(M->PTx1))+sin(dtor(M->PTx2)))/M->n;
	dL = 2.0/M->irho*sin(elem->L*M->irho/2.0)/M->n;
    }
  }

  switch (M->method) {

  case Meth_Linear:

  case Meth_First:
    cout << "Mpole_Pass: Meth_First not supported" << endl;
    exit_(1);
    break;

  case Meth_Second:
    cout << "Mpole_Pass: Meth_Second not supported" << endl;
    exit_(1);
    break;

  case Meth_Fourth:
    if (M->thick == thick) {
      dL1 = c_1*dL; dL2 = c_2*dL; dkL1 = d_1*dL; dkL2 = d_2*dL;

      dcurly_H = 0e0; dI4 = 0e0;
      for (seg = 1; seg <= M->n; seg++) {
	if (globval.emittance && (!globval.Cavity_on) && (M->irho != 0e0)) {
	  dcurly_H += is_tps<tps>::get_curly_H(x);
	  dI4 += is_tps<tps>::get_dI4(x);
	}

	Drift(dL1, x);
        thin_kick(M->order, M->bn, dkL1, M->irho, h_ref, x);
	Drift(dL2, x);
        thin_kick(M->order, M->bn, dkL2, M->irho, h_ref, x);

	if (globval.emittance && (!globval.Cavity_on) && (M->irho != 0e0)) {
	  dcurly_H += 4.0*is_tps<tps>::get_curly_H(x);
	  dI4 += 4.0*is_tps<tps>::get_dI4(x);
	}

	Drift(dL2, x);
        thin_kick(M->order, M->bn, dkL1, M->irho, h_ref, x);
	Drift(dL1, x);

	if (globval.emittance && (!globval.Cavity_on) && (M->irho != 0e0)) {
	  dcurly_H += is_tps<tps>::get_curly_H(x);
	  dI4 += is_tps<tps>::get_dI4(x);
	}
      }

      if (globval.emittance && (!globval.Cavity_on) && (M->irho != 0)) {
	dcurly_H /= 6.0*M->n;
	dI4 *= M->irho*(sqr(M->irho)+2.0*M->bnpar[Quad+HOMmax])/(6.0*M->n);

	I2 += elem->L*sqr(M->irho); I4 += elem->L*dI4;
	I5 += elem->L*fabs(cube(M->irho))*dcurly_H;
      }
    } else
      thin_kick(M->order, M->bn, 1e0, 0e0, 0e0, x);

    break;
  }

  if ((M->method == Meth_Second) || (M->method == Meth_Fourth)) {
    /* fringe fields */
    if (!globval.H_exact) {
      if (M->irho != 0e0) EdgeFocus(M->irho, M->tx2, M->gap, x);
    } else {
      bend_fringe(-M->irho, x); p_rot(M->tx2, x);
//      p_rot(elem->L*M->irho/2.0*180e0/M_PI, x);
    }
    if (globval.quad_fringe && (M->bn[Quad+HOMmax] != 0e0))
      quad_fringe(-M->bn[Quad+HOMmax], x);
  }

  /* Local -> Global */
  LtoG(x, Cell.dS, Cell.droll, M->c0, M->c1, M->s1);
}


template<typename T>
void Marker_Pass(CellType &Cell, ss_vect<T> &X)
{
  ElemType *elem;

  elem = &Cell.Elem;
  /* Global -> Local */
  GtoL(X, Cell.dS, Cell.droll, 0e0, 0e0, 0e0);
  /* Local -> Global */
  LtoG(X, Cell.dS, Cell.droll, 0e0, 0e0, 0e0);
}


template<typename T>
void Cav_Pass(CellType &Cell, ss_vect<T> &X)
{
  ElemType    *elem;
  CavityType  *C;
  T           delta;

  C = &Cell.Elem;
  if (globval.Cavity_on && C->volt != 0e0) {
    delta = -C->volt/(globval.Energy*1e9)
            *sin(2.0*M_PI*C->freq/c0*X[ct_]+C->phi);
    X[delta_] += delta;

    if (globval.radiation) globval.dE -= is_double<T>::cst(delta);

    if (globval.pathlength) X[ct_] -= C->h/C->freq*c0;
  }
}


template<typename T>
inline void get_Axy(const WigglerType *W, const double z,
		    ss_vect<T> &x, T AxoBrho[], T AyoBrho[])

{
  int     i;
  double  ky, kz_n;
  T       cx, cz, sx, sz, chy, shy;

  for (i = 0; i <= 3; ++i) {
    AxoBrho[i] = 0e0; AyoBrho[i] = 0e0;
  }

  for (i = 0; i < W->n_harm; i ++) {
    kz_n = W->harm[i]*2.0*M_PI/W->lambda; ky = sqrt(sqr(W->kxV[i])+sqr(kz_n));
    cx = cos(W->kxV[i]*x[x_]); sx = sin(W->kxV[i]*x[x_]);
    chy = cosh(ky*x[y_]); shy = sinh(ky*x[y_]); sz = sin(kz_n*z);

    AxoBrho[0] += W->BoBrhoV[i]/kz_n*cx*chy*sz;
    AyoBrho[0] += W->BoBrhoV[i]*W->kxV[i]/(ky*kz_n)*sx*shy*sz;

    // derivatives with respect to x
    AxoBrho[1] -= W->BoBrhoV[i]*W->kxV[i]/kz_n*sx*chy*sz;
    AyoBrho[1] += W->BoBrhoV[i]*sqr(W->kxV[i])/(ky*kz_n)*cx*shy*sz;

    // derivatives with respect to y
    AxoBrho[2] += W->BoBrhoV[i]*ky/kz_n*cx*shy*sz;
    AyoBrho[2] += W->BoBrhoV[i]*W->kxV[i]/kz_n*sx*chy*sz;

    if (globval.radiation) {
      cz = cos(kz_n*z);
      // derivatives with respect to z
      AxoBrho[3] += W->BoBrhoV[i]*cx*chy*cz;
      AyoBrho[3] += W->BoBrhoV[i]*W->kxV[i]/ky*sx*shy*cz;
    }
  }
}

/*
template<typename T>
inline void get_Axy_map(const FieldMapType *FM, const double z,
			const ss_vect<T> &x, T AxoBrho[], T AyoBrho[])
{
  float  y, ax0, ax1, ax2, ay0, ay1, ay2;

  const  float dy = 1e-3, dz = 1e-3;

  y = is_double<T>::cst(x[y_]);

  if ((z < FM->s_pos[1]) || (z > FM->s_pos[FM->n_s])) {
    cout << scientific << setprecision(3)
	 << "get_Axy_map: s out of range " << z << endl;
    exit_(1);
  }

  if ((y < FM->y_pos[1]) || (y > FM->y_pos[FM->m_y])) {
    cout << scientific << setprecision(3)
	 << "get_Axy_map: y out of range " << y << endl;
    exit_(1);
  }

  splin2(FM->y_pos, FM->s_pos, FM->AxoBrho, FM->AxoBrho2, FM->m_y, FM->n_s,
	 y, z, &ax1);
  AxoBrho[0] = FM->scl*ax1;

  splin2(FM->y_pos, FM->s_pos, FM->AyoBrho, FM->AyoBrho2, FM->m_y, FM->n_s,
	 y, z, &ay1);
  AyoBrho[0] = FM->scl*ay1;

  // derivatives with respect to x
  AxoBrho[1] = FM->scl*0e0; AyoBrho[1] = FM->scl*0e0;

  // derivatives with respect to y
  splin2(FM->y_pos, FM->s_pos, FM->AxoBrho, FM->AxoBrho2, FM->m_y, FM->n_s,
	 y+dy, z, &ax2);
  splin2(FM->y_pos, FM->s_pos, FM->AxoBrho, FM->AxoBrho2, FM->m_y, FM->n_s,
	 y-dy, z, &ax1);
  splin2(FM->y_pos, FM->s_pos, FM->AxoBrho, FM->AxoBrho2, FM->m_y, FM->n_s,
	 y, z, &ax0);
  AxoBrho[2] =
    (ax2-ax1)/(2.0*dy) + (ax2+ax1-2.0*ax0)/sqr(dy)*is_tps<T>::set_prm(y_+1);
  AxoBrho[2] *= FM->scl;

  splin2(FM->y_pos, FM->s_pos, FM->AyoBrho, FM->AyoBrho2, FM->m_y, FM->n_s,
	 y+dy, z, &ay2);
  splin2(FM->y_pos, FM->s_pos, FM->AyoBrho, FM->AyoBrho2, FM->m_y, FM->n_s,
	 y-dy, z, &ay1);
  splin2(FM->y_pos, FM->s_pos, FM->AyoBrho, FM->AyoBrho2, FM->m_y, FM->n_s,
	 y, z, &ay0);
  AyoBrho[2] =
    (ay2-ay1)/(2.0*dy) + (ay2+ay1-2.0*ay0)/sqr(dy)*is_tps<T>::set_prm(y_+1);
  AyoBrho[2] *= FM->scl;

  if (globval.radiation) {
    // derivatives with respect to z
    splin2(FM->y_pos, FM->s_pos, FM->AxoBrho, FM->AxoBrho2, FM->m_y, FM->n_s,
	   y, z+dz, &ax2);
    splin2(FM->y_pos, FM->s_pos, FM->AxoBrho, FM->AxoBrho2, FM->m_y, FM->n_s,
	   y, z-dz, &ax1);
    AxoBrho[3] = (ax2-ax1)/(2.0*dz); AxoBrho[3] *= FM->scl;

    splin2(FM->y_pos, FM->s_pos, FM->AyoBrho, FM->AyoBrho2, FM->m_y, FM->n_s,
	   y, z+dz, &ay2);
    splin2(FM->y_pos, FM->s_pos, FM->AyoBrho, FM->AyoBrho2, FM->m_y, FM->n_s,
	   y, z-dz, &ay1);
    AyoBrho[3] = (ay2-ay1)/(2.0*dz); AyoBrho[3] *= FM->scl;
    if (false)
      cout << fixed << setprecision(5)
	   << setw(8) << z << setw(9) << is_double<T>::cst(AxoBrho[3]) << endl;
  }
}
*/

template<typename T>
void Wiggler_pass_EF(const ElemType *elem, ss_vect<T> &x)
{
  // First order symplectic integrator for wiggler using expanded Hamiltonian

  int          i, nstep = 0;
  double       h, z;
  WigglerType  *W;
  FieldMapType *FM;
  T            AxoBrho[4], AyoBrho[4], psi, hodp, a12, a21, a22, det;
  T            d1, d2, a11, c11, c12, c21, c22, x2, B[3];

  switch (elem->kind) {
  case Wigl:
    W = (WigglerType*)elem;
    nstep = W->n;
    break;
  case FieldMap:
    FM = (FieldMapType*)elem;
    nstep = FM->n_step;
    break;
  default:
    cout << "Wiggler_pass_EF: unknown element type" << endl;
    exit_(1);
    break;
  }

  h = elem->L/nstep; z = 0e0;
  for (i = 1; i <= nstep; ++i) {
    switch (elem->kind) {
    case Wigl:
      get_Axy(W, z, x, AxoBrho, AyoBrho);
      break;
    case FieldMap:
//      get_Axy_map(FM, z, x, AxoBrho, AyoBrho);
      break;
    default:
      cout << "Wiggler_pass_EF: unknown element type" << endl;
      exit_(1);
      break;
    }

    psi = 1e0 + x[delta_]; hodp = h/psi;
    a11 = hodp*AxoBrho[1]; a12 = hodp*AyoBrho[1];
    a21 = hodp*AxoBrho[2]; a22 = hodp*AyoBrho[2];
    det = 1e0 - a11 - a22 + a11*a22 - a12*a21;
    d1 = hodp*AxoBrho[0]*AxoBrho[1]; d2 = hodp*AxoBrho[0]*AxoBrho[2];
    c11 = (1e0-a22)/det; c12 = a12/det; c21 = a21/det; c22 = (1e0-a11)/det;
    x2 = c11*(x[px_]-d1) + c12*(x[py_]-d2);

    x[py_] = c21*(x[px_]-d1) + c22*(x[py_]-d2); x[px_] = x2;
    x[x_] += hodp*(x[px_]-AxoBrho[0]); x[y_] += hodp*x[py_];
    x[ct_] += h*(sqr((x[px_]-AxoBrho[0])/psi)
	      + sqr((x[py_]-AyoBrho[0])/psi))/2.0;

    if (false)
      cout << scientific << setprecision(3)
	   << setw(8) << z
	   << setw(11) << is_double<T>::cst(x[x_])
	   << setw(11) << is_double<T>::cst(x[px_])
	   << setw(11) << is_double<T>::cst(x[y_])
	   << setw(11) << is_double<T>::cst(x[py_])
	   << endl;

    if (globval.pathlength) x[ct_] += h;

    if (globval.radiation || globval.emittance) {
      B[X_] = -AyoBrho[3]; B[Y_] = AxoBrho[3]; B[Z_] = AyoBrho[1] - AxoBrho[2];
      radiate(x, h, 0e0, B);
    }

    z += h;
  }
}


template<typename T>
inline void get_Axy2(const double z,
		     const double kxV, const double kxH, const double kz,
		     const double BoBrhoV, const double BoBrhoH,
		     const double phi,
		     ss_vect<T> &x, T AxoBrho[], T AyoBrho[])
{
  int  i;
  T    cx, sx, cz1, cz2, sz1, sz2, chy, shy, kyH, kyV, chx, shx, cy, sy;

  for (i = 0; i <= 3; ++i) {
    AxoBrho[i] = 0e0; AyoBrho[i] = 0e0;
  }

  kyV = sqrt(sqr(kz)+sqr(kxV)); kyH = sqrt(sqr(kz)+sqr(kxH));
  cx = cos(kxV*x[x_]); sx = sin(kxV*x[x_]);
  cy = cos(kxH*x[y_]); sy = sin(kxH*x[y_]);
  chx = cosh(kyH*x[x_]); shx = sinh(kyH*x[x_]);
  chy = cosh(kyV*x[y_]); shy = sinh(kyV*x[y_]);
  sz1 = sin(kz*z); sz2 = sin(kz*z+phi);

  AxoBrho[0] += BoBrhoV/kz*cx*chy*sz1;
  AxoBrho[0] -= BoBrhoH*kxH/(kyH*kz)*shx*sy*sz2;
  AyoBrho[0] += BoBrhoV*kxV/(kyV*kz)*sx*shy*sz1;
  AyoBrho[0] -= BoBrhoH/kz*chx*cy*sz2;

  /* derivatives with respect to x */
  AxoBrho[1] -= BoBrhoV*kxV/kz*sx*chy*sz1;
  AxoBrho[1] -= BoBrhoH*kxH/kz*chx*sy*sz2;
  AyoBrho[1] += BoBrhoV*sqr(kxV)/(kyV*kz)*cx*shy*sz1;
  AyoBrho[1] -= BoBrhoH*kyH/kz*shx*cy*sz2;

  /* derivatives with respect to y */
  AxoBrho[2] += BoBrhoV*kyV/kz*cx*shy*sz1;
  AxoBrho[2] -= BoBrhoH*sqr(kxH)/(kyH*kz)*shx*cy*sz2;
  AyoBrho[2] += BoBrhoV*kxV/kz*sx*chy*sz1;
  AyoBrho[2] += BoBrhoH*kxH/kz*chx*sy*sz2;

  if (globval.radiation) {
    cz1 = cos(kz*z); cz2=cos(kz*z+phi);
    /* derivatives with respect to z */
    AxoBrho[3] += BoBrhoV*cx*chy*cz1;
    AxoBrho[3] -= BoBrhoH*kxH/kyH*shx*sy*cz2;
    AyoBrho[3] += BoBrhoV*kxV/kyV*sx*shy*cz1;
    AyoBrho[3] -= BoBrhoH*chx*cy*cz2;
  }
}


template<typename T>
void Wiggler_pass_EF2(int nstep, double L, double kxV, double kxH, double kz,
		      double BoBrhoV, double BoBrhoH, double phi,
		      ss_vect<T> &x)
{
  // First order symplectic integrator for wiggler using expanded Hamiltonian

  int     i;
  double  h, z;
  T       hodp, B[3], px1, px2, px3, py1, py2, AxoBrho[4], AyoBrho[4], psi;
  T       px = 0e0, py = 0e0;

  h = L/nstep; z = 0e0;
  for (i = 1; i <= nstep; ++i) {
    get_Axy2(z, kxV, kxH, kz, BoBrhoV, BoBrhoH, phi, x, AxoBrho, AyoBrho);

    psi = 1e0 + x[delta_]; hodp = h/psi;

    px1 = (x[px_]-(AxoBrho[0]*AxoBrho[1]+AyoBrho[0]*AyoBrho[1])*hodp)
          *(1-AyoBrho[2]*hodp);
    px2 = (x[py_]-(AxoBrho[0]*AxoBrho[2]+AyoBrho[0]*AyoBrho[2])*hodp)
          *AyoBrho[1]*hodp;
    px3 = (1-AxoBrho[1]*hodp)*(1-AyoBrho[2]*hodp)
          - AxoBrho[2]*AyoBrho[1]*hodp*hodp;

    py1 = (x[py_]-(AxoBrho[0]*AxoBrho[2]+AyoBrho[0]*AyoBrho[2])*hodp)
          *(1-AxoBrho[1]*hodp);
    py2 = (x[px_]-(AxoBrho[0]*AxoBrho[1]+AyoBrho[0]*AyoBrho[1])*hodp)
          *AxoBrho[2]*hodp;

    py = (py1+py2)/px3; px = (px1+px2)/px3;
    x[x_] += hodp*(px-AxoBrho[0]); x[y_] += hodp*(py-AyoBrho[0]);
    x[ct_] += h*(sqr((px-AxoBrho[0])/psi) + sqr((py-AyoBrho[0])/psi))/2.0;

    if (globval.pathlength) x[ct_] += h;

    if (globval.radiation || globval.emittance) {
      B[X_] = -AyoBrho[3]; B[Y_] = AxoBrho[3]; B[Z_] = AyoBrho[1] - AxoBrho[2];
      radiate(x, h, 0e0, B);
    }

    z += h;
  }

  x[px_] = px; x[py_] = py;
}


template<typename T>
inline void get_Axy_EF3(const WigglerType *W, const double z,
		       const ss_vect<T> &x,
		       T &AoBrho, T dAoBrho[], T &dp, const bool hor)
{
  int     i;
  double  ky, kz_n;
  T       cx, sx, sz, chy, shy, cz;

  AoBrho = 0e0; dp = 0e0;

  for (i = 0; i < 3; i++)
    dAoBrho[i] = 0e0;

  for (i = 0; i < W->n_harm; i++) {
    kz_n = W->harm[i]*2.0*M_PI/W->lambda; ky = sqrt(sqr(W->kxV[i])+sqr(kz_n));

    cx  = cos(W->kxV[i]*x[x_]); sx = sin(W->kxV[i]*x[x_]);
    chy = cosh(ky*x[y_]); shy = sinh(ky*x[y_]); sz = sin(kz_n*z);

    if (hor) {
      // A_x/Brho
      AoBrho += W->BoBrhoV[i]/kz_n*cx*chy*sz;

      if (globval.radiation) {
	cz = cos(kz_n*z);
	dAoBrho[X_] -= W->BoBrhoV[i]*W->kxV[i]/kz_n*sx*chy*sz;
	dAoBrho[Y_] += W->BoBrhoV[i]*ky/kz_n*cx*shy*sz;
	dAoBrho[Z_] += W->BoBrhoV[i]*cx*chy*cz;
      }

      // dp_y
      if (W->kxV[i] == 0e0)
	dp += W->BoBrhoV[i]/kz_n*ky*x[x_]*shy*sz;
      else
	dp += W->BoBrhoV[i]/(W->kxV[i]*kz_n)*ky*sx*shy*sz;
    } else {
      // A_y/Brho
      AoBrho += W->BoBrhoV[i]*W->kxV[i]/(ky*kz_n)*sx*shy*sz;

      if (globval.radiation) {
	cz = cos(kz_n*z);
	dAoBrho[X_] += W->BoBrhoV[i]*sqr(W->kxV[i])/(ky*kz_n)*cx*shy*sz;
	dAoBrho[Y_] += W->BoBrhoV[i]*W->kxV[i]/kz_n*sx*chy*sz;
	dAoBrho[Z_] += W->BoBrhoV[i]*W->kxV[i]/ky*sx*shy*cz;
      }

      // dp_x
      dp += W->BoBrhoV[i]/kz_n*sqr(W->kxV[i]/ky)*cx*chy*sz;
    }
  }
}


template<typename T>
void Wiggler_pass_EF3(const ElemType *elem, ss_vect<T> &x)
{
  /* Second order symplectic integrator for insertion devices based on:

       E. Forest, et al "Explicit Symplectic Integrator for s-dependent
       Static Magnetic Field"                                                */

  int         i;
  double      h, z;
  WigglerType *W;
  T           hd, AxoBrho, AyoBrho, dAxoBrho[3], dAyoBrho[3], dpy, dpx, B[3];

  W = (WigglerType*)elem;

  h = elem->L/W->n; z = 0e0;
  for (i = 1; i <= W->n; i++) {
    hd = h/(1e0+x[delta_]);

    // 1: half step in z
    z += 0.5*h;

    // 2: half drift in y
    get_Axy_EF3(W, z, x, AyoBrho, dAyoBrho, dpx, false);

    x[px_] -= dpx; x[py_] -= AyoBrho;
    x[y_] += 0.5*hd*x[py_];
    x[ct_] += sqr(0.5)*hd*sqr(x[py_])/(1e0+x[delta_]);

    get_Axy_EF3(W, z, x, AyoBrho, dAyoBrho, dpx, false);

    x[px_] += dpx; x[py_] += AyoBrho;

    // 3: full drift in x
    get_Axy_EF3(W, z, x, AxoBrho, dAxoBrho, dpy, true);

    x[px_] -= AxoBrho; x[py_] -= dpy; x[x_] += hd*x[px_];
    x[ct_] += 0.5*hd*sqr(x[px_])/(1e0+x[delta_]);

    if (globval.pathlength) x[ct_] += h;

    get_Axy_EF3(W, z, x, AxoBrho, dAxoBrho, dpy, true);

    x[px_] += AxoBrho; x[py_] += dpy;

    // 4: a half drift in y
    get_Axy_EF3(W, z, x, AyoBrho, dAyoBrho, dpx, false);

    x[px_] -= dpx; x[py_] -= AyoBrho;
    x[y_] += 0.5*hd*x[py_];
    x[ct_] += sqr(0.5)*hd*sqr(x[py_])/(1e0+x[delta_]);

    get_Axy_EF3(W, z, x, AyoBrho, dAyoBrho, dpx, false);

    x[px_] += dpx; x[py_] += AyoBrho;

    // 5: half step in z
    z += 0.5*h;

    if (globval.radiation || globval.emittance) {
      get_Axy_EF3(W, z, x, AyoBrho, dAyoBrho, dpx, false);
      get_Axy_EF3(W, z, x, AxoBrho, dAxoBrho, dpy, true);
      B[X_] = -dAyoBrho[Z_]; B[Y_] = dAxoBrho[Z_];
      B[Z_] = dAyoBrho[X_] - dAxoBrho[Y_];
      radiate(x, h, 0e0, B);
    }
  }
}


template<typename T>
void Wiggler_Pass(CellType &Cell, ss_vect<T> &X)
{
  int         seg;
  double      L, L1, L2, K1, K2;
  ElemType    *elem;
  WigglerType *W;
  ss_vect<T>  X1;

  elem = &Cell.Elem; W = (WigglerType*)elem;
  GtoL(X, Cell.dS, Cell.droll, 0e0, 0e0, 0e0);
  switch (W->method) {

  case Meth_Linear:
    cout << "Wiggler_Pass: Meth_Linear not supported" << endl;
    exit_(1);
    break;

  case Meth_First:
    if ((W->BoBrhoV[0] != 0e0) || (W->BoBrhoH[0] != 0e0)) {
      if (!globval.EPU)
	Wiggler_pass_EF(Cell.Elem, X);
      else {
	Wiggler_pass_EF2(W->n, elem->L, W->kxV[0], W->kxH[0],
		2.0*M_PI/W->lambda, W->BoBrhoV[0], W->BoBrhoH[0],
		W->phi[0], X);
      }
    } else
      // drift if field = 0
      Drift(elem->L, X);
    break;

  case Meth_Second:
    if ((W->BoBrhoV[0] != 0e0) || (W->BoBrhoH[0] != 0e0)) {
      Wiggler_pass_EF3(Cell.Elem, X);
    } else
      // drift if field = 0
      Drift(elem->L, X);
    break;

  case Meth_Fourth:  /* 4-th order integrator */
    L = elem->L/W->n;
    L1 = c_1*L; L2 = c_2*L; K1 = d_1*L; K2 = d_2*L;
    for (seg = 1; seg <= W->n; seg++) {
      Drift(L1, X); X1 = X;
      thin_kick(W->order, W->bn, K1, 0e0, 0e0, X1);
      X[py_] = X1[py_]; Drift(L2, X); X1 = X;
      thin_kick(W->order, W->bn, K2, 0e0, 0e0, X1);
      X[py_] = X1[py_]; Drift(L2, X); X1 = X;
      thin_kick(W->order, W->bn, K1, 0e0, 0e0, X1);
      X[py_] = X1[py_]; Drift(L1, X);
    }
    break;
  }
  // Local -> Global
  LtoG(X, Cell.dS, Cell.droll, 0e0, 0e0, 0e0);
}

#undef eps
#undef kx


template<typename T>
void rk4_(const CellType &Cell, const ss_vect<T> &y, const ss_vect<T> dydx,
	  const double x, const double h, ss_vect<T> &yout,
	  void (*derivs)(const CellType &, const double, const ss_vect<T> &,
			 ss_vect<T> &))
{
  double      xh,hh,h6;
  ss_vect<T>  dym, dyt, yt;

  hh = h*0.5; h6 = h/6.0;
  xh = x + hh; yt = y + hh*dydx;
  (*derivs)(Cell, xh, yt, dyt); yt = y + hh*dyt;
  (*derivs)(Cell, xh, yt, dym); yt = y + h*dym; dym += dyt;
  (*derivs)(Cell, x+h, yt, dyt);
  yout = y + h6*(dydx+dyt+2.0*dym);
}


template<typename T>
inline T get_p_s_ps(const ss_vect<T> &x, const T &qop_Ax, const T &qop_Ay)
{
  // Compute p_s in phase space (with vector potential in axial gauge).

  return sqrt(sqr(1e0+x[delta_])-sqr(x[px_]-qop_Ax)-sqr(x[py_]-qop_Ay));
}


template<typename T>
void f_FM(const CellType &Cell, const double z, const ss_vect<T> &ps,
	  ss_vect<T> &Dps)
{
  // Coordinates are: [x, x', y, y', -ct, delta].

  int           j, kz;
  T             BoBrho[3], p_s;
  FieldMapType  *FM;

  const double  eps = 1e-5;


  FM = (FieldMapType*)Cell.Elem;

  kz = 0;
  for (j = 1; j <= FM->n[Z_]; j++)
    if (fabs(z-FM->x[Z_][j]) < eps) {
      kz = j;
      break;
    }

  if (kz == 0) {
    cout << fixed << setprecision(10)
	 << "z = " << setw(12) << z << " undefined" << endl;

    for (j = 0; j < ss_dim; j++)
      Dps[j] = NAN;

    return;
  }

  splin2_(FM->x[X_], FM->x[Y_], FM->BoBrho[X_][kz], FM->BoBrho2[X_][kz],
	   FM->n[X_], FM->n[Y_], ps[x_], ps[y_], BoBrho[X_]);

  if (BoBrho[X_] == NAN) {
    for (j = 0; j < ss_dim; j++)
      Dps[j] = NAN;
    return;
  }

  splin2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Y_][kz], FM->BoBrho2[Y_][kz],
	   FM->n[X_], FM->n[Y_], ps[x_], ps[y_], BoBrho[Y_]);

  if (BoBrho[Y_] == NAN) {
    for (j = 0; j < ss_dim; j++)
      Dps[j] = NAN;
    return;
  }

  splin2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Z_][kz], FM->BoBrho2[Z_][kz],
	   FM->n[X_], FM->n[Y_], ps[x_], ps[y_], BoBrho[Z_]);

  if (BoBrho[Z_] == NAN) {
    for (j = 0; j < ss_dim; j++)
      Dps[j] = NAN;
    return;
  }

  p_s = get_p_s_cs(ps);

  Dps[x_] = ps[px_];
  Dps[px_] =
    -(ps[px_]*ps[py_]*BoBrho[X_]-(1e0+sqr(ps[px_]))*BoBrho[Y_]
      +ps[py_]*BoBrho[Z_])/p_s;

  Dps[y_] = ps[py_];
  Dps[py_] =
    -((1e0+sqr(ps[py_]))*BoBrho[X_]-ps[px_]*ps[py_]*BoBrho[Y_]
      -ps[px_]*BoBrho[Z_])/p_s;

  Dps[ct_] = (1e0+ps[delta_])/p_s - 1e0;

  if (globval.pathlength) Dps[ct_] += 1e0;

  Dps[delta_] = 0e0;
}


template<typename T>
inline T get_p_s_cs(const ss_vect<T> &x)
{
  // Compute p_s in configuration space.

  return (1e0+x[delta_])/sqrt(1e0+sqr(x[px_])+sqr(x[py_]));
}


template<typename T>
void FieldMap_pass_RK(CellType &Cell, ss_vect<T> &ps, int k)
{
  int           i;
  double        h, z;
  T             p_s;
  ss_vect<T>    Dps;
  FieldMapType  *FM;

  const int     n_step = 2; // Each step needs: f(z_n), f(z_n+h), f(z_n+2h)

  FM = Cell.Elem.FM;

  switch (FieldMap_filetype) {
  case 2:
  case 3:
  case 4:
    // Transform to right handed system
    ps[x_] = -ps[x_]; ps[px_] = -ps[px_];
    break;
  }

  // transform to [x, x', y, y']
  p_s = get_p_s(ps); ps[px_] /= p_s; ps[py_] /= p_s;

  h = n_step*FM->dx[Z_]; z = FM->x[Z_][1]; FM->Lr = 0e0;
  if (trace)
    outf_ << scientific << setprecision(3)
	  << setw(11) << s_FM
	  << setw(11) << is_double<T>::cst(ps[x_])
	  << setw(11) << is_double<T>::cst(ps[px_])
	  << setw(11) << is_double<T>::cst(ps[y_])
	  << setw(11) << is_double<T>::cst(ps[py_])
	  << setw(11) << is_double<T>::cst(ps[delta_])
	  << setw(11) << is_double<T>::cst(ps[ct_])
	  << endl;
  for(i = 1+FM->cut; i < FM->n[Z_]-FM->cut; i += n_step) {
    if (i <= FM->n[Z_]-FM->cut-2) {
      f_FM(Cell, z, ps, Dps);

      if (Dps[x_] == NAN) {
	cout << "FieldMap_pass_RK: particle lost" << endl;
	cout << ps;
	return;
      }

      rk4_(Cell, ps, Dps, FM->x[Z_][i], h, ps, f_FM);

      z += h; FM->Lr += h; s_FM += h;
    } else {
      // Use 2nd order Runge-Kutta (aka Midpoint method)
      f_FM(Cell, z, ps, Dps);

      if (Dps[x_] == NAN) {
	cout << "FieldMap_pass_RK: particle lost" << endl;
	cout << ps;
	return;
      }

      ps += h/2.0*Dps;

      z += h/2.0; FM->Lr += h/2.0; s_FM += h/2.0;
    }

    if (trace)
      outf_ << scientific << setprecision(3)
	    << setw(11) << s_FM
	   << setw(11) << is_double<T>::cst(ps[x_])
	   << setw(11) << is_double<T>::cst(ps[px_])
	   << setw(11) << is_double<T>::cst(ps[y_])
	   << setw(11) << is_double<T>::cst(ps[py_])
	   << setw(11) << is_double<T>::cst(ps[delta_])
	   << setw(11) << is_double<T>::cst(ps[ct_])
	   << endl;
  }

  // transform back to [x, px, y, py] (A_x,y,z = 0)
  p_s = get_p_s_cs(ps); ps[px_] *= p_s; ps[py_] *= p_s;

  switch (FieldMap_filetype) {
  case 2:
  case 3:
  case 4:
    // Transform back to left handed system
    ps[x_] = -ps[x_]; ps[px_] = -ps[px_];
    break;
  }
}


template<typename T>
void FieldMap_pass_SI(CellType &Cell, ss_vect<T> &ps, int k)
{
  /* E. Chacon-Golcher, F. Neri "A Symplectic Integrator with Arbitrary
     Vector and Scalar Potentials" Phys. Lett. A 372 p. 4661-4666 (2008).    */

  int           i, j;
  double        h, z;
  T             hd, AoBrho[2], dAoBrho[2], AoBrho_int, ByoBrho;
  ss_vect<T>    ps1;
  FieldMapType  *FM;

  const int     n_step = 2;
  const double  d_diff = 1e0;


  FM = Cell.Elem.FM;

  switch (FieldMap_filetype) {
  case 2:
  case 3:
  case 4:
    // Transform to right handed system
    ps[x_] = -ps[x_]; ps[px_] = -ps[px_];
    break;
  }

  h = n_step*FM->dx[Z_]; z = 0e0; FM->Lr = 0e0;
  if (trace)
    outf_ << scientific << setprecision(3)
	  << setw(11) << s_FM
	  << setw(11) << is_double<T>::cst(ps[x_])
	  << setw(11) << is_double<T>::cst(ps[px_])
	  << setw(11) << is_double<T>::cst(ps[y_])
	  << setw(11) << is_double<T>::cst(ps[py_])
	  << setw(11) << is_double<T>::cst(ps[delta_])
	  << setw(11) << is_double<T>::cst(ps[ct_])
	  << endl;
  for (i = 1+FM->cut; i < FM->n[Z_]-FM->cut; i += n_step) {
    hd = h/(1e0+ps[delta_]);

    // 1. Half step in z.
    z += 0.5*h; j = i + 1; s_FM += 0.5*h;

    // 2. Half drift in y.
    ps1 = ps;

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps[x_], ps[y_], AoBrho[0]);

    if (AoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    ps1[y_] += (ps[py_]-AoBrho[0])*0.5*hd;
    ps1[ct_] += sqr(0.5)*hd*sqr(ps[py_]-AoBrho[0])/(1e0+ps[delta_]);

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps[x_]+d_diff*FM->dx[X_], ps[y_],
	    dAoBrho[1]);

    if (dAoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps[x_]-d_diff*FM->dx[X_], ps[y_],
	    dAoBrho[0]);

    if (dAoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    AoBrho_int = (dAoBrho[1]-dAoBrho[0])/(2.0*d_diff*FM->dx[X_]);

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps1[x_]+d_diff*FM->dx[X_], ps1[y_],
	    dAoBrho[1]);

    if (dAoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps1[x_]-d_diff*FM->dx[X_], ps1[y_],
	    dAoBrho[0]);

    if (dAoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    AoBrho_int += (dAoBrho[1]-dAoBrho[0])/(2.0*d_diff*FM->dx[X_]);

    // Trapezoidal rule
    ps1[px_] +=
      AoBrho_int*(is_double<T>::cst(ps1[y_])-is_double<T>::cst(ps[y_]))/2.0;

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps1[x_], ps1[y_], AoBrho[1]);

    if (AoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    ps1[py_] += AoBrho[1] - AoBrho[0];

    ps = ps1;

    // 3. Full drift in x.
    ps1 = ps;

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[X_][j], FM->AoBrho2[X_][j],
	    FM->n[X_], FM->n[Y_], ps[x_], ps[y_], AoBrho[0]);

    if (AoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    ps1[x_] += (ps[px_]-AoBrho[0])*hd;
    ps1[ct_] += 0.5*hd*sqr(ps[px_]-AoBrho[0])/(1e0+ps[delta_]);

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[X_][j], FM->AoBrho2[X_][j],
	    FM->n[X_], FM->n[Y_], ps1[x_], ps1[y_], AoBrho[1]);

    if (AoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    ps1[px_] += AoBrho[1] - AoBrho[0];

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[X_][j], FM->AoBrho2[X_][j],
	    FM->n[X_], FM->n[Y_], ps[x_], ps[y_]+d_diff*FM->dx[Y_],
	    dAoBrho[1]);

    if (dAoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[X_][j], FM->AoBrho2[X_][j],
	    FM->n[X_], FM->n[Y_], ps[x_], ps[y_]-d_diff*FM->dx[Y_],
	    dAoBrho[0]);

    if (dAoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    AoBrho_int = (dAoBrho[1]-dAoBrho[0])/(2.0*d_diff*FM->dx[Y_]);

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[X_][j], FM->AoBrho2[X_][j],
	    FM->n[X_], FM->n[Y_], ps1[x_], ps1[y_]+d_diff*FM->dx[Y_],
	    dAoBrho[1]);

    if (dAoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[X_][j], FM->AoBrho2[X_][j],
	    FM->n[X_], FM->n[Y_], ps1[x_], ps1[y_]-d_diff*FM->dx[Y_],
	    dAoBrho[0]);

    if (dAoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    AoBrho_int += (dAoBrho[1]-dAoBrho[0])/(2.0*d_diff*FM->dx[Y_]);

    // Trapezoidal rule
    ps1[py_] +=
      AoBrho_int*(is_double<T>::cst(ps1[x_])-is_double<T>::cst(ps[x_]))/2.0;

    ps = ps1;

    // 4. Half drift in y.
    ps1 = ps;

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps[x_], ps[y_], AoBrho[0]);

    if (AoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    ps1[y_] += (ps[py_]-AoBrho[0])*0.5*hd;
    ps1[ct_] += sqr(0.5)*hd*sqr(ps[py_]-AoBrho[0])/(1e0+ps[delta_]);

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps[x_]+d_diff*FM->dx[X_], ps[y_],
	    dAoBrho[1]);

    if (dAoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps[x_]-d_diff*FM->dx[X_], ps[y_],
	    dAoBrho[0]);

    if (dAoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    AoBrho_int = (dAoBrho[1]-dAoBrho[0])/(2.0*d_diff*FM->dx[X_]);

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps1[x_]+d_diff*FM->dx[X_], ps1[y_],
	    dAoBrho[1]);

    if (dAoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps1[x_]-d_diff*FM->dx[X_], ps1[y_],
	    dAoBrho[0]);

    if (dAoBrho[0] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    AoBrho_int += (dAoBrho[1]-dAoBrho[0])/(2.0*d_diff*FM->dx[X_]);

    // Trapezoidal rule
    ps1[px_] +=
      AoBrho_int*(is_double<T>::cst(ps1[y_])-is_double<T>::cst(ps[y_]))/2.0;

    splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	    FM->n[X_], FM->n[Y_], ps1[x_], ps1[y_], AoBrho[1]);

    if (AoBrho[1] == NAN) {
      for (j = 0; j < ss_dim; j++)
	ps[j] = NAN;
      return;
    }

    ps1[py_] += AoBrho[1] - AoBrho[0];

    ps = ps1;

    // 5. Half step in z.
    z += 0.5*h; j = i + 2; s_FM += 0.5*h;

    if (globval.pathlength) ps[ct_] += h;

    FM->Lr += h;

    if (globval.radiation || globval.emittance) {
//      B[X_] = -AoBrhoy[3]; B[Y_] = AoBrho[X_][3];
//      B[Z_] = AoBrhoy[1] - AoBrho[X_][2];
//      radiate(ps, h, 0e0, B);
    }

    if (trace) {
      splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[X_][j], FM->AoBrho2[X_][j],
	      FM->n[X_], FM->n[Y_], ps[x_], ps[y_], AoBrho[0]);
      splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	      FM->n[X_], FM->n[Y_], ps[x_], ps[y_], AoBrho[1]);
      splin2_(FM->x[X_], FM->x[Y_], FM->BoBrho[Y_][j], FM->BoBrho2[Y_][j],
	      FM->n[X_], FM->n[Y_], ps[x_], ps[y_], ByoBrho);

      outf_ << scientific << setprecision(3)
	   << setw(11) << s_FM
	   << setw(11) << is_double<T>::cst(ps[x_])
	   << setw(11) << is_double<T>::cst(ps[px_]-AoBrho[0])
	   << setw(11) << is_double<T>::cst(ps[y_])
	   << setw(11) << is_double<T>::cst(ps[py_]-AoBrho[1])
	   << setw(11) << is_double<T>::cst(ps[delta_])
	   << setw(11) << is_double<T>::cst(ps[ct_])
	   << setw(11) << is_double<T>::cst(AoBrho[0])
	   << setw(11) << is_double<T>::cst(AoBrho[1])
	   << setw(11) << is_double<T>::cst(ByoBrho)
	   << endl;
    }
  }

  // Change of gauge
  splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[X_][j], FM->AoBrho2[X_][j],
	  FM->n[X_], FM->n[Y_], ps[x_], ps[y_], AoBrho[0]);

  if (AoBrho[0] == NAN) {
    for (j = 0; j < ss_dim; j++)
      ps[j] = NAN;
    return;
  }

  ps[px_] -= AoBrho[0];

  splin2_(FM->x[X_], FM->x[Y_], FM->AoBrho[Y_][j], FM->AoBrho2[Y_][j],
	  FM->n[X_], FM->n[Y_], ps[x_], ps[y_], AoBrho[0]);

  if (AoBrho[0] == NAN) {
    for (j = 0; j < ss_dim; j++)
      ps[j] = NAN;
    return;
  }

  ps[py_] -= AoBrho[0];

  switch (FieldMap_filetype) {
  case 2:
  case 3:
  case 4:
    // Transform back to left handed system
    ps[x_] = -ps[x_]; ps[px_] = -ps[px_];
    break;
  }
}


template<typename T>
void FieldMap_Pass(CellType &Cell, ss_vect<T> &ps)
{
  int           k;
  FieldMapType  *FM;


  if (trace & first_FM) {
    file_wr(outf_, "FieldMap_pass.dat");
    s_FM = 0e0;
    first_FM = false;
  }

  FM = Cell.Elem.FM;

//  GtoL(ps, Cell.dS, Cell.droll, 0e0, 0e0, 0e0);

  ps[x_] += FM->x0;

  Drift(FM->L1/2.0, ps);

  p_rot(FM->phi/2.0*180e0/M_PI, ps);

  for (k = 1; k <= FM->n_step; k++) {
    if (sympl)
      FieldMap_pass_SI(Cell, ps, k);
    else
      FieldMap_pass_RK(Cell, ps, k);
  }

  Drift(FM->Ld, ps);

  p_rot(FM->phi/2.0*180e0/M_PI, ps);

  Drift(FM->L1/2.0, ps);

  ps[x_] -= FM->x0;

//  LtoG(ps, Cell.dS, Cell.droll, 0e0, 0e0, 0e0);

//  outf_.close();
}


template<typename T>
void Insertion_Pass(CellType &Cell, ss_vect<T> &x)
{
  /* Purpose:
       Track vector x through a insertion
       If radiation or cavity on insertion is like a drift

   Input:
       Cell element to track through
       x initial coordinates vector

   Output:
       x final coordinates vector

   Return:
       none

   Global variables:
       none

   Specific functions:
       LinearInterpolation2
       Drft
       CopyVec

   Comments:
       Outside of interpolation table simulated by putting 1 in x[4]
       01/07/03 6D tracking activated
       10/01/05 First order kick part added                                  */

  ElemType  *elem;
  double    LN = 0e0;
  T         tx2, tz2;      /* thetax and thetaz retrieved from
			      interpolation routine */
  T         d, B2_perp;
  double    alpha0 = 0e0;  // 1/ brh0
  double    alpha02= 0e0;  // alpha square
  int       Nslice = 0;
  int       i = 0;
  bool      outoftable = false;

  elem  = &Cell.Elem; Nslice = elem->ID->n;

  if (elem->ID->linear) {
    alpha0 = c0/globval.Energy*1E-9*elem->ID->scaling; alpha02 = alpha0*alpha0;
  } else
    alpha02 = 1e-6*elem->ID->scaling;

//  /* Global -> Local */
//  GtoL(X, Cell->dS, Cell->droll, 0e0, 0e0, 0e0);

  // (Nslice+1) drifts, nslice kicks
  // LN = elem->L/(Nslice+1);

  // Nslice drifts and kicks.
  LN = elem->L/Nslice;
  Drift(LN/2e0, x);

  for (i = 1; i <= Nslice; i++) {
    // Second order kick
    if (elem->ID->secondorder){
      // if (!elem->ID->linear)
      //   SplineInterpolation2(x[x_], x[y_], tx2, tz2, Cell, outoftable);
      // else {
        LinearInterpolation2(x[x_], x[y_], tx2, tz2, B2_perp, Cell,
			     outoftable, 2);

	// Scale locally with (Brho) (as above) instead of when the file
	// is read; since the beam energy might not be known at that time.
	if (globval.radiation || globval.emittance)
	  radiate_ID(x, LN, elem->ID->scaling*B2_perp);
      // }

      if (outoftable) {
	x[x_] = NAN;
        return;
      }

      d = alpha02/Nslice/(1e0+x[delta_]); x[px_] += d*tx2; x[py_] += d*tz2;
    }
    if (i != Nslice) Drift(LN, x);
  }

  Drift(LN/2e0, x);

//  CopyVec(6L, x, Cell->BeamPos);

//  /* Local -> Global */
//  LtoG(X, Cell->dS, Cell->droll, 0e0, 0e0, 0e0);
}


template<typename T>
void sol_pass(const ElemType &elem, ss_vect<T> &x)
{
  int     i;
  double  h, z;
  T       hd, AxoBrho, AyoBrho, dAxoBrho[3], dAyoBrho[3], dpy, dpx, B[3];

  h = elem.L/elem.Sol->n; z = 0e0;

  for (i = 1; i <= elem.Sol->n; i++) {
    hd = h/(1e0+x[delta_]);

    // 1: half step in z
    z += 0.5*h;

    // 2: half drift in y
    AyoBrho = elem.Sol->BoBrho*x[x_]/2.0; dpx = elem.Sol->BoBrho*x[y_]/2.0;
//    get_Axy_EF3(W, z, x, AyoBrho, dAyoBrho, dpx, false);

    x[px_] -= dpx; x[py_] -= AyoBrho;
    x[y_] += 0.5*hd*x[py_];
    x[ct_] += sqr(0.5)*hd*sqr(x[py_])/(1e0+x[delta_]);

    AyoBrho = elem.Sol->BoBrho*x[x_]/2.0; dpx = elem.Sol->BoBrho*x[y_]/2.0;
//    get_Axy_EF3(W, z, x, AyoBrho, dAyoBrho, dpx, false);

    x[px_] += dpx; x[py_] += AyoBrho;

    // 3: full drift in x
    AxoBrho = -elem.Sol->BoBrho*x[y_]/2.0; dpy = -elem.Sol->BoBrho*x[x_]/2.0;
//    get_Axy_EF3(W, z, x, AxoBrho, dAxoBrho, dpy, true);

    x[px_] -= AxoBrho; x[py_] -= dpy; x[x_] += hd*x[px_];
    x[ct_] += 0.5*hd*sqr(x[px_])/(1e0+x[delta_]);

    if (globval.pathlength) x[ct_] += h;

    AxoBrho = -elem.Sol->BoBrho*x[y_]/2.0; dpy = -elem.Sol->BoBrho*x[x_]/2.0;
//    get_Axy_EF3(W, z, x, AxoBrho, dAxoBrho, dpy, true);

    x[px_] += AxoBrho; x[py_] += dpy;

    // 4: a half drift in y
    AyoBrho = elem.Sol->BoBrho*x[x_]/2.0; dpx = elem.Sol->BoBrho*x[y_]/2.0;
//    get_Axy_EF3(W, z, x, AyoBrho, dAyoBrho, dpx, false);

    x[px_] -= dpx; x[py_] -= AyoBrho;
    x[y_] += 0.5*hd*x[py_];
    x[ct_] += sqr(0.5)*hd*sqr(x[py_])/(1e0+x[delta_]);

    AyoBrho = elem.Sol->BoBrho*x[x_]/2.0; dpx = elem.Sol->BoBrho*x[y_]/2.0;
//    get_Axy_EF3(W, z, x, AyoBrho, dAyoBrho, dpx, false);

    x[px_] += dpx; x[py_] += AyoBrho;

    // 5: half step in z
    z += 0.5*h;

    if (globval.radiation || globval.emittance) {
      dAxoBrho[X_] = 0e0;
      dAxoBrho[Y_] = -elem.Sol->BoBrho/2.0;
      dAxoBrho[Z_] = 0e0;
      dAyoBrho[X_] = elem.Sol->BoBrho/2.0;
      dAyoBrho[Y_] = 0e0;
      dAyoBrho[Z_] = 0e0;
//      get_Axy_EF3(W, z, x, AyoBrho, dAyoBrho, dpx, false);
//      get_Axy_EF3(W, z, x, AxoBrho, dAxoBrho, dpy, true);
      B[X_] = -dAyoBrho[Z_]; B[Y_] = dAxoBrho[Z_];
      B[Z_] = dAyoBrho[X_] - dAxoBrho[Y_];
      radiate(x, h, 0e0, B);
    }
  }
}


template<typename T>
void Solenoid_Pass(CellType &Cell, ss_vect<T> &ps)
{

  GtoL(ps, Cell.dS, Cell.droll, 0e0, 0e0, 0e0);

  sol_pass(Cell.Elem, ps);

  LtoG(ps, Cell.dS, Cell.droll, 0e0, 0e0, 0e0);
}


void getelem(long i, CellType *cellrec) { *cellrec = Cell[i]; }

void putelem(long i, CellType *cellrec) { Cell[i] = *cellrec; }


int GetnKid(const int Fnum1) { return (ElemFam[Fnum1-1].nKid); }


long Elem_GetPos(const int Fnum1, const int Knum1)
{
  long int  loc;

  if (ElemFam[Fnum1-1].nKid != 0)
    loc = ElemFam[Fnum1-1].KidList[Knum1-1];
  else {
    loc = -1;
    printf("Elem_GetPos: there are no kids in family %d (%s)\n",
	   Fnum1, ElemFam[Fnum1-1].Elem.name);
    exit_(0);
  }

  return loc;
}


static double thirdroot(double a)
{
  /* By substitution method */
  int    i;
  double x;

  x = 1e0; i = 0;
  do {
    i++; x = (x+a)/(x*x+1e0);
  } while (i != 250);
  return x;
}


void SI_init(void)
{
  // SI units are used internally
  // apart from globval.energy [GeV]
  /*  c_1 = 1/(2*(2-2^(1/3))),    c_2 = (1-2^(1/3))/(2*(2-2^(1/3)))
      d_1 = 1/(2-2^(1/3)),        d_2 = -2^(1/3)/(2-2^(1/3))                 */

  double C_gamma, C_u;

  c_1 = 1e0/(2e0*(2e0-thirdroot(2e0))); c_2 = 0.5e0 - c_1;
  d_1 = 2e0*c_1; d_2 = 1e0 - 2e0*d_1;

  // classical radiation
  // C_gamma = 4*pi*r_e [m]/(3*(m_e [GeV/c^2] *c^2)^3)
  C_gamma = 4.0*M_PI*r_e/(3.0*cube(1e-9*m_e));
  // P_gamma = e^2*c^3/(2*pi)*C_gamma*(E [GeV])^2*(B [T])^2
  // p_s = P_s/P, E = P*c, B/(Brho) = p/e
  cl_rad = C_gamma*cube(globval.Energy)/(2.0*M_PI);

  // eletron rest mass [GeV]: slightly off???
//  m_e_ = 0.5110034e-03;
  // quantum fluctuations
  C_u = 55.0/(24.0*sqrt(3.0));
  q_fluct =
    3.0*C_u*C_gamma*1e-9*h_bar*c0/(4.0*M_PI*cube(1e-9*m_e))
    *pow(globval.Energy, 5.0);
}


static void Mpole_Print(FILE *f, int Fnum1)
{
  ElemType  *elem;
  MpoleType *M;

  elem  = &ElemFam[Fnum1-1].Elem; M = elem->M;
  fprintf(f, "Element[%3d ] \n", Fnum1);
  fprintf(f, "   Name: %.*s,  kind:   mpole,  L=% .8E\n",
          SymbolLength, elem->name, elem->L);
  fprintf(f, "   Method: %d, N=%4d\n", M->method, M->n);
}


static void Drift_Print(FILE *f, int Fnum1)
{
  ElemFamType *elamfam;
  ElemType    *elem;

  elamfam = &ElemFam[Fnum1-1]; elem = &elamfam->Elem;
  fprintf(f, "Element[%3d ] \n", Fnum1);
  fprintf(f, "   Name: %.*s,  kind:   drift,  L=% .8E\n",
          SymbolLength, elem->name, elem->L);
  fprintf(f, "   nKid:%3d\n\n", elamfam->nKid);
}


static void Wiggler_Print(FILE *f, int Fnum1)
{
  ElemType *elem;

  elem = &ElemFam[Fnum1-1].Elem;
  fprintf(f, "Element[%3d ] \n", Fnum1);
  fprintf(f, "   Name: %.*s,  kind:   wiggler,  L=% .8E\n\n",
          NameLength, elem->name, elem->L);
}


static void Insertion_Print(FILE *f, int Fnum1)
{
  ElemType *elem;

  elem = &ElemFam[Fnum1-1].Elem;
  fprintf(f, "Element[%3d ] \n", Fnum1);
  fprintf(f, "   Name: %.*s,  kind:   wiggler,  L=% .8E\n\n",
          SymbolLength, elem->name, elem->L);
}


void Elem_Print(FILE *f, int Fnum1)
{
  int i;

  if (Fnum1 == 0) {
    // print all elements
    for (i = 1; i <= globval.Elem_nFam; i++)
      Elem_Print(f, i);
    return;
  }

  switch (ElemFam[Fnum1-1].Elem.kind) {
  case drift:
    Drift_Print(f, Fnum1);
    break;

  case Mpole:
    Mpole_Print(f, Fnum1);
    break;
  case Wigl:
    Wiggler_Print(f, Fnum1);
    break;
  case FieldMap:
    break;
  case Insertion:
    Insertion_Print(f, Fnum1);
    break;
  case Cavity:
    break;
  case marker:
    break;
  case Spreader:
    break;
  case Recombiner:
    break;
  case Solenoid:
    break;
  case undef:
    break;
  }
}


double Mpole_Getbn(int Fnum1, int Knum1, int Order);


double Elem_GetKval(int Fnum1, int Knum1, int Order)
{
  double        Result = 0e0;
  ElemType      *elem;
  MpoletypeType *M;

  if (Fnum1 > 0) {
    elem = &Cell[ElemFam[Fnum1-1].KidList[Knum1-1]].Elem;
    switch (elem->kind) {
    case drift:
      Result = 0e0;
      break;
    case marker:
      Result = 0e0;
      break;
    case Cavity:
      Result = 0e0;
      break;
    case Mpole: /* KL*/
      M = (MpoleType*)elem;
      if (elem->thick == thick)
	Result = elem->L*Mpole_Getbn(Fnum1, Knum1, Order);
      else
	Result = Mpole_Getbn(Fnum1, Knum1, Order);
      break;
    case Wigl:
      Result =
	elem->L*sqrt(2.0
	*Cell[ElemFam[Fnum1-1].KidList[Knum1-1]].Elem.W->bn[Order+HOMmax]);
      break;
    case FieldMap:
      Result = 0e0;
      break;
    case Insertion:
      Result = 0e0;
      break;
    case Spreader:
      Result = 0e0;
      break;
    case Recombiner:
      Result = 0e0;
      break;
    case Solenoid:
      Result = 0e0;
      break;
    case undef:
      break;
    }
  } else
    Result = 0e0;

  return Result;
}


#define n 4
void LinsTrans(Matrix &A, Vector &b)
{
  int    j;
  Vector c;

  CopyVec(n, b, c); /* c=b */
  LinTrans(n, A, c); /* c<-A*c */
  for (j = 0; j < n; j++)
    c[j] += A[j][n]*b[n] + A[n][j];
  CopyVec(n, c, b); /* b=c */
}
#undef n


#define n 4
void MulLsMat(Matrix &A, Matrix &B)
{
  int    i, k;
  Matrix C;

  CopyMat(n, B, C); /* C<-B */
  MulLMat(n, A, C); /* C<-A*C */
  for (i = 0; i < n; i++) {
    C[i][n] = A[i][n]; C[n][i] = 0e0;
    for (k = 0; k < n; k++) {
      C[i][n] += A[i][k]*B[k][n];
      C[n][i] += A[i][k]*B[n][k];
    }
  }
  C[n][n] = 1e0;
  CopyMat(n+1, C, B); /* B<-C */
}
#undef n


void Mpole_Setbn(int Fnum, int Knum, int Order)
{
  /*  called by Cell_SetdP
       Compute full multipole composent as sum of design, systematic
       and random part
       Set multipole order to Order if multipole (Order >2)                  */

  CellType *cell;  /* pointer on the Cell */
  ElemType *elem; /* pointer on the Elemetype */
  MpoleType *M;/* Pointer on the Multipole */

  M = (MpoleType*)Cell[ElemFam[Fnum-1].KidList[Knum-1]].Elem;
  M->bn[Order+HOMmax] =
    M->bnpar[Order+HOMmax] + M->bnsys[Order+HOMmax] +
    M->bnrms[Order+HOMmax]*M->bnrnd[Order+HOMmax];
  if (abs(Order) > M->order && M->bn[Order+HOMmax] != 0e0)
    M->order = abs(Order);
}


double Mpole_Getbn(int Fnum, int Knum, int Order)
{
  /*  Return multipole strength (of order Order) for Knum element of
      family Fnum
       Order =  2 for normal quadrupole
             = -2 for skew quadrupole                                        */

  MpoleType *M; /* Pointer on the multipole */

  M = (MpoleType*)Cell[ElemFam[Fnum-1].KidList[Knum-1]].Elem;
  return (M->bn[Order+HOMmax]);
}


void Mpole_Defbnpar(int Fnum, int Knum, int Order, double bnpar)
{
  ElemType   *elem;
  MpoleType  *M;

  M = (MpoleType*)Cell[ElemFam[Fnum-1].KidList[Knum-1]].Elem;
  M->bnpar[Order+HOMmax]=bnpar;
}


void Mpole_Defbnsys(int Fnum, int Knum, int Order, double bnsys)
{
  /*Fnum, Knum, Order : integer*/
  ElemType *elem;
  MpoleType *M;

  elem = (MpoleType*)Cell[ElemFam[Fnum-1].KidList[Knum-1]].Elem;
  M->bnsys[Order+HOMmax]=bnsys;
}


void Mpole_SetdS(int Fnum, int Knum)
{
  int       j;
  CellType  *cell;
  ElemType  *elem;
  MpoleType *M;

  cell = &Cell[ElemFam[Fnum-1].KidList[Knum-1]];
  M = (MpoleType*)cell->Elem;
  for (j = 0; j <= 1; j++)
    cell->dS[j] = M->dSsys[j] + M->dSrms[j]*M->dSrnd[j];
}

void Mpole_Setdroll(int Fnum, int Knum)
{
  CellType  *cell;
  MpoleType *M;

  cell = &Cell[ElemFam[Fnum-1].KidList[Knum-1]];
  M = (MpoleType*)cell->Elem;
  cell->droll[0] =
    cos(dtor(M->rollpar + M->rollsys + M->rollrms*M->rollrnd));
  cell->droll[1] = sin(
      dtor(M->rollpar + M->rollsys + M->rollrms*M->rollrnd));
  /* Calculate simplified p_rots */
  M->c0 = sin(elem->L*M->irho/2e0);
  M->c1 = cos(dtor(M->rollpar))*M->c0;
  M->s1 = sin(dtor(M->rollpar))*M->c0;
}


double Mpole_Getdroll(int Fnum, int Knum)
{
  MpoleType *M;

  M = (MpoleType*)Cell[ElemFam[Fnum-1].KidList[Knum-1]].Elem;

  return(M->rollpar + M->rollsys + M->rollrms*M->rollrnd);
}


void Mpole_Defdrollpar(int Fnum, int Knum, double rollpar)
{
  MpoleType *M;

  M = (MpoleType*)Cell[ElemFam[Fnum-1].KidList[Knum-1]].Elem;
  M->rollpar = rollpar;
}


void Mpole_Defdrollsys(int Fnum, int Knum, double rollsys)
{
  MpoleType *M;

  M = (MpoleType*)Cell[ElemFam[Fnum-1].KidList[Knum-1]].Elem; M = elem->M;

  M->rollsys=rollsys;
}


void Wiggler_Setbn(int Fnum, int Knum, int Order)
{
  CellType     *cell;
  WigglerType  *W;

  W = (WigglerType*)Cell[ElemFam[Fnum-1].KidList[Knum-1]].Elem;
  if (abs(Order) > W->order)
    W->order = abs(Order);
}


void Wiggler_SetdS(int Fnum, int Knum)
{
  int         j;
  CellType    *cell;
  ElemType    *elem;
  WigglerType *W;

  cell = &Cell[ElemFam[Fnum-1].KidList[Knum-1]];
  W = (WigglerType*)cell->Elem;
  for (j = 0; j <= 1; j++)
    cell->dS[j] = W->dSsys[j] + W->dSrms[j]*W->dSrnd[j];
}

void Wiggler_Setdroll(int Fnum, int Knum)
{
  CellType    *cell;
  WigglerType *W;

  cell = &Cell[ElemFam[Fnum-1].KidList[Knum-1]];
  W = (WigglerType*)cell->Elem;
  cell->droll[0] = cos(dtor(W->rollpar+W->rollsys+W->rollrms*W->rollrnd));
  cell->droll[1] = sin(dtor(W->rollpar+W->rollsys+W->rollrms*W->rollrnd));
}


// instantiate

template void spline_(const double [], const double [], const int,
		      const double, const double, double []);

template void spline_(const double [], const tps [], const int,
		      const double, const double, tps []);

template void splint_(const double [], const double [], const double [],
		      const int, const double &, double &);

template void splint_(const double [], const tps [], const tps [], const int,
		      const tps &, tps &);

template void splint_(const double [], const double [], const double [],
		      const int, const tps &, tps &);

template void splin2_(const double [], const double [], double **, double **,
		      const int, const int, const double &, const double &,
		      double &);

template void splin2_(const double [], const double [], double **, double **,
		      const int, const int, const tps &, const tps &, tps &);


template void f_FM(const CellType &, const double, const ss_vect<double> &,
		   ss_vect<double> &);
template void f_FM(const CellType &, const double, const ss_vect<tps> &,
		   ss_vect<tps> &);
template void rk4_(const CellType &, const ss_vect<double> &,
		   const ss_vect<double>, const double, const double,
		   ss_vect<double> &,
		   void (*derivs)(const CellType &, const double,
				  const ss_vect<double> &, ss_vect<double> &));
template void rk4_(const CellType &, const ss_vect<tps> &, const ss_vect<tps>,
		   const double, const double, ss_vect<tps> &,
		   void (*derivs)(const CellType &, const double,
				  const ss_vect<tps> &, ss_vect<tps> &));
template void FieldMap_pass_RK(CellType &, ss_vect<double> &, int k);
template void FieldMap_pass_RK(CellType &, ss_vect<tps> &, int k);
template void FieldMap_pass_SI(CellType &, ss_vect<double> &, int k);
template void FieldMap_pass_SI(CellType &, ss_vect<tps> &, int k);
