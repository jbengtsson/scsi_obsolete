/* Tracy-2

   J. Bengtsson, CBP, LBL      1990 - 1994   Pascal version
                 SLS, PSI      1995 - 1997
   M. Boege      SLS, PSI      1998          C translation
   L. Nadolski   SOLEIL        2002          Link to NAFF, Radia field maps
   J. Bengtsson  NSLS-II, BNL  2004 -


   To generate a lattice flat file.

   Type codes:

     marker     -1
     drift	 0
     multipole   1
     cavity      2
     thin kick   3
     wiggler     4

   Integration methods:

     linear, matrix style (obsolete)              0
     2nd order symplectic integrator (obsolete)   2
     4th order symplectic integrator              4

   Format:

     name, family no, kid no, element no
     type code, integration method, no of integration steps
     apertures: xmin, xmax, ymin, ymax

   The following lines follows depending on element type.

     type

     drift:	 L

     multipole:  hor., ver. displ., roll angle (design), roll angle (error)
                 L, 1/rho, entrance angle, exit angle
		 no of nonzero multipole coeff., n design
		 n, b , a
		     n   n
		     .
		     .
		     .

     wiggler:    L [m], lambda [m]
                 no of harmonics
                 harm no, kxV [1/m], BoBrhoV [1/m], kxH, BoBrhoH, phi
                    ...

     cavity:	 cavity voltage/beam energy [eV], omega/c, beam energy [eV]

     thin kick:	 hor., ver. displacement, roll angle (total)
		 no of nonzero multipole coeff.
		 n, b , a
		     n   n
		     .
		     .
		     .

     kick_map:   scale order <file name>

*/


#define snamelen   10

// numerical type codes
#define marker_   -1
#define drift_     0
#define mpole_     1
#define cavity_    2
#define thinkick_  3
#define wiggler_   4
#define insertion_ 6


ifstream  inf;


void get_kind(const int kind, ElemType *Elem)
{
  MarkerType    *Mrk;
  DriftType     *D;
  MpoleType     *M;
  CavityType    *C;
  WigglerType   *W;
  InsertionType *ID;

  switch (kind) {
  case marker_:
    Mrk = new MarkerType::MarkerType();
    Elem = dynamic_cast<ElemType*>(Mrk);
    Elem->kind = ElemKind(marker);
    break;
  case drift_:
    D = new DriftType::DriftType();
    Elem = dynamic_cast<ElemType*>(D);
    Elem->kind = ElemKind(drift);
    break;
  case mpole_:
    M = new MpoleType::MpoleType();
    Elem = dynamic_cast<ElemType*>(M);
    Elem->kind = ElemKind(Mpole); M->thick = thicktype(thick);
    break;
  case cavity_:
    C = new CavityType::CavityType();
    Elem = dynamic_cast<ElemType*>(C);
    Elem->kind = ElemKind(Cavity);
    break;
  case thinkick_:
    M = new MpoleType::MpoleType();
    Elem = dynamic_cast<ElemType*>(M);
   Elem->kind = ElemKind(Mpole); M->thick = thicktype(thin);
    break;
  case wiggler_:
    W = new WigglerType::WigglerType();
    Elem = dynamic_cast<ElemType*>(W);
    Elem->kind = ElemKind(Wigl);
    break;
  case insertion_:
    ID = new InsertionType::InsertionType();
    Elem = dynamic_cast<ElemType*>(ID);
    Elem->kind = ElemKind(Insertion);
    break;
  default:
    cout << "get_kind: unknown type " << kind << " " << Elem->name << endl;
    exit_(1);
    break;
  }
}


void rdmfile(const char *mfile_dat)
{
  char          line[max_str], file_name[max_str];
  int           j, k, nmpole, kind, method, n;
  long int      i;
  double        drollerr;
  MpoleType     *M;
  CavityType    *C;
  WigglerType   *W;
  InsertionType *ID;

  bool  prt = false;

  cout << endl;
  cout << "reading machine file: " << mfile_dat << endl;

  file_rd(inf, mfile_dat);

  while (inf.getline(line, max_str) != NULL) {
    if (prt) printf("%s\n", line);
    sscanf(line, "%*s %*d %*d %ld", &i);

    Cell[i].dS[X_] = 0e0; Cell[i].dS[Y_] = 0e0;
    Cell[i].droll[X_] = 1e0; Cell[i].droll[Y_] = 0e0;

    sscanf(line, "%s %d %d", Cell[i].Elem->name, &Cell[i].Fnum, &Cell[i].Knum);

    // For compability with lattice parser.
    k = 0;
    while (Cell[i].Elem->name[k] != '\0')
      k++;
    for (j = k; j < SymbolLength; j++)
      Cell[i].Elem->name[j] = ' ';

    if (Cell[i].Knum == 1) {
      strcpy(ElemFam[Cell[i].Fnum-1].Elem->name, Cell[i].Elem->name);
      globval.Elem_nFam = max(Cell[i].Fnum, globval.Elem_nFam);
    }

    if (i > 0) {
      ElemFam[Cell[i].Fnum-1].KidList[Cell[i].Knum-1] = i;
      ElemFam[Cell[i].Fnum-1].nKid =
	max(Cell[i].Knum, ElemFam[Cell[i].Fnum-1].nKid);
    }

    inf.getline(line, max_str);
    if (prt) printf("%s\n", line);
    sscanf(line, "%d %d %d", &kind, &method, &n);
    get_kind(kind, Cell[i].Elem);
    if (i > 0)
      ElemFam[Cell[i].Fnum-1].Elem->kind = Cell[i].Elem->kind;

    inf.getline(line, max_str);
    if (prt) printf("%s\n", line);
    sscanf(line, "%lf %lf %lf %lf",
	   &Cell[i].maxampl[X_][0], &Cell[i].maxampl[X_][1],
	   &Cell[i].maxampl[Y_][0], &Cell[i].maxampl[Y_][1]);

    Cell[i].Elem->L = 0e0;

    switch (Cell[i].Elem->kind) {
    case undef:
      cout << "rdmfile: unknown type " << i << endl;
      exit_(1);
      break;
    case marker:
      break;
    case drift:
      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf", &Cell[i].Elem->L);
      break;
    case Mpole:
      M = static_cast<MpoleType*>(Cell[i].Elem);
      M->method = method; M->n = n;

      if (M->thick == thick) {
	inf.getline(line, max_str);
	if (prt) printf("%s\n", line);
	sscanf(line, "%lf %lf %lf %lf",
	       &Cell[i].dS[X_], &Cell[i].dS[Y_], &M->rollpar, &drollerr);
	Cell[i].droll[X_] = cos(dtor(drollerr+M->rollpar));
	Cell[i].droll[Y_] = sin(dtor(drollerr+M->rollpar));
	M->rollrms = drollerr - M->rollpar;
	M->rollrnd = 1e0;

	inf.getline(line, max_str);
	if (prt) printf("%s\n", line);
	sscanf(line, "%lf %lf %lf %lf %lf",
	       &Cell[i].Elem->L, &M->irho, &M->tx1, &M->tx2, &M->gap);
	if (M->irho != 0e0) M->order = 1;
      } else {
	inf.getline(line, max_str);
	if (prt) printf("%s\n", line);
	sscanf(line, "%lf %lf %lf",
	       &Cell[i].dS[X_], &Cell[i].dS[Y_], &drollerr);
	Cell[i].droll[X_] = cos(dtor(drollerr));
	Cell[i].droll[Y_] = sin(dtor(drollerr));
	M->rollrms = drollerr; M->rollrnd = 1e0;
      }

      M->c0 = sin(Cell[i].Elem->L*M->irho/2.0);
      M->c1 = cos(dtor(M->rollpar))*M->c0;
      M->s1 = sin(dtor(M->rollpar))*M->c0;

      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%d %d", &nmpole, &M->n_design);
      for (j = 1; j <= nmpole; j++) {
	inf.getline(line, max_str);
	if (prt) printf("%s\n", line);
	sscanf(line, "%d", &n);
	sscanf(line, "%*d %lf %lf", &M->bn[HOMmax+n], &M->bn[HOMmax-n]);
	M->bnpar[HOMmax+n] = M->bn[HOMmax+n];
	M->bnpar[HOMmax-n] = M->bn[HOMmax-n];
	M->order = max(n, M->order);
      }
      break;
    case Cavity:
      C = static_cast<CavityType*>(Cell[i].Elem);
      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf %lf %d %lf",
	     &C->volt, &C->freq, &C->h, &globval.Energy);
      globval.Energy *= 1e-9;
      C->volt *= globval.Energy*1e9; C->freq *= c0/(2.0*M_PI);
     break;
    case Wigl:
      W = static_cast<WigglerType*>(Cell[i].Elem);
      W->method = method; W->n = n;

      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf %lf", &Cell[i].Elem->L, &W->lambda);

      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%d", &W->n_harm);

      if (Cell[i].Knum == 1)
	ElemFam[Cell[i].Fnum-1].Elem = new WigglerType::WigglerType();
      for (j = 0; j < W->n_harm; j++) {
	inf.getline(line, max_str);
	if (prt) printf("%s\n", line);
	sscanf(line, "%d %lf %lf %lf %lf %lf", &W->harm[j],
	       &W->kxV[j], &W->BoBrhoV[j], &W->kxH[j], &W->BoBrhoH[j],
	       &W->phi[j]);
	W->BoBrhoV[j] = W->BoBrhoV[j]; W->BoBrhoH[j] = W->BoBrhoH[j];
      }
      break;
    case Insertion:
      ID = static_cast<InsertionType*>(Cell[i].Elem);
      ID->method = method; ID->n = n;

      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf %d %s", &ID->scaling, &n, file_name);

      if (n == 1) {
	ID->firstorder = true;ID->secondorder = false;

	strcpy(ID->fname1, file_name);
	Read_IDfile(ID->fname1, Cell[i].Elem->L, ID->nx, ID->nz,
		    ID->tabx, ID->tabz, ID->thetax1, ID->thetaz1,
		    ID->long_comp, ID->B2);
      } else if (n == 2) {
	ID->firstorder = false;	ID->secondorder = true;

	strcpy(ID->fname2, file_name);
	Read_IDfile(ID->fname2, Cell[i].Elem->L, ID->nx, ID->nz,
		    ID->tabx, ID->tabz, ID->thetax, ID->thetaz,
		    ID->long_comp, ID->B2);
      } else {
	cout << "rdmfile: undef order " << n << endl;
	exit_(1);
      }

      if (ID->method == 1)
	ID->linear = true;
      else
	ID->linear = false;

      if (!ID->linear) {
	ID->mtx = gsl_matrix_alloc(ID->nz, ID->nx);
	GSL2NRDM2(pmtx, ID->mtx, ID->tx, 0);

	ID->mtz = gsl_matrix_alloc(ID->nz, ID->nx);
	GSL2NRDM2(pmtz, ID->mtz, ID->tz, 0);

	ID->tab1 = (double *)malloc((ID->nx)*sizeof(double));
	ID->tab2 = (double *)malloc((ID->nz)*sizeof(double));

	ID->mf2x = gsl_matrix_alloc(ID->nz, ID->nx);
	GSL2NRDM2(pmf2x, ID->mf2x, ID->f2x, 0);

	ID->mf2z = gsl_matrix_alloc(ID->nz, ID->nx);
	GSL2NRDM2(pmf2z, ID->mf2z, ID->f2z, 0);

	Matrices4Spline(ID);
      }

/*      free_matrix(tx, 1, nz, 1, nx); free_matrix(tz, 1, nz, 1, nx);
      free(tab1); free(tab2);
      free_matrix(f2x, 1, nz, 1, nx); free_matrix(f2z, 1, nz, 1, nx); */
      break;
    case FieldMap:
      break;
    default:
      cout << "rdmfile: unknown type" << endl;
      exit_(1);
      break;
    }

    if (i == 0)
      Cell[i].S = 0e0;
    else
      Cell[i].S = Cell[i-1].S + Cell[i].Elem->L;
  }

  globval.Cell_nLoc = i;

  globval.dPcommon = 1e-8; globval.CODeps = 1e-14; globval.CODimax = 40;

  SI_init();

  cout << endl;
  cout  << fixed << setprecision(5)
	<< "rdmfile: read " << globval.Cell_nLoc << " elements, C = "
	<< Cell[globval.Cell_nLoc].S << endl;

  inf.close();
}
