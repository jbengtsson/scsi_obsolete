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


void get_kind(const int kind, ElemType &Elem)
{

  switch (kind) {
  case marker_:
    Elem.kind = ElemKind(marker);
    break;
  case drift_:
    Elem.kind = ElemKind(drift);
    Drift_Alloc(&Elem);
    break;
  case mpole_:
    Elem.kind = ElemKind(Mpole);
    Mpole_Alloc(&Elem);
    Elem.M->thick = thicktype(thick);
    break;
  case cavity_:
    Elem.kind = ElemKind(Cavity);
    Cav_Alloc(&Elem);
    break;
  case thinkick_:
    Elem.kind = ElemKind(Mpole);
    Mpole_Alloc(&Elem);
    Elem.M->thick = thicktype(thin);
    break;
  case wiggler_:
    Elem.kind = ElemKind(Wigl);
    Wiggler_Alloc(&Elem);
    break;
  case insertion_:
    Elem.kind = ElemKind(Insertion);
    Insertion_Alloc(&Elem);
    break;
  default:
    cout << "get_kind: unknown type " << kind << " " << Elem.name << endl;
    exit_(1);
    break;
  }
}


void rdmfile(const char *mfile_dat)
{
  char      line[max_str], file_name[max_str];
  int       j, k, nmpole, kind, method, n;
  long int  i;
  double    drollerr;

  bool  prt = false;

  cout << endl;
  cout << "reading machine file: " << mfile_dat << endl;

  file_rd(inf, mfile_dat);

  while (inf.getline(line, max_str) != NULL) {
    if (prt) printf("%s\n", line);
    sscanf(line, "%*s %*d %*d %ld", &i);

    Cell[i].dS[X_] = 0.0; Cell[i].dS[Y_] = 0.0;
    Cell[i].droll[X_] = 1.0; Cell[i].droll[Y_] = 0.0;

    sscanf(line, "%s %d %d", Cell[i].Elem.name, &Cell[i].Fnum, &Cell[i].Knum);

    // For compability with lattice parser.
    k = 0;
    while (Cell[i].Elem.name[k] != '\0')
      k++;
    for (j = k; j < SymbolLength; j++)
      Cell[i].Elem.name[j] = ' ';

    if (Cell[i].Knum == 1) {
      strcpy(ElemFam[Cell[i].Fnum-1].ElemF.name, Cell[i].Elem.name);
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
      ElemFam[Cell[i].Fnum-1].ElemF.kind = Cell[i].Elem.kind;

    inf.getline(line, max_str);
    if (prt) printf("%s\n", line);
    sscanf(line, "%lf %lf %lf %lf",
	   &Cell[i].maxampl[X_][0], &Cell[i].maxampl[X_][1],
	   &Cell[i].maxampl[Y_][0], &Cell[i].maxampl[Y_][1]);

    Cell[i].Elem.L = 0.0;

    switch (Cell[i].Elem.kind) {
    case undef:
      cout << "rdmfile: unknown type " << i << endl;
      exit_(1);
      break;
    case marker:
      break;
    case drift:
      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf", &Cell[i].Elem.L);
      break;
    case Cavity:
      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf %lf %d %lf",
	     &Cell[i].Elem.C->volt, &Cell[i].Elem.C->freq,
	     &Cell[i].Elem.C->h, &globval.Energy);
      globval.Energy *= 1e-9;
      Cell[i].Elem.C->volt *= globval.Energy*1e9;
      Cell[i].Elem.C->freq *= c0/(2.0*M_PI);
     break;
    case Mpole:
      Cell[i].Elem.M->method = method; Cell[i].Elem.M->n = n;

      if (Cell[i].Elem.M->thick == thick) {
	inf.getline(line, max_str);
	if (prt) printf("%s\n", line);
	sscanf(line, "%lf %lf %lf %lf",
	       &Cell[i].dS[X_], &Cell[i].dS[Y_],
	       &Cell[i].Elem.M->rollpar, &drollerr);
	Cell[i].droll[X_] = cos(dtor(drollerr+Cell[i].Elem.M->rollpar));
	Cell[i].droll[Y_] = sin(dtor(drollerr+Cell[i].Elem.M->rollpar));
	Cell[i].Elem.M->rollrms = drollerr - Cell[i].Elem.M->rollpar;
	Cell[i].Elem.M->rollrnd = 1e0;

	inf.getline(line, max_str);
	if (prt) printf("%s\n", line);
	sscanf(line, "%lf %lf %lf %lf %lf",
	       &Cell[i].Elem.L, &Cell[i].Elem.M->irho,
	       &Cell[i].Elem.M->tx1, &Cell[i].Elem.M->tx2,
	       &Cell[i].Elem.M->gap);
	if (Cell[i].Elem.M->irho != 0.0) Cell[i].Elem.M->order = 1;
      } else {
	inf.getline(line, max_str);
	if (prt) printf("%s\n", line);
	sscanf(line, "%lf %lf %lf",
	       &Cell[i].dS[X_], &Cell[i].dS[Y_], &drollerr);
	Cell[i].droll[X_] = cos(dtor(drollerr));
	Cell[i].droll[Y_] = sin(dtor(drollerr));
	Cell[i].Elem.M->rollrms = drollerr; Cell[i].Elem.M->rollrnd = 1e0;
      }

      Cell[i].Elem.M->c0 = sin(Cell[i].Elem.L*Cell[i].Elem.M->irho/2.0);
      Cell[i].Elem.M->c1 = cos(dtor(Cell[i].Elem.M->rollpar))
	                   *Cell[i].Elem.M->c0;
      Cell[i].Elem.M->s1 = sin(dtor(Cell[i].Elem.M->rollpar))
	                   *Cell[i].Elem.M->c0;

      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%d %d", &nmpole, &Cell[i].Elem.M->n_design);
      for (j = 1; j <= nmpole; j++) {
	inf.getline(line, max_str);
	if (prt) printf("%s\n", line);
	sscanf(line, "%d", &n);
	sscanf(line, "%*d %lf %lf",
	       &Cell[i].Elem.M->bn[HOMmax+n], &Cell[i].Elem.M->bn[HOMmax-n]);
	Cell[i].Elem.M->bnpar[HOMmax+n] = Cell[i].Elem.M->bn[HOMmax+n];
	Cell[i].Elem.M->bnpar[HOMmax-n] = Cell[i].Elem.M->bn[HOMmax-n];
	Cell[i].Elem.M->order = max(n, Cell[i].Elem.M->order);
      }
      break;
    case Wigl:
      Cell[i].Elem.W->method = method; Cell[i].Elem.W->n = n;

      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf %lf", &Cell[i].Elem.L, &Cell[i].Elem.W->lambda);

      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%d", &Cell[i].Elem.W->n_harm);

      if (Cell[i].Knum == 1) Wiggler_Alloc(&ElemFam[Cell[i].Fnum-1].ElemF);
      for (j = 0; j < Cell[i].Elem.W->n_harm; j++) {
	inf.getline(line, max_str);
	if (prt) printf("%s\n", line);
	sscanf(line, "%d %lf %lf %lf %lf %lf", &Cell[i].Elem.W->harm[j],
	       &Cell[i].Elem.W->kxV[j], &Cell[i].Elem.W->BoBrhoV[j],
	       &Cell[i].Elem.W->kxH[j], &Cell[i].Elem.W->BoBrhoH[j],
	       &Cell[i].Elem.W->phi[j]);
	ElemFam[Cell[i].Fnum-1].ElemF.W->BoBrhoV[j]
	  = Cell[i].Elem.W->BoBrhoV[j];
	ElemFam[Cell[i].Fnum-1].ElemF.W->BoBrhoH[j]
	  = Cell[i].Elem.W->BoBrhoH[j];
      }
      break;
    case Insertion:
      Cell[i].Elem.ID->method = method; Cell[i].Elem.ID->n = n;

      inf.getline(line, max_str);
      if (prt) printf("%s\n", line);
      sscanf(line, "%lf %d %s", &Cell[i].Elem.ID->scaling, &n, file_name);

      if (n == 1) {
	Cell[i].Elem.ID->firstorder = true;
	Cell[i].Elem.ID->secondorder = false;

	strcpy(Cell[i].Elem.ID->fname1, file_name);
	Read_IDfile(Cell[i].Elem.ID->fname1, Cell[i].Elem.L,
		    Cell[i].Elem.ID->nx, Cell[i].Elem.ID->nz,
		    Cell[i].Elem.ID->tabx, Cell[i].Elem.ID->tabz,
		    Cell[i].Elem.ID->thetax1, Cell[i].Elem.ID->thetaz1,
		    Cell[i].Elem.ID->long_comp, Cell[i].Elem.ID->B2);
      } else if (n == 2) {
	Cell[i].Elem.ID->firstorder = false;
	Cell[i].Elem.ID->secondorder = true;

	strcpy(Cell[i].Elem.ID->fname2, file_name);
	Read_IDfile(Cell[i].Elem.ID->fname2, Cell[i].Elem.L,
		    Cell[i].Elem.ID->nx, Cell[i].Elem.ID->nz,
		    Cell[i].Elem.ID->tabx, Cell[i].Elem.ID->tabz,
		    Cell[i].Elem.ID->thetax, Cell[i].Elem.ID->thetaz,
		    Cell[i].Elem.ID->long_comp, Cell[i].Elem.ID->B2);
      } else {
	cout << "rdmfile: undef order " << n << endl;
	exit_(1);
      }

      if (Cell[i].Elem.ID->method == 1)
	Cell[i].Elem.ID->linear = true;
      else
	Cell[i].Elem.ID->linear = false;

      if (!Cell[i].Elem.ID->linear) {
	Cell[i].Elem.ID->mtx = gsl_matrix_alloc(Cell[i].Elem.ID->nz,
						Cell[i].Elem.ID->nx);
	GSL2NRDM2(pmtx, Cell[i].Elem.ID->mtx, Cell[i].Elem.ID->tx, 0);

	Cell[i].Elem.ID->mtz = gsl_matrix_alloc(Cell[i].Elem.ID->nz,
						Cell[i].Elem.ID->nx);
	GSL2NRDM2(pmtz, Cell[i].Elem.ID->mtz, Cell[i].Elem.ID->tz, 0);

	Cell[i].Elem.ID->tab1 = (double *)malloc((Cell[i].Elem.ID->nx)
						 *sizeof(double));
	Cell[i].Elem.ID->tab2 = (double *)malloc((Cell[i].Elem.ID->nz)
						 *sizeof(double));

	Cell[i].Elem.ID->mf2x = gsl_matrix_alloc(Cell[i].Elem.ID->nz,
						 Cell[i].Elem.ID->nx);
	GSL2NRDM2(pmf2x, Cell[i].Elem.ID->mf2x, Cell[i].Elem.ID->f2x, 0);

	Cell[i].Elem.ID->mf2z = gsl_matrix_alloc(Cell[i].Elem.ID->nz,
						 Cell[i].Elem.ID->nx);
	GSL2NRDM2(pmf2z, Cell[i].Elem.ID->mf2z, Cell[i].Elem.ID->f2z, 0);

	Matrices4Spline(Cell[i].Elem.ID);
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
      Cell[i].S = 0.0;
    else
      Cell[i].S = Cell[i-1].S + Cell[i].Elem.L;
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
