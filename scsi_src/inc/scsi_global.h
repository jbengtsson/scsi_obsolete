#define PLANES 2
#define IDXMAX 200
#define IDZMAX 100

const int Spreader_max = 10;

const int n_harm_max = 10;

class ElemType;

// LEGO block structure for the lattice cells.
struct CellType {
 public:
  int       Fnum;            // Element Family no.
  int       Knum;            // Element Kid no.
  double    S;               // Ring location.
  CellType  *next_ptr;       // Pointer to next cell (for ERLs, etc.)
  ElemType  *Elem;           // Pointer to the particular element.
  // Place holders for global properties.
  Vector2   dS,              // Transverse displacement
            droll;           // droll = (cos(roll), sin(roll))
  Vector2   Nu,              // Phase advances
            Alpha,           // Alpha functions (redundant)
            Beta,            // beta fonctions (redundant)
            Eta, Etap;       // dispersion and its derivative (redundant)
  Vector    BeamPos;         // Last position of the beam this cell
  Matrix    A,               // Floquet space to phase space transformation
            sigma;           // sigma matrix (redundant)
  Vector2   maxampl[PLANES]; // Horizontal and vertical physical apertures:
                             // maxampl[X_][0] < x < maxampl[X_][1]
                             // maxampl[Y_][0] < y < maxampl[Y_][1]
};


// Element base class.
class ElemType {
 public:
  partsName name;                             // Name.
  double    L;                                // Length[m]
  ElemKind  kind;                             // Type (enumerator).

  /* virtual bool Elem_Pass(const long i, ss_vect<double> &) = 0; */
  /* virtual bool Elem_Pass(const long i, ss_vect<tps> &) = 0; */
};


// Element families.
struct ElemFamType {
 public:
  ElemType   *Elem;            // Pointer to the particular element.
  int        nKid;             // No of kids.
  int        KidList[nKidMax];
  int        NoDBN;
  DBNameType DBNlist[nKidMax];
};


// Element derived classes.

class MarkerType : public ElemType {
 public:
  MarkerType();
};


class DriftType : public ElemType {
 public:
  DriftType();

  /* virtual bool Elem_Pass(const long i, ss_vect<double> &) { Drift(L, x); }; */
  /* virtual bool Elem_Pass(const long i, ss_vect<tps> &) { Drift(L, x); */
};


class MpoleType : public ElemType {
 public:
  int        method;     // Integration Method
  int        n;          // Number of integration steps
  // Displacement Errors
  Vector2    dSsys;      // systematic [m]
  Vector2    dSrms;      // rms [m]
  Vector2    dSrnd;      // random number
  // Roll angle
  double     rollpar;      // design [deg]
  double     rollsys;      // systematic [deg]
  double     rollrms;      // rms [deg]
  double     rollrnd;      // random number
  // Multipole strengths
  mpolArray  bnpar;      // design
  mpolArray  bnsys;      // systematic
  mpolArray  bnrms;      // rms
  mpolArray  bnrnd;      // random number
  mpolArray  bn;         // total
  int        order;      // The highest order in PB
  int        n_design;   // multipole order (design)
  thicktype  thick;
  // Bending Angles
  double     tx1;        // horizontal entrance angle [deg]
  double     tx2;        // horizontal exit angle [deg]
  double     gap;        // total magnet gap [m]
  double     irho;       // 1/rho [1/m]
  double     c0, c1, s1; // corrections for roll error of bend

  MpoleType();
};


class CavityType : public ElemType {
 public:
  double volt; // Vrf [V]
  double freq; // Vrf [Hz]
  double phi;   // RF phase
  int    h;    // Harmonic number

  CavityType();
};


class WigglerType : public ElemType {
 public:
  int       method;             // Integration Method
  int       n;                  // number of integration steps
  // Displacement Error
  Vector2   dSsys;              // systematic [m]
  Vector2   dSrms;              // rms [m]
  Vector2   dSrnd;              // random number
  // Roll angle
  double    rollpar;              // design [deg]
  double    rollsys;              // systematic [deg]
  double    rollrms;              // rms [deg]
  double    rollrnd;              // random number
  double    lambda;              // lambda
  int       n_harm;              // no of harmonics
  int       harm[n_harm_max];    // harmonic number
  double    BoBrhoV[n_harm_max]; // B/Brho vertical
  double    BoBrhoH[n_harm_max]; // B/Brho horizontal
  double    kxV[n_harm_max];     // kx
  double    kxH[n_harm_max];     // kx
  double    phi[n_harm_max];     // phi
  mpolArray bn;
  int       order;               // The highest order in PB

  WigglerType();
};


class InsertionType : public ElemType {
 public:
  int         method;       // Integration Method
  int         n;            // number of integration steps
  char        fname1[100];  // Filename for insertion description: first ordre
  char        fname2[100];  // Filename for insertion description: second ordre
  int         nx;           // Horizontal point number
  int         nz;           // Vertical point number
  double      scaling;      // static scaling factor as in BETA ESRF
  bool        linear;       // if true linear interpolation else spline
  bool        firstorder;   // true if first order kick map loaded
  bool        secondorder;  // true if second order kick map loaded
  double      tabx[IDXMAX]; // spacing in H-plane
  double      tabz[IDZMAX]; // spacing in V-plane
  // 1 for first order
  double      thetax[IDZMAX][IDXMAX], thetax1[IDZMAX][IDXMAX];
  double      thetaz[IDZMAX][IDXMAX], thetaz1[IDZMAX][IDXMAX];
  bool        long_comp;    // flag for longitudinal comp
  double      B2[IDZMAX][IDXMAX]; // B^2_perp

  //GSL add
  gsl_matrix* mtx;          // gsl matrix, coresponds to tx
  gsl_matrix* mtz;          // gsl matrix, coresponds to tz
  gsl_matrix* mf2x;         // gsl matrix, coresponds to f2x
  gsl_matrix* mf2z;         // gsl matrix, coresponds to f2z
  //end GSL add

  double      **tx, **tz, **f2x, **f2z;
  double      **tx1, **tz1, **f2x1, **f2z1; // a voir
  double      *tab1, *tab2; // tab of x and z meshes from Radia code

  // Displacement Error
  Vector2     dSsys;        // systematic [m]
  Vector2     dSrms;        // rms [m]
  Vector2     dSrnd;        // random number
  // Roll angle
  double      rollpar;      // design [deg]
  double      rollsys;      // systematic [deg]
  double      rollrms;      // rms [deg]
  double      rollrnd;      // random number
  // Strength
  /* double      lperiod;      // Length Period [m] */
  /* int         nperiod;      // Number of periods */
  /* double      BoBrho;       // B/Brho */
  /* double      Kx;           // kx */
  /* mpolArray   BW; */
  int         order;        // The highest order in PB

  InsertionType();
};


class FieldMapType : public ElemType {
 public:
  int     n_step;                       // number of integration steps
  int     n[3];                         // no of steps
  int     cut;                          // cut in z direction
  double  scl, phi, x0, Lr, Lm, Ld, L1;
  double  dx[3], *x[3];                 // [dx, dy, dz], [x, y, z]
  double  ***BoBrho[3], ***BoBrho2[3];  // [B_x, B_y, B_z]
  double  ***AoBrho[2], ***AoBrho2[2];  // [Ax(x, y, z), Ay(x, y, z)],
					// spline info
  //GSL add
  gsl_vector	*vx[3]; //coresponds to *x[3]
  //end GSL add

  FieldMapType();
};


class SpreaderType : public ElemType {
 public:
  double   E_max[Spreader_max];      // energy levels in increasing order
  CellType *Cell_ptrs[Spreader_max];

  SpreaderType();
};


class RecombinerType : public ElemType {
 public:
  double E_min;
  double E_max;

  RecombinerType();
};


class SolenoidType : public ElemType {
 public:
  int     n;        // Number of integration steps
  // Displacement Errors
  Vector2 dSsys;    // systematic [m]
  Vector2 dSrms;    // rms [m]
  Vector2 dSrnd;    // random number
  // Roll angle
  double  rollpar;  // design [deg]
  double  rollsys;  // systematic [deg]
  double  rollrms;  // rms [deg]
  double  rollrnd;  // random number
  double  BoBrho;   // normalized field strength

  SolenoidType();
};


typedef struct globvalrec {
  double  dPcommon,        // dp for numerical differentiation
          dPparticle;      // energy deviation
  double  delta_RF;        // RF acceptance
  Vector2 TotalTune;       // transverse tunes
  double  Omega,
          U0,              // energy lost per turn in keV
          Alphac;          // alphap
  Vector2 Chrom;           // chromaticities
  double  Energy;          // ring energy
  long    Cell_nLoc,       // number of elements
          Elem_nFam,       // number of families
          CODimax;         // maximum number of cod search before failing
  double  CODeps;          // precision for closed orbit finder
  Vector  CODvect;         // closed orbit
  int     bpm;             // bpm number
  int     hcorr;           // horizontal corrector number
  int     vcorr;           // vertical corrector number
  int     qt;              // vertical corrector number
  int     gs;              // girder start marker
  int     ge;              // girder end marker
  Matrix  OneTurnMat,      // oneturn matrix
          Ascr,
          Ascrinv,
          Vr,              // real part of the eigenvectors
          Vi;              // imaginal par of the eigenvectors

  bool    Cavity_on,       // if true, cavity turned on
          radiation,       // if true, radiation turned on
          emittance,
          quad_fringe,     // dipole- and quadrupole hard-edge fringe fields.
          H_exact,         // "small ring" Hamiltonian.
          pathlength,      // absolute path length
          stable,
          Aperture_on,
          EPU,
          wake_on;

  double  dE,              // energy loss
          alpha_rad[DOF],  // damping coeffs.
          D_rad[DOF],      // diffusion coeffs (Floquet space)
          J[DOF],          // partition numbers
          tau[DOF];        // damping times
  bool    IBS;             // intrabeam scattering
  double  Qb,              // bunch charge
          D_IBS[DOF];      // diffusion matrix (Floquet space)
  Vector  wr, wi;          // real and imaginary part of eigenvalues
  double  eps[DOF],        // 3 motion invariants
       	  epsp[DOF],       // transverse and longitudinal projected emittances
          alpha_z, beta_z; // longitudinal alpha and beta
  int     RingType;        // 1 if a ring (0 if transfer line)
} globvalrec;
