
// Beam line class.

class ElemType;

// LEGO block structure for the lattice cells.
struct CellType {
 public:
  int
    Fnum,                // Element Family no.
    Knum;                // Element Kid no.
  double
    S;                   // Ring location.
  CellType
    *next_ptr;           // Pointer to next cell (for ERLs, etc.)
  ElemType
    *Elem;               // Pointer to the particular element.
  // Place holders for global properties.
  double
    dS[2],               // Transverse displacement
    droll[2],            // droll = (cos(roll), sin(roll))
    Nu[2],               // Phase advances
    Alpha[2],            // Alpha functions (redundant)
    Beta[2],             // beta fonctions (redundant)
    Eta[2],
    Etap[2];             // dispersion and its derivative (redundant)
  std::vector<double>
    maxampl[2];          // Horizontal and vertical physical apertures:
                         // maxampl[X_][0] < x < maxampl[X_][1]
                         // maxampl[Y_][0] < y < maxampl[Y_][1]
  Vector
    BeamPos;             // Last position of the beam this cell
  Matrix
    A,                   // Floquet space to phase space transformation
    sigma;               // sigma matrix (redundant)
};


// Element base class.
class ElemType {
 public:
  partsName
    name;                // Name.
  double
    L;                   // Length[m]
  ElemKind
    kind;                // Type (enumerator).

  /* virtual bool Elem_Pass(const long i, ss_vect<double> &) = 0; */
  /* virtual bool Elem_Pass(const long i, ss_vect<tps> &) = 0; */
};


// Element families.
struct ElemFamType {
 public:
  ElemType
    *Elem;               // Pointer to the particular element.
  int
    nKid,                // No of kids.
    KidList[nKidMax],
    NoDBN;
  std::vector<DBNameType>
    DBNlist[nKidMax];
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
  int        method,     // Integration Method
  int        n;          // Number of integration steps
  // Displacement Errors
  Vector2    dSsys,      // systematic [m]
  Vector2    dSrms,      // rms [m]
  Vector2    dSrnd;      // random number
  // Roll angle
  double     rollpar,    // design [deg]
  double     rollsys,    // systematic [deg]
  double     rollrms,    // rms [deg]
  double     rollrnd;    // random number
  // Multipole strengths
  mpolArray  bnpar,      // design
  mpolArray  bnsys,      // systematic
  mpolArray  bnrms,      // rms
  mpolArray  bnrnd,      // random number
  mpolArray  bn;         // total
  int        order,      // The highest order in PB
  int        n_design;   // multipole order (design)
  thicktype  thick;
  // Bending Angles
  double     tx1,        // horizontal entrance angle [deg]
  double     tx2,        // horizontal exit angle [deg]
  double     gap,        // total magnet gap [m]
  double     irho,       // 1/rho [1/m]
  double     c0, c1, s1; // corrections for roll error of bend

  MpoleType();
};


class CavityType : public ElemType {
 public:
  double volt,           // Vrf [V]
  double freq,           // Vrf [Hz]
  double phi;            // RF phase
  int    h;              // Harmonic number

  CavityType();
};


class WigglerType : public ElemType {
 public:
  int
    method,              // Integration Method
    n,                   // number of integration steps
    order,               // The highest order in PB
    n_harm;              // no of harmonics
  // Displacement Error
  double
    dSsys[2],            // systematic [m]
    dSrms[2],            // rms [m]
    dSrnd[2],            // random number
  // Roll angle
    rollpar,             // design [deg]
    rollsys,             // systematic [deg]
    rollrms,             // rms [deg]
    rollrnd,             // random number
    lambda;              // lambda
  std::vector<int>
    harm            ;    // harmonic number
  std::vector<double>
    BoBrhoV,             // B/Brho vertical
    BoBrhoH,             // B/Brho horizontal
    kxV,                 // kx
    kxH,                 // kx
    phi[n_harm_max];     // phi
  mpolArray
    bn;

  WigglerType();
};


class InsertionType : public ElemType {
 public:
  int
    method,              // Integration Method
    n,                   // number of integration steps
    order;               // The highest order in PB
  string
    fname1,              // Filename for insertion description: first ordre
    fname2;              // Filename for insertion description: second ordre
  int
    nx,                  // Horizontal point number
    nz;                  // Vertical point number
  double
    scaling;             // static scaling factor as in BETA ESRF
  bool
    linear,              // if true linear interpolation else spline
    firstorder,          // true if first order kick map loaded
    secondorder;         // true if second order kick map loaded
  std::vector<double>
    tabx,                // spacing in H-plane
    tabz;                // spacing in V-plane
  // 1 for first order
  std::vector<std::vector<double>>
    thetax,
    thetax1,
    thetaz,
    thetaz1;
  bool
    long_comp;           // flag for longitudinal comp
  std::vector<std::vector<double>>
    B2[IDZMAX][IDXMAX];  // B^2_perp

  //GSL add
  gsl_matrix*
    mtx,                 // gsl matrix, coresponds to tx
    mtz,                 // gsl matrix, coresponds to tz
    mf2x,                // gsl matrix, coresponds to f2x
    mf2z;                // gsl matrix, coresponds to f2z
  //end GSL add

  double
    **tx,
    **tz,
    **f2x,
    **f2z,
    **tx1,
    **tz1,
    **f2x1,
    **f2z1,              // a voir
    *tab1,
    *tab2;               // tab of x and z meshes from Radia code

  // Displacement Error
  Vector2
    dSsys,               // systematic [m]
    dSrms,               // rms [m]
    dSrnd;               // random number
  // Roll angle
  double
    rollpar,             // design [deg]
    rollsys,             // systematic [deg]
    rollrms,             // rms [deg]
    rollrnd;             // random number
  // Strength
  // int
  //   nperiod;             // Number of periods 
  // double
  //   lperiod,             // Length Period [m]
  //   BoBrho,              // B/Brho 
  //   Kx;                  // kx 
  //  mpolArray
  //    BW; 

  InsertionType();
};


class FieldMapType : public ElemType {
 public:
  int
    n_step,              // number of integration steps
    n[3],                // no of steps
    cut;                 // cut in z direction
  double
    scl,
    phi,
    x0,
    Lr,
    Lm,
    Ld,
    L1,
    dx[3],
    *x[3];               // [dx, dy, dz], [x, y, z]
  double
    ***BoBrho[3],
    ***BoBrho2[3],       // [B_x, B_y, B_z]
    ***AoBrho[2],
    ***AoBrho2[2];       // [Ax(x, y, z), Ay(x, y, z)], spline info
  //GSL add
  gsl_vector
    *vx[3];              // coresponds to *x[3]
  //end GSL add

  FieldMapType();
};


class SpreaderType : public ElemType {
 public:
  std::vector<double>
    E_max;               // energy levels in increasing order.
  std::vector<CellType>
    *Cell_ptrs;

  SpreaderType();
};


class RecombinerType : public ElemType {
 public:
  double
    E_min,
    E_max;

  RecombinerType();
};


class SolenoidType : public ElemType {
 public:
  int
    n;                   // Number of integration steps
  // Displacement Errors
  double dSsys[2],       // systematic [m]
    dSrms[2],            // rms [m]
    dSrnd[2],            // random number
  // Roll angle
    rollpar,             // design [deg]
    rollsys,             // systematic [deg]
    rollrms,             // rms [deg]
    rollrnd,             // random number
    BoBrho;              // normalized field strength

  SolenoidType();
};


typedef struct globvalrec {
  bool
    Cavity_on,           // if true, cavity turned on
    radiation,           // if true, radiation turned on
    emittance,
    quad_fringe,         // dipole- and quadrupole hard-edge fringe fields.
    H_exact,             // "small ring" Hamiltonian.
    pathlength,          // absolute path length
    stable,
    Aperture_on,
    EPU,
    wake_on,
    IBS;                 // intrabeam scattering
  long int
    Cell_nLoc,           // number of elements
    Elem_nFam,           // number of families
    CODimax;             // maximum number of cod search before failing
  int
    RingType,            // 1 if a ring (0 if transfer line)
    bpm,                 // bpm number
    hcorr,               // horizontal corrector number
    vcorr,               // vertical corrector number
    qt,                  // vertical corrector number
    gs,                  // girder start marker
    ge;                  // girder end marker
  double
    dPcommon,            // dp for numerical differentiation
    dPparticle,          // energy deviation
    delta_RF,            // RF acceptance
    Omega,
    U0,                  // energy lost per turn in keV
    Alphac,              // alphap
    Energy,              // ring energy
    CODeps,              // precision for closed orbit finder
    TotalTune[2],        // transverse tunes
    Chrom[2],            // chromaticities
    dE,                  // energy loss
    alpha_rad[DOF],      // damping coeffs.
    D_rad[DOF],          // diffusion coeffs (Floquet space)
    J[DOF],              // partition numbers
    tau[DOF],            // damping times
    Qb,                  // bunch charge
    D_IBS[DOF],          // diffusion matrix (Floquet space)
    eps[DOF],            // 3 motion invariants
    epsp[DOF],           // transverse and longitudinal projected emittances
    alpha_z, beta_z;     // longitudinal alpha and beta
  Vector
    CODvect,             // closed orbit
    wr,
    wi;                  // real and imaginary part of eigenvalues
  Matrix
    OneTurnMat,          // oneturn matrix
    Ascr,
    Ascrinv,
    Vr,                  // real part of the eigenvectors
    Vi;                  // imaginal par of the eigenvectors
} globvalrec;
