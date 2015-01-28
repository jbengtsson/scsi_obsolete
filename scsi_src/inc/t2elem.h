/* Tracy-3

   J. Bengtsson, BNL 2007

*/

extern bool    sympl;
extern int     FieldMap_filetype;
extern double  cl_rad, q_fluct, I2, I4, I5;

void Drift_Init(int Fnum);

void Mpole_Init(int Fnum);

void Wiggler_Init(int Fnum);

void FieldMap_Init(int Fnum);

void Cav_Init(int Fnum);

void Marker_Init(int Fnum);

void Insertion_Init(int Fnum);

void Spreader_Init(int Fnum);

void Recombiner_Init(int Fnum);

void Solenoid_Init(int Fnum);


template<typename T>
T get_p_s(const ss_vect<T> &);

void getelem(long i, CellType *cellrec);

void putelem(long i, CellType *cellrec);


int GetnKid(const int Fnum);

long Elem_GetPos(const int Fnum, const int Knum);


template<typename T>
void GtoL(ss_vect<T> &X, Vector2 &S, Vector2 &R,
	  const double c0, const double c1, const double s1);

template<typename T>
void LtoG(ss_vect<T> &X, Vector2 &S, Vector2 &R,
	  double c0, double c1, double s1);

template<typename T>
void p_rot(double phi, ss_vect<T> &x);


template<typename T>
void get_B2(const double h_ref, const T B[], const ss_vect<T> &xp,
	    T &B2_perp, T &B2_par);

template<typename T>
void radiate(ss_vect<T> &x, const double L, const double h_ref, const T B[]);

template<typename T>
void Drift(double L, ss_vect<T> &x);

template<typename T>
void bend_fringe(double hb, ss_vect<T> &x);

template<typename T>
void EdgeFocus(double irho, double phi, double gap, ss_vect<T> &x);

template<typename T>
void quad_fringe(double b2, ss_vect<T> &x);


template<typename T>
void Drift_Pass(CellType &Cell, ss_vect<T> &x);

template<typename T>
void thin_kick(int Order, double MB[], double L, double h_bend, double h_ref,
	       ss_vect<T> &x);

template<typename T>
void Mpole_Pass(CellType &Cell, ss_vect<T> &x);

template<typename T>
void Marker_Pass(CellType &Cell, ss_vect<T> &X);

template<typename T>
void Cav_Pass(CellType &Cell, ss_vect<T> &X);

template<typename T>
void Wiggler_pass_EF(ElemType *elem, ss_vect<T> &x);

template<typename T>
void Wiggler_pass_EF2(int nstep, double L, double kxV, double kxH, double kz,
		      double BoBrhoV, double BoBrhoH, double phi,
		      ss_vect<T> &x);

template<typename T>
void Wiggler_pass_EF3(ElemType *elem, ss_vect<T> &x);

template<typename T>
void Wiggler_Pass(CellType &Cell, ss_vect<T> &X);

template<typename T>
void FieldMap_Pass(CellType &Cell, ss_vect<T> &X);

template<typename T>
void Insertion_Pass(CellType &Cell, ss_vect<T> &x);

template<typename T>
void sol_pass(const ElemType *elem, ss_vect<T> &x);

template<typename T>
void Solenoid_Pass(CellType &Cell, ss_vect<T> &x);


void Mpole_Setbn(int Fnum, int Knum, int Order);

void Wiggler_Setbn(int Fnum, int Knum, int Order);


void MulLsMat(Matrix &A, Matrix &B);

void LinsTrans(Matrix &A, Vector &b);


void SI_init(void);

void get_B(const char *file_name, FieldMapType *FM);

double Elem_GetKval(int Fnum, int Knum, int Order);

void Mpole_SetdS(int Fnum, int Knum);

void Mpole_Setdroll(int Fnum, int Knum);
