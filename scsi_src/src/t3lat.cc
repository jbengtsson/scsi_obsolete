/*
 * Copyright (C) 2008 Lingyun Yang.
 *
 * For more information please contact lingyun.yang@gmail.com
 *
 */

/**
 * This is an example of how to use the general lattice parser, i.e. glps
 */

#include "../glps/glps.h"
#include "../glps/glpsutil.h"
#include "../glps/glpserror.h"
extern int yy_flex_debug;
extern int yydebug;

// using namespace std;
extern GLPS_SYMB* symbtab;

jmp_buf env0;

void GetRingType();/* define whether a ring or a transfer line */
void GetEnergy();  /* define particle energy */
void GetCODEPS();  /* define COD precision */
void GetDP();      /* define energy offset */
bool CheckWiggler(long i);
void ClearHOM(double *B, bool *BA);
void AssignHOM(long elem,  double *B, bool * BA );
long CheckElementtable(const string name);
// void AssignHarm(long elem, int *harm, double * kxV, double * BoBrhoV, double * kxH, double * BoBrhoH, double * phi);

void RegisterKids(); /* Check wether too many elements */

long ElemIndex(const string &name)
{
  long    i, j;
  string  name1;

  name1 = name;
  //  j = (signed)name.length();
  //  for (i = 0; i < j; i++)
  //    name1[i] = tolower(name1[i]);
  //  for (i = j; i < SymbolLength; i++)
  //    name1 += '\0';

  trace=0;

  if (trace) {
    cout << endl;
    cout << "ElemIndex: " << name << " (";
    for (i = 0; i < (signed)name1.length(); i++)
      cout << setw(4) << (int)name1[i];
    cout << setw(4) << (int)name1[name1.length()] << " )" << endl;
    cout << endl;
  }

  if (globval.Elem_nFam > Elem_nFamMax) {
    printf("ElemIndex: Elem_nFamMax exceeded: %ld(%d)\n",
           globval.Elem_nFam, Elem_nFamMax);
    exit_(1);
  }

  i = 1;
  while (i <= globval.Elem_nFam) {
    if (trace) {
      cout << setw(2) << (name1 == ElemFam[i-1].ElemF.PName)
           << " " << name1 << " " << ElemFam[i-1].ElemF.PName << " (";
      for (j = 0; j < SymbolLength; j++)
        cout << setw(4) << (int)ElemFam[i-1].ElemF.PName[j];
      cout  << " )" << endl;
    }

    if (name1 == ElemFam[i-1].ElemF.PName) break;

    i++;
  }

  if (name1 != ElemFam[i-1].ElemF.PName) {
    cout << "ElemIndex: undefined element " << name << endl;
    exit_(1);
  }

  return i;
}


bool Lattice_Read(const char *fi_)
{
  printf(" Lattice Read (t3lat):");
  printf(" NO %d\n",NO); 
  if (setjmp(env0)) return(false);
    
  globval.Cell_nLoc = 0; globval.Elem_nFam = 0; 

  globval.CODeps   = 0.0; globval.dPcommon = 0.0; globval.Energy   = 0.0;


  // get the version
  int vmaj, vmin, vpatch, vrevision;
  glps_version(&vmaj, &vmin, &vpatch, &vrevision);

  // parse the input lattice from command line.
  int ret = parse_lattice(fi_);
  if ( ret ) {
    cerr << "ERROR: " << ret << endl;
    free_lattice();
    longjmp(env0, 1);
  }

  cout << "# glps version " << vmaj << '.' << vmin << '.' << vpatch << endl;

  GetEnergy(); /* define particle energy */
  /* Energy must be defined before ID kickmap is parsed as energy is needed
     to scale B2_perp for ID radiation calculation */

  string use;

  std::vector<string> elemsdefined;
  elemsdefined.push_back("drift");
  elemsdefined.push_back("bending");
  elemsdefined.push_back("quadrupole");
  elemsdefined.push_back("sextupole");
  elemsdefined.push_back("cavity");
  elemsdefined.push_back("corrector");
  elemsdefined.push_back("bpm");
  elemsdefined.push_back("marker");
  elemsdefined.push_back("multipole");
  elemsdefined.push_back("wiggler");
  elemsdefined.push_back("fieldmap");
  elemsdefined.push_back("insertion");
  elemsdefined.push_back("spreader");
  elemsdefined.push_back("recombiner");
  elemsdefined.push_back("solenoid");

  std::list<string> elems = glps_elements();
  string str;
  int pos;
  for (std::list<string>::iterator i = elems.begin(); i != elems.end(); ++i) {

    glps_read(*i, str);
    for (int j=0; j<(int)str.length();++j)
      str[j]=tolower(str[j]); //str.replace(i,1,1,tolower(str[i]));
    for (pos=0; pos<(int)elemsdefined.size();++pos)
      if (!str.compare(elemsdefined[pos])) break;
    if (pos == (int)elemsdefined.size()) {
      string unknown=*i; //strdup(*i);
      for(int k=0; k<(int)unknown.length();++k)
	unknown[k]=tolower(unknown[k]);
      if (!unknown.compare("use")){ 
	use = str;
	continue;
      }
      else {
	cerr << "ERROR: undefined element - " << str << endl;
	longjmp(env0, 1);
      }
    }

    mpolArray B;
    bool BA[2*HOMmax+1];
 
    double val;
    std::vector<double> vec;

    double         t, t1, t2, gap, QL = 0.0, QK;
    double         QKV, QKH, QKxV, QKxH, QPhi, QKS;
    double         dt, Frf, Vrf;
    long           k1, k2, harnum;
    //      Lat_symbol     sym1;
    //      symset         mysys, SET;
    //      long           SET1[(long)lsym / 32 + 2];
    ElemFamType    *WITH;
    elemtype       *WITH1;
    MpoleType      *WITH2;
    //      symset         SET2;
    CavityType     *WITH3;
    //      long           SET3[(long)possym / 32 + 2];
    //      long           SET4[(long)dbnsym / 32 + 2];
    WigglerType    *WITH4;
    FieldMapType   *WITH6;
    /* ID Laurent */
    InsertionType  *WITH5;
    SolenoidType   *WITH7;
    //      char str1[100] = "";
    //      char str2[100] = "";
    bool firstflag  = true; // flag for first kick input
    bool secondflag = true; // flag for second kick input
    //      long           i;
    int            kx, kz;
    double         scaling;
                
    int harm[n_harm_max];
    double kxV[n_harm_max],kxH[n_harm_max];
    double BoBrhoV[n_harm_max],BoBrhoH[n_harm_max];
    double phi[n_harm_max];
    int n_harm=0;
    string str1="";
    string str2="";

    switch(pos)
      {
      case 0: //drift

	glps_read(*i,"l",val);

	globval.Elem_nFam++;
	if (globval.Elem_nFam <= Elem_nFamMax) {
	  WITH = &ElemFam[globval.Elem_nFam-1];
	  WITH1 = &WITH->ElemF;
	  memset(WITH1->PName,0,sizeof(partsName));
	  memcpy(WITH1->PName, (*i).c_str(), (*i).length());
	  WITH1->PL = val;
	  WITH1->Pkind = PartsKind(drift);
	  Drift_Alloc(&WITH->ElemF);
	} else {
	  cerr << "Elem_nFamMax exceeded: " << globval.Elem_nFam
	       << "(" << (long)Elem_nFamMax << ")" << endl;
	  longjmp(env0, 1);
	}
    
	break;

      case 1: //bending
    
	QL = 0.0;   /* L */
	QK = 0.0;   /* K */
	k1 = 0;     /* N */
	t  = 0.0;   /* T */
	t1 = 0.0;   /* T1 */
	t2 = 0.0;   /* T2 */
	gap = 0.0;  /* gap */
	k2 = Meth_Linear;   /* method */
	dt = 0.0;
	ClearHOM(B, BA);

	if (GLPS_SUCCESS==glps_read(*i,"l",val)) QL = val;
	if (GLPS_SUCCESS==glps_read(*i,"k",val)) QK = val;
	if (GLPS_SUCCESS==glps_read(*i,"n",val)) k1 = (long) floor(val+0.5);
	if (GLPS_SUCCESS==glps_read(*i,"t",val)) t = val;
	if (GLPS_SUCCESS==glps_read(*i,"roll",val)) dt = val;
	if (GLPS_SUCCESS==glps_read(*i,"t1",val)) t1 = val;
	if (GLPS_SUCCESS==glps_read(*i,"t2",val)) t2 = val;
	if (GLPS_SUCCESS==glps_read(*i,"gap",val)) gap = val;
	if (GLPS_SUCCESS==glps_read(*i,"method",val))
	  k2 = (long) floor(val+0.5);
	if (k2 != Meth_Linear) globval.MatMeth = false;
	if ((unsigned int)k2 >= 32 ||
	    ((1 << k2) &
	     ((1 << Meth_Linear) | (1 << Meth_Second) |
	      (1 << Meth_Fourth))) == 0)
	  longjmp(env0,1);
	//              +   +
	if (GLPS_SUCCESS==glps_read(*i,"hom",vec)) 
	  for (int k=0; k<(int)vec.size()/3;k++) {
	    long order = (long) floor(vec[3*k] + 0.5);
	    if (order < 1 || order > HOMmax) {
	      cerr << "invalide value detected" << endl;
	      longjmp(env0, 1);
	    }
	    BA[order+HOMmax] = true; BA[HOMmax-order] = true;
	    B[order+HOMmax] = vec[3*k+1]; B[HOMmax-order] = vec[3*k+2];
	  }

	//    case dbnsym:
	//      GetDBN_(&V);
	//      break;
      
	//              +   +
	globval.Elem_nFam++;
	if (globval.Elem_nFam <= Elem_nFamMax) {
	  WITH = &ElemFam[globval.Elem_nFam-1];
	  WITH1 = &WITH->ElemF;
	  memset(WITH1->PName,0,sizeof(partsName));
	  memcpy(WITH1->PName, (*i).c_str(), (*i).length());
	  WITH1->PL = QL;
	  WITH1->Pkind = Mpole;
	  Mpole_Alloc(&WITH->ElemF);
	  WITH2 = WITH1->M;
	  WITH2->Pmethod = k2;
	  WITH2->PN = k1;
	  if (WITH1->PL != 0.0)
	    WITH2->Pirho = t * M_PI / 180.0 / WITH1->PL;
	  else
	    WITH2->Pirho = t * M_PI / 180.0;
	  WITH2->PTx1 = t1; WITH2->PTx2 = t2; WITH2->Pgap = gap;
	  WITH2->n_design = Dip;
	  AssignHOM(globval.Elem_nFam, B, BA);
	  //      SetDBN(&V);
	  //          +   +   +   +
	  WITH2->PBpar[HOMmax+2] = QK; WITH2->PdTpar = dt;
	} else {
	  cerr << "Elem_nFamMax exceeded: " << globval.Elem_nFam
	       << "(" << (long)Elem_nFamMax << ")" << endl;
	  longjmp(env0, 1);
	}
	break;

      case 2: //quadrupole

	QL = 0.0;   /* L */
	QK = 0.0;   /* K */
	k1 = 0;   /* N */
	k2 = Meth_Fourth;   /* method */
	dt = 0.0;
	ClearHOM(B, BA);

	if (GLPS_SUCCESS==glps_read(*i,"l",val)) QL = val;
	if (GLPS_SUCCESS==glps_read(*i,"k",val)) QK = val;
	if (GLPS_SUCCESS==glps_read(*i,"n",val)) k1 = (long) floor(val+0.5);
	if (GLPS_SUCCESS==glps_read(*i,"t",val)) t = val;
	if (GLPS_SUCCESS==glps_read(*i,"roll",val)) dt = val;
	if (GLPS_SUCCESS==glps_read(*i,"method",val))
	  k2 = (long) floor(val+0.5);
	if (k2 != Meth_Linear) globval.MatMeth = false;
	if ((unsigned int)k2 >= 32 ||
	    ((1 << k2) &
	     ((1 << Meth_Linear) | (1 << Meth_First) |
	      (1 << Meth_Second) | (1 << Meth_Fourth))) == 0)
	  longjmp(env0, 1);

	if (GLPS_SUCCESS==glps_read(*i,"hom",vec))
	  for (int k=0; k<(int)vec.size()/3;k++)
	    {
	      long order = (long) floor(vec[3*k] + 0.5);
	      if (order < 1 || order > HOMmax) {
		cerr << "invalide value detected" << endl;
		longjmp(env0, 1);
	      }
	      BA[order+HOMmax] = true; BA[HOMmax-order] = true;
	      B[order+HOMmax] = vec[3*k+1]; B[HOMmax-order] = vec[3*k+2];
	    }


	//    case dbnsym:
	//      GetDBN_(&V);
      

	globval.Elem_nFam++;
	if (globval.Elem_nFam <= Elem_nFamMax) {
	  WITH = &ElemFam[globval.Elem_nFam-1];
	  WITH1 = &WITH->ElemF;
	  memset(WITH1->PName,0,sizeof(partsName));
	  memcpy(WITH1->PName, (*i).c_str(), (*i).length());
	  WITH1->PL = QL; WITH1->Pkind = Mpole;
	  Mpole_Alloc(&WITH->ElemF);
	  WITH2 = WITH1->M;
	  WITH2->Pmethod = k2; WITH2->PN = k1; WITH2->PdTpar = dt;
	  AssignHOM(globval.Elem_nFam, B,BA);
	  //    SetDBN(&V);
	  WITH2->n_design = Quad; WITH2->PBpar[HOMmax+2] = QK;
	} else {
	  cerr << "Elem_nFamMax exceeded: " << globval.Elem_nFam
	       << "(" << (long)Elem_nFamMax << ")" << endl;
	  longjmp(env0, 1);
	}
	break;

      case 3: //sextupole
	QL = 0.0;   /* L */
	QK = 0.0;   /* K */
	k1 = 0;   /* N */
	k2 = Meth_Fourth;   /* method */
	dt = 0.0;
	ClearHOM(B,BA);

	if (GLPS_SUCCESS==glps_read(*i,"l",val)) QL = val;
	if (GLPS_SUCCESS==glps_read(*i,"k",val)) QK = val;
	if (GLPS_SUCCESS==glps_read(*i,"n",val)) k1 = (long) floor(val+0.5);
	if (GLPS_SUCCESS==glps_read(*i,"roll",val)) dt = val;
	if (GLPS_SUCCESS==glps_read(*i,"method",val)) k2 = (long) floor(val+0.5);
	if (k2 != Meth_Linear) globval.MatMeth = false;
	if ((unsigned int)k2 >= 32 ||
	    ((1 << k2) &
	     ((1 << Meth_Linear) | (1 << Meth_Second) |
	      (1 << Meth_Fourth))) == 0)
	  longjmp(env0, 1);

	if (GLPS_SUCCESS==glps_read(*i,"hom",vec))
	  for (int k=0; k<(int)vec.size()/3;k++)
	    {
	      long order = (long) floor(vec[3*k] + 0.5);
	      if (order < 1 || order > HOMmax) {
		cerr << "invalide value detected" << endl;
		longjmp(env0, 1);
	      }
	      BA[order+HOMmax] = true; BA[HOMmax-order] = true;
	      B[order+HOMmax] = vec[3*k+1]; B[HOMmax-order] = vec[3*k+2];
	    }

	//        case dbnsym:
	//          GetDBN_(&V);
	//          break;

	globval.Elem_nFam++;
	if (globval.Elem_nFam <= Elem_nFamMax) {
	  WITH = &ElemFam[globval.Elem_nFam-1];
	  WITH1 = &WITH->ElemF;
	  memset(WITH1->PName,0,sizeof(partsName));
	  memcpy(WITH1->PName, (*i).c_str(), (*i).length());
	  //    memcpy(WITH1->PName, ElementName, sizeof(partsName));
	  WITH1->PL = QL;
	  WITH1->Pkind = Mpole;
	  Mpole_Alloc(&WITH->ElemF);
	  WITH2 = WITH1->M;
	  WITH2->Pmethod = k2;
	  WITH2->PN = k1;
	  if (WITH1->PL != 0.0)
	    WITH2->Pthick = pthicktype(thick);
	  else
	    WITH2->Pthick = pthicktype(thin);
	  WITH2->PdTpar = dt; WITH2->n_design = Sext;
	  AssignHOM(globval.Elem_nFam, B, BA);
	  //    SetDBN(&V);
	  WITH2->PBpar[HOMmax + 3] = QK;
	} else {
	  cerr << "Elem_nFamMax exceeded: " << globval.Elem_nFam
	       << "(" << (long)Elem_nFamMax << ")" << endl;
	  longjmp(env0, 1);
	}
	break;

      case 4: //Cavity
	//              ClearHOM(B,BA);
	Frf = 0.0;   /* Frf */
	Vrf = 0.0;   /* Vrf */
	QPhi = 0.0;
	harnum = 0;   /* Voff */

	if (GLPS_SUCCESS==glps_read(*i,"frequency",val)) Frf = val;
	if (GLPS_SUCCESS==glps_read(*i,"voltage",val)) Vrf = val;
	if (GLPS_SUCCESS==glps_read(*i,"phi",val)) QPhi = val;
	if (GLPS_SUCCESS==glps_read(*i,"harnum",val))
	  harnum = (long) floor(val+0.5);

	//    case dbnsym:
	//      GetDBN_(&V);
	//      break;
      
	globval.Elem_nFam++;
	if (globval.Elem_nFam <= Elem_nFamMax) {
	  WITH = &ElemFam[globval.Elem_nFam-1];
	  WITH1 = &WITH->ElemF;
	  memset(WITH1->PName,0,sizeof(partsName));
	  memcpy(WITH1->PName, (*i).c_str(), (*i).length());
	  //    memcpy(WITH1->PName, ElementName, sizeof(partsName));
	  WITH1->Pkind = Cavity;
	  Cav_Alloc(&WITH->ElemF);
	  WITH3 = WITH1->C;
	  WITH3->Pvolt = Vrf;   /* Voltage [V] */
	  WITH3->Pfreq = Frf;   /* Frequency in Hz */
	  WITH3->phi = QPhi*M_PI/180.0;
	  WITH3->Ph = harnum;
	  //    SetDBN(&V);
	} else {
	  cerr << "Elem_nFamMax exceeded: " << globval.Elem_nFam
	       << "(" << (long)Elem_nFamMax << ")" << endl;
	  longjmp(env0, 1);
	}
	break;


      case  5:  //Corrector

	QL = 0.0;   /* L */
	QK = 1.0;   /* K */
	k1 = 0;     /* N */
	k2 = Meth_Linear;   /* method */
	dt = 0.0;
	//              ClearHOMandDBN(&V);
	if (GLPS_SUCCESS==glps_read(*i,"l",val)) QL = val;
	if (GLPS_SUCCESS==glps_read(*i,"n",val)) k1 = (long) floor(val+0.5);
	if (GLPS_SUCCESS==glps_read(*i,"method",val))
	  k2 = (long) floor(val+0.5);
	if (GLPS_SUCCESS==glps_read(*i,"plane",val)) QK = val;
	if (GLPS_SUCCESS==glps_read(*i,"roll",val)) dt = val;
	if (k2 != Meth_Linear) globval.MatMeth = false; 
	if ((unsigned int)k2 >= 32 ||
	    ((1 << k2) &
	     ((1 << Meth_Linear) | (1 << Meth_Second) |
	      (1 << Meth_Fourth))) == 0)
	  longjmp(env0, 1);


	//            case dbnsym:
	//              GetDBN_(&V);
	//              break;

	globval.Elem_nFam++;
	if (globval.Elem_nFam <= Elem_nFamMax) {
	  WITH = &ElemFam[globval.Elem_nFam-1];
	  WITH1 = &WITH->ElemF;
	  memset(WITH1->PName,0,sizeof(partsName));
	  memcpy(WITH1->PName, (*i).c_str(), (*i).length());
	  //    memcpy(WITH1->PName, ElementName, sizeof(partsName));
	  WITH1->PL = QL;
	  WITH1->Pkind = Mpole;
	  Mpole_Alloc(&WITH->ElemF);
	  WITH2 = WITH1->M;
	  //    SetDBN(&V);
	  if (WITH1->PL != 0.0)
	    WITH2->Pthick = pthicktype(thick);
	  else
	    WITH2->Pthick = pthicktype(thin);
	  WITH2->Pmethod = k2;
	  WITH2->PN = k1;
	  WITH2->PdTpar = dt;
	} else {
	  cerr << "Elem_nFamMax exceeded: " << globval.Elem_nFam
	       << "(" << (long)Elem_nFamMax << ")" << endl;
	  longjmp(env0, 1);
	}
	break;

      case   6: //BPM


	//  ClearHOMandDBN(&V);
	//    if (sym1 == dbnsym)
	//      GetDBN_(&V);
    
	globval.Elem_nFam++;
	if (globval.Elem_nFam <= Elem_nFamMax) {
	  WITH = &ElemFam[globval.Elem_nFam-1];
	  WITH1 = &WITH->ElemF;
	  memset(WITH1->PName,0,sizeof(partsName));
	  memcpy(WITH1->PName, (*i).c_str(), (*i).length());
	  //    memcpy(WITH1->PName, ElementName, sizeof(partsName));
	  WITH1->Pkind = Mpole;
	  Mpole_Alloc(&WITH->ElemF);
	  WITH2 = WITH1->M;
	  WITH2->Pthick = pthicktype(thin);
	  //    SetDBN(&V);
	} else {
	  cerr << "Elem_nFamMax exceeded: " << globval.Elem_nFam
	       << "(" << (long)Elem_nFamMax << ")" << endl;
	  longjmp(env0, 1);
	}
	break;

      case    7: //Marker

	//  ClearHOMandDBN(&V);
	//      GetDBN_(&V);

	globval.Elem_nFam++;
	if (globval.Elem_nFam <= Elem_nFamMax) {
	  WITH = &ElemFam[globval.Elem_nFam-1];
	  WITH1 = &WITH->ElemF;
	  memset(WITH1->PName,0,sizeof(partsName));
	  memcpy(WITH1->PName, (*i).c_str(), (*i).length());
	  //    memcpy(WITH1->PName, ElementName, sizeof(partsName));
	  WITH1->PL = 0.0;
	  WITH1->Pkind = PartsKind(marker);
	  //    SetDBN(&V);
	} else {
	  cerr << "Elem_nFamMax exceeded: " << globval.Elem_nFam
	       << "(" << (long)Elem_nFamMax << ")" << endl;
	  longjmp(env0, 1);
	}
	break;

      case    8 : // Multipole

	QL = 0.0;   /* L */
	QK = 0.0;   /* K */
	k1 = 0;   /* N */
	t = 0.0;   /* T */
	t1 = 0.0;   /* T1 */
	t2 = 0.0;   /* T2 */
	gap = 0.0;   /* gap */
	k2 = Meth_Linear;   /* method */
	dt = 0.0;
	ClearHOM(B, BA);

	if (GLPS_SUCCESS==glps_read(*i,"l",val)) QL = val;
	if (GLPS_SUCCESS==glps_read(*i,"n",val)) k1 = (long) floor(val+0.5);
	if (GLPS_SUCCESS==glps_read(*i,"t",val)) t = val;
	if (GLPS_SUCCESS==glps_read(*i,"t1",val)) t1 = val;
	if (GLPS_SUCCESS==glps_read(*i,"t2",val)) t2 = val;
	if (GLPS_SUCCESS==glps_read(*i,"gap",val)) gap = val;
	if (GLPS_SUCCESS==glps_read(*i,"roll",val)) dt = val;
	if (GLPS_SUCCESS==glps_read(*i,"method",val))
	  k2 = (long) floor(val+0.5);
	if (k2 != Meth_Linear) globval.MatMeth = false;
	if ((unsigned int)k2 >= 32 ||
	    ((1 << k2) &
	     ((1 << Meth_Linear) | (1 << Meth_Second) |
	      (1 << Meth_Fourth))) == 0)
	  longjmp(env0, 1);

	if (GLPS_SUCCESS==glps_read(*i,"hom",vec))
	  for (int k=0; k<(int)vec.size()/3;k++)
	    {
	      long order = (long) floor(vec[3*k] + 0.5);
	      if (order < 1 || order > HOMmax) {
		cerr << "invalide value detected" << endl;
		longjmp(env0, 1);
	      }
	      BA[order+HOMmax] = true; BA[HOMmax-order] = true;
	      B[order+HOMmax] = vec[3*k+1]; B[HOMmax-order] = vec[3*k+2];
	    }


	//    case dbnsym:
	//      GetDBN_(&V);
	//      break;
      
	globval.Elem_nFam++;
	if (globval.Elem_nFam <= Elem_nFamMax) {
	  WITH = &ElemFam[globval.Elem_nFam-1];
	  WITH1 = &WITH->ElemF;
	  Mpole_Alloc(&WITH->ElemF);
	  WITH2 = WITH1->M;
	  memset(WITH1->PName,0,sizeof(partsName));
	  memcpy(WITH1->PName, (*i).c_str(), (*i).length());
	  //    memcpy(WITH1->PName, ElementName, sizeof(partsName));
	  WITH1->Pkind = Mpole;
	  WITH1->PL = QL;
	  if (WITH1->PL != 0e0) {
	    WITH2->Pthick = pthicktype(thick);
	    WITH2->Pirho = t * M_PI / 180.0 / WITH1->PL;
	  } else {
	    WITH2->Pthick = pthicktype(thin);
	    WITH2->Pirho = t * M_PI / 180.0;
	  }
	  WITH2->PN = k1; WITH2->Pmethod = k2;
	  WITH2->PTx1 = t1; WITH2->PTx2 = t2; WITH2->Pgap = gap;
	  WITH2->PdTpar = dt;
	  AssignHOM(globval.Elem_nFam, B,BA);
	  WITH2->n_design = WITH2->Porder;
	  //    SetDBN(&V);
	} else {
	  cerr << "Elem_nFamMax exceeded: " << globval.Elem_nFam
	       << "(" << (long)Elem_nFamMax << ")" << endl;
	  longjmp(env0, 1);
	}
	break;
 
      case  9 :// Wiggler

	QL = 0e0; QK = 0e0; QKV = 0e0; QKH = 0e0; QKxV = 0e0; QKxH = 0e0;
	QPhi = 0e0; QKS = 0e0; k1 = 0; k2 = Meth_Linear; dt = 0e0;
	//                    ClearHOMandDBN(&V);
                
	if (GLPS_SUCCESS==glps_read(*i,"l",val)) QL = val;
	if (GLPS_SUCCESS==glps_read(*i,"bobrhov",val)) QKV = val;
	if (GLPS_SUCCESS==glps_read(*i,"bobrhoh",val)) QKH = val;
	if (GLPS_SUCCESS==glps_read(*i,"kxv",val)) QKxV = val;
	if (GLPS_SUCCESS==glps_read(*i,"kxh",val)) QKxH = val;
	if (GLPS_SUCCESS==glps_read(*i,"phi",val)) QPhi = val;
	if (GLPS_SUCCESS==glps_read(*i,"n",val)) k1 = (long) floor(val+0.5);
	if (GLPS_SUCCESS==glps_read(*i,"lambda",val)) QKS = val;
	if (GLPS_SUCCESS==glps_read(*i,"roll",val)) dt = val;
	if (GLPS_SUCCESS==glps_read(*i,"method",val)) k2 = (long) floor(val+0.5);
	if (k2 != Meth_Linear) globval.MatMeth = false;
	if ((unsigned int)k2 >= 32 ||
	    ((1 << k2) &
	     ((1 << Meth_Linear) | (1 << Meth_First) | (1 << Meth_Second) |
	      (1 << Meth_Fourth) | (1 << Meth_genfun))) == 0)
	  longjmp(env0, 1);
	if (GLPS_SUCCESS==glps_read(*i,"harm",vec))
	  for(int n=0; n<(int)(vec.size())/6;++n)
	    { 
	      n_harm++;harm[n] = (long) floor(vec[6*n] + 0.5);
	      if (harm[n] < 1 || n_harm_max < n + 2)
		{
		  cerr << "invalid value detected" << endl;
		  longjmp(env0, 1);
		}
	      kxV[n] = vec[6*n+1];
	      BoBrhoV[n] = vec[6*n+2];
	      kxH[n] = vec[6*n+3];
	      BoBrhoH[n] = vec[6*n+4];
	      phi[n] = vec[6*n+5];
	    }



	//   case dbnsym:
	//     GetDBN_(&V);
	//     break;

	globval.Elem_nFam++;
	if (globval.Elem_nFam <= Elem_nFamMax) {
	  WITH = &ElemFam[globval.Elem_nFam-1]; WITH1 = &WITH->ElemF;
	  memset(WITH1->PName,0,sizeof(partsName));
	  memcpy(WITH1->PName, (*i).c_str(), (*i).length());
	  //    memcpy(WITH1->PName, ElementName, sizeof(partsName));
	  WITH1->PL = QL; WITH1->Pkind = Wigl;
	  Wiggler_Alloc(&WITH->ElemF); WITH4 = WITH1->W;
	  WITH4->Pmethod = k2; WITH4->PN = k1;
	  WITH4->PdTpar = dt;
	  //    SetDBN(&V);
	  WITH4->lambda = QKS; WITH4->n_harm = 1; WITH4->harm[0] = 1;
	  WITH4->kxV[0] = QKxV; WITH4->BoBrhoV[0] = QKV;
	  WITH4->kxH[0] = QKxH; WITH4->BoBrhoH[0] = QKH;
	  WITH4->phi[0] = QPhi;
	  //                    AssignHarm(globval.Elem_nFam, &V);
	  //   W = ElemFam[elem-1].ElemF.W; 
	  WITH4->n_harm += n_harm;
	  // the fundamental is stored in harm[0], etc.
	  for (int k = 1; k < WITH4->n_harm; k++) {
	    WITH4->harm[k] = harm[k-1];
	    WITH4->kxV[k] = kxV[k-1]; WITH4->BoBrhoV[k] = BoBrhoV[k-1];
	    WITH4->kxH[k] = kxH[k-1]; WITH4->BoBrhoH[k] = BoBrhoH[k-1];
	    WITH4->phi[k] = phi[k-1];
	  }

	  /* Equivalent vertically focusing gradient */
	  WITH4->PBW[HOMmax+2] = -QK*QK/2e0;
	  if (!CheckWiggler(globval.Elem_nFam))
	    longjmp(env0, 1); //Always returns true????
	} else {
	  cerr << "Elem_nFamMax exceeded: " << globval.Elem_nFam
	       << "(" << (long)Elem_nFamMax << ")" << endl;
	  longjmp(env0, 1);
	}
	break;

      case  10: //Fieldmap
	QL = 0.0; k1 = 0;
	if (GLPS_SUCCESS==glps_read(*i,"l",val)) QL = val;
	if (GLPS_SUCCESS==glps_read(*i,"t",val)) t = val;
	if (GLPS_SUCCESS==glps_read(*i,"n",val)) k1 = (long) floor(val+0.5);
	if (GLPS_SUCCESS!=glps_read(*i,"file1",str1)) str1="";
	if (GLPS_SUCCESS==glps_read(*i,"scaling",val)) scaling = val;
	globval.Elem_nFam++;
	if (globval.Elem_nFam <= Elem_nFamMax) {
	  WITH = &ElemFam[globval.Elem_nFam-1]; WITH1 = &WITH->ElemF;
	  memset(WITH1->PName,0,sizeof(partsName));
	  memcpy(WITH1->PName, (*i).c_str(), (*i).length());
	  //    memcpy(WITH1->PName, ElementName, sizeof(partsName));
	  WITH1->PL = QL; WITH1->Pkind = FieldMap;
	  FieldMap_Alloc(WITH1);
	  //,true);
	  WITH6 = WITH1->FM; 
	  WITH6->phi = t*M_PI/180.0;WITH6->n_step = k1; WITH6->scl = scaling;

	  if (glps_read("energy",val)==GLPS_SUCCESS) {
	    globval.Energy=val;
	    if (str1.length() > 0) get_B(str1.c_str(), WITH6);
	  } else {
	    cerr << "Fieldmap: energy not defined" << endl;
	    longjmp(env0, 1);
	  }
	} else {
	  cerr << "Elem_nFamMax exceeded: " << globval.Elem_nFam
	       << "(" << (long)Elem_nFamMax << ")" << endl;
	  longjmp(env0, 1);
	}
	break;

      case   11: // insertion
	QK  = 0e0; QKxV = 0e0; QKS = 0e0;
	k1  = 1;
	k2  = 1;       // linear interpolation by default
	dt  = 0e0;
	scaling = 1.0; // scaling factor
	//                 string str1,str2;
	//     bool firstflag=true, secondflag=true;
	if (GLPS_SUCCESS==glps_read(*i,"n",val)) k1 = (long) floor(val+0.5);
	if (GLPS_SUCCESS!=glps_read(*i,"file1",str1))
	  {str1=""; firstflag=false; }
	if (GLPS_SUCCESS!=glps_read(*i,"file2",str2))
	  {str2=""; secondflag=false; }
	if (GLPS_SUCCESS==glps_read(*i,"scaling",val))
	  scaling = abs((long)floor(val));      
	if (GLPS_SUCCESS==glps_read(*i,"method",val))
	  // method for interpolation: 1 means linear 2 spline
	  k2 = (long) floor(val+0.5);
	if (k2 != Meth_Linear) globval.MatMeth = false;

	globval.Elem_nFam++;
	/* Fills up the ID */
	if (globval.Elem_nFam <= Elem_nFamMax) {
	  WITH  = &ElemFam[globval.Elem_nFam-1];
	  WITH1 = &WITH->ElemF;
	  memset(WITH1->PName,0,sizeof(partsName));
	  memcpy(WITH1->PName, (*i).c_str(), (*i).length());
	  //    memcpy(WITH1->PName, ElementName, sizeof(partsName));
	  WITH1->Pkind = Insertion;
	  Insertion_Alloc(&WITH->ElemF);
	  WITH5 = WITH1->ID;
	  WITH5->Pmethod = k2;
	  WITH5->PN = k1;
	  WITH5->scaling = scaling;

	  // if (glps_read("energy",val)==GLPS_SUCCESS) {
	  //   globval.Energy=val;
	  //   if (str1.length() > 0) get_B(str1.c_str(), WITH6);
	  // } else {
	  //   cerr << "Insertion_Alloc: energy not defined" << endl;
	  //   longjmp(env0, 1);
	  // }

	  // Check if filename given for first order kicks
	  //    if (firstflag) 
	  if (str1.length() > 0) {
	    strcpy(WITH5->fname1,str1.c_str());
	    WITH5->firstorder = true;
	    Read_IDfile(WITH5->fname1, WITH1->PL, WITH5->nx, WITH5->nz,
			WITH5->tabx, WITH5->tabz, WITH5->thetax1,
			WITH5->thetaz1, WITH5->long_comp, WITH5->B2);
	    // scale factor from Radia: Tmm to get Tm.
	    for (kx = 0; kx < WITH5->nx; kx++) {
	      for (kz = 0; kz < WITH5->nz; kz++) {
		WITH5->thetax1[kz][kx] = 1e-3*WITH5->thetax1[kz][kx];
		WITH5->thetaz1[kz][kx] = 1e-3*WITH5->thetaz1[kz][kx];
	      }
	    }
	  } else {
	    strcpy(WITH5->fname1,"/*No_Filename1_Given*/");
	  }

	  // Check if filename given for Second order kicks
	  //    if (secondflag) 
	  if (str2.length() >0) {
	    //        strcpy(WITH5->fname2,"/*No_Filename2_Given*/");
	    strcpy(WITH5->fname2,str2.c_str());
	    WITH5->secondorder = true;;
	    // Read Id file for second order kicks
	    Read_IDfile(WITH5->fname2, WITH1->PL, WITH5->nx, WITH5->nz,
			WITH5->tabx, WITH5->tabz, WITH5->thetax, WITH5->thetaz,
			WITH5->long_comp, WITH5->B2);
	  } else {
	    strcpy(WITH5->fname2,"/*No_Filename2_Given*/");
	  }

	  // check whether no Radia filename read: something is wrong
	  //    if (!firstflag && !secondflag) 
	  if (str1.length()<1 && str2.length()<1) {
	    cerr << "Error: no Insertion filename found as"
		 << " an input in lattice file" << endl;
	    longjmp(env0, 1);
	  }

	  if (k2 != 1) { // linear interpolation
	    WITH5->linear = false;
	  } else { // cubic interpolation
	    WITH5->linear = true;
	  }
	  // stuff for spline interpolation
	  if (!WITH5->linear) {
	    WITH5->mtx = gsl_matrix_alloc(WITH5->nz,WITH5->nx);
	    GSL2NRDM2(pmtx,WITH5->mtx,WITH5->tx,0);
	    WITH5->mtz = gsl_matrix_alloc(WITH5->nz,WITH5->nx);
	    GSL2NRDM2(pmtz,WITH5->mtz,WITH5->tz,0);
	    WITH5->tab1 = (double *)malloc((WITH5->nx)*sizeof(double));
	    WITH5->tab2 = (double *)malloc((WITH5->nz)*sizeof(double));
	    WITH5->mf2x = gsl_matrix_alloc(WITH5->nz,WITH5->nx);
	    GSL2NRDM2(pmf2x,WITH5->mf2x,WITH5->f2x,0);
	    WITH5->mf2z = gsl_matrix_alloc(WITH5->nz,WITH5->nx);
	    GSL2NRDM2(pmf2z,WITH5->mf2z,WITH5->f2z,0);
	    Matrices4Spline(WITH5);
	  }

	  // to put somewhere
	  //      /** Free memory **/
	  //      free(tab1);
	  //      free(tab2);
	  //
	  //      free_matrix(tx,1,nz,1,nx);
	  //      free_matrix(tz,1,nz,1,nx);
	  //      free_matrix(f2x,1,nz,1,nx);
	  //      free_matrix(f2z,1,nz,1,nx);
	} else {
	  cerr << "Elem_nFamMax exceeded: " << globval.Elem_nFam
	       << "(" << (long)Elem_nFamMax << ")" << endl;
	  longjmp(env0, 1);
	}
	break;

      case   12: //Spreader
	//    if (sym1 == dbnsym)
	//        GetDBN_(&V);
	globval.Elem_nFam++;
	if (globval.Elem_nFam <= Elem_nFamMax) {
	  WITH = &ElemFam[globval.Elem_nFam-1];
	  WITH1 = &WITH->ElemF;
	  memset(WITH1->PName,0,sizeof(partsName));
	  memcpy(WITH1->PName, (*i).c_str(), (*i).length());
	  //    memcpy(WITH1->PName, ElementName, sizeof(partsName));
	  if (GLPS_SUCCESS==glps_read(*i,"l",val))  WITH1->PL = val;
	  //    WITH1->PL = *V.rnum;
	  WITH1->Pkind = PartsKind(Spreader);
	  Spreader_Alloc(&WITH->ElemF);
	} else {
	  cerr << "Elem_nFamMax exceeded: " << globval.Elem_nFam
	       << "(" << (long)Elem_nFamMax << ")" << endl;
	  longjmp(env0, 1);
	}
	break;

      case   13: // recombiner

	//      GetDBN_(&V);
	globval.Elem_nFam++;
	if (globval.Elem_nFam <= Elem_nFamMax) {
	  WITH = &ElemFam[globval.Elem_nFam-1];
	  WITH1 = &WITH->ElemF;
	  memset(WITH1->PName,0,sizeof(partsName));
	  memcpy(WITH1->PName, (*i).c_str(), (*i).length());
	  if (GLPS_SUCCESS==glps_read(*i,"l",val))  WITH1->PL = val;
	  //    WITH1->PL = *V.rnum;
	  WITH1->Pkind = PartsKind(Recombiner);
	  Recombiner_Alloc(&WITH->ElemF);
	} else {
	  cerr << "Elem_nFamMax exceeded: " << globval.Elem_nFam
	       << "(" << (long)Elem_nFamMax << ")" << endl;
	  longjmp(env0, 1);
	}
	break;

      case    14: // Solenoid
	QL = 0.0; /* L */
	QK = 0.0; /* K */
	k1 = 0;   /* N */
	if (GLPS_SUCCESS==glps_read(*i,"l",val)) QL = val;
	if (GLPS_SUCCESS==glps_read(*i,"bobrho",val)) QK = val;
	if (GLPS_SUCCESS==glps_read(*i,"n",val)) k1 = (long) floor(val+0.5);
	globval.Elem_nFam++;
	if (globval.Elem_nFam <= Elem_nFamMax) {
	  WITH = &ElemFam[globval.Elem_nFam-1];
	  WITH1 = &WITH->ElemF;
	  Solenoid_Alloc(&WITH->ElemF);
	  WITH7 = WITH1->Sol;
	  memset(WITH1->PName,0,sizeof(partsName));
	  memcpy(WITH1->PName, (*i).c_str(), (*i).length());
	  WITH1->Pkind = Solenoid;
	  WITH1->PL = QL; WITH7->N = k1; WITH7->BoBrho = QK;
	} else {
	  cerr << "Elem_nFamMax exceeded: " << globval.Elem_nFam
	       << "(" << (long)Elem_nFamMax << ")" << endl;
	  longjmp(env0, 1);
	}
	break;

      default:

	cout << " Not Defined " << str << endl;
      }  //switch


  }   // elem
            

  cout << "Beamline <" << use << "> will be used" << endl;
  //    return 0;
    
  std::list<string> lines = glps_lines();
  //  cout << "Number of lines: " << lines.size() << endl;
  string beamline;
  for (std::list<string>::iterator i = lines.begin(); i != lines.end(); ++i) {
    beamline=*i;
    for (int k=0; k<(int)beamline.length();k++)
      beamline[k]=tolower(beamline[k]);
    if (!use.compare(beamline)) break;
  }
  //  cout << *i << ": " << endl;
  std::list<string> blelems = glps_elements(beamline, true);
  int cnum = 0;
  static int fnum = 0;
  for (std::list<string>::iterator j = blelems.begin(); j != blelems.end(); ++j) {
    size_t pt = (*j).rfind('.');
    string tmp = (*j).substr(pt+1);
    cnum++;
    fnum = CheckElementtable(tmp);
    if (cnum <= Cell_nLocMax)
      Cell[cnum].Fnum = fnum;
    else {
      cerr << "** Cell_nLocMax exhausted: " << cnum
	   << "(" << Cell_nLocMax << ")" << endl;
      longjmp(env0, 1);
    }
    //    cout << cnum << " elements in Cell" << endl;
  }

  if (cnum <= Cell_nLocMax)
    globval.Cell_nLoc = cnum;   /*number of Elements in a cell*/
  else {
    cerr << "** Cell_nLocMax exhausted: " << cnum
	 << "(" << Cell_nLocMax << ")" << endl;
    longjmp(env0, 1);
  }
  GetRingType();/* define whether a ring or a transfer line */
  GetCODEPS();  /* define COD precision */
  GetDP();      /* define energy offset */
  free_lattice();


  RegisterKids(); /* Check wether too many elements */
  //  PrintResult();  /* Print lattice statistics */

  return (true);
}



void ClearHOM(double *B, bool *BA)
{
  long i;

  for (i = -HOMmax; i <= HOMmax; i++) {
    B[i + HOMmax] = 0.0;
    BA[i + HOMmax] = false;
  }
  //  LINK->DBNsavemax = 0;
}


void AssignHOM(long elem,  double *B, bool * BA )
{
  long       i;
  MpoleType  *M;

  M = ElemFam[elem-1].ElemF.M;
  for (i = -HOMmax; i <= HOMmax; i++) {
    if (BA[i+HOMmax]) {
      M->PBpar[i+HOMmax] = B[i+HOMmax];
      M->Porder = max(abs(i), M->Porder);
    }
  }
}

/*
  void AssignHarm(long elem, int *harm, double * kxV, double * BoBrhoV, double * kxH, double * BoBrhoH, double * phi)
  {
  long         i;
  WigglerType  *W;

  W = ElemFam[elem-1].ElemF.W; W->n_harm += LINK->n_harm;
  // the fundamental is stored in harm[0], etc.
  for (i = 1; i < W->n_harm; i++) {
  W->harm[i] = LINK->harm[i-1];
  W->kxV[i] = LINK->kxV[i-1]; W->BoBrhoV[i] = LINK->BoBrhoV[i-1];
  W->kxH[i] = LINK->kxH[i-1]; W->BoBrhoH[i] = LINK->BoBrhoH[i-1];
  W->phi[i] = LINK->phi[i-1];
  }
  }
*/

bool CheckWiggler( long i)
{
  bool         Result;
  double       a, Lambda, L, diff;
  long         NN;
  ElemFamType  *WITH;
  elemtype     *WITH1;
  WigglerType  *WITH2;

  Result = false;
  WITH = &ElemFam[i-1]; WITH1 = &WITH->ElemF; WITH2 = WITH1->W;
  Lambda = WITH2->lambda;
  L = WITH1->PL; a = L/Lambda;
  NN = (long)floor(a+0.01+0.5);
  diff = fabs((L-NN*Lambda)/L);
  if (diff < 1e-5) return true;
  cerr << endl << ">>> Incorrect definition of " << WITH1->PName << endl;
  cerr         << "    L      ( total length ) = " << L << endl;           //%20.12f [m]\n", L);
  cerr         << "    Lambda ( wave  length ) = " << Lambda << endl;      //%20.12f [m]\n", L);
  cerr         << "    # of Period = L /Lambda = " << L/Lambda << "?????" << endl;     //     
  return false;
}
/*
  static void CheckWiggler(long i, struct LOC_Lat_DealElement *LINK)
  {
  if (!Lat_CheckWiggler(LINK->fo, i, LINK->LINK))
  longjmp(LINK->_JL9999, 1);
  }
*/


long CheckElementtable(const string name)
{
  /* globval.Elem_nFam = Number of parts in a Element */
  long  i, j, FORLIM;

  j = 0;
  if (globval.Elem_nFam > Elem_nFamMax) {
    cerr << "Elem_nFamMax exceeded: " << globval.Elem_nFam
	 << "(" << (long)Elem_nFamMax << ")" << endl;
    longjmp(env0, 1);
  }

  //  if (strstr(LINK->line,"insertion") != NULL) return 0;

  FORLIM = globval.Elem_nFam;
  for (i = 1; i <= FORLIM; i++) {
    //  if (!strncmp(ElemFam[i-1].ElemF.PName, name.c_str(), name.length()))
    if (!strcmp(ElemFam[i-1].ElemF.PName, name.c_str()))
      j = i;
  }
  return j;
}

void GetEnergy()
{
  double val;

  if (glps_read("energy", val)!=GLPS_SUCCESS )
    {
      cout << "> Energy is not defined, default is 3.0 GeV" << endl;
      globval.Energy = 3.0;
    } else 
    globval.Energy = val;
}


void GetRingType()
{
  double val;

  if (glps_read("ringtype", val)!=GLPS_SUCCESS || floor(val+0.5) < 0 ||
      floor(val+0.5) > 1)
    {
      cout << "> Ring type is not defined or not correctly, default is ring."
	   << endl;
      globval.RingType = 1;
    } else 
    globval.RingType = floor(val+0.5);
}


void GetDP()
{
  double val;

  if (glps_read("dp", val)!=GLPS_SUCCESS )
    {
      cout << "> dP/P is not defined, default is 1.0e-8." << endl;
      globval.dPcommon = 1e-8;
    } else 
    globval.dPcommon = val;
}

void GetCODEPS()
{
  double val;

  if (glps_read("codeps", val)!=GLPS_SUCCESS )
    {
      cout << "> CODEPS is not defined, default is 1.0e-12." << endl;
      globval.CODeps = 1e-12;
    } else 
    globval.CODeps = val;
}

double Circumference()
{
  long i;
  double S;
  long FORLIM;

  S = 0.0;
  FORLIM = globval.Cell_nLoc;
  for (i = 1; i <= FORLIM; i++)
    S += ElemFam[Cell[i].Fnum - 1].ElemF.PL;
  return S;
}



void RegisterKids()
{
  long i, FORLIM;
  ElemFamType *WITH;

  if (setjmp(env0)) return; 
  if (globval.Elem_nFam <= Elem_nFamMax) {
    FORLIM = globval.Elem_nFam;
    for (i = 0; i < FORLIM; i++)
      ElemFam[i].nKid = 0;
  } else {
    cerr << "Elem_nFamMax exceeded: " << globval.Elem_nFam
	 << "(" << (long)Elem_nFamMax << ")" << endl;
    longjmp(env0, 1);
  }

  FORLIM = globval.Cell_nLoc;
  for (i = 1; i <= FORLIM; i++) {
    WITH = &ElemFam[Cell[i].Fnum - 1];
    WITH->nKid++;
    if (WITH->nKid <= nKidMax) {
      WITH->KidList[WITH->nKid - 1] = i;
      Cell[i].Knum = WITH->nKid;
    } else
      {
	cerr << "Elem_nFamMax exceeded: "
	     << WITH->nKid << "(" << nKidMax << ")" << endl;
	longjmp(env0, 1);
      }
  }
}


/*
  void PrintResult()
  {
  long j, nKid, FORLIM;
  struct tm *newtime;

  // Get time and date 
  newtime = GetTime();

  printf("\n");
  printf("  TRACY III v. 3.5 compiled on %s\n",__DATE__);
  printf("\n");
  printf("  LATTICE Statistics for today %s \n\n", asctime2(newtime));
  printf("  Number of constants: UDIC                 =%5ld"
  ", UDImax          =%5d\n",
  LINK->UDIC, UDImax);
  printf("  Number of keywords : nkw                  =%5ld"
  ", Lat_nkw_max     =%5d\n",
  LINK->nkw, Lat_nkw_max);
  printf("  Number of Families : globval.Elem_nFam    =%5ld"
  ", Elem_nFamMax    =%5d\n",
  globval.Elem_nFam, Elem_nFamMax);
  nKid = 0L;
  FORLIM = globval.Elem_nFam;
  for (j = 0L; j < FORLIM; j++) {
  if (ElemFam[j].nKid > nKid)
  nKid = ElemFam[j].nKid;
  }
  printf("  Max number of Kids : nKidMax              =%5ld"
  ", nKidMax         =%5d\n",
  nKid, nKidMax);
  printf("  Number of Blocks   : NoB                  =%5ld"
  ", NoBmax          =%5d\n",
  LINK->NoB, NoBmax);
  printf("  Max Block size     : NoBE                 =%5ld"
  ", NoBEmax         =%5d\n",
  LINK->Bpointer, NoBEmax);
  printf("  Number of Elements : globval.Cell_nLoc    =%5ld"
  ", Cell_nLocmax    =%5d\n",
  globval.Cell_nLoc, Cell_nLocMax);
  printf("  Circumference      : %12.7f [m]\n", Circumference(LINK));
  printf("\n");
  printf("\n");
  }

*/
