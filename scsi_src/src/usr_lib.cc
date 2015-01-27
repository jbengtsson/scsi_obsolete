

void get_bn(const string type, const int Fnum, const int Knum,
	    const int n, double &bn, double &an)
{
  MpoleType *M;

  if ((n < 1) || (n > HOMmax)) {
    cout << "get_bn: n < 1 (" << n << ")" << endl;
    exit(1);
  }

  M = static_cast<MpoleType*>(Cell[Elem_GetPos(Fnum, Knum)].Elem);

  if (type.compare("dsgn") == 0) {
    bn = M->bnpar[HOMmax+n]; an = M->bnpar[HOMmax-n];
  } else if (type.compare("sys") == 0) {
    bn = M->bnsys[HOMmax+n]; an = M->bnsys[HOMmax-n];
  } else if (type.compare("rms") == 0) {
    bn = M->bnrms[HOMmax+n]; an = M->bnrms[HOMmax-n];
  } else {
    cout << "get_bn: undef. type" << type << endl;
    exit(1);
  }
}


void set_bn(const string type, const int Fnum, const int Knum,
	    const int n, const double bn, const double an)
{
  MpoleType *M;

  if ((n < 1) || (n > HOMmax)) {
    cout << "set_bn: n < 1 (" << n << ")" << endl;
    exit(1);
  }

  M = static_cast<MpoleType*>(Cell[Elem_GetPos(Fnum, Knum)].Elem);

  if (type.compare("dsgn") == 0) {
    M->bnpar[HOMmax+n] = bn; M->bnpar[HOMmax-n] = an;
  } else if (type.compare("sys") == 0) {
    M->bnsys[HOMmax+n] = bn; M->bnsys[HOMmax-n] = an;
  } else if (type.compare("rms") == 0) {
    M->bnrms[HOMmax+n] = bn; M->bnrms[HOMmax-n] = an;
  } else {
    cout << "get_bn: undef. type" << type << endl;
    exit(1);
  }

  Mpole_Setbn(Fnum, Knum, n); Mpole_Setbn(Fnum, Knum, -n);
}


void set_bn(const string type, const int Fnum,
	    const int n, const double bn, const double an)
{
  int k;

  if ((n < 1) || (n > HOMmax)) {
    cout << "set_bn: n < 1 (" << n << ")" << endl;
    exit(1);
  }

  for (k = 1; k <= GetnKid(Fnum); k++)
    set_bn(type, Fnum, k, n, bn, an);
}
