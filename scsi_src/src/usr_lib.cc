

void get_bn(const string type, const int Fnum, const int Knum,
	    const int n, double &bn, double &an)
{
  ElemType  elem;

  if ((n < 1) || (n > HOMmax)) {
    cout << "get_bn: n < 1 (" << n << ")" << endl;
    exit(1);
  }

  elem = Cell[Elem_GetPos(Fnum, Knum)].Elem;

  if (type.compare("dsgn") == 0) {
    bn = elem.M->PBpar[HOMmax+n]; an = elem.M->PBpar[HOMmax-n];
  } else if (type.compare("sys") == 0) {
    bn = elem.M->PBsys[HOMmax+n]; an = elem.M->PBsys[HOMmax-n];
  } else if (type.compare("rms") == 0) {
    bn = elem.M->PBrms[HOMmax+n]; an = elem.M->PBrms[HOMmax-n];
  } else {
    cout << "get_bn: undef. type" << type << endl;
    exit(1);
  }
}


void set_bn(const string type, const int Fnum, const int Knum,
	    const int n, const double bn, const double an)
{
  ElemType  elem;

  if ((n < 1) || (n > HOMmax)) {
    cout << "set_bn: n < 1 (" << n << ")" << endl;
    exit(1);
  }

  elem = Cell[Elem_GetPos(Fnum, Knum)].Elem;

  if (type.compare("dsgn") == 0) {
    elem.M->PBpar[HOMmax+n] = bn; elem.M->PBpar[HOMmax-n] = an;
  } else if (type.compare("sys") == 0) {
    elem.M->PBsys[HOMmax+n] = bn; elem.M->PBsys[HOMmax-n] = an;
  } else if (type.compare("rms") == 0) {
    elem.M->PBrms[HOMmax+n] = bn; elem.M->PBrms[HOMmax-n] = an;
  } else {
    cout << "get_bn: undef. type" << type << endl;
    exit(1);
  }

  Mpole_SetPB(Fnum, Knum, n); Mpole_SetPB(Fnum, Knum, -n);
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
