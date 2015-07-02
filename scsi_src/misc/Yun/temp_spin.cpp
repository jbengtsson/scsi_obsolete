
/* My code is only valid for protons, not for electrons.

   Coordinates are x[9]: x[0] to x[5] as yours.
   x[6] to x[8] are spin vectors for x, y, s direction.


Yun                                                                           */



//=================================
//
//      spin precession
//
//=================================


template <class T> void cal_Bfield(T x[], int Nint, int Norder, double Angle, double KNL[11], double KNSL[11],T BLbrho[])
{
  int i;
  int fac=1;
  T  Xn, Yn, Xn0, Yn0;
  T  By, Bx;

  By=KNL[0];
  Bx=KNSL[0];
  Xn=1.;
  Yn=0.;

  for(i=1;i<Norder+1;i++){
    Xn0=Xn;
    Yn0=Yn;
    Xn=Xn0*x[x_]-Yn0*x[y_];
    Yn=Xn0*x[y_]+Yn0*x[x_];
    fac=fac*i;
    if ( KNL[i] != 0. || KNSL[i] !=0. ) {
      By=By+(KNL[i]*Xn-KNSL[i]*Yn)/fac;
      Bx=Bx+(KNL[i]*Yn+KNSL[i]*Xn)/fac;
    }
  }

  BLbrho[0]=Bx/Nint;
  BLbrho[1]=(Angle + By)/Nint;
  BLbrho[2]=0.;
}


void spin_kick_pass(double x[], double angle, double href, double BLbrho[])
{
  int i,j;
  double vp[3], v, u[3];
  double Bpar[3], udotB;
  double temp, pz;
  double delta1=sqrt( 1.0 + 2.*x[pt_]/GP.beta+x[pt_]*x[pt_])-1.0;
  double gamma1=GP.gamma+sqrt(GP.gamma*GP.gamma-1.)*x[pt_];
  double Gr   = GP.G * gamma1;
  double F[3], Fabs, uF[3];
  double theta, costheta, sintheta;
  double rot[9];
  double spin0[3], spin[3];

  if(  BLbrho[0]* BLbrho[0] +  BLbrho[1]* BLbrho[1] +  BLbrho[2]* BLbrho[2] !=0. ){

    for(i=0;i<3;i++) spin0[i]=x[6+i];

    temp=sqrt( (1+delta1) * (1+delta1)-x[px_]*x[px_]-x[py_]*x[py_] );
    u[0] = x[px_] / temp;
    u[1] = x[py_] / temp;
    u[2] = 1.0 ;
    temp=sqrt( u[0]*u[0] +  u[1]*u[1] +  u[2]*u[2] );
    for(i=0;i<3;i++) u[i]=u[i]/temp;
    temp=(1 + href * x[0])* temp;

    udotB =  u[0] * BLbrho[0] +  u[1] * BLbrho[1] +  u[2] * BLbrho[2];
    for(i=0;i<3;i++) Bpar[i] =  udotB * u[i];
    for(i=0;i<3;i++) F[i] = ( 1+ Gr ) * BLbrho[i]  + GP.G*(1-gamma1)*Bpar[i] ;
    for(i=0;i<3;i++) F[i] = F[i] * temp / ( 1+ delta1) ;
    F[1]=F[1] - angle;
    Fabs=sqrt(F[0]* F[0] + F[1]*F[1]  +  F[2]* F[2] );
    for(i=0;i<3;i++) uF[i] = F[i]/Fabs;
    theta=-Fabs; costheta = cos(theta); sintheta= sin(theta);

    rot[0] = uF[0]*uF[0]*(1-costheta) + costheta ;
    rot[1] = uF[0]*uF[1]*(1-costheta) - uF[2]*sintheta;
    rot[2] = uF[0]*uF[2]*(1-costheta) + uF[1]*sintheta;
    rot[3] = uF[0]*uF[1]*(1-costheta) + uF[2]*sintheta;
    rot[4] = uF[1]*uF[1]*(1-costheta) + costheta ;
    rot[5] = uF[1]*uF[2]*(1-costheta) - uF[0]*sintheta;
    rot[6] = uF[0]*uF[2]*(1-costheta) - uF[1]*sintheta;
    rot[7] = uF[1]*uF[2]*(1-costheta) + uF[0]*sintheta;
    rot[8] = uF[2]*uF[2]*(1-costheta) + costheta;
    for(i=0;i<3;i++){
      spin[i]=0.;
      for(j=0;j<3;j++) spin[i] = spin[i] + rot[i*3 + j] * spin0[j];
    }

    for(i=0;i<3;i++) x[6+i]= spin[i];
  }
}


void SBEND_sPass(double x[], double L, int Nint, double Angle, double E1, double E2)
{
  int i,j;
  double href=Angle/L;
  double Lint=L/Nint;
  double edgefocus;
  double BLbrho[3];
  double x1[9];

  if(Angle==0.){
    DRIFT_Pass(x,L); return;
  }

  edgefocus = tan(E1)*href;
  x[1] = x[1]+edgefocus*x[0];
  x[3] = x[3]-edgefocus*x[2];

  BLbrho[0] = 0.;
  BLbrho[1] = Angle/Nint;
  BLbrho[2] = 0.;
  if(GP.H_expand == true){
    for(i=0;i<Nint;i++){
      for(j=0;j<9;j++) x1[j]= x[j];
      DRIFT_Pass(x,Fdrift1*Lint);
      bend_kick_pass(x, Fkick1*Lint, href);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_kick_pass(x, Fkick2*Lint, href);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_kick_pass(x, Fkick1*Lint, href);
      DRIFT_Pass(x,Fdrift1*Lint);
      for(j=0;j<9;j++) x1[j]=(x[j]+x1[j])/2.;
      spin_kick_pass(x1, Angle/Nint, href, BLbrho);
      x[6]= x1[6];  x[7]= x1[7];    x[8]= x1[8];
    }
  }
  else{
    for(i=0;i<Nint;i++){
      for(j=0;j<9;j++) x1[j]= x[j];
      sbend_exact_pass(x,Lint,Angle/Nint,cos(Angle/Nint),sin(Angle/Nint));
      for(j=0;j<9;j++) x1[j]=(x[j]+x1[j])/2.;
      spin_kick_pass(x1, Angle/Nint, href, BLbrho);
      x[6]= x1[6];  x[7]= x1[7];    x[8]= x1[8];
    }
  }

  edgefocus = tan(E2)*href;
  x[1] = x[1]+edgefocus*x[0];
  x[3] = x[3]-edgefocus*x[2];
}



void MULT_sPass(double x[], double L, int Nint, int Norder, double KNL[11], double KNSL[11])
{
  int i,j;
  double Lint=L/Nint;
  double preal, pimag;
  double BLbrho[3];
  double x1[9];

  if(L==0.)
    {
      cal_Bfield(x, 1, Norder, 0, KNL, KNSL, BLbrho);
      spin_kick_pass(x,0,0,BLbrho);
      mult_kick_pass(x, Norder, KNL, KNSL);
    }
  else
    {
      double knl_kick1[11],knsl_kick1[11];
      double knl_kick2[11],knsl_kick2[11];
      for(i=0;i<11;i++) {
	knl_kick1[i] =Fkick1*KNL[i]/Nint;
	knsl_kick1[i]=Fkick1*KNSL[i]/Nint;
      }
      for(i=0;i<11;i++) {
	knl_kick2[i] =Fkick2*KNL[i]/Nint;
	knsl_kick2[i]=Fkick2*KNSL[i]/Nint;
      }
      for(i=0;i<Nint;i++){
        for(j=0;j<9;j++) x1[j]=x[j];
	DRIFT_Pass(x,Fdrift1*Lint);
	mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
	DRIFT_Pass(x,Fdrift2*Lint);
	mult_kick_pass(x, Norder,knl_kick2, knsl_kick2);
	DRIFT_Pass(x,Fdrift2*Lint);
	mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
	DRIFT_Pass(x,Fdrift1*Lint);
	for(j=0;j<9;j++) x1[j]=(x[j]+x1[j])/2.;
	cal_Bfield(x1, Nint, Norder, 0, KNL, KNSL, BLbrho);
	spin_kick_pass(x1,0,0,BLbrho);
        x[6]= x1[6];  x[7]= x1[7];    x[8]= x1[8];
      }
    }
}

void SNAKE_sPass(double x[9], double L, double n[3], double angle)
{
  int i,j;
  double spin0[3], spin[3];
  double F, uF[3], costheta=cos(-angle), sintheta=sin(-angle), rot[9];

  DRIFT_Pass(x, L/2.);

  for(i=0;i<3;i++) spin0[i]=x[6+i];

  F=sqrt( n[0]*n[0] +  n[1]*n[1] +  n[2]*n[2] );
  for(i=0;i<3;i++) uF[i]= n[i]/F;

  rot[0] = uF[0]*uF[0]*(1-costheta) + costheta ;
  rot[1] = uF[0]*uF[1]*(1-costheta) - uF[2]*sintheta;
  rot[2] = uF[0]*uF[2]*(1-costheta) + uF[1]*sintheta;
  rot[3] = uF[0]*uF[1]*(1-costheta) + uF[2]*sintheta;
  rot[4] = uF[1]*uF[1]*(1-costheta) + costheta ;
  rot[5] = uF[1]*uF[2]*(1-costheta) - uF[0]*sintheta;
  rot[6] = uF[0]*uF[2]*(1-costheta) - uF[1]*sintheta;
  rot[7] = uF[1]*uF[2]*(1-costheta) + uF[0]*sintheta;
  rot[8] = uF[2]*uF[2]*(1-costheta) + costheta;

  for(i=0;i<3;i++){
    spin[i]=0.;
    for(j=0;j<3;j++) spin[i] = spin[i] + rot[i*3 + j] * spin0[j];
  }
  for(i=0;i<3;i++) x[6+i]= spin[i];

  DRIFT_Pass(x,L/2.);
}

