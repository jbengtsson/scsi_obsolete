
canonical variables

enum ps_index { x_ = 0, px_ = 1, y_ = 2, py_ = 3,  z_ = 4, pt_= 5 }; 
//  coordinate system: right handed, x^ cross y^ =  s^, bending angle >0, bend particle to -x direction.
//  canonical variables: (x, px, y, py, -c*dt, pt); 
//  z= -cdt, difference in time of flight divided by c,  z > 0 arrive earlier than reference particle.
//  pt_= dE/(P_0 c)= (E- E_0)/ (P_0 c), P_0 on momentum, E_0, on-momentum total enrgy



template <class T>  void DRIFT_Pass(T x[], double L)
{
  if( GP.ultra_rel == true){ 
    T u;
    u= L/(1.0+x[pt_]); 
    x[x_]=x[x_]+x[px_]*u;
    x[y_]=x[y_]+x[py_]*u;
    x[z_]=x[z_]-(x[px_]*x[px_]+x[py_]*x[py_])*u/2.0/(1+x[pt_]);
  }
  else{
    T pz, betaz;
    T delta1,gamma1,beta1;

    delta1= sqrt(1.0 + 2*x[pt_]/GP.beta+x[pt_]*x[pt_]) -1.0;
    gamma1 =GP.gamma+sqrt(GP.gamma*GP.gamma-1)*x[pt_];
    beta1 = sqrt(1.0-1.0/gamma1/gamma1);

    pz = sqrt( (1.0+delta1)*(1.0+delta1)-x[px_]*x[px_]-x[py_]*x[py_]);
    betaz=pz*beta1/(1.0+delta1);
    x[x_]=x[x_]+x[px_]*L/pz; 
    x[y_]=x[y_]+x[py_]*L/pz;
    x[z_]=x[z_]-(L/betaz-L/GP.beta);
  }
}


template <class T>  void bend_exact_pass(T x[], double L, double Angle)
{
  int i;
  T x1[6];
  T angle1, rho1, rho;
  T xo1, yo1, xp0, yp0, vx0, vy0, vt, vz0, pz0;
  T pos_x0, pos_y0, pos_x1, pos_y1, pos_o1x, pos_o1y, theta;
  T dt1, dt0;
  T vx1, vz1, xp1,yp1,temp1,temp2;
  
  T pz, betaz;
  T delta1,gamma1,beta1;
  delta1= sqrt(1.0 + 2*x[pt_]/GP.beta+x[pt_]*x[pt_]) -1.0;
  gamma1 =GP.gamma+sqrt(GP.gamma*GP.gamma-1)*x[pt_];
  beta1 = sqrt(1.0-1.0/gamma1/gamma1);
  
  if(Angle > 0.) {
    pz0=sqrt( (1+delta1)*(1+delta1)-x[1]*x[1]-x[3]*x[3]);
    xp0=x[1]/pz0;  yp0=x[3]/pz0; 
    pz0=sqrt(xp0*xp0+yp0*yp0+1.0);
    vx0=beta1*xp0/pz0;  vy0=beta1*yp0/pz0; 
    vz0=beta1/pz0;      vt=sqrt(vx0*vx0+vz0*vz0);  
    
    rho=L/Angle;
    rho1=sqrt(xp0*xp0+1)*(1.+delta1)*rho/pz0;
    temp1=sin(Angle);
    temp2=sqrt(1.0- temp1*temp1);
    pos_x0=(rho + x[0])*temp2;
    pos_y0=(rho + x[0])*temp1;
    
    theta=atan(vx0/vz0);
    theta=theta+Angle;
    pos_o1x=pos_x0+rho1*cos(PI+theta);
    pos_o1y=pos_y0+rho1*sin(PI+theta);
    pos_x1=sqrt(rho1*rho1-pos_o1y*pos_o1y)+ pos_o1x;
    pos_y1=0.;
    
    temp1=sqrt((pos_x1-pos_x0)*(pos_x1-pos_x0) +(pos_y1-pos_y0)*(pos_y1-pos_y0));
    angle1=2.*asin(temp1/rho1/2.);
    dt1= angle1*rho1/vt;
    dt0= L/GP.beta;
    
    x[0]=pos_x1-rho;
    x[2]=x[2]+vy0*dt1;
    x[4]=x[4]-(dt1-dt0);
    
    theta=asin(pos_o1y/rho1);
    temp1=sin(theta); temp2 =sqrt(1.-temp1*temp1);
    vx1=-vt*temp1;    vz1=vt*temp2; 
    xp1=-temp1/temp2; yp1=vy0/vz1;
    temp1=(1+delta1)*(1+delta1)/(1. + 1./(xp1*xp1+yp1*yp1));
    x[1]=xp1*sqrt((1+delta1)*(1+delta1)-temp1);
    x[3]=yp1*sqrt((1+delta1)*(1+delta1)-temp1); 
  }
  else if(Angle<0.){
    Angle=-Angle;
    pz0=sqrt( (1+delta1)*(1+delta1)-x[px_]*x[px_]-x[py_]*x[py_]);
    xp0=x[px_]/pz0;  yp0=x[py_]/pz0; 
    pz0=sqrt(xp0*xp0+yp0*yp0+1.0);
    vx0=beta1*xp0/pz0;  vy0=beta1*yp0/pz0; 
    vz0=beta1/pz0;      vt=sqrt(vx0*vx0+vz0*vz0);  
    
    rho=L/Angle;
    rho1=sqrt(xp0*xp0+1)*(1.+delta1)*rho/pz0;
    temp1=sin(Angle);
    temp2=sqrt(1.0- temp1*temp1);
    pos_x0=(rho - x[0])*temp2;
    pos_y0=(rho - x[0])*temp1;
    
    theta=atan(vx0/vz0);
    theta=Angle-theta;
    pos_o1x=pos_x0+rho1*cos(PI+theta);
    pos_o1y=pos_y0+rho1*sin(PI+theta);
    pos_x1=sqrt(rho1*rho1-pos_o1y*pos_o1y)+ pos_o1x;
    pos_y1=0.;
    
    temp1=sqrt((pos_x1-pos_x0)*(pos_x1-pos_x0) + pos_y0*pos_y0);
    angle1=2.*asin(temp1/rho1/2.);
    dt1= angle1*rho1/vt;
    dt0= L/GP.beta;
    
    x[x_]=-pos_x1+rho;
    x[y_]=x[y_]+vy0*dt1;
    x[z_]=x[z_]-(dt1-dt0);
    
    theta=asin(pos_o1y/rho1);
    temp1=sin(theta); temp2 =sqrt(1.-temp1*temp1);
    vx1=vt*temp1;     vz1=vt*temp2; 
    xp1=temp1/temp2;  yp1=vy0/vz1;
    temp1=(1+delta1)*(1+delta1)/(1. + 1./(xp1*xp1+yp1*yp1));
    x[px_]=xp1*sqrt((1+delta1)*(1+delta1)-temp1);
    x[py_]=yp1*sqrt((1+delta1)*(1+delta1)-temp1); 
  }
  else{
    DRIFT_Pass(x,L);
  }
}


template <class T>  void BEND_Pass(T x[], double L, int Nint, double Angle, double E1, double E2) 
{
  int i;
  double href=Angle/L;
  double Lint=L/Nint;
  
  x[1] = x[1]+ tan(E1)*x[0]*href;   
  x[3] = x[3]- tan(E1)*x[2]*href;
  
  if(GP.ultra_rel == true ){
    for(i=0;i<Nint;i++){
      DRIFT_Pass(x,Fdrift1*Lint);
      bend_kick_pass(x, Fkick1*Lint, href);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_kick_pass(x, Fkick2*Lint, href);
      DRIFT_Pass(x,Fdrift2*Lint);
      bend_kick_pass(x, Fkick1*Lint, href);
      DRIFT_Pass(x,Fdrift1*Lint); 
    }
  }
  else{
    bend_exact_pass(x,L,Angle);
  }
  
  x[1] = x[1]+ tan(E2)*x[0]*href;   
  x[3] = x[3]- tan(E2)*x[2]*href;
}


template <class T> void MULT_Pass(T x[], double L, int Nint, int Norder, double KNL[11], double KNSL[11]) 
{
  int i;
  double Lint=L/Nint;
  
  if(L==0.) 
    {
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
	DRIFT_Pass(x,Fdrift1*Lint);
	mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
	DRIFT_Pass(x,Fdrift2*Lint);
	mult_kick_pass(x, Norder,knl_kick2, knsl_kick2);
	DRIFT_Pass(x,Fdrift2*Lint);
	mult_kick_pass(x, Norder,knl_kick1, knsl_kick1);
	DRIFT_Pass(x,Fdrift1*Lint);
      }
    }
} 


template <class T> void GMULT_Pass(T x[], double L, int Nint,  int Norder, double Angle, double E1, double E2, double KNL[11], double KNSL[11])
{
  if(L==0.) 
    {
      mult_kick_pass(x, Norder, KNL, KNSL);
    }
  else 
    {
      int i,j;
      double href=Angle/L;
      double Lint=L/Nint;
      double knl_kick1[11],knsl_kick1[11];
      double knl_kick2[11],knsl_kick2[11];
    
      x[1] = x[1]+ tan(E1)*x[0]*href;   
      x[3] = x[3]- tan(E1)*x[2]*href;  
      
      if(GP.ultra_rel == true){
	
	for(j=0;j<11;j++) {
	  knl_kick1[j] =Fkick1*KNL[j]/Nint;
	  knsl_kick1[j]=Fkick1*KNSL[j]/Nint;
	}
	for(j=0;j<11;j++) {
	  knl_kick2[j] =Fkick2*KNL[j]/Nint;
	  knsl_kick2[j]=Fkick2*KNSL[j]/Nint;
	}

	for(i=0;i<Nint;i++){
	  DRIFT_Pass(x,Fdrift1*Lint);
	  gmult_kick_pass(x, Fkick1*Lint, href, Norder,knl_kick1, knsl_kick1);
	  DRIFT_Pass(x,Fdrift2*Lint);
	  gmult_kick_pass(x, Fkick2*Lint, href, Norder,knl_kick2, knsl_kick2);
	  DRIFT_Pass(x,Fdrift2*Lint);
	  gmult_kick_pass(x, Fkick1*Lint, href, Norder,knl_kick1, knsl_kick1);
	  DRIFT_Pass(x,Fdrift1*Lint);
	}
      }
      else{  
	
	for(j=0;j<11;j++) {
	  knl_kick1[j] = KNL[j]/Nint/4.;
	  knsl_kick1[j]=KNSL[j]/Nint/4.;
	}

	for(i=0;i<Nint;i++){

  	  bend_exact_pass(x, Lint/8, Angle/Nint/8.);
	  mult_kick_pass(x , Norder,knl_kick1, knsl_kick1);
	  bend_exact_pass(x, Lint/4, Angle/Nint/4.);
	  mult_kick_pass(x , Norder,knl_kick1, knsl_kick1);
	  bend_exact_pass(x, Lint/4, Angle/Nint/4.);
	  mult_kick_pass(x , Norder,knl_kick1, knsl_kick1);
	  bend_exact_pass(x, Lint/4, Angle/Nint/4.);
	  mult_kick_pass(x , Norder,knl_kick1, knsl_kick1);
	  bend_exact_pass(x, Lint/8, Angle/Nint/8.);
	}
      }
      
      x[1] = x[1]+ tan(E2)*x[0]*href;   
      x[3] = x[3]- tan(E2)*x[2]*href; 
    }
} 
