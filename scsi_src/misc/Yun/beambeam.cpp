//----------------------BEAMBEAM------------------------------
void cerrf( double xx, double yy, double & wx, double & wy )
{
  int n,nc,nu ;
  double x,y,q,h,xl,xh,yh,tx,ty,tn,sx,sy,saux;
  double rx[33],ry[33];
  double cc=1.12837916709551, zero=0, one=1, two=2, half=.5, 
         xlim=5.33, ylim=4.29, fac1=3.2, fac2=23, fac3=21;

  x = abs(xx);
  y = abs(yy);
  if (y < ylim  and  x  < xlim) 
    {
      q  = (one - y / ylim) * sqrt(one - (x/xlim)*(x/xlim) );
      h  = one / (fac1 * q);
      nc = 7 + int(fac2*q);
      xl = pow(h,(1 - nc));
      xh = y + half/h;
      yh = x;
      nu = 10 + int(fac3*q);
      rx[nu+1] = zero;
      ry[nu+1] = zero;
      for ( n=nu; n>=1; n--) {
	tx = xh + n * rx[n+1];
	ty = yh - n * ry[n+1];
	tn = tx*tx + ty*ty;
	rx[n] = half * tx / tn;
	ry[n] = half * ty / tn;
      }
      sx = zero;
      sy = zero;
      for( n=nc;n>=1; n--) {
	saux = sx + xl;
	sx = rx[n] * saux - ry[n] * sy;
	sy = rx[n] * sy + ry[n] * saux;
	xl = h * xl;
      }
      wx = cc * sx;
      wy = cc * sy;
    }
  else 
    {
      xh = y;
      yh = x;
      rx[1] = zero;
      ry[1] = zero;
      for ( n = 9; n>= 1; n--) {
	tx = xh + n * rx[1];
	ty = yh - n * ry[1];
	tn = tx*tx + ty*ty;
	rx[1] = half * tx / tn;
	ry[1] = half * ty / tn;
      }
      wx = cc * rx[1];
      wy = cc * ry[1];
    }
  if(yy < zero) 
    {
      wx =   two * exp(y*y-x*x) * cos(two*x*y) - wx;
      wy = - two * exp(y*y-x*x) * sin(two*x*y) - wy;
      if(xx >  zero) 	wy = -wy;
    }
  else
    {
      if(xx <  zero) wy = -wy;
    }	
}

void BB4D(double x, double y, double gamma, double N, double sigmax, double sigmay, double & Dpx, double & Dpy)
{
  double rp=1.534698e-18;
  double x1, y1, x2, y2, signx, signy;

  signx = ((x>=0.0)?(1):(-1));
  signy = ((y>=0.0)?(1):(-1));

  if (x==0. && y==0. )
    {
      Dpx=0.;
      Dpy=0.;
    }
  else
    {
      if ( abs(sigmax -sigmay)/sigmax < 1.0e-6 )
	{  
	  double r2= x*x+ y*y;
	  double temp1= 2*N*rp/gamma/r2;
	  double temp2= 1-exp(-r2/2/sigmax/sigmax);
	  Dpx=temp1 * x * temp2;
	  Dpy=temp1 * y * temp2;
	}
      else if ( sigmax > sigmay) 
	{
	  double temp1=(2.*N*rp/gamma)*sqrt(PI/2./(sigmax*sigmax-sigmay*sigmay));
	  double temp2=exp(-x*x/2./sigmax/sigmax - y*y/2./sigmay/sigmay);
	  double temp3=sqrt(2.*(sigmax*sigmax-sigmay*sigmay));
	  double w1_real, w1_imag;
	  double w2_real, w2_imag;
	  x1 = fabs(x)/temp3;
	  y1 = fabs(y)/temp3;
	  x2 = fabs(x)*sigmay/sigmax/temp3;
	  y2 = fabs(y)*sigmax/sigmay/temp3;
	  cerrf( x1, y1, w1_real, w1_imag);
	  cerrf( x2, y2, w2_real, w2_imag);
	  Dpy=temp1 *( w1_real-temp2* w2_real)*signy;
	  Dpx=temp1 *( w1_imag-temp2* w2_imag)*signx;  
	}
      else if ( sigmax < sigmay) 
	{
	  double temp1=(2*N*rp/gamma)*sqrt(PI/2/(sigmay*sigmay-sigmax*sigmax));      
	  double temp2=exp(-x*x/2/sigmax/sigmax - y*y/2/sigmay/sigmay);
	  double temp3=sqrt(2*(sigmay*sigmay-sigmax*sigmax));
	  double w1_real, w1_imag;
	  double w2_real, w2_imag;
	  x1 = fabs(x)/temp3;
	  y1 = fabs(y)/temp3;
	  x2 = fabs(x)*sigmay/sigmax/temp3;
	  y2 = fabs(y)*sigmax/sigmay/temp3;
	  cerrf( y1, x1, w1_real, w1_imag);
	  cerrf( y2, x2, w2_real, w2_imag);
	  Dpy=temp1 *( w1_imag-temp2* w2_imag)*signy; 
	  Dpx=temp1 *( w1_real-temp2* w2_real)*signx;
	}
    }
}

void BB6D(double x[6], double gamma, double Np, double sigma_l, int N_slice, 
                       double emitx_rms,  double betax_star, double alfx_star,
	               double emity_rms,  double betay_star, double alfy_star)
 // emitx_rms, emity_rms:  un-normalized rms emittance,  sigma=SQRT[ emitx_rms * betax ]  
{
  int i,j; 
  double rp= 1.534698e-18;
  double x0[6];
  double Np_slice[N_slice], z_star[N_slice];
  double S, gx_star, gy_star, betax, alfx, betay, alfy;
  double sigmax, sigmay, dsigmax2ds, dsigmay2ds, dUdsigmax2, dUdsigmay2; 
  double X, Y,  Dpx, Dpy, temp1, temp2, temp3, temp4;
  
  //----center of each slice of strong bunch, each slice has not same particle population.

  for(i=0;i<N_slice;i++)  Np_slice[i] =  1.0* Np /  N_slice;

  if(N_slice == 11 ) {
    z_star[0]=  1.7996765362153655;
    z_star[1]=  1.1049618149811977;
    z_star[2]=  0.75072324235797372;
    z_star[3]=  0.47407703558819131;
    z_star[4]=  0.23041113693387344;
    z_star[5]=  0.;
    z_star[6]= -0.23041113693387344;
    z_star[7]= -0.47407703558819131; 
    z_star[8]= -0.75072324235797372;
    z_star[9]= -1.1049618149811977;
    z_star[10]=-1.7996765362153655; 
    for(i=0;i<11;i++) z_star[i]=z_star[i] * sigma_l;
  }
  else{
    double y[N_slice-1];
    for(i=0;i<N_slice-1;i++){
      y[i]= gsl_cdf_ugaussian_Pinv(  (i+1)*1.0/N_slice   );
    }
    z_star[0]= N_slice *( exp(-y[0]*y[0]/2 ) ) / sqrt(2.0*3.14159265);
    z_star[N_slice-1]= N_slice *( 0-exp(-y[0]*y[0]/2 ) )/ sqrt(2.0*3.14159265);
    for(i=1;i<N_slice-1;i++)  z_star[i]= N_slice *( exp(-y[i]*y[i]/2) - exp(-y[i-1]*y[i-1]/2 )  ) / sqrt(2.0*3.14159265);
    for(i=0;i<N_slice;i++) z_star[i]=  z_star[i] * sigma_l;
    //for(i=0;i<N_slice;i++) cout<<  z_star[i] /  sigma_l<<endl; 
  }

  //-----calculate the changes in x[6]
  for(i=0; i<N_slice;i++) {

    for (j=0;j<6;j++) x0[j]=x[j];
    S=(x0[4]-z_star[i])/2.0;

    gx_star= (1.0 + alfx_star * alfx_star )/ betax_star;
    gy_star= (1.0 + alfy_star * alfy_star )/ betay_star;
    betax= betax_star - 2* S * alfx_star + S * S *  gx_star; 
    alfx =(alfx_star -gx_star * S );
    betay=  betay_star - 2* S * alfy_star + S * S *  gy_star; 
    alfy =(alfy_star -gy_star * S );

    sigmax  =  sqrt( emitx_rms * betax ) ;
    sigmay  =  sqrt( emity_rms * betay ) ;
    dsigmax2ds =-2.0*alfx* ( emitx_rms );
    dsigmay2ds =-2.0*alfy* ( emity_rms );
    
    //----calculate the kicks for each slice
    X=x0[0] + x0[1]*S;
    Y=x0[2] + x0[3]*S;
    BB4D(X, Y, gamma, Np_slice[i], sigmax, sigmay, Dpx, Dpy);
    temp1 = 1.0/2./(sigmax*sigmax-sigmay*sigmay) ;
    temp2 = X *  Dpx + Y *  Dpy;
    temp3 = 2.* Np_slice[i] *rp /gamma ;
    temp4 = exp ( - X * X / 2. / sigmax /sigmax -  Y * Y / 2. / sigmay /sigmay ) ;
    
    x[0]=x0[0] - S * Dpx * GP.bbscale ;
    x[1]=x0[1] + Dpx  * GP.bbscale;   
    x[2]=x0[2] - S * Dpy * GP.bbscale ;
    x[3]=x0[3] + Dpy * GP.bbscale;   
    x[4]=x0[4];
    
    if ( abs(sigmax- sigmay)/sigmax < 1.0e-6 ) 
      {
      	x[5]=x0[5]+ 0.5 * Dpx * GP.bbscale * ( x0[1] + 0.5* Dpx * GP.bbscale ) 
                  + 0.5*  Dpy * GP.bbscale * ( x0[3] + 0.5* Dpy * GP.bbscale ) 
                  + (1.0/sigmax/sigmax ) * dsigmax2ds  * ( temp3 * GP.bbscale/ 4. ) * temp4;
      }
    else if (sigmax > sigmay ) 
      {
	dUdsigmax2 =  temp1 * GP.bbscale * ( temp2  + temp3 *( sigmay/ sigmax * temp4 -1. ) );
	dUdsigmay2 = -temp1 * GP.bbscale * ( temp2  + temp3 *( sigmax/ sigmay * temp4 -1. ) );
	x[5]= x0[5] + 0.5 * Dpx * GP.bbscale * ( x0[1] + 0.5* Dpx * GP.bbscale) 
                    + 0.5 * Dpy * GP.bbscale * ( x0[3] + 0.5* Dpy * GP.bbscale) 
                    - ( 0.5 * dsigmax2ds * dUdsigmax2 + 0.5 * dsigmay2ds *  dUdsigmay2 ) ;
      }
    else
      {
	dUdsigmay2 = -temp1 *  GP.bbscale * ( temp2  + temp3 *( sigmax/ sigmay * temp4 -1. ) );
	dUdsigmax2 =  temp1 *  GP.bbscale * ( temp2  + temp3 *( sigmay/ sigmax * temp4 -1. ) );
	x[5]= x0[5] + 0.5 * Dpx *GP.bbscale * ( x0[1] + 0.5* Dpx* GP.bbscale) 
                    + 0.5 * Dpy *GP.bbscale *  ( x0[3] + 0.5* Dpy * GP.bbscale) 
                    - ( 0.5 * dsigmax2ds * dUdsigmax2 + 0.5 * dsigmay2ds *  dUdsigmay2 ) ;
      }
  }
}

void BB6D_Angle(double x[6], double gamma, double Np, double sigma_l, double theta, int N_slice, 
                       double emitx_rms,  double betax_star, double alfx_star,
	               double emity_rms,  double betay_star, double alfy_star)
 // emitx_rms, emity_rms:  un-normalized rms emittance,  sigma=SQRT[ emitx_rms * betax ]  
{
  int i; 
  double rp= 1.534698e-18;
  double Np_slice[N_slice], Z0[N_slice], Z0_star[N_slice];
  double S, gx_star, gy_star, betax, alfx, betay, alfy;
  double sigmax, sigmay, dsigmax2ds, dsigmay2ds, dUdsigmax2, dUdsigmay2; 
  double X, Y,  Dpx, Dpy, temp1, temp2, temp3, temp4;
  
  double  x0, y0, z0, px0, py0, pz0, ps0, h0;

  double  x_star, y_star, z_star, px_star, py_star, pz_star;
  double  hx_star, hy_star, hz_star, ps_star, h_star;

  double  x_star_0, y_star_0, z_star_0, px_star_0, py_star_0, pz_star_0;
  double a1, b1, a2, b2;

  //----center of each slice of strong bunch, here each slice has not same particle population.
  for(i=0;i<N_slice;i++)  Np_slice[i] =  1.0* Np /  N_slice;

  if(N_slice == 11 ) {
    Z0[0]=  1.7996765362153655;
    Z0[1]=  1.1049618149811977;
    Z0[2]=  0.75072324235797372;
    Z0[3]=  0.47407703558819131;
    Z0[4]=  0.23041113693387344;
    Z0[5]=  0.;
    Z0[6]= -0.23041113693387344;
    Z0[7]= -0.47407703558819131; 
    Z0[8]= -0.75072324235797372;
    Z0[9]= -1.1049618149811977;
    Z0[10]=-1.7996765362153655; 
    for(i=0;i<11;i++) Z0[i]=Z0[i] * sigma_l;
  }
  else{
    double y[N_slice-1];
    for(i=0;i<N_slice-1;i++){
      y[i]= gsl_cdf_ugaussian_Pinv(  (i+1)*1.0/N_slice   );
    }
    Z0[0]= N_slice *( exp(-y[0]*y[0]/2 ) ) / sqrt(2.0*3.14159265);
    Z0[N_slice-1]= N_slice *( 0-exp(-y[0]*y[0]/2 ) )/ sqrt(2.0*3.14159265);
    for(i=1;i<N_slice-1;i++)  Z0[i]= N_slice *( exp(-y[i]*y[i]/2) - exp(-y[i-1]*y[i-1]/2 )  ) / sqrt(2.0*3.14159265);
    for(i=0;i<N_slice;i++) Z0[i]=  Z0[i] * sigma_l;
    //for(i=0;i<N_slice;i++) cout<<  Z0[i] /  sigma_l<<endl; 
  }

  //----Perform Lorentz booster to the test particle's (x,px,y,py,z,delta)
  x0= x[0];
  px0=x[1];
  y0= x[2];
  py0=x[3];
  z0= x[4];
  pz0=x[5];
  ps0= sqrt( (1+pz0 )*(1+pz0)- px0*px0 - py0*py0 );  
  h0= (1+pz0)- ps0;

  px_star = ( px0 - tan(theta) * h0 )  / cos(theta );
  py_star = py0 / cos(theta);
  pz_star = pz0 - tan(theta) * px0 + tan(theta) * tan(theta) * h0;

  ps_star = sqrt( (1+pz_star )*(1+pz_star)- px_star * px_star - py_star * py_star ); 
  h_star  = (1+pz_star) - ps_star;
  hx_star = px_star / ps_star;  
  hy_star = py_star / ps_star; 
  hz_star =-h_star  / ps_star;

  x_star  = tan(theta)*z0 + ( 1 + hx_star * sin(theta) ) * x0 ;
  y_star  = y0 + sin(theta) * hy_star * x0;
  z_star  = z0/cos(theta) + hz_star * sin(theta) * x0;

  //----Perform Lorentz booster to the center of  sliced strong beam
  for(i=0;i<N_slice;i++)
    Z0_star[i]=Z0[i]/cos(theta); 

  //----Twiss parmeters and emittances in the new frame
  gx_star= ((1.0 + alfx_star * alfx_star )/ betax_star ) ;
  gy_star= ((1.0 + alfy_star * alfy_star )/ betay_star ) ;
  
  gx_star=gx_star/cos(theta);
  gy_star=gy_star/cos(theta);
  betax_star= betax_star *cos(theta);
  betay_star= betay_star *cos(theta);
  emitx_rms = emitx_rms / cos(theta);
  emity_rms = emity_rms / cos(theta);

  gx_star=gx_star;
  gy_star=gy_star;
  betax_star= betax_star;
  betay_star= betay_star;
  emitx_rms =  emitx_rms;
  emity_rms =  emity_rms;

  //-----head-on collision in head-on frame
  for(i=0; i<N_slice;i++) {

    x_star_0 =x_star; 
    px_star_0=px_star;
    y_star_0 =y_star; 
    py_star_0=py_star;
    z_star_0 =z_star; 
    pz_star_0=pz_star;

    S=(z_star_0 - Z0_star[i])/2.0;

    betax= betax_star - 2* S * alfx_star + S * S *  gx_star; 
    alfx =(alfx_star -gx_star * S );
    betay= betay_star - 2* S * alfy_star + S * S *  gy_star; 
    alfy =(alfy_star -gy_star * S );

    sigmax  =  sqrt( emitx_rms * betax ) ;
    sigmay  =  sqrt( emity_rms * betay ) ;
    dsigmax2ds =-2.0*alfx* ( emitx_rms );
    dsigmay2ds =-2.0*alfy* ( emity_rms );

    X= x_star_0 + px_star_0 * S - Z0_star[i] * sin(theta);
    Y= y_star_0 + py_star_0 * S - 0 ;
    BB4D(X, Y, gamma, Np_slice[i], sigmax, sigmay, Dpx, Dpy);
    temp1 = 1.0/2./(sigmax*sigmax-sigmay*sigmay) ;
    temp2 = X *  Dpx + Y *  Dpy;
    temp3 = 2.* Np_slice[i] *rp /gamma ;
    temp4 = exp ( - X * X / 2. / sigmax /sigmax -  Y * Y / 2. / sigmay /sigmay ) ;
    
    x_star =x_star_0  -  S * Dpx * GP.bbscale ;
    px_star=px_star_0 +  Dpx * GP.bbscale;   
    y_star =y_star_0  -  S * Dpy * GP.bbscale ;
    py_star=py_star_0 +  Dpy * GP.bbscale;   
    z_star =z_star_0;

    if ( abs(sigmax- sigmay)/sigmax < 1.0e-6 ) 
      {
      	pz_star=pz_star_0 + 0.5 * Dpx * GP.bbscale * ( px_star_0 + 0.5* Dpx * GP.bbscale ) 
                          + 0.5*  Dpy * GP.bbscale * ( py_star_0 + 0.5* Dpy * GP.bbscale ) 
                          + (1.0/sigmax/sigmax ) * dsigmax2ds  * ( temp3 * GP.bbscale/ 4. ) * temp4;
      }
    else if (sigmax > sigmay ) 
      {
	dUdsigmax2 =  temp1 * GP.bbscale * ( temp2  + temp3 *( sigmay/ sigmax * temp4 -1. ) );
	dUdsigmay2 = -temp1 * GP.bbscale * ( temp2  + temp3 *( sigmax/ sigmay * temp4 -1. ) );
	pz_star= pz_star_0 + 0.5 * Dpx * GP.bbscale * ( px_star_0 + 0.5* Dpx * GP.bbscale) 
                           + 0.5 * Dpy *GP.bbscale  * ( py_star_0 + 0.5* Dpy *  GP.bbscale) 
                           - ( 0.5 * dsigmax2ds * dUdsigmax2 + 0.5 * dsigmay2ds *  dUdsigmay2 ) ;
      }
    else
      {
	dUdsigmay2 = -temp1 *  GP.bbscale * ( temp2  + temp3 *( sigmax/ sigmay * temp4 -1. ) );
	dUdsigmax2 =  temp1 *  GP.bbscale * ( temp2  + temp3 *( sigmay/ sigmax * temp4 -1. ) );
	pz_star    =  pz_star_0 + 0.5 * Dpx *GP.bbscale * ( px_star_0 + 0.5* Dpx* GP.bbscale) 
                                + 0.5 * Dpy *GP.bbscale * ( py_star_0 + 0.5* Dpy * GP.bbscale) 
                                - ( 0.5 * dsigmax2ds * dUdsigmax2 + 0.5 * dsigmay2ds *  dUdsigmay2 ) ;
      }
   }

  //----perform inverse Lorentz booster to the test particle's coordinates
  ps_star = sqrt( (1+pz_star )*(1+pz_star)- px_star * px_star - py_star * py_star );
  h_star  = 1+ pz_star -  ps_star;
  hx_star = px_star / ps_star ;
  hy_star = py_star / ps_star ;
  hz_star = pz_star / ps_star ;

  px0 = (px_star + h_star * sin(theta)) * cos(theta)  ;
  py0 = py_star * cos(theta);
  pz0 = px_star * sin(theta) + pz_star;

  a1=1+hx_star*sin(theta);
  b1=tan(theta);
  a2=hz_star*sin(theta);
  b2=1/cos(theta);
  LinearEquations( a1, b1, x_star, a2, b2, z_star, x0, z0);
  y0=y_star - sin(theta)*hy_star * x0;
  
  x[0]=x0 ;
  x[1]=px0;
  x[2]=y0 ;
  x[3]=py0;
  x[4]=z0 ;
  x[5]=pz0;
}

void BEAMBEAM_Pass(double x[6], int TREATMENT, double NP, double SIGMAL, int NSLICE, double EMITX,  double BETAX, double ALFAX, double EMITY, double BETAY, double ALFAY)
{
  if ( NP != 0. ) {
      if(int(TREATMENT) == 6){
	BB6D(x, GP.gamma, NP, SIGMAL, NSLICE, EMITX,  BETAX, ALFAX, EMITY, BETAY, ALFAY); 
      }
      else{
	double Dpx, Dpy;
	BB4D(x[0], x[2], GP.gamma, NP, sqrt( EMITX * BETAX), sqrt( EMITY * BETAY ), Dpx,  Dpy);
	x[1]=x[1]+Dpx * GP.bbscale ;
	x[3]=x[3]+Dpy * GP.bbscale;
      }
  } 
}
