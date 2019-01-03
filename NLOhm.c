  //
  //  Created by mori shoji on 2013/07/15.
  //
#include "Headers.h"    
void n_calc(double T_n,double zeta,double rho_gas,double f_dg,double E,double r_d,double *n_e,double *n_i){
  double E_crit;//電子の電場加熱が効いてくる電場の値。=sqrt(6me/mn)kT/el 引用。Landau&lifshitz 1993
  double T_e;//電子の運動エネルギーを温度に換算したもの。1/2mv^2=3/2kT
  double T_i;//陽子の運動エネルギーを温度に換算したもの。
  double K_de;//吸着速度係数[cm^3/s]。電子がダストに吸着する速度のようなもの。ダストに単位時間に吸着する電荷量I_a=n_a*q_a*K_da
  double K_di;//吸着速度係数[cm^3/s]。陽イオンがダストに吸着する速度のようなもの。
  double K_rec;//再結合速度定数[cm^3/s]。電子と陽イオンが再結合する速度のようなもの。引用したものを使う。
  double n_n,n_d;
  double Z;//ダストに吸着している電荷数。ダストが正の帯電なら正の整数、負の帯電なら負の整数になる。ダスト一個？全体？空間あたり？
  double mean_free_time_e,mean_free_time_i;//電子・陽イオンの平均自由時間。
  double v_drift_e,v_drift_i;
	//double delta;//ニュートン法の際に、傾きを求めるために使われるZの微小幅。
	//double Z_0;//ニュートン法の初めのZの位置。
	//E=1e-50;
  if(E<1.0e-50)E=1.0e-50;
  n_n=rho_gas/m_n;
  n_d=0.75*m_n*f_dg*n_n/(pi*rho_d*r_d*r_d*r_d);
  E_crit=sqrt(6.0*m_e/m_n)*k_Boltzmann*n_n*sigma_en*T_n/e_charge;
  T_e=T_e_func(T_n,E,E_crit);
  T_i=T_i_func(T_n,E,E_crit);
  
	//////////////Solve by Newton's method////////////////////
  Z=Z_calc_by_NR(E,T_e,T_i,m_n,n_n,T_n,m_i,q_i,n_d,r_d,zeta);
  
  if(isnan(Z)){
	*n_i=-1.0;
	*n_e=-1.0;
	return;
  }
	////////////////E  N  D////////////////////
	//printf("b\n");
	//Z_mem=Z;
	////////////////E  N  D////////////////////
  mean_free_time_i = MeanFreeTime_i_func(n_n);
  mean_free_time_e = MeanFreeTime_e_func(n_n, T_e);
  v_drift_e = drift_velocity_e_func(n_n, T_e, E);///
  v_drift_i = (m_n+m_i)/(m_i*m_n)*q_i*E*mean_free_time_i;//
  K_de = K_de_func(r_d, T_e, e_charge*Z/r_d);
  K_di = K_di_func(r_d,T_n,T_i,e_charge*Z/r_d,v_drift_i);
  K_rec = K_rec_func(T_e);
	// A_value=K_de*K_di*n_d*n_d/zeta/n_n/K_rec;
  *n_i=n_i_func(K_de,K_di,K_rec,n_n,n_d,zeta);
  *n_e=n_e_func(K_de,K_di,K_rec,n_n,n_d,zeta);
	//printf("%5.3e %5.3e %5.3e  ",T_n,*n_e,*n_i);
  
  
	//exit(1);
}
void n_calc_woEH(double T,double zeta,double rho_gas,double f_dg, double r_d,double *n_e,double *n_i, double *S_grain){
  double K_de, K_di, K_rec;
  double n_n, n_d, Z;
  n_n = rho_gas/m_n;
  //n_H2 = rho_gas/(1.425*one_amu);
  //n_n=n_H2;
  n_d = 0.75*rho_gas*f_dg/( pi*rho_d*r_d*r_d*r_d );
	//////////////Solve by Newton's method////////////////////
  Z = Z_calc_by_NR_woE(T,T,n_n,T,n_d,r_d,zeta);
  if( isnan(Z) ){ printf("Z error \n"); exit(0); }
	////////////////E  N  D////////////////////
  K_de  = K_de_func(r_d, T, e_charge*Z/r_d);
  K_di  = K_di_func_wo_drift(r_d, T, e_charge*Z/r_d);
  K_rec = K_rec_func(T);
  *n_i  = n_i_func(K_de,K_di,K_rec,n_n,n_d,zeta);
  *n_e  = n_e_func(K_de,K_di,K_rec,n_n,n_d,zeta);
  *S_grain = K_de*K_di*n_d*n_d/( zeta*n_n*K_rec );
}










double sigma_conductivity_func(double T_n,double zeta,double rho_gas,double f_dg,double E,double r_d){
  double E_crit;//電子の電場加熱が効いてくる電場の値。=sqrt(6me/mn)kT/el 引用。Landau&lifshitz 1993
  double T_e;//電子の運動エネルギーを温度に換算したもの。1/2mv^2=3/2kT
  double T_i;//陽子の運動エネルギーを温度に換算したもの。
  double K_de;//吸着速度係数[cm^3/s]。電子がダストに吸着する速度のようなもの。ダストに単位時間に吸着する電荷量I_a=n_a*q_a*K_da
  double K_di;//吸着速度係数[cm^3/s]。陽イオンがダストに吸着する速度のようなもの。
  double K_rec;//再結合速度定数[cm^3/s]。電子と陽イオンが再結合する速度のようなもの。引用したものを使う。
  double n_n,n_e,n_i,n_d;
  double Z;//ダストに吸着している電荷数。ダストが正の帯電なら正の整数、負の帯電なら負の整数になる。ダスト一個？全体？空間あたり？
  double mean_free_time_e,mean_free_time_i;//電子・陽イオンの平均自由時間。
  double v_drift_e,v_drift_i;
  double J_tot;//電流密度
			   //double delta;//ニュートン法の際に、傾きを求めるために使われるZの微小幅。
			   //double Z_0;//ニュートン法の初めのZの位置。
  if(E<1.0e-50)E=1.0e-50;
  n_n=rho_gas/m_n;
  n_d=0.75*m_n*f_dg*n_n/(pi*rho_d*r_d*r_d*r_d);
  E_crit=sqrt(6.0*m_e/m_n)*k_Boltzmann*n_n*sigma_en*T_n/e_charge;
  T_e=T_e_func(T_n,E,E_crit);
	//T_e = T_n;
  
  T_i=T_i_func(T_n,E,E_crit);
	//if(T_i>10.0*T_n)printf("Ion heating occurs!!T= %5.3e rho =%5.3e \n",T_n,rho_gas);
	//in this paper, we neglect the ion heating.
	//T_i = T_n;
  
	//printf("E=%e Ecrit=%e\n",E,E_crit);
	//printf("Te=%e Ti=%e\n",T_e,T_i);
	//////////////Solve by Newton's method////////////////////
	//Z_0=-0.01;
	//delta=1e-3*Z_0;
	//////////////Solve by Newton's method////////////////////
	//printf("a\n");
  
  Z=Z_calc_by_NR(E,T_e,T_i,m_n,n_n,T_n,m_i,q_i,n_d,r_d,zeta);
  if(isnan(Z))return NAN;
	////////////////E  N  D////////////////////
	//printf("b\n");
	//Z_mem=Z;
	////////////////E  N  D////////////////////
	//printf("E=%e \n",E);
	//printf("Z=%e \n",Z);
  mean_free_time_i = MeanFreeTime_i_func(n_n);
  mean_free_time_e = MeanFreeTime_e_func(n_n, T_e);
  v_drift_e = drift_velocity_e_func(n_n, T_e, E);///
  v_drift_i = (m_n+m_i)/(m_i*m_n)*q_i*E*mean_free_time_i;//
  K_de = K_de_func(r_d, T_e, e_charge*Z/r_d);
  K_di = K_di_func(r_d,T_n,T_i,e_charge*Z/r_d,v_drift_i);
  K_rec = K_rec_func(T_e);
  n_i = n_i_func(K_de,K_di,K_rec,n_n,n_d,zeta);
  n_e = n_e_func(K_de,K_di,K_rec,n_n,n_d,zeta);
	// A_value=K_de*K_di*n_d*n_d/zeta/n_n/K_rec;
	////////////Calculate current density/////////////
  J_tot = -e_charge * n_e * v_drift_e + q_i * n_i * v_drift_i ;
  
	//if(T_n<40.0)printf("T = %8.4e , E/Ecrit = %8.4e :sig_e/sig_i=%8.4e \n x_e/x_i = %8.4e , t_e/t_i = %8.4e\n",T_n,E/E_crit,n_e * v_drift_e/(n_i * v_drift_i) , n_e/n_i , mean_free_time_e/mean_free_time_i );
  
  return J_tot/E;
}
double MeanFreeTime_e_func(double n_n,double T_e){
	//return 3.0/(16.0 * sigma_en * n_n * sqrt( k_Boltzmann * T_e / ( 2.0 * pi * m_e ) ) );
	// return 1.0/(sigma_en * n_n * sqrt(8.0*k_Boltzman*T_e/m_e /pi));//in cgsG //changed 20150105
  return 1.0/(sigma_en * n_n * sqrt(3.0*k_Boltzmann*T_e/m_e ));//in cgsG //changed 20150105
}

double MeanFreeTime_i_func(double n_n){
  return 1.0/(K_in*n_n);//in cgsG
}
double T_e_func(double T_n, double E ,double E_crit){
  double x = E/E_crit;
	// return T_n*(0.5 + sqrt(0.25 + 9.0*pi/64.0 *x*x ));
	// printf("E = %e , Ecr = %e , Te  %e \n",E,E_crit,T_n*(0.5 + sqrt(0.25 + 2.0/3.0 *E*E/E_crit/E_crit)));
  return T_n*(0.5 + sqrt(0.25 + 2.0/3.0 *x*x));
}
double T_i_func(double T_n, double E , double E_crit){
  double x = E/E_crit;
	//double mu_in = m_i*m_n/(m_i + m_n);
	//double kappa_in = 2.0* m_i*m_n/(m_i+m_n)/(m_i+m_n);
	//return T_n;
	//return T_n * ( 1.0 + 4.0*m_e*k_Boltzmann*T_n*sigma_en*sigma_en/(mu_in*kappa_in*m_n*K_in*K_in) *x*x);
  return T_n * ( 1.0 + 7.6e-7*(T_n/100.0)*x*x);
  
	// double T_i = T_n*(1.0 + 2.0*m_e*k_Boltzman*100*pow((m_i+m_n)/m_n,3.0)*pow(sigma_en/K_in/m_i,2.0)*(100.0/100.0)*E*E/E_crit/E_crit);
}
double drift_velocity_e_func(double n_n,double T_e,double E){
	//double mfp = 1.0/n_n/sigma_en;
	//return - 3.0 * sqrt( 3.0 * pi ) * e_charge * mfp * E /( 16.0 * sqrt( m_e * 1.5*k_Boltzmann*T_e ) );
  return (m_e+m_n)/(m_e*m_n)*(-e_charge)*E*MeanFreeTime_e_func(n_n, T_e);
  
}
double K_de_func(double r_d,double T_e,double phi_d){
  double K_de = pi * r_d*r_d * sqrt( 8.0 *k_Boltzmann*T_e / ( pi * m_e ) );
  double Psi_e = -e_charge*phi_d/( k_Boltzmann * T_e );
  if(Psi_e < 0.0){
	K_de *= ( 1.0 - Psi_e );
  }else if(Psi_e >= 0.0){
	K_de *= exp( -Psi_e );
	  //        if(Psi_e<700)K_de*=exp( -Psi_e );
	  //        else         K_de*=expl( -Psi_e );
  }
  if(K_de==0.0){
	  //printf("!!K_de error!! \n psi_e=%3.1e , exp psi_e=%3.1e\n",Psi_e,exp( -Psi_e));
	return NAN;
  }
  if(isnan(K_de)) printf("[K_de_func]: NAN \n");
  return K_de;
}
double K_di_func(double r_d,double T_n,double T_i,double phi_d,double v_drift){
  double K_di=0.0;
  double theta=T_i/T_n;
  double u=v_drift/sqrt(k_Boltzmann*T_n/m_i);
  double Psi=-q_i*phi_d/( k_Boltzmann*T_n );
  double u_thermal = sqrt(8.0*k_Boltzmann*T_n/( pi*m_i ) );
  if(phi_d < 0.0 ){
	
	K_di   =  0.5*sqrt(theta)*exp(-0.5*u*u/theta);
	if(isnan(K_di)) printf("[K_di_func]: NAN 1, phid<0\n");
	K_di  +=  sqrt(pi/8.0)*(theta+2.0*Psi+u*u)/u*erf( u/sqrt(2.0*theta) );
	if(isnan(K_di)) printf("[K_di_func]: NAN 2, phid<0\n");
	K_di *=  pi*r_d*r_d*u_thermal;
	if(isnan(K_di)) printf("[K_di_func]: NAN 3, phid<0\n");
	
  }else if(phi_d >= 0.0 ){
	  //printf("NOT CLEAR Kdi\n");
	  // printf("ph:%8.2e\n",phi_d);
	  //T_i=T_i-m_i*v_drift*v_drift/(3.0*k_Boltzman);
	K_di = 4.0*pi*r_d*r_d*sqrt(k_Boltzmann*T_i/(2.0*pi*m_i))*exp( -q_i*phi_d/k_Boltzmann/T_i );
	if(isnan(K_di)) printf("[K_di_func]: NAN , phid>0\n");
  }
  if(K_di==0.0){
	printf("[K_di_func]: !!error!! \n");
	return NAN;exit(1);
  }
  if(isnan(K_di)) printf("[K_di_func]: NAN \n");
  return K_di;
}

double K_di_func_wo_drift(double r_d,double T, double phi_d){
  double K_di  = pi * r_d*r_d * sqrt( 8.0 *k_Boltzmann* T / ( pi * m_i ) );
  double Psi_i = q_i*phi_d/( k_Boltzmann * T );
  if(Psi_i< 0.0){
	K_di *= ( 1.0 - Psi_i );
  }else if(Psi_i >= 0.0){
	K_di *= exp( -Psi_i );
  }
  if(K_di==0.0){
	  //printf("!!K_de error!! \n psi_e=%3.1e , exp psi_e=%3.1e\n",Psi_e,exp( -Psi_e));
	return NAN;
  }
  if(K_di<0) printf("[K_di_func]: negative \n");
  if(isnan(K_di)) printf("[K_di_func]: NAN \n");
  return K_di;
}



double K_rec_func(double T_e){
  return 2.4e-7*pow(T_e/300.0,-0.69);
	//再結合速度定数[cm^3/s]。引用したものを使う。Ganguli et al. (1988)
}

double n_e_func(double K_de,double K_di,double K_rec,double n_n,double n_d,double zeta){
  double A_i =  K_rec *zeta *n_n /( K_di*K_de*n_d*n_d );
  return zeta * n_n/( K_de*n_d * ( 0.5+sqrt( 0.25 + A_i )) );
}

double n_i_func(double K_de,double K_di,double K_rec,double n_n,double n_d,double zeta){
  double A_i =  K_rec * zeta *n_n /( K_di*K_de*n_d*n_d );
  return zeta * n_n/( K_di*n_d * ( 0.5+sqrt( 0.25 + A_i )) );
}

double Z_new_func(double n_i,double n_e,double n_d){
  return ( n_e-n_i )/n_d;
}

double E_crit_func(double T,double n_n){
	//return 1.0e-9*(T/100.0)*(n_n/1.0e12);// approximated equation
  return 1.078e-9*(T/100.0)*(n_n/1.0e12);// approximated equation
}

double Z_calc_by_NR(double E,double T_e,double T_i,double m_n,double n_n,double T_n,double m_i,double q_i,double n_d,double r_d,double zeta){
  double v_drift_i;
  double K_de,K_di,K_rec;
  double n_e,n_i;
  double F_0,F_d;
  double Z=0.0,Z_start;
  double mean_free_time_i;
  int i=0;
  double Z_0=-2.0*(T_n/100.0);
  double delta=1.0e-3*Z_0;
  int EH_off = 0;
  if(E<1e-20) EH_off = 1;
  Z_start=Z_0;
  mean_free_time_i=MeanFreeTime_i_func(n_n);
  v_drift_i=(m_n+m_i)/(m_i*m_n)*q_i*E*mean_free_time_i;
  K_rec=K_rec_func(T_e);
  i=0;
  while (1){
	i+=1;
	K_de=K_de_func(r_d, T_e, e_charge*Z_0/r_d);
	  ////
	if(isnan(K_de)){printf("Kde error \n"); return NAN;}
	if(i>1000){printf("No convergence \n"); return NAN;}
	  //////
	
	if(EH_off == 0) K_di = K_di_func(r_d, T_n,T_i, e_charge*Z_0/r_d,v_drift_i);
	if(EH_off == 1) K_di = K_di_func_wo_drift(r_d, T_i, e_charge*Z_0/r_d);
	n_e=n_e_func(K_de,K_di,K_rec,n_n,n_d,zeta);
	n_i=n_i_func(K_de,K_di,K_rec,n_n,n_d,zeta);
	F_0=Z_0*n_d+n_i-n_e;
	
	K_de=K_de_func(r_d, T_e, e_charge*(Z_0+delta)/r_d);
	if(EH_off == 0) K_di = K_di_func(r_d, T_n,T_i, e_charge*(Z_0+delta)/r_d,v_drift_i);
	if(EH_off == 1) K_di = K_di_func_wo_drift(r_d, T_i, e_charge*(Z_0+delta)/r_d);
	n_e=n_e_func(K_de,K_di,K_rec,n_n,n_d,zeta);
	n_i=n_i_func(K_de,K_di,K_rec,n_n,n_d,zeta);
	F_d=(Z_0+delta)*n_d+n_i-n_e;
	  //printf("Kde:%Le \n",K_de);
	  //printf("F1:%e DF:%e \n",F_0,F_d-F_0);
	  //printf("Z:%e Z0:%e Zk:%e\n",Z,Z_0,Z_k);
	Z = Z_0 - delta*F_0/(F_d - F_0);
	if(isnan(Z)){printf("Z error in Zcalc 1 \n"); return NAN;}
	if(fabs(Z-Z_0)<1e-8*fabs(Z_0)){break;}
	if(fabs(Z-Z_0)>1e20*fabs(Z_0)){printf("???\n");return NAN; exit(1);}
	  //Z_k=Z_0;
	Z_0=Z;
	if(fabs(Z)<1e-100){Z_start*=1.1;Z=Z_start;Z_0=Z*0.1;}
  }
  return Z;
}

double Z_calc_by_NR_woE(double T_e,double T_i,double n_n,double T_n,double n_d,double r_d,double zeta){
	// K_di is neglecting dust-ion reation by drift since there is no electric field.
	// K_rec in Ganguli 1988 is used
  
  double K_de,K_di,K_rec, n_e,n_i, F_0,F_d;
  double Z = 0.0 , Z_0 = -20.0*(T_n/100.0) , Z_start = Z_0 ;
  double delta = 1.0e-3*Z_0;
  int i = 0;
  
  K_rec = K_rec_func(T_e);
  
  while (1){
	i += 1;
	K_de = K_de_func(r_d, T_e, e_charge*Z_0/r_d);
	  ////
	if(isnan(K_de)){printf("Kde error \n"); return NAN;}
	if(i>1000){printf("No convergence \n"); return NAN;}
	  //////
	
	K_di = K_di_func_wo_drift(r_d, T_i, e_charge*Z_0/r_d);
	n_e  = n_e_func(K_de,K_di,K_rec,n_n,n_d,zeta);
	n_i  = n_i_func(K_de,K_di,K_rec,n_n,n_d,zeta);
	F_0  = Z_0*n_d + n_i - n_e;
	
	K_de = K_de_func(r_d, T_e, e_charge*(Z_0+delta)/r_d);
	K_di = K_di_func_wo_drift(r_d, T_i, e_charge*(Z_0+delta)/r_d);
	n_e  = n_e_func(K_de,K_di,K_rec,n_n,n_d,zeta);
	n_i  = n_i_func(K_de,K_di,K_rec,n_n,n_d,zeta);
	F_d  = (Z_0+delta)*n_d + n_i - n_e;
	
	Z    = Z_0 - delta*F_0/(F_d - F_0);
	if(isnan(Z)){printf("Z error in Zcalc 1 \n"); return NAN;}
	if(fabs(Z-Z_0)<1e-8*fabs(Z_0)){break;}
	if(fabs(Z-Z_0)>1e20*fabs(Z_0)){printf("???\n");return NAN; exit(1);}
	  //Z_k=Z_0;
	Z_0  = Z;
	printf("%d",i);
	if(fabs(Z)<1e-100){Z_start*=1.1; Z=Z_start ; Z_0=Z*0.1;}
  }
  printf("\n --> %g %g\n",Z_start,Z);
  return Z;
}

///////////////////////////////////////////////////////////
// For calc eta_O eta_H eta_A

//void calc_sigma_OHP_wograin_woEH(double rho, double T, double n_e, double n_i, double B, double *sig_O, double *sig_H, double *sig_P){
void calc_sigma_OHP_wograin_woEH(double rho, double T, double n_e, double n_i, double B, double *sig_O, double *sig_H, double *sig_P, double *sig_O__P){
	// e.g.  see  Wardle 2007
  double be, bi, ec_B = e_charge*c_light/B;
  Hall_param_ei( rho , T , B, &be, &bi );
  *sig_O =  n_e*be             + n_i*bi            ;
  //  *sig_H = -n_e   /(1.+be*be)  + n_i   /(1.+bi*bi) ; This is using assumption of Sum n* Z* = 0
  *sig_H = -n_e*be*be/(1.+be*be)  + n_i*bi*bi/(1.+bi*bi) ;
  *sig_P =  n_e*be/(1.+be*be)  + n_i*bi/(1.+bi*bi) ;
  *sig_O__P = n_e*be*be*be/(1.+be*be)  + n_i*bi*bi*bi/(1.+bi*bi) ;

  *sig_O    *= ec_B;
  *sig_H    *= -ec_B;
  *sig_P    *= ec_B;
  *sig_O__P *= ec_B;
#ifdef DEBUG
  printf("be bi = %g %g \n", be, bi);
  printf("sig_OHP = %g %g %g\n", *sig_O, *sig_H, *sig_O__P);
  printf("eta_A = %g \n",(n_e*n_i*be*bi*(be-bi)/(1+be*be)/(1+bi*bi)));
#endif
}

//void calc_sigma_OHP_wgrain_woEH(double rho, double T, double f_dg, double r_d, double n_e, double n_i, double B, double *sig_O, double *sig_H, double *sig_P){
void calc_sigma_OHP_wgrain_woEH(double rho, double T, double f_dg, double r_d, double n_e, double n_i, double B, double *sig_O, double *sig_H,  double *sig_P, double *sig_O__P){
	// e.g.  see  Wardle 2007
	// Here I use n_g intead of n_d, but these are same.
	// mg nd = rho fdg --> nd = rho fdg / (4/3 pi rho_d rd**3  ) 
  double n_g      = 0.75*rho*f_dg /( pi*rho_d*r_d*r_d*r_d );
  double Z        = (n_e  - n_i)/n_g;
  double bg = Hall_param_g( rho, T , B , r_d , rho_d , Z );
  double ec_B = e_charge * c_light / B; 
  calc_sigma_OHP_wograin_woEH( rho, T, n_e, n_i, B, &*sig_O, &*sig_H, &*sig_P, &*sig_O__P);

  *sig_O += n_g * fabs(Z) *bg             * ec_B;
//  *sig_H += n_g * Z           /(1.+bg*bg) * ec_B;
  *sig_H += n_g * Z *bg*bg/(1.+bg*bg) *ec_B ;
  *sig_P += n_g * fabs(Z) *bg /(1.+bg*bg) * ec_B;
  *sig_P += n_g * fabs(Z) *bg*bg*bg       * ec_B;
  
#ifdef DEBUG
  printf("sig_O sig_H sig_P = %g %g %g ,  %g  %g \n",*sig_O,*sig_H,*sig_P, *sig_O * *sig_P , *sig_H* *sig_H + *sig_P* *sig_P);
  printf( "sig_O_g / sig_O_ei = %g \n" , n_g * fabs(Z) *bg             * ec_B/ *sig_O );
  printf( "sig_H_g / sig_H_ei = %g \n" , n_g * Z           /(1.+bg*bg) * ec_B/ *sig_H );
  printf( "sig_P_g / sig_P_ei = %g \n" , n_g * fabs(Z) *bg /(1.+bg*bg) * ec_B/ *sig_P );
  printf("sig_O sig_H sig_P = %g %g %g ,  %g  %g \n",*sig_O,*sig_H,*sig_P, *sig_O * *sig_P , *sig_H* *sig_H + *sig_P* *sig_P);
#endif
}

void Hall_param_ei( double rho , double T, double B,  double *be, double *bi ){
	int mode = 1;

	if(mode==1){ // Thuis mode is for X.Bai
	double K_en, K_in, gamma_e, gamma_i ;
	// caution : momentum transfer rate is not same as the collisonal frequency. You must see (14) in Draine 1983.
	// Here I copy the values from Bai2011b ,(14)--(16)
	K_en = 8.3e-9 * MAX(  1. , sqrt(T/100.) )              ;
 	K_in = 2.0e-9 * sqrt( one_amu / (m_i*m_n)*(m_i+m_n) )  ;
	gamma_e = K_en / ( m_n + m_e ) ;
	gamma_i = K_in / ( m_n + m_i ) ;
	*be = ( e_charge * B )/( m_e * c_light * gamma_e * rho );
	*bi = ( q_i      * B )/( m_i * c_light * gamma_i * rho );
#ifdef DEBUG
	  printf("be bi = %g %g \n", *be, *bi);
#endif
	}

	if(mode==2){  // This mode is for S.Okuzumi
	double nu_en, nu_in , n_n = rho/m_n ;	
	nu_en = 16./3.* n_n *sigma_en * sqrt( k_Boltzmann * T / ( 2*pi*m_e ) );
	nu_in = K_in * n_n;// K_in is 1.6e-9 cm3/s (Nakano & Umebayashi 1986)
	*be = ( e_charge * B )/( m_e * c_light ) * ( m_e + m_n )/( m_n * nu_en );
	*bi = ( q_i      * B )/( m_i * c_light ) * ( m_i + m_n )/( m_n * nu_in );
	}
}
double Hall_param_g( double rho, double T , double B , double r_d , double rho_d , double Z  ){
//  double n_g      = 0.75*rho*f_dg /( pi*rho_d*r_d*r_d*r_d );
  double m_g     = 1.3333333333333333*pi*r_d*r_d*r_d*rho_d;
  double gamma_g = MAX( 1.3e-9*fabs(Z) , 1.6e-7*(r_d*r_d/1e-8)*sqrt(T/100.) )/(m_g+m_n);
  return (  e_charge * fabs(Z) * B   )/( m_g * c_light * gamma_g * rho ); // see Bai2011b
}
