#include "Headers.h"
  // See Bai's papers in 2011 and 2013:
  // Title : THE ROLE OF TINY GRAINS ON THE ACCRETION PROCESS IN PROTOPLANETARY DISKS
  // Title : WIND-DRIVEN ACCRETION IN PROTOPLANETARY DISKS. I
#define SQR(x) ((x)*(x))
#define WITH_GRAIN
//#define DEBUG
int main(int argc, char *argv[]){
  ////////////////////////////////////////////////////
  double zeta      = atof(argv[1]);
  double rho       = atof(argv[2]);
  double T         = atof(argv[3]);
  double f_dg      = atof(argv[4]);
  double r_d       = atof(argv[5]);
  ////////////////////////////////////////////////////
  double n_n       = rho/m_n ;
  double n_e, n_i, S_grain, sig_O, sig_H, sig_P, sig_O__P;
  // Assume small B strength for
  double B = 1e-3 * 1./2.1 * (n_n/1.e15)/MAX(1,sqrt(T/100.0));
  double cc_4p = c_light*c_light/(4.*pi);  

  n_calc_woEH(T, zeta, rho, f_dg, r_d, &n_e, &n_i, &S_grain);
  
#ifdef WITH_GRAIN
	calc_sigma_OHP_wgrain_woEH(rho, T, f_dg, r_d, n_e, n_i, B, &sig_O, &sig_H, &sig_P, &sig_O__P);
#elif
	calc_sigma_OHP_wograin_woEH(rho, T, n_e, n_i, B, &sig_O, &sig_H, &sig_P, &sig_O__P);
#endif

  double sig_perp2 = SQR(sig_H) + SQR(sig_P) ;
  double eta_O     = cc_4p / sig_O;
  double eta_H     = cc_4p * sig_H/sig_perp2 ;
  double eta_A     = cc_4p * ( sig_P*sig_O__P - sig_H*sig_H ) /( sig_perp2*sig_O );
  double Q_H = eta_H/B   ;
  double Q_A = eta_A/B/B ;
  
#ifdef DEBUG
  printf("B = %g , n = %g \n",B,rho/m_n);
  printf("n_e n_i Znd = %g %g %g \n", n_e , n_i , n_e - n_i);
  printf("x_e x_i x_d = %g %g %g \n", n_e / n_n, n_i/n_n , (n_e - n_i)/n_n);
#endif
  printf("%e %e 0.0 %e %e 0.0 %e 0.0\n", zeta , rho , eta_O , Q_H , Q_A );

  return 1;
}
