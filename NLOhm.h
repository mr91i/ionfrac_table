






#ifndef EHzone_NLOhm_h
#define EHzone_NLOhm_h
double MeanFreeTime_e_func(double n_n,double T_e);
double MeanFreeTime_i_func(double n_n);
double T_e_func(double T_n, double E ,double E_crit);
double T_i_func(double T_n, double E ,double E_crit);
double K_de_func(double r_d,double T_e,double phi_d);
double K_di_func(double r_d,double T_n,double T_i,double phi_d,double v_drift);
double K_di_func_wo_drift(double r_d,double T, double phi_d);
double K_rec_func(double T_e);
double n_e_func(double K_de,double K_di,double K_rec,double n_n,double n_d,double zeta);
double n_i_func(double K_de,double K_di,double K_rec,double n_n,double n_d,double zeta);
double Z_new_func(double n_i,double n_e,double n_d);
double Z_calc_by_NR(double E,double T_e,double T_i,double m_n,double n_n,double T_n,double m_i,double q_i,double n_d,double r_d,double zeta);
double A_func(double K_di,double K_de,double K_rec,double zeta,double n_d,double n_n);
double F_calc(double Z,double T_e,double T_i,double zeta,double n_n,double T_n,double r_d);
double sigma_conductivity_func(double T_n,double zeta,double rho_gas,double f_dg,double E,double r_d);
void n_calc(double T_n,double zeta,double rho_gas,double f_dg,double E,double r_d,double *n_e,double *n_i);
double E_crit_func(double T,double n_n);
double drift_velocity_e_func(double n_n,double T_e,double E);
double MeanFreeTime_e_func(double n_n,double T_e);
void n_calc_woEH(double T,double zeta,double rho_gas,double f_dg, double r_d,double *n_e,double *n_i, double *S_grain);
double Z_calc_by_NR_woE(double T_e,double T_i,double n_n,double T_n,double n_d,double r_d,double zeta);

void calc_sigma_OHP_wograin_woEH(double rho, double T, double n_e, double n_i, double B, double *sig_O, double *sig_H, double *sig_P,  double *sig_O__P);
void calc_sigma_OHP_wgrain_woEH(double rho, double T, double f_dg, double r_d, double n_e, double n_i, double B, double *sig_O, double *sig_H, double *sig_P,  double *sig_O__P);
void Hall_param_ei( double rho , double T, double B,  double *be, double *bi );
double Hall_param_g( double rho, double T , double B , double r_d , double rho_d , double Z  );
#endif

