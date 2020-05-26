#ifndef _foreground_params_h_
#define _foreground_params_h_


//This header contains parameters for foreground to 21cm signal

#define freq_ast (double) 150 //in MHz

//********************Synchrotron
#define a_syn (double) 335.4 //in K
#define alpha_syn (double) 2.8
#define dalpha_syn (double) 0.16
//*******************************


//**********************free-free
#define a_ff (double) 33.5 //in K
#define alpha_ff (double) 2.15
#define dalpha_ff (double) 0.01
//*******************************

//******************Point Sources
#define a_ps (double) 90.45 //in K (Contributes about 27%)
#define alpha_ps (double) 2.81
#define dalpha_ps (double) 0.25
#define beta_ps (double) 0.125
#define scut_ps (double) 0.1 //in mJy
//*******************************


//*****************Detector Noise
#define t_sys (double) 200 //in K
//#define tau (double) 3600 //1hr=3600sec
//*******************************

double synchro(double);
double ffree(double);
double pointy(double);
double detector_sigma();

#endif