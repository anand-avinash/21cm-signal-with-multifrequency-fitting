#ifndef _power_params_h_
#define _power_params_h_

//****Variables and functions for computing the power spectrum

//Cosmological Parameters
#define o_m (double) 0.3
#define o_r (double) 0.0
#define o_l (double) 0.7
#define h (double) 0.7
//**********************

#define k_factor (double) 1.5 //For linear binning in log scale
//k_factor=1.5 for ska_data
//k_factor=1.3 for sim_data

#define step_size (double) 0.001 //Step size for trapezoidal integration
int iterate;

//Functions declaration:
double func_Dc(double); //Integrand
double compute_box_size(); //To compute the box size in comoving coordinates
double k_array_bins();
double compute_ps(double *, char *);

//Others:
double data_cube_average, norm, kf_z, kf_xy;
double *k_xy, *k_z;
unsigned long long *bin_weight_3d, *bin_weight_2d;
double *k_values, *ps_3d, *ps_2d, *ps_aux, *ps_2d_final;
int bins_2d;

#endif