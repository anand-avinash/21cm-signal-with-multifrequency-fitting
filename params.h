#ifndef _params_h_
#define _params_h_

#include <stdio.h>
#include <time.h>
#include <unistd.h> //For getpid() function
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>


//****Simulation Parameters:
#define N_ps (int) 64 //Number of points in PS file
#define N_gen (int) 4 //To generate the 21cm field
#define LEN (double) 8.0 //Length of simulation box (in Mpc h^{-1})
#define freq_width (double) 2.0 //Width of frequency channel (in MHz)

//****Foreground Removal Parameters:
double *expo;
int expo_array;
int count(int);

//****Function declarations:
double interp_ps_lin(double);
double read_21cm_cube(char *);
double sim_21cm_cube();
double generate_full_signal();
double count_freq_channels();
double noise_removal();

//****Multifrequency fitting declarations:
double *y, *x, *N, *Ni, *xtNi, *xtNix, *xtNix_i, *big_matrix, *a, *xa;

//****GSL rng declarations:
const gsl_rng_type *T;
gsl_rng *rand1, *rand2;
long seed;
double a_q, b_q;

//****Others:
int SimON, XY_pix, Z_pix, Npol; //XY_pix and Z_pix are the box dimensions in pixels
double frequency, original_21cm_mean, total_21cm_mean;
double *original_21cm_field, *total_21cm_field, *detector_noise_field, *total_noise_field , *recovered_21cm_field, *parameter_a_cube;
double k[N_ps], Pk[N_ps];
double *freq_array;
double sigma, nu, tau, offset;
double XY_LEN, Z_LEN; //Length of the box

FILE *file;
char cmnd[1000], header[100], polystr[100];
//****polystr is to store the exponents of variables of the fitting polynomial

#endif