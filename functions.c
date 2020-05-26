#include "params.h"

int count(int i){
  int var=0;
  while (i>0){
    i/=10;
    var++;
  }//while
  return var;
}//count


double interp_ps_lin(double kval){
  double ans;
  int i;
  if(kval<k[0] && kval>k[N_ps-1]){
    printf("k should lie between k[0] and k[N_ps-1]\n");
    return -1;
  }//if

  for(i=0;i<N_ps;i++){
    if(kval==k[i]){
      ans=Pk[i];
      break;
    }//if
    
    if(kval>k[i] && kval<k[i+1]){
      ans=Pk[i] + ((Pk[i+1]-Pk[i]) *(kval-k[i]) / (k[i+1]-k[i]));
      break;
    }//if
  }//for
  return ans;
}//interp_ps_lin


double sim_21cm_cube(){
  double a, b;
  FILE *F;
  int i;

  sprintf(cmnd, "ps_table/freq_%d.dat", (int)frequency);
  F = fopen(cmnd,"r");

  for(i=0;i<N_ps;i++){
    fscanf(F, "%lE\t%lE", &a, &b);
    k[i]=a;
    Pk[i]=b*2*M_PI*M_PI/(a*a*a); //b was dimensionless PS
  }//for
  fclose(F);

  //************************Declarations for RNG********************
  T=gsl_rng_default;
  rand1=gsl_rng_alloc(T);
  seed=time(NULL)*getpid();
  gsl_rng_set(rand1, seed);
  //****************************************************************

  //***********************Generating 21cm field********************
  for(int i=0;i<XY_pix;i++){
    for(int j=0;j<XY_pix;j++){
      for(int k=0;k<Z_pix;k++){
        double x=i*LEN/(double)XY_pix;
        double y=j*LEN/(double)XY_pix;
        original_21cm_field[k+Z_pix*(j+XY_pix*i)]=0.0;

        for(int n=1;n<=N_gen;n++){
          a_q=gsl_ran_gaussian_ziggurat(rand1, sqrt(interp_ps_lin(2*M_PI*n/LEN)/2));
          b_q=gsl_ran_gaussian_ziggurat(rand1, sqrt(interp_ps_lin(2*M_PI*n/LEN)/2));

          original_21cm_field[k+Z_pix*(j+XY_pix*i)]+=a_q*cos(M_PI*n*sqrt(x*x+y*y)/LEN) + b_q*sin((M_PI*n*sqrt(x*x+y*y)/LEN));  
        }//for n
        original_21cm_field[k+Z_pix*(j+XY_pix*i)]/=sqrt(N_gen);
      }//for k
    }//for j
  }//for i
  //****************************************************************

  gsl_rng_free(rand1);
  return 0.0;
}//sim_21cm_cube


double count_freq_channels(){

  int ch;
  //****Reading the frequencies started
  sprintf(cmnd, "./ska_data/frequency.dat");
  file=fopen(cmnd, "r"); 

  Z_pix=0;
  while(!feof(file)){
    ch = fgetc(file);
    if(ch == '\n'){
      Z_pix++;
    }//if
  }//while
  fclose(file); 

  freq_array=(double *)malloc(sizeof(double)*Z_pix);

  printf("Number of frequency channels: %d\nFrequencies (in MHz):\n", Z_pix);
  file=fopen(cmnd, "r"); 
  for(int i=0;i<Z_pix;i++){
    fscanf(file, "%d", &ch);
    freq_array[i]=ch;
    printf("%f\n", freq_array[i]);
  }//for
  fclose(file);
  //****Reading the frequencies finished
  //Note that synchro, free-free and pointy varies only along the los direction
  
  return 0.0;
}//count_freq_channels


double read_21cm_cube(char *header_1){

  sprintf(cmnd, "ska_data/%s_pre_total_field_%d_%1.1f_hr-image.bin", header_1, XY_pix, tau);
  printf("Reading total_21cm_field from %s\n", cmnd);
  file=fopen(cmnd, "rb");
  fread(total_21cm_field, sizeof(double), XY_pix*XY_pix*Z_pix, file);
  fclose(file);

  return 0.0;
}//read_21cm_cube


double noise_removal(){

  //****Notations are directly taken from https://doi.org/10.1088/1674-4527/9/6/004
  //****Notation table for this code:
  //    Matrix referred in paper                c array       gsl_array                 size of matrix
  //    [y] (log(intensity))                    y             matrix_y                  Z_pix x 1
  //    [x] ((log\nu)^n...)                     x             matrix_x                  Z_pix x Npol
  //    [N] (Covariance matrix)                 N             matrix_N                  Z_pix x Z_pix
  //    [N^{-1}] (inverse [N])                  Ni            matrix_Ni                 Z_pix x Z_pix
  //    [x^t N^{-1}]                            xtNi          matrix_xtNi               Npol x Z_pix
  //    [x^t N^{-1} x]                          xtNix         matrix_xtNix              Npol x Npol
  //    [x^t N^{-1} x]^{-1}                     xtNix_i       matrix_xtNix_i            Npol x Npol
  //    [x^t N^{-1} x]^{-1}.[x^t N^{-1}]        big_matrix    matrix_big_matrix         Npol x Z_pix
  //    [x^t N^{-1} x]^{-1}.[x^t N^{-1}].[y]    a             matrix_a                  Npol x 1
  //    [xa] (fitted log polynomial)            xa            matrix_xa                 Z_pix x 1
  //*********************************************************************************************

  printf("Preparing for the multifrequency fitting\n");

  //****Allocating the memories for the required matrices
  x=(double *)malloc(sizeof(double)*Z_pix*Npol);

  //****Computing the matrix element of x matrix
  for(int i=0;i<Z_pix;i++){
    for(int j=0;j<Npol;j++){
      x[j+Npol*i]=pow(freq_array[i],expo[j]);
    }//for_j
  }//for_i

  //****Creating gsl matrix x
  gsl_matrix_view matrix_x = gsl_matrix_view_array(x,Z_pix,Npol);

  //****Constructing the N matrix
  //****Using the prescription from 2006, Wang, Tegmark, Santos, Knox; setting N matrix to Idenetity
  N=(double *)malloc(sizeof(double)*Z_pix*Z_pix);
  Ni=(double *)malloc(sizeof(double)*Z_pix*Z_pix);

  //****Computing the inverse of matrix N
  //****Warning****//
  //****Matrix N and Ni can be really large, and it may take a large time to compute the inverse of N. So, better set N=Ni if N is idenetity OR if N is not identity, compute the inverse of N in another external code, store the Ni in a bin file and take it as input here.

  //Setting N and Ni to unity
  for(int i=0;i<Z_pix;i++){
    for(int j=0;j<Z_pix;j++){
      (i==j)?(N[j+Z_pix*i]=1.0,Ni[j+Z_pix*i]=1.0):(N[j+Z_pix*i]=0.0,Ni[j+Z_pix*i]=0.0);
    }//for_j
  }//for_i
  
  gsl_matrix_view matrix_N = gsl_matrix_view_array(N,Z_pix,Z_pix);
  gsl_matrix_view matrix_Ni = gsl_matrix_view_array(Ni,Z_pix,Z_pix);

  int sign; //****sign is used twice in this code
  /*
  //****Computing the inverse of N
  gsl_permutation *permute = gsl_permutation_alloc(Z_pix);
  gsl_linalg_LU_decomp(&matrix_N.matrix, permute, &sign);
  gsl_linalg_LU_invert(&matrix_N.matrix, permute, &matrix_Ni.matrix);
  gsl_permutation_free(permute);
  */

  //****Printing matrix N and Ni (N-inverse) for debugging
  /*
  printf("Matrix N:\n");
  for (int i=0;i<Z_pix;i++){
    for (int j=0;j<Z_pix;j++){
      printf(j!=Z_pix-1 ? "%f\t":"%f\n", gsl_matrix_get(&matrix_N.matrix,i,j));
    }//for_j
  }//for_i

  printf("Matrix N-inverse:\n");
  for (int i=0;i<Z_pix;i++){
    for (int j=0;j<Z_pix;j++){
      printf(j!=Z_pix-1 ? "%f\t":"%f\n", gsl_matrix_get(&matrix_Ni.matrix,i,j));
    }//for_j
  }//for_i
  */
  
  free(N); //No use of matrix N now

  //****Creating xtNi by multiplying transverse of x with N inverse
  xtNi=(double *)malloc(sizeof(double)*Npol*Z_pix);
  gsl_matrix_view matrix_xtNi = gsl_matrix_view_array(xtNi,Npol,Z_pix);
  gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,&matrix_x.matrix,&matrix_Ni.matrix,0.0,&matrix_xtNi.matrix); //Note that CblasTrans will multiply the transverse of matrix_x with matrix_Ni (CblasNoTrans: No transpose or hermitian)
  //Refer: GSL manual and https://www.ibm.com/support/knowledgecenter/en/SSFHY8_5.5.0/com.ibm.cluster.essl.v5r5.essl100.doc/am5gr_hsgemm.htm

  free(Ni); //No use of Ni now

  //****Creating xtNix by multiplying xtNi with x
  xtNix=(double *)malloc(sizeof(double)*Npol*Npol);
  gsl_matrix_view matrix_xtNix = gsl_matrix_view_array(xtNix,Npol,Npol);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,&matrix_xtNi.matrix,&matrix_x.matrix,0.0,&matrix_xtNix.matrix);

  //****Taking inverse of xtNix
  xtNix_i=(double *)malloc(sizeof(double)*Npol*Npol);
  gsl_matrix_view matrix_xtNix_i = gsl_matrix_view_array(xtNix_i,Npol,Npol);
  gsl_permutation *permute = gsl_permutation_alloc(Npol);
  gsl_linalg_LU_decomp(&matrix_xtNix.matrix, permute, &sign);
  gsl_linalg_LU_invert(&matrix_xtNix.matrix, permute, &matrix_xtNix_i.matrix);
  gsl_permutation_free(permute);

  //****Printing matrix xtNix and xtNix_i (xtNix-inverse) for debugging
  /*
  printf("Matrix xtNix:\n");
  for (int i=0;i<Npol;i++){
    for (int j=0;j<Npol;j++){
      printf(j!=Npol-1 ? "%f\t":"%f\n", gsl_matrix_get(&matrix_xtNix.matrix,i,j));
    }//for_j
  }//for_i

  printf("Matrix xtNix-inverse:\n");
  for (int i=0;i<Npol;i++){
    for (int j=0;j<Npol;j++){
      printf(j!=Npol-1 ? "%f\t":"%f\n", gsl_matrix_get(&matrix_xtNix_i.matrix,i,j));
    }//for_j
  }//for_i
  */

  free(xtNix); //No use of xtNix now

  //****Creating big_matrix by multiplying xtNix_i with xtNi
  big_matrix=(double *)malloc(sizeof(double)*Npol*Z_pix);
  gsl_matrix_view matrix_big_matrix = gsl_matrix_view_array(big_matrix,Npol,Z_pix);
  gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,&matrix_xtNix_i.matrix,&matrix_xtNi.matrix,0.0,&matrix_big_matrix.matrix);

  free(xtNi); //No use of xtNi now
  free(xtNix_i); //No use of xtNix_i now

  //****Declarations for multifrequency fitting
  recovered_21cm_field=(double *)malloc(sizeof(double)*XY_pix*XY_pix*Z_pix);
  parameter_a_cube=(double *)malloc(sizeof(double)*XY_pix*XY_pix*Npol); //This is to store the fitted parameter 'a' for each los

  y=(double *)malloc(sizeof(double)*Z_pix);
  a=(double *)malloc(sizeof(double)*Npol);
  xa=(double *)malloc(sizeof(double)*Z_pix);

  printf("Started the multifrequency fitting\n");

  gsl_matrix_view matrix_y = gsl_matrix_view_array(y,Z_pix,1);
  gsl_matrix_view matrix_a = gsl_matrix_view_array(a,Npol,1);
  gsl_matrix_view matrix_xa = gsl_matrix_view_array(xa,Z_pix,1);

  for(int i=0;i<XY_pix;i++){
    for(int j=0;j<XY_pix;j++){
      
      for(int k=0;k<Z_pix;k++){
        y[k]=log(total_21cm_field[k+Z_pix*(j+XY_pix*i)]+offset);
      }//for_k
      
      gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,&matrix_big_matrix.matrix,&matrix_y.matrix,0.0,&matrix_a.matrix);
    
      gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,&matrix_x.matrix,&matrix_a.matrix,0.0,&matrix_xa.matrix);

      for(int k=0;k<Npol;k++){
        parameter_a_cube[k+Npol*(j+XY_pix*i)]=a[k];
      }//for_k

      for(int k=0;k<Z_pix;k++){
        recovered_21cm_field[k+Z_pix*(j+XY_pix*i)]=(y[k]-xa[k])*exp(xa[k]);
      }//for_k

    }//for_j
  }//for_i

  printf("Finished the multifrequency fitting\n");

  free(x);
  free(big_matrix);
  free(y);
  free(a);
  free(xa);

  return 0.0;
}//noise_removal