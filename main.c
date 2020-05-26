#include "params.h"
#include "foreg_params.h"
#include "power_params.h"


int main(int argc, char **argv){
  
  if(argc==6){

    //SimON=1 (To simulate the field)
    frequency=atof(argv[1]); //Observation frequency
    XY_pix=atoi(argv[2]);
    Z_pix=atoi(argv[3]);
    tau=atof(argv[4]); //Integration time
    expo_array=atoi(argv[5]);
    offset=0.0;
    SimON=1;
    freq_array=(double *)malloc(sizeof(double)*Z_pix);

  }//if
  else if(argc==7){

    //SimON=0 (To read the field from SKA data)
    count_freq_channels();
    strcpy(header,argv[1]);
    XY_pix=atoi(argv[2]);
    tau=atof(argv[3]); //Integration time
    offset=atof(argv[5]);
    expo_array=atoi(argv[6]);
    SimON=0;

  }//else if
  else{

    printf("Usage:\n./execute [frequency (in MHz)] [XY_pix] [Z_pix] [Int_time (hrs)] [expo_array] (To simulate the fields)\nOR\n"); //argc=6 case
    //**Available frequency is given in folder ps_table

    printf("./execute [header] [XY_pix] [Int_time (hrs)] -o [offset] [expo_array] (To read SKA challenge data; Z_pix=21 here)\n"); //argc=7 case
    //**Available headers are: eor_p_noise, eor_p_fg, eor_p_noise_p_fg, noise (Refer to the folder ska_data)
    //**Z_pix is the number of data points in ska_data/frequency.dat
    //**XY_pix=300 (Refer to the variable N in SKA_Blind_Challenge_Pilot/images/code_1.py)
    //**Available choices for Int_time are 10.0, 100.0, 1000.0 hrs
    return -1;

  }//else


  //************************************************
  //****To create an array with the exponents of the terms in the fitting polynomial
  Npol=count(expo_array); //Npol is the number of terms in the fitting polynomial
  //The function count(int) is defined in functions.c

  expo=(double *)malloc(sizeof(double)*Npol); //expo stores the exponent of the terms of the fitting polynomial
  sprintf(polystr,"%d",expo_array);

  //Storing the exponent values in the expo array
  int m=0;
  while(expo_array>0){
    expo[m]=(double)(expo_array%10);
    expo_array/=10;
    m++;
  }//while

  printf("Terms in the fitting polynomial have following exponents:\n");
  for(int i=0;i<Npol;i++){
    printf("%1.1f\t",expo[i]);
  }//for
  printf("\n");
  //************************************************

  

  total_21cm_field=(double *)malloc(sizeof(double)*XY_pix*XY_pix*Z_pix);


  //****************************************************************
  //****Reading/Simulating the 21cm field
  //****For SimON=0, the program will read the SKA data challenge's data
  //****For SimON=1, the program will simulate the 21cm field from using the available PS.


  if(SimON==0){

    printf("Reading 21cm field.\n");
    read_21cm_cube(header);
    printf("Reading the 21cm field is completed.\n");

  }//if
  else if(SimON==1){

    //****Started simulating 21cm field
    printf("Simulating 21cm field.\n");
    original_21cm_field=(double *)malloc(sizeof(double)*XY_pix*XY_pix*Z_pix);
    sim_21cm_cube();
    //****Simulating the 21cm field finished

    sprintf(cmnd, "results/sim_cube/cube/simulated_original_21cm_field_nu_%d_XY_%d_Z_%d_%1.1f_hrs_Npol_%s.bin", (int)frequency, XY_pix, Z_pix, tau, polystr);
    printf("Writing simulated original_21cm_field to the file %s\n", cmnd);
    file=fopen(cmnd, "wb");
    fwrite(original_21cm_field, sizeof(double), XY_pix*XY_pix*Z_pix, file);
    fclose(file);

    detector_noise_field=(double *)malloc(sizeof(double)*XY_pix*XY_pix*Z_pix);
    total_noise_field=(double *)malloc(sizeof(double)*XY_pix*XY_pix*Z_pix);

    printf("Started generating the complete signal.\n");
    generate_full_signal();
    free(detector_noise_field);

    sprintf(cmnd, "results/sim_cube/cube/simulated_total_21cm_field_nu_%d_XY_%d_Z_%d_%1.1f_hrs_Npol_%s.bin", (int)frequency, XY_pix, Z_pix, tau, polystr);
    printf("Writing simulated total_21cm_field to the file %s\n", cmnd);
    file=fopen(cmnd, "wb");
    fwrite(total_21cm_field, sizeof(double), XY_pix*XY_pix*Z_pix, file);
    fclose(file);

    sprintf(cmnd, "results/sim_cube/cube/simulated_total_noise_field_nu_%d_XY_%d_Z_%d_%1.1f_hrs_Npol_%s.bin", (int)frequency, XY_pix, Z_pix, tau, polystr);
    printf("Writing simulated total_noise_field to the file %s\n", cmnd);
    file=fopen(cmnd, "wb");
    fwrite(total_noise_field, sizeof(double), XY_pix*XY_pix*Z_pix, file);
    fclose(file);

    //****Calculating the mean of 21cm field
    original_21cm_mean=0.0;
    for(int i=0;i<XY_pix;i++){
      for(int j=0;j<XY_pix;j++){
        for(int k=0;k<Z_pix;k++){
          original_21cm_mean+=original_21cm_field[k+Z_pix*(j+XY_pix*i)];
        }//k
      }//j
    }//i
    original_21cm_mean/=(XY_pix*XY_pix*Z_pix);
    printf("Mean of original_21cm_field = %f\n", original_21cm_mean);

  }//else if
  else{
    printf("Invalid value of SimON was provided. Code Terminated.\n");
    return -1;
  }//else

  //****Reading/Simulating the 21cm field finished
  //****************************************************************

  //****************************************************************
  //****Noise removal started
  noise_removal(); //Noise removal
  //****Noise removal finished
  //****************************************************************

  //****Writing parameter_a_cube to file
  if(SimON==0){
    sprintf(cmnd, "results/ska_cube/cube/%s_parameter_a_cube_XY_%d_%1.1f_hrs_Npol_%s_offset_%f.bin", header, XY_pix, tau, polystr, offset);
  }//if
  else{
    sprintf(cmnd, "results/sim_cube/cube/parameter_a_cube_nu_%d_XY_%d_Z_%d_%1.1f_hrs_Npol_%s.bin", (int)frequency, XY_pix, Z_pix, tau, polystr);
  }//else

  printf("Writing parameter_a_cube to the file %s\n", cmnd);
  file=fopen(cmnd, "wb");
  fwrite(parameter_a_cube, sizeof(double), XY_pix*XY_pix*Npol, file);
  fclose(file);

  free(parameter_a_cube);

  //****Writing recovered 21cm field to file
  if(SimON==0){
    sprintf(cmnd, "results/ska_cube/cube/%s_recovered_21cm_signal_XY_%d_%1.1f_hrs_Npol_%s_offset_%f.bin", header, XY_pix, tau, polystr, offset);
  }//if
  else{
    sprintf(cmnd, "results/sim_cube/cube/recovered_21cm_signal_nu_%d_XY_%d_Z_%d_%1.1f_hrs_Npol_%s.bin", (int)frequency, XY_pix, Z_pix, tau, polystr);

    //****Setting the mean of the recovered 21cm field to the mean of simulated one
    //****For debugging
    /*
    double recovered_21cm_mean=0.0;
    double offset_mean;
    for(int i=0;i<XY_pix;i++){
      for(int j=0;j<XY_pix;j++){
        for(int k=0;k<Z_pix;k++){
          recovered_21cm_mean+=recovered_21cm_field[k+Z_pix*(j+XY_pix*i)];
        }//k
      }//j
    }//i
    offset_mean=recovered_21cm_mean-original_21cm_mean;
    for(int i=0;i<XY_pix;i++){
      for(int j=0;j<XY_pix;j++){
        for(int k=0;k<Z_pix;k++){
          recovered_21cm_field[k+Z_pix*(j+XY_pix*i)]-=offset_mean;
        }//k
      }//j
    }//i
    */
  }//else
  
  printf("Writing recovered 21cm field to the file %s\n", cmnd);
  file=fopen(cmnd, "wb");
  fwrite(recovered_21cm_field, sizeof(double), XY_pix*XY_pix*Z_pix, file);
  fclose(file);

  //****Preparing to compute the power spectra
  compute_box_size();

  //****To compute the k_x, k_y, k_z (Fourier space coordinates) and to calculate the weights of bins in linear scale
  k_array_bins();


  //****To compute the PS of various fields
  if(SimON==0){
    sprintf(cmnd,"%s_total_21cm_field_XY_%d_%1.1f_hrs_Npol_%s_offset_%f", header, XY_pix, tau, polystr, offset);
    compute_ps(total_21cm_field,cmnd);
    printf("PS calc done for total field\n");

    sprintf(cmnd,"%s_recovered_21cm_field_XY_%d_%1.1f_hrs_Npol_%s_offset_%f", header, XY_pix, tau, polystr, offset);
    compute_ps(recovered_21cm_field,cmnd);
    printf("PS calc done for recovered field\n");
  }//if
  else{
    sprintf(cmnd,"original_21cm_field_nu_%d_XY_%d_Z_%d_%1.1f_hrs_Npol_%s", (int)frequency, XY_pix, Z_pix, tau, polystr);
    compute_ps(original_21cm_field,cmnd);
    printf("PS calc done for original_21cm_field\n");

    //**********************Additional***
    //***********************************
    double *new_array;
    new_array=(double *)malloc(sizeof(double)*XY_pix*XY_pix*Z_pix);
    //double sigma=detector_sigma();
    for(int i=0;i<XY_pix;i++){
      for(int j=0;j<XY_pix;j++){
        for(int k=0;k<Z_pix;k++){
          nu=freq_array[k];
          new_array[k+Z_pix*(j+XY_pix*i)]=0.0;
          new_array[k+Z_pix*(j+XY_pix*i)]+=(original_21cm_field[k+Z_pix*(j+XY_pix*i)]+ffree(nu));
        }//for k
      }//for j
    }//for i

    sprintf(cmnd,"original_p_ff_nu_%d_XY_%d_Z_%d_%1.1f_hrs_Npol_%s", (int)frequency, XY_pix, Z_pix, tau, polystr);
    compute_ps(new_array,cmnd);
    printf("PS calc done for original+ffree\n");


    for(int i=0;i<XY_pix;i++){
      for(int j=0;j<XY_pix;j++){
        for(int k=0;k<Z_pix;k++){
          nu=freq_array[k];
          new_array[k+Z_pix*(j+XY_pix*i)]=0.0;
          new_array[k+Z_pix*(j+XY_pix*i)]+=(original_21cm_field[k+Z_pix*(j+XY_pix*i)]+ffree(nu)+pointy(nu));
        }//for k
      }//for j
    }//for i

    sprintf(cmnd,"original_p_ff_p_point_nu_%d_XY_%d_Z_%d_%1.1f_hrs_Npol_%s", (int)frequency, XY_pix, Z_pix, tau, polystr);
    compute_ps(new_array,cmnd);
    printf("PS calc done for original+ffree+pointy\n");

    for(int i=0;i<XY_pix;i++){
      for(int j=0;j<XY_pix;j++){
        for(int k=0;k<Z_pix;k++){
          nu=freq_array[k];
          new_array[k+Z_pix*(j+XY_pix*i)]=0.0;
          new_array[k+Z_pix*(j+XY_pix*i)]+=(original_21cm_field[k+Z_pix*(j+XY_pix*i)]+synchro(nu)+ffree(nu)+pointy(nu));
        }//for k
      }//for j
    }//for i

    sprintf(cmnd,"original_p_fg_nu_%d_XY_%d_Z_%d_%1.1f_hrs_Npol_%s", (int)frequency, XY_pix, Z_pix, tau, polystr);
    compute_ps(new_array,cmnd);
    printf("PS calc done for original+foreground\n");
    //***********************Additional part over*****
    //************************************************

    

    sprintf(cmnd,"total_noise_field_nu_%d_XY_%d_Z_%d_%1.1f_hrs_Npol_%s", (int)frequency, XY_pix, Z_pix, tau, polystr);
    compute_ps(total_noise_field,cmnd);
    printf("PS calc done for total_noise_field\n");

    free(original_21cm_field);
    free(total_noise_field);

    sprintf(cmnd,"total_21cm_field_nu_%d_XY_%d_Z_%d_%1.1f_hrs_Npol_%s", (int)frequency, XY_pix, Z_pix, tau, polystr);
    compute_ps(total_21cm_field,cmnd);
    printf("PS calc done for total field\n");

    sprintf(cmnd,"recovered_21cm_field_nu_%d_XY_%d_Z_%d_%1.1f_hrs_Npol_%s", (int)frequency, XY_pix, Z_pix, tau, polystr);
    compute_ps(recovered_21cm_field,cmnd);
    printf("PS calc done for recovered field\n");

    free(new_array);
  }//else

  free(k_xy);
  free(k_z);
  free(bin_weight_2d);

  free(freq_array);
  free(total_21cm_field);
  free(recovered_21cm_field);

  free(expo);

  return 0;
}//main