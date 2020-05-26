#include "params.h"
#include "power_params.h"

//*********************************************************************
//****Going to calculate the distance between two given frequencies (redshift)
//****Comoving distance between redshifts z_1 and z_2 is
//D_c(z_1,z_2)=(c/H_0)*\int_{z_1}^{z_2}\frac{dz'}{E(z')}
//E(z)=\sqrt{\frac{1}{\Omega_m*(1+z)^3+\Omega_r*(1+z)^4+\Omega_\Lambda}}

double func_Dc(double z){
  z+=1;
  return sqrt(1/( o_m*pow(z,3)+o_r*pow(z,4)+o_l ));
}//func_Dc

double compute_box_size(){
  double freq_low, freq_high, z_low, z_high;
  unsigned long long iterate;

  printf("Computing the box sizes\n");

  Z_LEN=0.0;

  if(SimON==0){//ska
    freq_low=freq_array[0]/1e6;
    freq_high=freq_array[Z_pix-1]/1e6;
  }//if
  else{//SimON==1
    freq_low=freq_array[0];
    freq_high=freq_array[Z_pix-1];
  }//else

  z_low=(1420.0/freq_high)-1.0;
  z_high=(1420.0/freq_low)-1.0;
  iterate=(int)((z_high-z_low)/step_size);

  //Integration by trapezoidal rule
  for(int i=1;i<=iterate;i++){
    Z_LEN+=(func_Dc(z_low+step_size*(i-1))+func_Dc(z_low+step_size*i))*step_size/2;
  }//for_i

  Z_LEN*=(3.0e8/100.0e3); //Unit is Mpc h^{-1}
  //c=3\times 10^8, H_0=100h \times 10^3 ms^{-1} Mpc^{-1}

  if(SimON==0){//ska
    XY_LEN=830.0; //Mpc h^{-1}
  }//if
  else{//SimON==1
    XY_LEN=LEN;
  }//else

  printf("Z_LEN: %f Mpc h^{-1}\n",Z_LEN);

  return 0.0;
}
//*********************************************************************


double k_array_bins(){

  //Computing index for k_xy
  //Linear binning, needed to calculate k_perp and k_para
  k_xy=(double *)malloc(sizeof(double)*XY_pix);
  kf_xy=2.0*M_PI/(double)XY_LEN; //Fundamental mode
  for(int i=0;i<XY_pix;i++){
    k_xy[i]=i;
    if(i>XY_pix/2){
      k_xy[i]=i-XY_pix;
    }//if
    k_xy[i]*=kf_xy;
  }//for

  //Computing index for k_z
  //Linear binning, needed to calculate k_perp and k_para
  k_z=(double *)malloc(sizeof(double)*Z_pix);
  kf_z=2.0*M_PI/(double)Z_LEN; //Fundamental mode
  for(int i=0;i<Z_pix;i++){
    k_z[i]=i;
    if(i>Z_pix/2){
      k_z[i]=i-Z_pix;
    }//if
    //k_z[i]*=kf_z;  //Will calculate k later
  }//for

  /*
  //Ska asks for linear binning in logscale for 3d PS. So I am gonna bin the data points while calculating 3d PS itself
  //Assigning weights to each bins for 3d ps
  diagonal_len=sqrt( 2*XY_pix*XY_pix + Z_pix*Z_pix );
  bin_weight_3d=(int *)malloc(sizeof(int)*diagonal_len); 
  for(int i=0;i<diagonal_len;i++){
    bin_weight_3d[i]=0.0;
  }//for
  //Note that maximum iteration below will be sqrt( 2*(XY_pix-1)*(XY_pix-1) + (Z_pix-1)*(Z_pix-1) ) + 1 which is less than diagonal_len
  for(int i=0;i<XY_pix;i++){
    for(int j=0;j<XY_pix;j++){
      for(int k=0;k<Z_pix;k++){
        bin=0.5+sqrt( k_xy[i]*k_xy[i] + k_xy[j]*k_xy[j] + k_z[k]*k_z[k] );
        bin_weight_3d[bin]+=1;
      }//for k
    }//for j
  }//for i
  */

  int diagonal_len, bin, zeros;

  //Assigning weights to each bins for k_perp (linear scale binning)
  diagonal_len=sqrt(2*XY_pix*XY_pix);
  bin_weight_2d=(unsigned long long *)malloc(sizeof(unsigned long long)*diagonal_len);
  for(int i=0;i<diagonal_len;i++){
    bin_weight_2d[i]=0.0;
  }//for

  for(int i=0;i<XY_pix;i++){
    for(int j=0;j<XY_pix;j++){
      bin=0.5+sqrt( (k_xy[i]*k_xy[i] + k_xy[j]*k_xy[j])/(kf_xy*kf_xy) );
      bin_weight_2d[bin]+=1;
    }//for j
  }//for i

  zeros=0;
  for(int i=0;i<diagonal_len;i++){
    if(bin_weight_2d[i]==0)
      zeros+=1;
  }//for

  bins_2d=diagonal_len-zeros;
  

  //Finally computing k_z in Fourier space
  for(int i=0;i<Z_pix;i++){
    k_z[i]*=kf_z;
  }//for

  //****Printing the k_perp and k_z arrays in a file

  //****Computing k_perp
  double *k_perp;
  double k_val;
  k_perp=(double *)malloc(sizeof(double)*bins_2d);

  for(int i=0;i<bins_2d;i++){
    k_perp[i]=0.0;
  }//for i

  bin=0;
  for(int i=0;i<XY_pix;i++){
    for(int j=0;j<XY_pix;j++){
      k_val=sqrt( (k_xy[i]*k_xy[i] + k_xy[j]*k_xy[j])/(kf_xy*kf_xy) );

      if(k_val>0.0){
        bin=0.5+k_val;
        k_perp[bin]+=k_val;
      }//if
    }//for j
  }//for i

  for(int i=0;i<bins_2d;i++){
    k_perp[i]/=(double)bin_weight_2d[i];
  }//for i
  //****Computed k_perp

  //****Now will print k_perp and k_z
  if(SimON==0){
    sprintf(cmnd,"results/ska_cube/ps/%s_k_perp_XY_%d_%1.1f_hrs_Npol_%s_offset_%f.dat", header, XY_pix, tau, polystr, offset);
  }//if
  else{
    sprintf(cmnd, "results/additional/ps/k_perp_nu_%d_XY_%d_Z_%d_%1.1f_hrs_Npol_%s.dat", (int)frequency, XY_pix, Z_pix, tau, polystr);
  }//else

  file=fopen(cmnd,"w");
  for(int i=1;i<bins_2d;i++){
    fprintf(file,"%f\n",k_perp[i]);
  }//for i

  fclose(file);
  free(k_perp);

  if(SimON==0){
    sprintf(cmnd,"results/ska_cube/ps/%s_k_z_XY_%d_%1.1f_hrs_Npol_%s_offset_%f.dat", header, XY_pix, tau, polystr, offset);
  }//if
  else{
    sprintf(cmnd, "results/additional/ps/k_z_nu_%d_XY_%d_Z_%d_%1.1f_hrs_Npol_%s.dat", (int)frequency, XY_pix, Z_pix, tau, polystr);
  }//else

  file=fopen(cmnd,"w");
  for(int i=1;i<(Z_pix/2);i++){
    fprintf(file,"%f\n",k_z[i]);
  }//for i
  
  printf("Computed the k-arrays and bin weights\n");
  fclose(file);

  return 0.0;
}//k_array_bins


double compute_ps(double *data_cube, char *filename){

  char outfname[1000];

  //****To cast the data_cube into a fftw_complex array
  fftw_complex *r2c_data_cube;
  r2c_data_cube=(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*XY_pix*XY_pix*Z_pix);

  //****Setting the data_cube to zero mean delta
  data_cube_average=0.0;

  for(int i=0;i<XY_pix;i++){
    for(int j=0;j<XY_pix;j++){
      for(int k=0;k<Z_pix;k++){
        data_cube_average+=data_cube[k+Z_pix*(j+XY_pix*i)];
      }//for_k
    }//for_j
  }//for_i

  data_cube_average/=((double)(Z_pix*XY_pix*XY_pix));

  for(int i=0;i<XY_pix;i++){
    for(int j=0;j<XY_pix;j++){
      for(int k=0;k<Z_pix;k++){
        data_cube[k+Z_pix*(j+XY_pix*i)]/=data_cube_average;
        data_cube[k+Z_pix*(j+XY_pix*i)]-=1.0;

        //****Casting the data_cube into a fftw_complex array
        r2c_data_cube[k+Z_pix*(j+XY_pix*i)][0]=data_cube[k+Z_pix*(j+XY_pix*i)];
        r2c_data_cube[k+Z_pix*(j+XY_pix*i)][1]=0.0;
      }//for_k
    }//for_j
  }//for_i

  //****Computing the Fourier transform of data_cube
  fftw_complex *fftw_cube;
  fftw_plan plan;

  fftw_cube=(fftw_complex *)fftw_malloc(sizeof(fftw_complex)*XY_pix*XY_pix*Z_pix);

  plan = fftw_plan_dft_3d(XY_pix, XY_pix, Z_pix, r2c_data_cube, fftw_cube, FFTW_FORWARD, FFTW_ESTIMATE);
  printf("Doing FFT\n");
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  fftw_free(r2c_data_cube);

  //Normalization of FFT
  norm=sqrt(XY_pix*XY_pix*Z_pix);
  for(int i=0;i<XY_pix;i++){
    for(int j=0;j<XY_pix;j++){
      for(int k=0;k<Z_pix;k++){
        fftw_cube[k+Z_pix*(j+XY_pix*i)][0]/=(double)norm;
        fftw_cube[k+Z_pix*(j+XY_pix*i)][1]/=(double)norm;
      }//for_k
    }//for_j
  }//for_i

  //*******************************************************
  //****Calculating 3d PS

  //Declaration for binning and calculating 3d ps
  double k_max, k_min, k_low, k_high, k_mag;
  int bins;

  k_max=sqrt(2*k_xy[XY_pix/2]*k_xy[XY_pix/2]+k_z[Z_pix/2]*k_z[Z_pix/2]);
  k_min=( (k_xy[1]<=k_z[1]) ? k_xy[1]:k_z[1]); //min of k_xy and k_z

  //Counting the number of bins using the k_factor
  bins=0; //Will store the number of bins
  k_low=0.0;
  k_high=k_min;
  while(k_high<k_max){
    bins+=1;
    k_low=k_high;
    k_high*=k_factor;
  }//while

  printf("Computing 3D power spectrum\n");

  k_values=(double *)malloc(sizeof(double)*bins);
  ps_3d=(double *)malloc(sizeof(double)*bins);
  bin_weight_3d=(unsigned long long *)malloc(sizeof(unsigned long long)*bins);

  for(int i=0;i<bins;i++){
    k_values[i]=0.0;
    ps_3d[i]=0.0;
    bin_weight_3d[i]=0;
  }//for

  int index;
  double real, imag;
  //Note that due to hermiticity, half of the points in the fft_cube are redundant. So a factor 1/2 in k iteration
  for(int i=0;i<XY_pix;i++){
    for(int j=0;j<XY_pix;j++){
      for(int k=0;k<(Z_pix/2);k++){

        k_mag=sqrt( k_xy[i]*k_xy[i] + k_xy[j]*k_xy[j] + k_z[k]*k_z[k] );
        index=0;
        k_low=0.0;
        k_high=k_min;

        while(k_high<k_max){

          if( (k_mag>=k_low) && (k_mag<k_high)){
            real=fftw_cube[k+(Z_pix/2)*(j+XY_pix*i)][0];
            imag=fftw_cube[k+(Z_pix/2)*(j+XY_pix*i)][1];
            bin_weight_3d[index]+=1;
            k_values[index]+=k_mag;
            ps_3d[index]+= pow(k_mag,3)*(real*real+imag*imag)/(2.0*M_PI*M_PI); //Dimensionless PS
            break;
          }//if

          index+=1;
          k_low=k_high;
          k_high*=k_factor;

        }//while
        
      }//for k
    }//for j
  }//for i

  if(SimON==0){
    sprintf(outfname,"results/ska_cube/ps/3d_ps_%s.dat",filename);
  }//if
  else{
    sprintf(outfname,"results/additional/ps/3d_ps_%s.dat",filename);
  }//else

  file=fopen(outfname, "w");
  //Starting from i=1 to avoid DC mode
  for(int i=1;i<bins;i++){
    fprintf(file,"%d\t%f\t%f\n", i, k_values[i]/(double)bin_weight_3d[i], ps_3d[i]/(double)bin_weight_3d[i]);
    //printf("%lld\n",bin_weight_3d[i]);
  }//for
  fclose(file);

  free(k_values);
  free(ps_3d);
  free(bin_weight_3d);
  
  //*******************************************************
  
  //*******************************************************
  //****Calculating k_perp, k_para PS
  //k_values=(double *)malloc(sizeof(double)*bins_2d);
  printf("Computing 2D power spectra\n");
  ps_2d=(double *)malloc(sizeof(double)*bins_2d*(Z_pix/2));
  ps_aux=(double *)malloc(sizeof(double)*bins_2d);

  for(int i=0;i<bins_2d*(Z_pix/2);i++){
    //k_values[i]=0.0;
    ps_2d[i]=0.0;
  }//for

  index=0;
  for(int k=0;k<(Z_pix/2);k++){

    for(int j=0;j<bins_2d;j++){
      ps_aux[j]=0.0;
    }//for j

    for(int i=0;i<XY_pix;i++){
      for(int j=0;j<XY_pix;j++){

        k_mag=sqrt( (k_xy[i]*k_xy[i] + k_xy[j]*k_xy[j])/(kf_xy*kf_xy) );

        if(k_mag>0.0){
          index=0.5+k_mag;
          real=fftw_cube[k+(Z_pix/2)*(j+XY_pix*i)][0];
          imag=fftw_cube[k+(Z_pix/2)*(j+XY_pix*i)][1];
          ps_aux[index]+=(real*real+imag*imag);//Not dimensionless

        }//if
        
      }//for j
    }//for i

    for(int j=0;j<bins_2d;j++){
      ps_2d[j+bins_2d*k]=ps_aux[j]/(double)bin_weight_2d[j];
    }//for j

  }//for k
  free(ps_aux);

  //Removing the first row and first column of ps_2d to remove the modes that corresponds to k_z=0 and k_perp=0.
  int i2, j2;
  ps_2d_final=(double *)malloc(sizeof(double)*(bins_2d-1)*((Z_pix/2)-1));

  for(int i=0;i<(bins_2d-1)*((Z_pix/2)-1);i++){
    ps_2d_final[i]=0.0;
  }//for i

  for(int i=1;i<(Z_pix/2);i++){
    for(int j=1;j<bins_2d;j++){
      i2=i-1;
      j2=j-1;
      ps_2d_final[j2+(bins_2d-1)*i2] = ps_2d[j+bins_2d*i];
    }//for j
  }//for i

  //Printing ps_2d in file
  if(SimON==0){
    sprintf(outfname,"results/ska_cube/ps/2d_ps_%s.bin",filename);
  }//if
  else{
    sprintf(outfname,"results/additional/ps/2d_ps_%s.bin",filename);
  }//else
  file=fopen(outfname, "wb");
  fwrite(ps_2d_final, sizeof(double), (bins_2d-1)*((Z_pix/2)-1),file);
  fclose(file);

  free(ps_2d);
  free(ps_2d_final);

  

  //*******************************************************


  fftw_free(fftw_cube);

  return 0.0;
}//compute_ps