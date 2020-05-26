#include "params.h"
#include "foreg_params.h"


double synchro(double freq){
  double expo;
  expo=-1*(alpha_syn+dalpha_syn*log(freq/freq_ast));
  return a_syn*pow(freq/freq_ast,expo);
}//synchro


double ffree(double freq){
  double expo;
  expo=-1*(alpha_ff+dalpha_ff*log(freq/freq_ast));
  return a_ff*pow(freq/freq_ast,expo);
}//ffree


double pointy(double freq){
  double expo;
  expo=-1*(alpha_ps+dalpha_ps*log(freq/freq_ast));
  return pow(scut_ps,beta_ps)*a_ps*pow(freq/freq_ast,expo);
}//pointy


//*****Standard deviation of the detector noise
double detector_sigma(){
  double bandwd;
  bandwd=(freq_width/(double)Z_pix)*pow(10,6);
  tau*=3600;  //Conversion from hrs to sec
  return t_sys/sqrt(tau*bandwd);
}//detector


double generate_full_signal(){
  T=gsl_rng_default;
  rand2=gsl_rng_alloc(T);
  seed=time(NULL)*getpid();
  gsl_rng_set(rand2, seed);

  double sigma=detector_sigma();
  for(int i=0;i<XY_pix;i++){
    for(int j=0;j<XY_pix;j++){
      for(int k=0;k<Z_pix;k++){
        nu=frequency+k*freq_width/(double)Z_pix;
        freq_array[k]=nu;
        total_21cm_field[k+Z_pix*(j+XY_pix*i)]=0.0;
        detector_noise_field[k+Z_pix*(j+XY_pix*i)]=0.0;
        total_noise_field[k+Z_pix*(j+XY_pix*i)]=0.0;
        detector_noise_field[k+Z_pix*(j+XY_pix*i)]+=gsl_ran_gaussian_ziggurat(rand2,sigma);
        total_noise_field[k+Z_pix*(j+XY_pix*i)]+=(synchro(nu)+ffree(nu)+pointy(nu)+detector_noise_field[k+Z_pix*(j+XY_pix*i)]);
        total_21cm_field[k+Z_pix*(j+XY_pix*i)]+=(original_21cm_field[k+Z_pix*(j+XY_pix*i)]+total_noise_field[k+Z_pix*(j+XY_pix*i)]);
      }//for k
    }//for j
  }//for i

  tau/=3600;
  
  sprintf(cmnd, "results/sim_cube/cube/simulated_synchro_nu_%d_XY_%d_Z_%d_%1.1f_hrs_Npol_%s.bin", (int)frequency, XY_pix, Z_pix, tau, polystr);
  file=fopen(cmnd,"w");
  for(int i=0;i<Z_pix;i++){
    fprintf(file,"%d\t%f\t%f\n",i,freq_array[i],synchro(freq_array[i]));
  }//for i
  fclose(file);
  
  sprintf(cmnd, "results/sim_cube/cube/simulated_ffree_nu_%d_XY_%d_Z_%d_%1.1f_hrs_Npol_%s.bin", (int)frequency, XY_pix, Z_pix, tau, polystr);
  file=fopen(cmnd,"w");
  for(int i=0;i<Z_pix;i++){
    fprintf(file,"%d\t%f\t%f\n",i,freq_array[i],ffree(freq_array[i]));
  }//for i
  fclose(file);

  sprintf(cmnd, "results/sim_cube/cube/simulated_pointy_nu_%d_XY_%d_Z_%d_%1.1f_hrs_Npol_%s.bin", (int)frequency, XY_pix, Z_pix, tau, polystr);
  file=fopen(cmnd,"w");
  for(int i=0;i<Z_pix;i++){
    fprintf(file,"%d\t%f\t%f\n",i,freq_array[i],pointy(freq_array[i]));
  }//for i
  fclose(file);

  gsl_rng_free(rand2);
  return 0.0;
}//generate_full_signal