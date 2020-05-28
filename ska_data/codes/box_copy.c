//Joins all the dat created by code_1.py and copy it to a bin file.

#include <stdio.h>
#include <stdlib.h>

int N_pix, ch;
double int_time, var;
int *frequency;
double *total_signal;
char *fname1;
FILE *F;

int main(int argc,char **argv){

  if(argc != 4){
    printf("Usage: ./execute [header] [N_pix] [int_time]\n");
    //N_pix=300
    return -1;
  }//if

  int_time=atof(argv[3]);
  N_pix=atoi(argv[2]);

  //****Reading the frequencies started
  fname1=(char *)malloc(sizeof(char)*400);
  sprintf(fname1,"../frequency.dat");
  F=fopen(fname1, "r"); 

  int count=0;
  while(!feof(F)){
    ch = fgetc(F);
    if(ch == '\n'){
      count++;
    }//if
  }//while
  fclose(F); 
  frequency=(int *)malloc(sizeof(int)*count);

  printf("Number of frequency channels: %d\nFrequencies (in MHz):\n", count);
  F=fopen(fname1, "r"); 
  for(int i=0;i<count;i++){
    fscanf(F, "%d", &ch);
    frequency[i]=ch;
    printf("%d\n", frequency[i]);
  }//for
  fclose(F);
  //****Reading the frequencies finished

  //****Reading the datacube started
  total_signal=(double *)malloc(sizeof(double)*count*N_pix*N_pix);
  for(int k=0; k<count; k++){
    sprintf(fname1,"../%s_%d_Hz_t_obsv_%1.1f_hr-image.fits.txt",argv[1],frequency[k],int_time);
    printf("Reading: %s\n",fname1);
    F=fopen(fname1, "r");
    if(F==NULL){
      perror("Not able to open file\n");
      return -1;
    }//if

    for(int i=0;i<N_pix;i++){
      for(int j=0;j<N_pix;j++){
        fscanf(F, "%lE", &var);
        total_signal[k+count*(j+N_pix*i)]=var;
      }//for_j
    }//for_i
    fclose(F);
  }//for_k
  //****Reading the datacube finished
  

  /* //****Distribution output to file started
  sprintf(fname1,"./%s_dat/pre_total_field_%d_%1.1f_hr-image.dat",argv[1],N_pix,int_time);
  printf("Writing to %s\n",fname1);
  F=fopen(fname1, "w");
  for(int i=0;i<N_pix;i++){
		for(int j=0;j<N_pix;j++){
			for(int k=0;k<count;k++){
				fprintf(F,"%0.8lE\n",total_signal[k+count*(j+N_pix*i)]);
			}//for k
		}//for j
	}//for i

  fclose(F);
  //****Distribution output to file finished
   */

  //****Writing the distribution to binary file started
  sprintf(fname1,"../%s_pre_total_field_%d_%1.1f_hr-image.bin",argv[1],N_pix,int_time);
  printf("Writing to %s\n",fname1);
  F=fopen(fname1, "wb");
  fwrite(total_signal, sizeof(double), count*N_pix*N_pix,F);
  fclose(F);
  //****Writing the distribution to binary file finished

  free(frequency);
  free(fname1);
  free(total_signal);

  return(0);
}