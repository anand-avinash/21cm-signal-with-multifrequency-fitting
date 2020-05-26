//An example to demonstrate the plotting of 2D pixel plot from array saved as bin file in row major order
//Refer to bin_1_plot.py

#include <stdlib.h>
#include <stdio.h>

double *array;
int main(){
  double parray[10]={-5,-4,-3,-2,-1,0,1.0,2,3,4};
  array=(double *)malloc(sizeof(double)*5*2); //let it to be x*y, x=2, y=5
  for(int i=0; i<2; i++){
    for(int j=0; j<5; j++){
      array[j+5*i]=parray[j+5*i];
      printf("%f\t",array[j+5*i]);
    }//for j
    printf("\n");
  }//for i
  FILE *file;
  file=fopen("bindata_1.bin","wb");
  fwrite(array,sizeof(double),10,file);
  fclose(file);
  free(array);
}//main