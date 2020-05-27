#****Compilation begins

gcc -c -g functions.c -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas -lm
gcc -c -g power_spec.c -I/usr/local/include -L/usr/local/lib -lfftw3 -lm
gcc -c -g foreground.c -lm
gcc -g functions.o foreground.o power_spec.o main.c -I/usr/local/include -L/usr/local/lib -lfftw3 -lgsl -lgslcblas -lm -o execute
rm *.o

#****Compilation done

#****Making the directory tree

function dir {
  if [ ! -d "$1" ]; then 
    mkdir "$1"
  fi
}

dir results

cd results
dir additional
dir sim_cube
dir ska_cube

cd additional
dir ps

cd ../sim_cube
dir cube
dir ps

cd ../ska_cube
dir cube
dir ps

cd ../../plots
dir ska
dir add_ps

cd ska
dir los_1
dir los_2
dir los_3
dir 2dps_2
dir 3dps_2

cd ../add_ps
dir 2d_100
dir 2d_150
dir 2d_200
dir 2d_50
dir 3d

cd ../../
