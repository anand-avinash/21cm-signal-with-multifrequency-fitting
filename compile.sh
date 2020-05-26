gcc -c -g functions.c -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas -lm
gcc -c -g power_spec.c -I/usr/local/include -L/usr/local/lib -lfftw3 -lm
gcc -c -g foreground.c -lm
gcc -g functions.o foreground.o power_spec.o main.c -I/usr/local/include -L/usr/local/lib -lfftw3 -lgsl -lgslcblas -lm -o execute
rm *.o
