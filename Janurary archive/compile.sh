gcc -Wall -c testrepwvl.c -I .
gcc -Wall -c ascii.c -I .
gcc -Wall -c repwvl_thermal.c
gcc -o testrepwvl testrepwvl.o repwvl_thermal.o ascii.o -lnetcdf -lm
