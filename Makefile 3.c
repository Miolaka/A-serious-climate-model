repwvl_thermal.o: repwvl_thermal.c
    gcc -c $< -o $@ -I /opt/homebrew/include

    	gcc -I/opt/homebrew/include -c repwvl_thermal.c -L/opt/homebrew/lib -lnetcdf