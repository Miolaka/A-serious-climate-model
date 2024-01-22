MyModel: repwvl_thermal.o ascii.o radiation_model.o functions.c
	gcc -o MyModel repwvl_thermal.o ascii.o radiation_model_rep_wvl.o -lm -L/opt/homebrew/lib -lnetcdf

repwvl_thermal.o: repwvl_thermal.c
	gcc -c $< -o $@ -I /opt/homebrew/include

ascii.o: ascii.c
	gcc -c ascii.c

radiation_model.o: radiation_model_rep_wvl.c
	gcc -I/opt/homebrew/include -c radiation_model_rep_wvl.c
