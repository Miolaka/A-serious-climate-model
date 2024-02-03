MyModel: repwvl_thermal.o repwvl_solar.o ascii.o radiation_model.o functions.o
	gcc -o MyModel repwvl_thermal.o repwvl_solar.o ascii.o radiation_model_rep_wvl.o functions.o -lm -L/opt/homebrew/lib -lnetcdf

repwvl_thermal.o: repwvl_thermal.c
	gcc -c $< -o $@ -I /opt/homebrew/include

repwvl_solar.o: repwvl_solar.c
	gcc -c $< -o $@ -I /opt/homebrew/include

ascii.o: ascii.c
	gcc -c ascii.c

radiation_model.o: radiation_model_rep_wvl.c
	gcc -I/opt/homebrew/include -c radiation_model_rep_wvl.c functions.c

functions.o: functions.c
	gcc -c functions.c -o functions.o