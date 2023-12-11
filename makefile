MyModel: main.o functions.o
        gcc -o MyModel main.o ascii.o -lm

functions.o: ascii.c
        gcc -c ascii.c

main.o: main.c
        gcc -c main.c