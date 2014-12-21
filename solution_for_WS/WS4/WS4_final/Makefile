CC = mpicc
CFLAGS = -Wall -pedantic

SRC_SIM 	= 	helper 		\
  		    	init 		\
      			boundary_val 	\
		      	uvp 		\
			visual 		\
			parallel 	\
			sor		\
		      	main


all:
	@make sim

sim: $(SRC_SIM:%=%.o)
	$(CC) $(CFLAGS) -o sim $(SRC_SIM:%=%.o)  -lm


%.o : %.c
	$(CC) -c $(CFLAGS) $*.c -o $*.o

clean:
	/bin/rm -f $(SRC_SIM:%=%.o) sim

helper.o     	: helper.h 
init.o       	: helper.h init.h 
boundary_val.o  : helper.h boundary_val.h 
uvp.o        	: helper.h uvp.h
visual.o     	: helper.h
parallel.o 	: parallel.h
sor.o		: sor.h

$main.o       	: helper.h init.h boundary_val.h uvp.h visual.h parallel.h sor.h
