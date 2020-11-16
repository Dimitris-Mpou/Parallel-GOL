OBJS = main.o functions.o
SOURCE = main.c functions.c
HEADER = functions.h
OUT = gol
CC = mpicc
FLAG = -c

all: $(OBJS)
	$(CC) -O3 $(OBJS) -fopenmp -o $(OUT) -lm

main.o: main.c
	$(CC) $(FLAG) main.c

functions.o: functions.c
	$(CC) $(FLAG) functions.c


clean:
	rm -f $(OBJS) $(OUT)
