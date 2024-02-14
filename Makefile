CC=gcc
CFLAGS=-ggdb -fopenmp -O3
LIBS=-lm

all: Troll

Troll: matriz.o teste.o
	$(CC) $(CFLAGS) matriz.o teste.o -o Troll $(LIBS)

matriz.o: matriz.c
	$(CC) $(CFLAGS) -c matriz.c

teste.o: teste.c
	$(CC) $(CFLAGS) -c teste.c

clean:
	rm -f *.o Troll
