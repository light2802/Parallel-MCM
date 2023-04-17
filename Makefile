CC = g++ #-std=c99
cc = gcc
CFLAGS = -fopenmp -Drestrict=__restrict__ -O2 -DNDEBUG -ffast-math # -g -pg
LDFLAGS = -O2

all: msBFSGraft

graphgenBP.o: graphgenBP.h graphgenBP.cpp ThreadedMMReader.h
	$(CC) $(CFLAGS) -c -o graphgenBP.o graphgenBP.cpp 

msBFSGraft.o: msBFSGraft.cpp graphgenBP.h graphgenBP.o maximalMatching.o
	$(CC) $(CFLAGS) -c -o msBFSGraft.o msBFSGraft.cpp 

msBFSGraft: msBFSGraft.o maximalMatching.o
	$(CC) $(CFLAGS) $(LDFLAGS)  -o msBFSGraft msBFSGraft.o graphgenBP.o maximalMatching.o -lm -lstdc++

maximalMatching.o: maximalMatching.cpp maximalMatching.h
	$(CC) $(CFLAGS) -c -o maximalMatching.o maximalMatching.cpp 

clean:
	-rm -f msBFSGraft *.o
