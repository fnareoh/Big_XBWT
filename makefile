# compilation flags
CXX_FLAGS=-std=c++11 -O3 -Wall -Wextra -g
CFLAGS=-O3 -Wall -std=c99 -g
CC=gcc
BAM_FLAGS=-I $(HOME)/include/bamtools/ -L $(HOME)/lib/

# executables not using threads (and therefore not needing the thread library)
EXECS_NT=newscanNT.x

# targets not producing a file declared phony
.PHONY: all clean tarfile

newscanNT.x: newscan.cpp utils.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl

newscanNT_BAM_READER.x: newscan.cpp malloc_count.o utils.o
	$(CXX) $(CXX_FLAGS) $(BAM_FLAGS) -o $@ $^ -ldl -DBAM_READER -lbamtools -lz


%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

tarfile:
		tar -zcf bigbwt.tgz bigbwt newscan.[ch]pp pscan.[ch]pp pfbwt.cpp pfthreads.hpp simplebwt.c bwtparse.c unparse.c remap.c makefile utils.[ch] xerrors.[ch] f2s.py gsa/gsacak.[ch] gsa/LICENSE gsa/README.md malloc_count.[ch]

clean:
	rm -f $(EXECS) $(EXECS_NT) *.o gsa/*.o
