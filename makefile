# compilation flags
CXX_FLAGS=-std=c++17 -O3 -Wall -Wextra -g -lstdc++
CFLAGS=-O3 -Wall -std=c99 -g
CC=gcc
BAM_FLAGS=-I $(HOME)/include/bamtools/ -L $(HOME)/lib/
SDSL_FLAGS= -I ~/include -L ~/lib
EXECS=scan.x scan_BAM_READER.x scan_extended.x bwtparse.x pfbwt.x pfbwt64.x

# targets not producing a file declared phony
.PHONY: all clean tarfile

all: scan.x scan_BAM_READER.x bwtparse.x pfbwt.x pfbwt64.x

external/gsa/gsacak.o: external/gsa/gsacak.c external/gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $<

external/gsa/gsacak64.o: external/gsa/gsacak.c external/gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $< -DM64

scan.x: src/scan.cpp src/parameters.cpp external/utils.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl

scan_extended.x: src/scan.cpp src/parameters.cpp external/utils.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl -DOUTPUT_EXTENDED_READ

scan_BAM_READER.x: src/scan.cpp src/parameters.cpp external/malloc_count.o external/utils.o
	$(CXX) $(CXX_FLAGS) $(BAM_FLAGS) -o $@ $^ -ldl -DBAM_READER -lbamtools -lz

bwtparse.x: src/bwtparse.cpp external/malloc_count.o external/utils.o
	$(CC) $(CXX_FLAGS) -o $@ $^ -ldl -lstdc++

# prefix free BWT construction
pfbwt.x: src/pfbwt.cpp src/parameters.cpp external/gsa/gsacak.o external/utils.o external/malloc_count.o
	$(CXX) $(CXX_FLAGS) $(SDSL_FLAGS) -o $@ $^ -ldl -lsdsl

pfbwt64.x: src/pfbwt.cpp src/parameters.cpp external/gsa/gsacak64.o external/utils.o external/malloc_count.o
	$(CXX) $(CXX_FLAGS) $(SDSL_FLAGS) -o $@ src/pfbwt.cpp src/parameters.cpp external/gsa/gsacak64.o external/utils.o external/malloc_count.o -ldl -DM64 -lsdsl


%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(EXECS)  external/*.o external/gsa/*.o

clean_all:
	rm -f $(EXECS)  external/*.o external/gsa/*.o test/xbwt_of_reference.x data/*.occ data/*.dict data/*.parse data/*.parse_old data/*.ilist data/*.extended_input data/*.full_children data/*.limits data/*.bwt_limits data/*.bwt
