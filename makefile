# compilation flags
CXX_FLAGS=-std=c++17 -O3 -Wall -Wextra -g -lstdc++
CFLAGS=-O3 -Wall -std=c99 -g
CC=gcc
BAM_FLAGS=-I $(HOME)/include/bamtools/ -L $(HOME)/lib/
EXECS=scan.x scan_BAM_READER.x scan_extended.x bwtparse.x pfbwt.x

# targets not producing a file declared phony
.PHONY: all clean tarfile

all: scan.x bwtparse.x pfbwt.x

scan.x: src/scan.cpp external/utils.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl

scan_extended.x: src/scan.cpp external/utils.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl -DOUTPUT_EXTENDED_READ

scan_BAM_READER.x: src/scan.cpp external/malloc_count.o external/utils.o
	$(CXX) $(CXX_FLAGS) $(BAM_FLAGS) -o $@ $^ -ldl -DBAM_READER -lbamtools -lz

bwtparse.x: src/bwtparse.cpp external/malloc_count.o external/utils.o
	$(CC) $(CXX_FLAGS) -o $@ $^ -ldl

# prefix free BWT construction
pfbwt.x: src/pfbwt.cpp external/gsa/gsacak.o external/utils.o external/malloc_count.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl

%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -f $(EXECS)  external/*.o external/gsa/*.o

clean_all:
	rm -f $(EXECS)  external/*.o external/gsa/*.o test/xbwt_of_reference.x data/*.occ data/*.dict data/*.parse data/*.parse_old data/*.ilist data/*.extended_input data/*.full_children data/*.bwt
