# compilation flags
CXX_FLAGS=-std=c++11 -O3 -Wall -Wextra -g
CFLAGS=-O3 -Wall -std=c99 -g
CC=gcc

EXECS=bwtparse bwtparse64 simplebwt simplebwt64 newscan.x pfbwt.x pfbwt64.x pfbwtNT.x pfbwtNT64.x

# targets not producing a file declared phony
.PHONY: all clean tarfile

all: $(EXECS) 


gsa/gsacak.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $<

gsa/gsacak64.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $< -DM64

bwtparse: bwtparse.c gsa/gsacak.o utils.o malloc_count.o
	$(CC) $(CFLAGS) -o $@ $^ -ldl

bwtparse64: bwtparse.c gsa/gsacak64.o utils.o malloc_count.o
	$(CC) $(CFLAGS) -o $@ $^ -ldl -DM64

simplebwt: simplebwt.c gsa/gsacak.o
	$(CC) $(CFLAGS) -o $@ $^

simplebwt64: simplebwt.c gsa/gsacak64.o
	$(CC) $(CFLAGS) -o $@ $^ -DM64

# cnewscan executable to scan gzipped files (currently not active) 
cnewscan.x: newscan.cpp malloc_count.o utils.o
	$(CXX) $(CXX_FLAGS) -DGZSTREAM -o $@ $^ -lgzstream -lz -ldl 

newscan.x: newscan.cpp malloc_count.o utils.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl


# prefix free BWT construction. malloc_count not used since not compatible with -pthread
pfbwt.x: pfbwt.cpp pfthreads.hpp gsa/gsacak.o utils.o xerrors.o malloc_count.o
	$(CXX) $(CXX_FLAGS) -o $@ pfbwt.cpp gsa/gsacak.o utils.o xerrors.o malloc_count.o -pthread -ldl

pfbwt64.x: pfbwt.cpp pfthreads.hpp gsa/gsacak64.o utils.o xerrors.o malloc_count.o
	$(CXX) $(CXX_FLAGS) -o $@ pfbwt.cpp gsa/gsacak64.o utils.o xerrors.o malloc_count.o -pthread -ldl -DM64

# TO BE REMOVED? (now pfbwt*.x works with malloc_count)
# prefix free BWT construction without threads: useful since supports malloc_count
pfbwtNT.x: pfbwt.cpp gsa/gsacak.o utils.o malloc_count.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl -DNOTHREADS

# as above for large files
pfbwtNT64.x: pfbwt.cpp gsa/gsacak64.o utils.o malloc_count.o
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl -DNOTHREADS -DM64

%.o: %.c %.h
	$(CC) $(CFLAGS) -c -o $@ $<

tarfile:
		tar -zcf bigbwt.tgz bigbwt newscan.cpp pfbwt.cpp pfthreads.hpp simplebwt.c bwtparse.c makefile utils.[ch] xerrors.[ch] f2s.py gsa/gsacak.[ch] gsa/LICENSE gsa/README.md malloc_count.[ch]

clean:
	rm -f $(EXECS) *.o gsa/*.o
