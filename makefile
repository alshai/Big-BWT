# compilation flag
CXX_FLAGS=-std=c++11 -O2 -Wall -Wextra -g
CFLAGS=-O3 -Wall -std=c99 -g
CC=gcc

EXECS=bwtparse bwtparse64 simplebwt simplebwt64 newscan.x pfbwt.x pfbwt64.x 

all: $(EXECS)

gsa/gsacak.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $<

gsa/gsacak64.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $< -DM64

utils.o: utils.c utils.h
	$(CC) $(CFLAGS) -c -o $@ $<

xerrors.o: xerrors.c xerrors.h
	$(CC) $(CFLAGS) -c -o $@ $<

bwtparse: bwtparse.c gsa/gsacak.o utils.o malloc_count.o
	$(CC) $(CFLAGS) -o $@ $^ -ldl

bwtparse64: bwtparse.c gsa/gsacak64.o utils.o malloc_count.o
	$(CC) $(CFLAGS) -o $@ $^ -ldl -DM64

simplebwt: simplebwt.c gsa/gsacak.o
	$(CC) $(CFLAGS) -o $@ $^

simplebwt64: simplebwt.c gsa/gsacak64.o
	$(CC) $(CFLAGS) -o $@ $^ -DM64

# cnewscan executable to scan gzipped files (currently not active) 
cnewscan.x: newscan.cpp malloc_count.o 
	$(CXX) $(CXX_FLAGS) -DGZSTREAM -o $@ $^ -lgzstream -lz -ldl 

newscan.x: newscan.cpp malloc_count.o  
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl

# prefix free BWT construction. malloc_count not used since not compatible with -pthread
pfbwt.x: pfbwt.cpp pfthreads.hpp gsa/gsacak.o utils.o xerrors.o
	$(CXX) $(CXX_FLAGS) -o $@ pfbwt.cpp gsa/gsacak.o utils.o xerrors.o -pthread 

pfbwt64.x: pfbwt.cpp pfthreads.hpp gsa/gsacak64.o utils.o xerrors.o
	$(CXX) $(CXX_FLAGS) -o $@ pfbwt.cpp gsa/gsacak64.o utils.o xerrors.o -pthread -DM64

tarfile:
		tar -zcf bigbwt.tgz bigbwt newscan.cpp pfbwt.cpp pfthreads.hpp simplebwt.c bwtparse.c makefile gsa/gsacak.[ch] utils.[ch] xerrors.[ch] gsa/LICENSE gsa/README.md malloc_count.[ch]

clean:
	rm -f $(EXECS) *.o gsa/*.o
