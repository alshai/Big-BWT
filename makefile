# compilation flag
CXX_FLAGS=-std=c++11 -O2 -Wall -Wextra -g
CFLAGS=-O2 -Wall -std=c99 -g
CC=gcc

EXECS=bwtparse bwtparse64 simplebwt simplebwt64 newscan.x pfbwt.x pfbwt64.x 

all: $(EXECS)

gsa/gsacak.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $<

gsa/gsacak64.o: gsa/gsacak.c gsa/gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $< -DM64

utils.o: utils.c utils.h
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
c%.x: %.cpp malloc_count.o 
	$(CXX) $(CXX_FLAGS) -DGZSTREAM -o $@ $^ -lgzstream -lz -ldl 

newscan.x: newscan.cpp malloc_count.o  
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl

pfbwt.x: pfbwt.cpp pfthreads.hpp gsa/gsacak.o utils.o malloc_count.o 
	$(CXX) $(CXX_FLAGS) -o $@ pfbwt.cpp gsa/gsacak.o utils.o -lpthread

pfbwt64.x: pfbwt.cpp pfthreads.hpp gsa/gsacak64.o utils.o malloc_count.o 
	$(CXX) $(CXX_FLAGS) -o $@ pfbwt.cpp gsa/gsacak64.o utils.o malloc_count.o -lpthread  -ldl -DM64

tarfile:
		tar -zcf bigbwt.tgz bigbwt newscan.cpp pfbwt.cpp simplebwt.c bwtparse.c makefile gsa/gsacak.[ch] utils.[ch] gsa/LICENSE gsa/README.md malloc_count.[ch]

clean:
	rm -f $(EXECS) *.o gsa/*.o
