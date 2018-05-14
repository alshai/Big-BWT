# compilation flag
CXX_FLAGS=-std=c++11 -O2 -Wall -Wextra -g 
CFLAGS=-O2 -Wall -std=c99 -g
CC=gcc

EXECS=bwtparse simplebwt simplebwt64 newscan.x pfbwt.x pfbwt64.x 

all: $(EXECS)

gsacak.o: gsacak.c gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $<

gsacak64.o: gsacak.c gsacak.h
	$(CC) $(CFLAGS) -c -o $@ $< -DM64

bwtparse: bwtparse.c gsacak.o malloc_count.o
	$(CC) $(CFLAGS) -o $@ $^ -ldl

simplebwt: simplebwt.c gsacak.o
	$(CC) $(CFLAGS) -o $@ $^

simplebwt64: simplebwt.c gsacak64.o
	$(CC) $(CFLAGS) -o $@ $^ -DM64

# (new)scan executable to scan gzipped files 
c%.x: %.cpp malloc_count.o 
	$(CXX) $(CXX_FLAGS) -DGZSTREAM -o $@ $^ -lgzstream -lz -ldl 

newscan.x: newscan.cpp malloc_count.o  
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl

pfbwt.x: pfbwt.cpp gsacak.o malloc_count.o 
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl

pfbwt64.x: pfbwt.cpp gsacak64.o malloc_count.o 
	$(CXX) $(CXX_FLAGS) -o $@ $^ -ldl -DM64

# scanfq executable to process gzipped fastq files  
scanfq.x: scanfq.cpp
	$(CXX) $(CXX_FLAGS) -o $@ $^  -lgzstream -lz
	
tarfile:
		tar -zcf bigbwt.tgz bigbwt newscan.cpp pfbwt.cpp simplebwt.c bwtparse.c gsacak.[ch] makefile malloc_count.[ch]

clean:
	rm -f $(EXECS) gsacak.o
