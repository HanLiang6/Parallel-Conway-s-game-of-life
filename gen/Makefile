CXX=mpic++
# CCFLAGS=-Wall -g
# activate for compiler optimizations:
CCFLAGS=-Wall --std=c++17 -O3 -g -Wall
LDFLAGS=

all: generate

generate: generate.cpp
	$(CXX) $(CCFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm -f *.o life
	
new:
	mpirun -np 1 ./generate 1000 1000 new_grid.txt
