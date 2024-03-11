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
	mpirun -np 1 ./generate 5000 5000 new_grid_5K.txt
	mpirun -np 1 ./generate 2000 2000 new_grid_2K.txt
	mpirun -np 1 ./generate 1000 1000 new_grid_1K.txt
	mpirun -np 1 ./generate 500 500 new_grid_500.txt
	mpirun -np 1 ./generate 250 250 new_grid_250.txt
	mpirun -np 1 ./generate 100 100 new_grid_100.txt
