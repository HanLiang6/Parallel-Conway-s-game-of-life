# Parallel implementation of Conway's game of life
This is a parallel implementation of Conway's game of life.

## Conway's game of life
The game takes place in a grid (or board) of cells (or squares). Each cell can be empty or occupied. The occupied cells change from one iteration (or generation) to the next. To go from one generation to the next, each cell in the grid is examined to see if it will be occupied or not in the next generation. This determination is made according to the following rules:
- If an occupied cell has less than 2 neighbors, it will become empty.
* If an occupied cell has more than 3 neighbors, it will become empty.
+ If an empty cell has exactly 3 neighbors, it will become occupied.

The neighbors of a cell are the 8 cells adjacent to it along the directions N, NW, W, SW, S, SE, E, NE. Cells along an edge or at a corner have fewer neighbors.

The arguments to the program consists of two integers m and n specifying the size of the life grid (_m_ × _n_ grid) and an integer _gen_ specifying the number of generations.

The program take input from a file containing an initial state of each grid cell given one row of the grid per line. An empty cell is indicated with a 0 and an occupied cell is indicated with a 1. Consecutive entries in a row are separated by a single space. The program will run the specified number of generations and print the resulting grid in the same format.

The grid of life is partitioned on to the virtual grid of processors such that each processor is responsible for a subgrid of the grid of life. The subgrids assigned to processors is as close in size as possible.

## Interface
The code is runnable with the CLI interface:

`mpirun -np P ./life m n gen input output`

where _P_ is the number of processors, _m x n_ is the size of the grid, _gen_ is the number of generations and _input_ and _output_ are the files that record the initial and end state of the grid.

## Run the code
There are a few examples and reference that can be used to test the correctness of the code.

To generate a random initial state of a grid with given size, 

`mpirun -np 1 gen/generate m n output`

## Acknowledgement
This work is based on a course project of High Performance Computing, from the School of Computational Science and Engineering, Georgia Institute of Technology.  
