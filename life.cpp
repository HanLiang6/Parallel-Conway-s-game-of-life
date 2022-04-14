#include <mpi.h>
#include <omp.h>

#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#define NDIM 2

/* Note: All code is provided as a potential reference serial implementation.
 *
 * You are not required to use said code
 *
 * Said code is not necessarily efficient, but is provided so that you can focus on the parallel portion
 *
 * Said code is also not necessarily representative of good C/C++ programming practices
 *
 * Said code may not necessarily lead to a nice parallel implementation
 *
 * You are able to make modifications to the code if required for your parallel implementation
 *
 * If you do find any bugs in the code, please report them on piazza
 *
 */

/** @brief Provide a default parameter for vec_to_str **/
template <typename T>
std::string vec_to_str(const std::vector<T>& vec,
                       const std::string& delim = ",");

/** @brief A function which converts a vector to a string **/
template <typename T>
std::string vec_to_str(const std::vector<T>& vec, const std::string& delim) {
  std::ostringstream oss;
  if (!vec.empty()) {
    std::copy(vec.begin(), vec.end() - 1,
              std::ostream_iterator<T>(oss, delim.c_str()));
    oss << vec.back();
  }
  return oss.str();
}

int linear_index(int m, int n, int row, int col) {
  return row * n + col;
}

/* A sequential function to update a grid,
 * uses LDA (leading dimension) to allow for subgrid considerations */
void update_state(int m, int n, const int* in_grid, int* out_grid) {
  
  for (int i = 0; i < m; i++) { // For each row
    for (int j = 0; j < n; j++) { // For each column
      //Consider a single element
      int lin_loc = linear_index(m, n, i, j); //This is the linear index of the element

      int alive = 0; 

      /* Look at each neighbor */
      for (int k = -1; k < 2; k++) {
        for (int l = -1; l < 2; l++) {
          /* Figure out the index associated with each neighbor */
          int y_loc = i + k;
          int x_loc = j + l;
          int neighbor_lin_loc = linear_index(m, n, y_loc, x_loc);

          /* Ensure the considered neighbor is in bounds */
          if ((x_loc >= 0) && (y_loc >= 0) && (y_loc < m) && (x_loc < n)) {
            /* Check that the neighbor is actually a neighbor */
            if (!(k == 0 && l == 0)) {
              /* If it is alive, count it as alive */
              if (in_grid[neighbor_lin_loc]) {
                alive++;
              }
            }
          }
        }
      }

      /* Based on the number of alive neighbors, update the output accordingly */
      if (in_grid[lin_loc]) {
        if (alive < 2) {
          out_grid[lin_loc] = 0;
        } else if (alive > 3) {
          out_grid[lin_loc] = 0;
        } else {
          out_grid[lin_loc] = 1;
        }
      } else {
        if (alive == 3) {
          out_grid[lin_loc] = 1;
        } else {
          out_grid[lin_loc] = 0;
        }
      }
    }
  }
}

/* Read in the data from an input file */
void read_data(const std::string& input_filename, int m, int n,
               std::vector<int>& output_data) {
  output_data.reserve(m * n);
  std::ifstream input_file(input_filename, std::ios::in);
  for (int i = 0; i < m * n; i++) {
    std::string mystring;
    input_file >> mystring;
    output_data.push_back(std::stoi(mystring));
  }
}

/* Write output data to a file */
void write_data(const std::string& output_filename, int m, int n,
                const std::vector<int>& output_data) {
  std::ofstream output_file(output_filename, std::ios::out);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      output_file << output_data[linear_index(m, n, i, j)];
      if (j < n - 1) {
        output_file << " ";
      }
    }
    output_file << "\n";
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);

  const int m = std::stoi(argv[1]);
  const int n = std::stoi(argv[2]);
  const int gen = std::stoi(argv[3]);
  const std::string input_file(argv[4]);
  const std::string output_file(argv[5]);

  int rank;
  int size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  std::vector<int> global_data;

  /* On root, read in the data */
  if (rank == 0) {
    read_data(input_file, m, n, global_data);
  }

  std::vector<int> output_data;


  if (size == 1) { /* If serial, use the serial code */
    /* Allocate the output data buffer */
    output_data.reserve(global_data.size());

    /* For each generation update the state */
    for (int i = 0; i < gen; i++) {
      update_state(m, n, global_data.data(), output_data.data());

      /* Swap the input and output */
      if (i < gen - 1) {
        std::swap(global_data, output_data);
      }
    }
  } else {
    /* You implement this */
    //assert(0);
      int crank;
      int dims[NDIM] = {0,0};
      int period[NDIM] = {0,0};
      int coords[NDIM];
      MPI_Comm comm;
      
      MPI_Dims_create(size, NDIM, dims);
      MPI_Cart_create(MPI_COMM_WORLD, NDIM, dims, period, 1, &comm);
      
      MPI_Comm_rank(comm, &crank);
      MPI_Cart_coords(comm, crank, NDIM, coords);
      
      MPI_Comm row_comm, column_comm;
      MPI_Comm_split(comm, coords[0], coords[1], &row_comm);
      MPI_Comm_split(comm, coords[1], coords[0], &column_comm);
      
      int croot;
      if(rank==0){
          croot = crank;
      }
      MPI_Bcast(&croot, 1, MPI_INT, crank, MPI_COMM_WORLD);
      //Broadcast the rank of the root in the Cartesian network
      
      int croot_coords[2];
      MPI_Cart_coords(comm, croot, NDIM, croot_coords);
      
      
  }

  /* On root, output the data */
  if (rank == 0) {
    write_data(output_file, m, n, output_data);
  }

  MPI_Finalize();
}

//This is on a new branch
