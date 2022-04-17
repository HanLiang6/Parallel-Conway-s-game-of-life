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

//update the internal subgrid
void update_innerstate(int m, int n, const int* in_grid, int* out_grid) {
  
  for (int i = 1; i < m-1; i++) { // For each row
    for (int j = 1; j < n-1; j++) { // For each column
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

            /* Check that the neighbor is actually a neighbor */
            if (!(k == 0 && l == 0)) {
              /* If it is alive, count it as alive */
                if (in_grid[neighbor_lin_loc]) {
                alive++;
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
      MPI_Bcast(&croot, 1, MPI_INT, 0, MPI_COMM_WORLD);
      //Broadcast the rank of the root in the Cartesian network
      
      int croot_coords[2];
      MPI_Cart_coords(comm, croot, NDIM, croot_coords);
      
      int sub_m, sub_n, sub_m_last, sub_n_last;
      sub_m = (m % dims[0])? (m / dims[0] + 1) : (m / dims[0]);
      sub_n = (n % dims[1])? (n / dims[1] + 1) : (n / dims[1]);
      sub_m_last = m - sub_m * (dims[0] - 1);
      sub_n_last = n - sub_n * (dims[1] - 1);
      
      std::vector<int> local_data;
      
      if(coords[0] == dims[0] - 1){
          if(coords[1] == dims[1] - 1){
              local_data.reserve(sub_m_last * sub_n_last);
          }else{
              local_data.reserve(sub_m_last * sub_n);
          }
      }else{
          if(coords[1] == dims[1] - 1){
              local_data.reserve(sub_m * sub_n_last);
          }else{
              local_data.reserve(sub_m * sub_n);
          }
      }
      
      MPI_Datatype tmp, col_type, col_type_last;
      MPI_Type_vector(sub_m, 1, n, MPI_INT, &tmp);
      MPI_Type_create_resized(tmp, 0, sizeof(int), &col_type);
      MPI_Type_commit(&col_type);
      MPI_Type_vector(sub_m_last, 1, n, MPI_INT, &tmp);
      MPI_Type_create_resized(tmp, 0, sizeof(int), &col_type_last);
      MPI_Type_commit(&col_type_last);
      
      int srank;
      int scoords[NDIM];
      MPI_Request req;
      MPI_Status status;
      
      if(rank==0){
          int displs;
      for(int i=0; i<dims[0]; i++){
          for(int j=0; j<dims[1]; j++){
              displs = i * n * sub_m + j * sub_n;
              scoords[0] = i;
              scoords[1] = j;
              MPI_Cart_rank(comm, scoords, &srank);
              if(i<dims[0]-1){
                  if(j<dims[1]-1){
                      MPI_Isend(&global_data[displs], sub_n, col_type,
                                srank, 1, comm, &req);
                  }else{
                      MPI_Isend(&global_data[displs], sub_n_last, col_type,
                                srank, 1, comm, &req);
                  }
              }else{
                  if(j<dims[1]-1){
                      MPI_Isend(&global_data[displs], sub_n, col_type_last,
                                srank, 1, comm, &req);
                  }else{
                      MPI_Isend(&global_data[displs], sub_n_last, col_type_last,
                                srank, 1, comm, &req);
                  }
              }
          }
      }
      }
      
      if(coords[0]<dims[0]-1){
          if(coords[1]<dims[1]-1){
              MPI_Recv(&local_data[0], sub_m * sub_n, MPI_INT,
                        croot, 1, comm, &status);
          }else{
              MPI_Recv(&local_data[0], sub_m * sub_n_last, MPI_INT,
                        croot, 1, comm, &status);
          }
      }else{
          if(coords[1]<dims[1]-1){
              MPI_Recv(&local_data[0], sub_m_last * sub_n, MPI_INT,
                        croot, 1, comm, &status);
          }else{
              MPI_Recv(&local_data[0], sub_m_last * sub_n_last, MPI_INT,
                        croot, 1, comm, &status);
          }
      }
      //Finish distributa data
      
      //persistent communication
      std::vector<int> edge_up, edge_down, edge_left, edge_right;
      int height, width;
      if(coords[0]<dims[0]-1){
          if(coords[1]<dims[1]-1){
              height = sub_m;
              width = sub_n
          }else{
              height = sub_m;
              width = sub_n_last;
          }
      }else{
          if(coords[1]<dims[1]-1){
              height = sub_m_last;
              width = sub_n
          }else{
              height = sub_m_last;
              width = sub_n_last;
          }
      }
      edge_up.resize(width+2,0);
      edge_down.resize(width+2,0);
      edge_left.resize(height,0);
      edge_right.resize(height,0);
      
      MPI_Request request_send_up, request_send_down, request_send_left, request_send_right;
      MPI_Request request_recv_up, request_recv_down, request_recv_left, request_recv_right;
      
      MPI_Datatype local_col_type;
      MPI_Type_vector(height, 1, width, MPI_INT, &tmp);
      MPI_Type_create_resized(tmp, 0, sizeof(int), &local_col_type);
      MPI_Type_commit(&local_col_type);
      
      if(coords[1]!=0){
          MPI_Send_init(&local_data[0], 1, local_col_type, coords[1]-1, 1, row_comm, &request_send_left);
      }
      if(coords[1]!=0){
          MPI_Recv_init(&edge_left[0], height, MPI_INT, coords[1]-1, 1, row_comm, &request_recv_left);
      }
      
      if(coords[1]!=dims[1]-1){
          MPI_Send_init(&local_data[0], 1, local_col_type, coords[1]+1, 1, row_comm, &request_send_right);
      }
      if(coords[1]!=dims[1]-1){
          MPI_Recv_init(&edge_left[0], height, MPI_INT, coords[1]+1, 1, row_comm, &request_recv_right);
      }
      
      //Collect data to rank 0 processor
      output_data.reserve(m * n);

      if(coords[0]<dims[0]-1){
          if(coords[1]<dims[1]-1){
              MPI_Isend(&local_data[0], sub_m * sub_n, MPI_INT,
                        croot, 1, comm, &req);
          }else{
              MPI_Isend(&local_data[0], sub_m * sub_n_last, MPI_INT,
                        croot, 1, comm, &req);
          }
      }else{
          if(coords[1]<dims[1]-1){
              MPI_Isend(&local_data[0], sub_m_last * sub_n, MPI_INT,
                        croot, 1, comm, &req);
          }else{
              MPI_Isend(&local_data[0], sub_m_last * sub_n_last, MPI_INT,
                        croot, 1, comm, &req);
          }
      }
      
      if(rank==0){
          int displs;
      for(int i=0; i<dims[0]; i++){
          for(int j=0; j<dims[1]; j++){
              displs = i * n * sub_m + j * sub_n;
              scoords[0] = i;
              scoords[1] = j;
              MPI_Cart_rank(comm, scoords, &srank);
              if(i<dims[0]-1){
                  if(j<dims[1]-1){
                      MPI_Recv(&output_data[displs], sub_n, col_type,
                                srank, 1, comm, &status);
                  }else{
                      MPI_Recv(&output_data[displs], sub_n_last, col_type,
                                srank, 1, comm, &status);
                  }
              }else{
                  if(j<dims[1]-1){
                      MPI_Recv(&output_data[displs], sub_n, col_type_last,
                                srank, 1, comm, &status);
                  }else{
                      MPI_Recv(&output_data[displs], sub_n_last, col_type_last,
                                srank, 1, comm, &status);
                  }
              }
          }
      }
      }
    
      MPI_Comm_free(&column_comm);
      MPI_Comm_free(&row_comm);
      MPI_Comm_free(&comm);
  }

  /* On root, output the data */
  if (rank == 0) {
    write_data(output_file, m, n, output_data);
  }

  MPI_Finalize();
}

//This is on a new branch
