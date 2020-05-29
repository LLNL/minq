#include <iostream>
#include <cstdlib>
#include <cassert>
#include <mpi.h>

int main(int argc, char ** argv){

  int provided = 0;
  auto err = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided );
  assert( err == 0 );
  assert( provided == MPI_THREAD_MULTIPLE );

  int size = 0;
  int rank = 0;

  err = MPI_Comm_size( MPI_COMM_WORLD, &size );
  assert( err == 0 );
  err = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  assert( err == 0 );

  if(rank == 0){  
    if(argc != 3){
      std::cerr << "Please use '" << argv[0] << " <nstates> <npoints>'" << std::endl; 
      exit(1);
    }
  }
    
  long nstates = atoi(argv[1]);
  long npoints = atoi(argv[2]);

  if(rank == 0){  
    if(nstates < 0) {
      std::cerr << "The number of states must be positive. The value received was " << nstates << std::endl; 
      exit(1);
    }
    
    if(npoints < 0) {
      std::cerr << "The number of points must be positive. The value received was " << nstates << std::endl; 
      exit(1);
    }
  }
  
  int dims[2] = {0, 0};
  
  err = MPI_Dims_create(size, 2, dims);
  assert(err == 0);
  assert(dims[0]*dims[1] == size);

  int periods[2] = {1, 1};
  MPI_Comm cart_comm;
  err = MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart_comm);
  assert(err == 0);

  if(rank == 0){

    std::cout << "minq input:" << std::endl;
    std::cout << "  total procs =  " <<  size   << std::endl;
    std::cout << "  states procs = " << dims[0] << std::endl;
    std::cout << "  points procs = " << dims[1] << std::endl;
    std::cout << "  nstates  =     " << nstates << std::endl;
    std::cout << "  npoints  =     " << npoints << std::endl;

  }
  
  
  err = MPI_Finalize();
  assert( err == 0 );

}
