/* 
Copyright 2020 Xavier Andrade <xavier@llnl.gov>

Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions
   and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions
   and the following disclaimer in the documentation and/or other materials provided with the
   distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse
   or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "run.hpp"

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <mpi.h>
#include <complex>

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

  minq::run<std::complex<double>>(nstates, npoints, cart_comm);
  
  err = MPI_Finalize();
  assert( err == 0 );

}
