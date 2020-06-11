#ifndef MINQ__RUN
#define MINQ__RUN

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

#include "randomize.hpp"
#include "orthogonalize.hpp"
#include "subspace_diagonalization.hpp"

#include <mpi.h>
#include <slate/slate.hh>

namespace minq {

template <class Type>
void run(long nstates, long npoints, MPI_Comm comm){

  int rank;
  auto ierr = MPI_Comm_rank(comm, &rank);
  assert(ierr == 0);
  
  using matrix = slate::Matrix<Type>;
  
  int nprocs[2];
  int periods[2];
  int coords[2];

  auto err = MPI_Cart_get(comm, 2, nprocs, periods, coords);
  assert(err == 0);

  //These are the blocksizes, essentially we need one block per
  //process because of the limitation of other operations like a 3D
  //FFT
  auto nbs = (nstates + nprocs[0] - 1)/nprocs[0];
  auto nbp = (npoints + nprocs[1] - 1)/nprocs[1];

  if(rank == 0) std::cout << "Allocating wave functions       :";
  std::cout.flush();
  matrix wavefunction(nstates, npoints, nbs, nbp, nprocs[0], nprocs[1], comm);
  wavefunction.insertLocalTiles();
  matrix hwavefunction(nstates, npoints, nbs, nbp, nprocs[0], nprocs[1], comm);
  hwavefunction.insertLocalTiles();
  MPI_Barrier(comm);
  if(rank == 0) std::cout << "    [  DONE  ]" << std::endl;
  
  if(rank == 0) std::cout << "Randomizing wave functions      :";
  std::cout.flush();
  {
    aux::randomize(wavefunction);
    aux::randomize(wavefunction);
  }
  MPI_Barrier(comm);
  if(rank == 0) std::cout << "    [  DONE  ]" << std::endl;

  if(rank == 0) std::cout << "Orthogonalizing wave functions  :";
  std::cout.flush();
  {
    orthogonalize(wavefunction);
  }
  MPI_Barrier(comm);
  if(rank == 0) std::cout << "    [  DONE  ]" << std::endl;

  aux::check_orthogonalization(wavefunction);

  if(rank == 0) std::cout << "Subspace diagonalizations       :";
  std::cout.flush();
  {
    subspace_diagonalization(hwavefunction, wavefunction);
  }
  MPI_Barrier(comm);
  if(rank == 0) std::cout << "    [  DONE  ]" << std::endl;
  
}

}

#endif

