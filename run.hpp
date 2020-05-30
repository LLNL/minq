#include "randomize.hpp"
#include "orthogonalize.hpp"

#include <mpi.h>
#include <slate/slate.hh>

namespace minq {

template <class Type>
void run(long nstates, long npoints, MPI_Comm comm){

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
  
  matrix wavefunction(nstates, npoints, nbs, nbp, nprocs[0], nprocs[1], comm);
  wavefunction.insertLocalTiles();

  //Randomize the wavefunction. This is just to initialize the
  //matrix and it is not really representative of an intensive
  //operation in a DFT code.
  randomize(wavefunction);

  orthogonalize(wavefunction);

  
  
}

}
