#include "overlap.hpp"

#include <slate/slate.hh>

#include <iostream>
namespace minq {

template <class Type>
auto orthogonalize(slate::Matrix<Type>  & wavefunction){

  int nprocs[2];
  int periods[2];
  int coords[2];
  
  auto err = MPI_Cart_get(wavefunction.mpiComm(), 2, nprocs, periods, coords);
  assert(err == 0);  
  
  auto nstates = wavefunction.m();
  auto nbs = (nstates + nprocs[0] - 1)/nprocs[0];
  
  slate::HermitianMatrix<Type> overlap(slate::Uplo::Lower, nstates, nbs, nprocs[0], nprocs[1], wavefunction.mpiComm());
  overlap.insertLocalTiles();

  minq::overlap(wavefunction, overlap);

  slate::potrf(overlap);

  auto cholesky = slate::TriangularMatrix<Type>(slate::Diag::NonUnit, overlap);

  slate::trsm(slate::Side::Left, Type(1.0), cholesky, wavefunction);

}

}
