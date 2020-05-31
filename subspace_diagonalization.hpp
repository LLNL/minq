#ifndef MINQ__SUBSPACE_DIAGONALIZATION
#define MINQ__SUBSPACE_DIAGONALIZATION

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

#include "overlap.hpp"

#include <slate/slate.hh>

#include <iostream>
namespace minq {

template <class Type>
auto subspace_diagonalization(slate::Matrix<Type>  & hwavefunction, slate::Matrix<Type>  & wavefunction){

  int nprocs[2];
  int periods[2];
  int coords[2];
  
  auto err = MPI_Cart_get(wavefunction.mpiComm(), 2, nprocs, periods, coords);
  assert(err == 0);  
  
  auto nstates = wavefunction.m();
  auto nbs = (nstates + nprocs[0] - 1)/nprocs[0];

  slate::Matrix<Type> hamiltonian(nstates, nstates, nbs, nbs, nprocs[0], nprocs[1], wavefunction.mpiComm());
  hamiltonian.insertLocalTiles();

  auto trans = conj_transpose(hwavefunction);
  
  slate::gemm(Type(0.1), wavefunction, trans, Type(0.0), hamiltonian);

  slate::HermitianMatrix<Type> hermitian_hamiltonian(slate::Uplo::Upper, hamiltonian);

  std::vector<double> eigenvalues(nstates);
  slate::heev(hermitian_hamiltonian, eigenvalues);

}

}


#endif
