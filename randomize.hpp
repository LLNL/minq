#ifndef MINQ__RANDOMIZE
#define MINQ__RANDOMIZE

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

#include <slate/slate.hh>

#include <cstdlib>
#include <complex>

namespace minq {
namespace aux {

template <class Type>
Type random();

template <>
double random(){
  return 2.0*(drand48() - 0.5);
}

template <>
std::complex<double> random(){
  return std::complex<double>(random<double>(), random<double>());
}

//Randomize the wavefunction. This is just to initialize the
//matrix and it is not really representative of an intensive
//operation in a DFT code.

template <class Type>
void randomize(slate::Matrix<Type> & matrix){

  int rank;
  auto ierr = MPI_Comm_rank(matrix.mpiComm(), &rank);
  assert(ierr == 0);

  // we don't care about the quality of the numbers, just that they are different in every processor
  srand48(rank);
  
  for (long jj = 0; jj < matrix.nt(); ++jj) {
    for (long ii = 0; ii < matrix.mt(); ++ii) {
      if(not matrix.tileExists(ii, jj)) continue;

      auto tile = matrix(ii, jj);
      
      for(long jtile = 0; jtile < tile.nb(); jtile++){
        for(long itile = 0; itile < tile.mb(); itile++){
          tile.data()[itile + tile.stride()*jtile] = random<Type>();
        }
      }
      
    }
  }
  
}

}
}

#endif
