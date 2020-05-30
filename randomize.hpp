#include <cstdlib>

namespace minq {

template <class MatrixType>
void randomize(MatrixType & matrix){

  for (long jj = 0; jj < matrix.nt(); ++jj) {
    for (long ii = 0; ii < matrix.mt(); ++ii) {
      if (matrix.tileIsLocal( ii, jj )) {
        auto tile = matrix(ii, jj);

        for(long jtile = 0; jtile < tile.nb(); jtile++){
          for(long itile = 0; itile < tile.mb(); itile++){
            tile.data()[itile + tile.stride()*jtile] = drand48();
          }
        }

      }
    }
  }
  
}

}
