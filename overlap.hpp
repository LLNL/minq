#include <slate/slate.hh>

namespace minq {

template <class Type>
auto overlap(slate::Matrix<Type> & wavefunction, slate::HermitianMatrix<Type> & olap){
  slate::herk(1.0, wavefunction, 0.0, olap);
}

}
