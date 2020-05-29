#include <iostream>
#include <cstdlib>

int main(int argc, char ** argv){
  
  if(argc != 3){
    std::cerr << "Please use '" << argv[0] << " <nstates> <npoints>'" << std::endl; 
    exit(1);
  }

  long nstates = atoi(argv[1]);
  long npoints = atoi(argv[2]);

  if(nstates < 0) {
    std::cerr << "The number of states must be positive. The value received was " << nstates << std::endl; 
    exit(1);
  }

  if(npoints < 0) {
    std::cerr << "The number of points must be positive. The value received was " << nstates << std::endl; 
    exit(1);
  }

  std::cout << "Minq input:" << std::endl;
  std::cout << "  nstates = " << nstates << std::endl;
  std::cout << "  npoints = " << npoints << std::endl;

}
