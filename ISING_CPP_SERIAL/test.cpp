#include "spin.hpp"

int main(){
  unsigned int N = 20;
  Spin lattice(N);
  lattice.configuration_update(0.0, 1.0, 0.0, 100000);

  std::cout << lattice.getMagnetization()/(2*100000.0)  << std::endl;
  //lattice.lattice_cmap();



  return 0;
}

