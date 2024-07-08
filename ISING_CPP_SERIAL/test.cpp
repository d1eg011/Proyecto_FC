#include "spin.hpp"

int main(){
  unsigned int N = 500;
  Spin lattice(N);
  lattice.configuration_update(0.0, 1.0, 0.0, 5000000);

  lattice.lattice_cmap();

  //lattice.configuration_update(2.0, 1.0, 0.0, 2000000);

  std::cout << "Energy: " << lattice.calcEnergy() << " " <<
               "Magnetization: " << lattice.calcMagnetization() << std::endl;

  return 0;
}

