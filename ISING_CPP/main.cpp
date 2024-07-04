#include "spin.hpp"

int main(){
  unsigned int N = 20;
  Spin lattice(N);
  lattice.configuration_update(0.0, 1.0, 0.0, 10000);

  lattice.configuration_update(2.0, 1.0, 0.0, 2000000);

  std::cout << "Energy: " << lattice.calcEnergy() << " " <<
               "Magnetization: " << lattice.calcMagnetization() << std::endl;

  return 0;
}

