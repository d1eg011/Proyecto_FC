#include "spin.hpp"

int main(){
  unsigned int N = 20;
  Spin lattice(N);
  lattice.configuration_update(0.5, 1.0, 0.0, 100000);


  std::cout << lattice.getSqenergy()/(4.0*100000.0) << "  " << std::pow(lattice.getEnergy()/(100000.0),2)/4  << std::endl;
  //lattice.lattice_cmap();



  return 0;
}

