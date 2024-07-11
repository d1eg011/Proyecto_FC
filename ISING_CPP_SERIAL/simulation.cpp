#include "ising.hpp"

int main(){

  unsigned int size = 500;
  unsigned int max_iter = 2000000;
  unsigned int step_size = 10000;
  double J = 1.0;
  double H = 0.0;
  double beta_hot = 0.2;
  double beta_cold = 2.0;
  
  unsigned int steps = (max_iter/step_size);

  Spin lattice(size);

  for (unsigned int i = 1; i < steps+1; ++i){
    lattice.configuration_update(beta_hot , J, H, i*step_size);
    lattice.lattice_cmap(i);
  }
  
  for (unsigned int i = 1; i < steps+1; ++i){
    lattice.configuration_update(beta_cold, J, H, i*step_size);
    lattice.lattice_cmap( i + steps);
  }


  return 0;
}
