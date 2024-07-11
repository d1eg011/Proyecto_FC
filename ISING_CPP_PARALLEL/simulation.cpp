#include "ising.hpp"

int main(){

  unsigned int size = 500;
  unsigned int max_iter = 1000000;
  unsigned int step_size = 10000;
  double J = 1.0;
  double H = 0.0;
  double beta_hot = 0.2;
  double beta_cold = 2.0;
  
  unsigned int steps = (max_iter/step_size);

  Ising ising(size);

  for (unsigned int i = 1; i < steps+1; ++i){
    ising.configuration_update(beta_hot , J, H, i*step_size);
    ising.lattice_cmap(i, 1/beta_hot, i*step_size);
  }
  
  for (unsigned int i = 1; i < steps+1; ++i){
    ising.configuration_update(beta_cold, J, H, i*step_size*10);
    ising.lattice_cmap( i + steps, 1/beta_cold, i*step_size*10);
  }


  return 0;
}
