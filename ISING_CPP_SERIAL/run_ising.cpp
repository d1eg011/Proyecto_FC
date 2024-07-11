#include "ising.hpp"

int main(){
  //needed data
  float J = 1.0;
  float H = 0.0;
  unsigned int size = 50;
  unsigned int its = 500000;
  unsigned int steps = 50.0;
  double T_min = 0.01;
  double T_max = 4.0;

  Spin lattice(size);

  lattice.simulTemp(T_max, T_min, steps, its, J, H);

  return 0;
}
