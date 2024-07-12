#include "ising.hpp"
#include <sys/time.h>

double seconds()
{
  struct timeval tmp;
  double sec;
  gettimeofday( &tmp, (struct timezone *)0 );
  sec = tmp.tv_sec + ((double)tmp.tv_usec)/1000000.0;

  return sec;
}

int main(){
  unsigned int size = 1000;
  unsigned int its = 10000000;
  double beta = 0.1;
  double J = 1.0;
  double H = 0.0;
  
  Ising ising(size);

  double time_1 = seconds();
  ising.configuration_update(beta, J, H, its);  
  double time_2 = seconds();
  
  std::cout << "Time to complete loop: " << time_2 - time_1 << std::endl;
  return 0;
}
