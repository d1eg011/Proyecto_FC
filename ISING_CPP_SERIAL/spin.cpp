#include "spin.hpp"

Spin::Spin(const unsigned int &size){
  N = size;
  lattice = std::vector<int>(size*size, 1);
}

Spin::Spin(const Spin &obj){
  N = obj.N;
  lattice = obj.lattice;
}

Spin::~Spin(){

}

Spin &Spin::operator=(const Spin &obj){
  N = obj.N;
  lattice = obj.lattice;
  
  return *this;
}

int Spin::get_N(){
  return N;
}

void Spin::print_lattice(){
  for (unsigned int i = 0; i < N; ++i){
    for (unsigned int j = 0; j < N; ++j){
      std::cout << lattice[ (i*N) + j ] << " ";
    }
    std::cout << std::endl;
  }
}

int Spin::close_neighbord_energy(const unsigned int &row, const unsigned int &col){
  int dE = 0;
  dE = 2*lattice[ (row*N) + col ]*( lattice[ ((row-1)*N)%N + col ] + 
                                    lattice[ ((row+1)*N)%N + col ] + 
                                    lattice[ (row*N) + (col-1)%N ] + 
                                    lattice[ (row*N) + (col+1)%N ] );
  return dE;
}

int Spin::calcEnergy(){
  int energy = 0;
  for (unsigned int row = 0; row < N; ++row){
    for (unsigned int col = 0; col < N; ++col){
      energy += 0.5*lattice[ (row*N) + col ]*(lattice[ ((row-1)*N)%N + col ] +
                                              lattice[ ((row+1)*N)%N + col ] +
                                              lattice[ (row*N) + (col-1)%N ] +
                                              lattice[ (row*N) + (col+1)%N ] );
    }
  }

  return energy;
}
  
int Spin::calcMagnetization(){
  int magnetization = 0;
  for (unsigned int i = 0; i < N*N; ++i){
    magnetization += lattice[i];
  }

  return magnetization;
}

void Spin::configuration_update(const double &beta, const float &J, const float &H, const unsigned int &max_iter){

  int dE = 0;
  int dM = 0;
  
  std::random_device rd;
  std::mt19937 rng(rd());
  std::uniform_real_distribution<double> rndb(0, 1);

  std::unordered_map<int, double> probE = {
    {-8, exp(8 * beta * J)},
    {-4, exp(4 * beta * J)},
    {0, 1},
    {4, exp(-4 * beta * J)},
    {8, exp(-8 * beta * J)},
  };

  std::unordered_map<int, double> probM{
    {-2, exp(-2*beta *H)},
    {2, exp(2*beta*H)},
  };

  for (unsigned int i = 0; i < max_iter; ++i){
	  unsigned int row = static_cast<unsigned int>( floor(rndb(rng) * N));
    unsigned int col = static_cast<unsigned int>( floor(rndb(rng) * N));

    dE = close_neighbord_energy(row, col);
    dM = 2*lattice[ (row*N) + col ];
    if (rndb(rng) < probM[dM]*probE[dE]){
      this->lattice[ (row*N) + col ] = (this->lattice[ (row*      N) + col ] > 0) ? -1 : 1;

    }
  }

}
