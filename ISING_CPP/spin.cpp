#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>
#include <cmath>
#include <omp.h>
#include "spin.hpp"

Spin::Spin(const unsigned int &size){
  N = size;
  lattice = std::vector<int>(size*size, 1);
  
  #pragma omp parallel reduction(+ : energy, magnetization)
  {
    for (unsigned int row = 0; row < size; ++row){
      #pragma omp for schedule(static)
      for (unsigned int col = 0; col < size; ++col){
        #pragma omp atomic
        energy += lattice[ (row*N) + col ]*( lattice[ ((row-1)*N)%N + col ] + 
                                      lattice[ ((row+1)*N)%N + col ] + 
                                      lattice[ (row*N) + (col-1)%N ] + 
                                      lattice[ (row*N) + (col+1)%N ] );
        magnetization += lattice[ (row*N) + col];
      }
    }
  }

}

Spin::Spin(const Spin &obj){
  N = obj.N;
  lattice = obj.lattice;
  energy = obj.energy;
  magnetization = obj.magnetization;
}

Spin::~Spin(){

}

Spin &Spin::operator=(const Spin &obj){
  N = obj.N;
  lattice = obj.lattice;
  energy = obj.energy;
  
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

int Spin::get_energy(){
  return energy;
  }

int Spin::get_magnetization(){
  return magnetization;
}

int Spin::close_neighbord_energy(const unsigned int &row, const unsigned int &col){
  int dE = 0;
  dE = 2*lattice[ (row*N) + col ]*( lattice[ ((row-1)*N)%N + col ] + 
                                    lattice[ ((row+1)*N)%N + col ] + 
                                    lattice[ (row*N) + (col-1)%N ] + 
                                    lattice[ (row*N) + (col+1)%N ] );
  return dE;
}

void Spin::calcEM(){
  energy = 0;
  magnetization = 0;
  #pragma omp parallel reduction(+ : energy, magnetization)
  {
    for (unsigned int row = 0; row < N; ++row){
      #pragma omp for schedule(static)
      for (unsigned int col = 0; col < N; ++col){
        #pragma omp atomic
        energy += lattice[ (row*N) + col ]*(lattice[ ((row-1)*N)%N + col ] +
                                      lattice[ ((row+1)*N)%N + col ] +
                                      lattice[ (row*N) + (col-1)%N ] +
                                      lattice[ (row*N) + (col+1)%N ] );
        magnetization += lattice[ (row*N) + col];
      }
    }
  }

}
  

Spin Spin::configuration_update(const double &beta, const float &J, const float &H, const unsigned int &max_iter){
  int dE = 0;
  int dM = 0;
  
  std::random_device rd;
  std::mt19937 rng(rd());
  //std::uniform_int_distribution<unsigned int> rndint(0, N - 1);
  std::uniform_real_distribution<double> rndb(0, 1);

  std::unordered_map<int, double> probE = {
    {-8, exp(-8 * beta * J)},
    {-4, exp(-4 * beta * J)},
    {0, 1},
    {4, exp(4 * beta * J)},
    {8, exp(8 * beta * J)},
  };

  std::unordered_map<int, double> probM{
    {-2, exp(-2*beta *H)},
    {2, exp(2*beta*H)},
  };

  for (unsigned int i = 0; i < max_iter; ++i){
	  unsigned int row = static_cast<unsigned int>( floor(rndb(rng) * N));
    unsigned int col = static_cast<unsigned int>( floor(rndb(rng) * N));
    Spin lattice_copy = *this;   
    lattice_copy.lattice[ (row*N) + col ] = (lattice_copy.lattice[ (row*N) + col ] > 0) ? -1 : 1;

    dE = lattice_copy.close_neighbord_energy(row, col);
    dM = 2*lattice_copy.lattice[ (row*N) + col ];
    if (rndb(rng) < probM[dM]*probE[dE]){
      *this = lattice_copy;
    }
  }

  calcEM();

  return *this;
}
