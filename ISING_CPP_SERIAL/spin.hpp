#ifndef SPIN_HPP 
#define SPIN_HPP

#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>
#include <cmath>
#include <cstdio>
#include "gnuplot-iostream.h"

class Spin
{
  private:
    Spin(); //Default constructor
    unsigned int N; //Dimension of the lattice
    int energy;
    int magnetization;
    std::vector<int> lattice;
    //private methods access only by configuration update
    int close_neighbord_energy(const unsigned int &row, const unsigned int &col);

  public:

    //Constructors
    Spin(const unsigned int &size);
    ~Spin();
    Spin(const Spin &obj);
    Spin &operator=(const Spin &obj);
    
    //Methods
    unsigned int get_N();
    int get_energy();
    int get_magnetization();
    std::vector<int> get_lattice();
    int calcEnergy();
    int calcMagnetization();
    void print_lattice();
    void lattice_cmap();
    void configuration_update(const double &beta, const float &J, const float &H, const unsigned int &max_iter); 
};
#endif

