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
    int sqenergy;
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
    int getEnergy();
    int getMagnetization();
    int getSqenergy();
    std::vector<int> get_lattice();
    double calcEnergy();
    double calcMagnetization();
    void print_lattice();
    void lattice_cmap();
    void configuration_reset();
    void configuration_update(const double &beta, const double &J, const double &H, const unsigned int &max_iter); 
};
#endif

