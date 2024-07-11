#ifndef ISING_HPP 
#define ISING_HPP

#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include "gnuplot-iostream.h"

class Ising
{
  private:
    Ising(); //Default constructor
    unsigned int N; //Dimension of the lattice
    int energy;
    int magnetization;
    int sqenergy;
    std::vector<int> lattice;
    //private methods access only by configuration update
    int close_neighbord_energy(const unsigned int &row, const unsigned int &col);
    void plot(std::vector<double>& x, std::vector<double>& y, std::string file_name, std::string Xlabel, std::string Ylabel, std::string color, std::string Legend);
    std::vector<double> sampling(const double &init, const double &end, const unsigned int &steps);

  public:

    //Constructors
    Ising(const unsigned int &size);
    ~Ising();
    Ising(const Ising &obj);
    Ising &operator=(const Ising &obj);
    
    //Methods
    unsigned int get_N();
    int getEnergy();
    int getMagnetization();
    int getSqenergy();
    std::vector<int> get_lattice();
    void lattice_cmap(const unsigned int& filenumber);
    void configuration_reset();
    void configuration_update(const double &beta, const double &J, const double &H, const unsigned int &max_iter);
    void simulTemp(const double& T_max, const double& T_min, const unsigned int& steps, const unsigned int& iterations, const double& J, const double& H);

};

#endif

