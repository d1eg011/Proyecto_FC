#ifndef SPIN_HPP 
#define SPIN_HPP

#include <iostream>
#include <vector>
#include <unordered_map>
#include <random>
#include <cmath>

class Spin
{
  private:
    Spin(); //Default constructor
    unsigned int N; //Dimension of the lattice
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
    void print_lattice();
    int get_N();
    int calcEnergy();
    int calcMagnetization();
    void configuration_update(const double &beta, const float &J, const float &H, const unsigned int &max_iter); 
};
#endif

