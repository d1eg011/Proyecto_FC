#ifndef SPIN_HPP 
#define SPIN_HPP

class Spin
{
  private:
    Spin(); //Default constructor
    unsigned int N; //Dimension of the lattice
    int energy = 0;
    int magnetization = 0;
    std::vector<int> lattice;
    //private methods access only by configuration update
    int close_neighbord_energy(const unsigned int &row, const unsigned int &col);
    void calcEM();

  public:

    //Constructors
    Spin(const unsigned int &size);
    ~Spin();
    Spin(const Spin &obj);
    Spin &operator=(const Spin &obj);
    
    //Methods
    void print_lattice();
    int get_N();
    int get_energy();
    int get_magnetization();
    Spin configuration_update(const double &beta, const float &J, const float &H, const unsigned int &max_iter); 
};
#endif

