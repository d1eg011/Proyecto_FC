#include "ising.hpp"

Spin::Spin(const unsigned int &size){
  N = size;
  energy = 0;
  sqenergy = 0;
  magnetization = size*size;
  lattice = std::vector<int>(size*size, 1);
}

Spin::Spin(const Spin &obj){
  N = obj.N;
  energy = obj.energy;
  sqenergy = obj.sqenergy;
  magnetization = obj.magnetization;
  lattice = obj.lattice;
}

Spin::~Spin(){

}

Spin &Spin::operator=(const Spin &obj){
  N = obj.N;
  energy = obj.energy;
  sqenergy = obj.sqenergy;
  magnetization = obj.magnetization;
  lattice = obj.lattice;
  
  return *this;
}

unsigned int Spin::get_N(){
  return N;
}

int Spin::getEnergy(){
  return energy;
}

int Spin::getMagnetization(){
  return magnetization;
}

int Spin::getSqenergy(){
  return sqenergy;
}

std::vector<int> Spin::get_lattice(){
  return lattice;
}

void Spin::lattice_cmap(const unsigned int& filenumber){
   
    Gnuplot gp;

    std::string filename = "lattice_" + std::to_string(filenumber);

    gp << "set terminal epslatex color size 6.0in,5.0in standalone font \" 14\"\n";
    gp << "set output '" << filename << ".tex'\n";

    gp << "set border linewidth 6\n";

    gp << "set title '2D Lattice'\n";
  
    gp << "set palette defined (-1 'white', 0 'gray', 1 'black')\n";

    // Set matrix dimensions
    gp << "set size ratio -1\n";
    gp << "set xrange [0:" << N - 1 << "]\n";
    gp << "set yrange [0:" << N - 1 << "]\n";
    gp << "set cbrange [-1:1]\n";
    gp << "set pm3d map\n";

    // Plot the matrix as an image
    gp << "plot '-' matrix with image title ' '\n";

    // Send matrix data to Gnuplot
    for (unsigned int row = 0; row < N; ++row) {
        for (unsigned int col = 0; col < N; ++col) {
            gp << lattice[(row * N) + col] << " ";
        }
        gp << "\n";
    }   
    gp << "e\n";
    gp << "set output\n";

    // Flush the Gnuplot commands to ensure the file is written
    gp.flush();

   // std::string cmdlatex = "latex " + filename + ".tex";
   // std::string cmddvips = "dvips " + filename + ".dvi";
   // std::string cmdps2pdf = "ps2pdf " + filename + ".ps";

    std::system( ("latex " + filename + ".tex").c_str() );
    std::system( ("dvips " + filename + ".dvi" ).c_str() );
    std::system( ("ps2pdf " + filename + ".ps" ).c_str() );

    std::remove( (filename + ".tex").c_str() );
    std::remove( (filename + ".log").c_str() );
    std::remove( (filename + ".aux").c_str() );
    std::remove( (filename + "-inc.eps").c_str() );
    std::remove( (filename + ".dvi").c_str() );
    std::remove( (filename + ".ps").c_str() );
}

int Spin::close_neighbord_energy(const unsigned int &row, const unsigned int &col){
  int dE = 0;
  dE = 2*lattice[ (row*N) + col ]*( lattice[ ((row-1)*N)%N + col ] + 
                                    lattice[ ((row+1)*N)%N + col ] + 
                                    lattice[ (row*N) + (col-1)%N ] + 
                                    lattice[ (row*N) + (col+1)%N ] );
  return dE;
}

void Spin::configuration_reset(){
  energy = 0;
  magnetization = 0;
  sqenergy = 0;
  lattice = std::vector<int>(N*N, 1);
}

void Spin::configuration_update(const double &beta, const double &J, const double &H, const unsigned int &max_iter){
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
      this->lattice[ (row*N) + col ] = (this->lattice[ (row*N) + col ] > 0) ? -1 : 1;
    }
    this->energy -= dE;
    this->magnetization += dM;
    this->sqenergy += dE*dE;
  }

}
