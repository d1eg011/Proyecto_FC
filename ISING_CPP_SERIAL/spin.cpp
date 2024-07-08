#include "spin.hpp"

Spin::Spin(const unsigned int &size){
  N = size;
  energy = 2*size*size;
  magnetization = size*size;
  lattice = std::vector<int>(size*size, 1);
}

Spin::Spin(const Spin &obj){
  N = obj.N;
  energy = obj.energy;
  magnetization = obj.magnetization;
  lattice = obj.lattice;
}

Spin::~Spin(){

}

Spin &Spin::operator=(const Spin &obj){
  N = obj.N;
  energy = obj.energy;
  magnetization = obj.magnetization;
  lattice = obj.lattice;
  
  return *this;
}

unsigned int Spin::get_N(){
  return N;
}

int Spin::get_energy(){
  return energy;
}

int Spin::get_magnetization(){
  return magnetization;
}

std::vector<int> Spin::get_lattice(){
  return lattice;
}

void Spin::lattice_cmap(){
// Initialize Gnuplot object
    Gnuplot gp;
  
    // Set output terminal to JPEG format
    gp << "set terminal epslatex color size 6.0in,5.0in standalone font "" 14\n";

    // Set output file name
    gp << "set output 'Figure.tex'\n";
  
    gp << "set border linewidth 6\n";
  
    // Set plot title and labels 
    gp << "set title '2D Lattice'\n";
  
    // Define color palette
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
            gp << lattice[ (row * N) + col] << " ";
        }
        gp << "\n";
    }
    gp << "e\n";
// Close the output file
    gp << "set output\n";

    // Flush the Gnuplot commands to ensure the file is written
    gp.flush();

    std::system("latex Figure.tex");
    std::system("dvips Figure.dvi");
    std::system("ps2pdf Figure.ps");

    // Remove unnecessary files
    std::remove("Figure.tex");
    std::remove("Figure.log");
    std::remove("Figure.aux");
    std::remove("Figure-inc.eps");
    std::remove("Figure.dvi");
    std::remove("Figure.ps");
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
  energy = 0;
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
  magnetization = 0;
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
      energy += dE;
      magnetization += dM;
    }
  }

}
