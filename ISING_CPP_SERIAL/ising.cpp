#include "ising.hpp"

Ising::Ising(const unsigned int &size){
  N = size;
  energy = 0;
  sqenergy = 0;
  magnetization = size*size;
  lattice = std::vector<int>(size*size, 1);
}

Ising::Ising(const Ising &obj){
  N = obj.N;
  energy = obj.energy;
  sqenergy = obj.sqenergy;
  magnetization = obj.magnetization;
  lattice = obj.lattice;
}

Ising::~Ising(){

}

Ising &Ising::operator=(const Ising &obj){
  N = obj.N;
  energy = obj.energy;
  sqenergy = obj.sqenergy;
  magnetization = obj.magnetization;
  lattice = obj.lattice;
  
  return *this;
}

unsigned int Ising::get_N(){
  return N;
}

int Ising::getEnergy(){
  return energy;
}

int Ising::getMagnetization(){
  return magnetization;
}

int Ising::getSqenergy(){
  return sqenergy;
}

std::vector<int> Ising::get_lattice(){
  return lattice;
}

void Ising::lattice_cmap(const unsigned int& filenumber){
   
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

    gp.flush();

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

void Ising::plot(std::vector<double>& x, std::vector<double>& y, std::string file_name, std::string Xlabel, std::string Ylabel, std::string color, std::string Legend){
   Gnuplot gp;

    // Set the terminal to pdfcairo for direct PDF output
    gp << "set terminal pdfcairo enhanced color font 'CMU Serif,10' linewidth 2 size 6.0in, 5.0in\n";
    gp << "set output '" << file_name << ".pdf'\n";
  
    gp << "set border linewidth 2\n";
    gp << "set xlabel '" << Xlabel << "'\n";
    gp << "set ylabel '" << Ylabel << "'\n";
    gp << "set ytics scale 2\n";
    gp << "set xtics scale 2\n";

    gp << "plot '-' with points pointtype 7 pointsize 0.75 linecolor rgb" << color << "title '" << Legend << "'\n";
    gp.send1d(boost::make_tuple(x, y));

    std::cout << "Plot saved to " << file_name << ".pdf" << std::endl;
}

std::vector<double> Ising::sampling(const double &init, const double &end, const unsigned int &steps){
  
  std::vector<double> array(steps);
  double step_size = (end-init)/(steps);

  for (unsigned int i = 1 ; i < steps+1 ; ++i){
    array[i-1] = end - step_size*i;
  }

  return array;
}

void Ising::simulTemp(const double& T_max, const double& T_min, const unsigned int& steps, const unsigned int& iterations, const double& J, const double& H){

  std::vector<double> T(steps , 0.0);
  std::vector<double> mean_energy( steps, 0.0);
  std::vector<double> mean_magnetization(steps, 0.0);
  std::vector<double> specificHeat( steps, 0.0 );

  T = sampling(T_min, T_max, steps);

  for (unsigned i = 0; i < steps; ++i){
    configuration_update(1/T[i], J, H, iterations);

    mean_energy[i] = energy / (4.0*iterations);
    mean_magnetization[i] = magnetization / (2.0*iterations);
    specificHeat[i] = ( 1/std::pow(T[i],2) ) * ( sqenergy /(8.0*iterations) - std::pow( energy /(1.0*iterations), 2)/8.0 );

    configuration_reset();
  }

  plot(T, mean_energy, "Energy_plot", "Temperature (kT)", "Energy", "'red'" ,"mean energy per spin");

  plot(T, mean_magnetization, "Magnetization_plot", "Temperature (kT)", "Magnetization", "'blue'" ,"mean magnetization per spin");

  plot(T, specificHeat, "SpecHeat_plot", "Temperature (kT)", "Specific Heat","'#00B1A0'"  , "specific heat");

}

int Ising::close_neighbord_energy(const unsigned int &row, const unsigned int &col){
  int dE = 0;
  dE = 2*lattice[ (row*N) + col ]*( lattice[ ((row + N - 1) % N) * N + col ] + 
                                    lattice[((row + 1) % N) * N + col] + 
                                    lattice[row * N + (col + N - 1) % N] +
                                    lattice[row * N + (col + 1) % N]
                                   );
  return dE; 
}

void Ising::configuration_reset(){
  energy = 0;
  sqenergy = 0;
  magnetization = 0;
  lattice = std::vector<int>(N*N, 1);
}

void Ising::configuration_update(const double &beta, const double &J, const double &H, const unsigned int &max_iter){
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
      lattice[ (row*N) + col ] = -lattice[ (row*N) + col ];
      energy += dE;
      sqenergy += dE*dE;
    }
    magnetization += dM;
  }

}
