#include "spin.hpp"

void plot(std::vector<double>& x, std::vector<double>& y, std::string file_name, std::string Title, std::string Xlabel, std::string Ylabel, std::string color, std::string Legend) {
    Gnuplot gp;

    // Set the terminal to pdfcairo for direct PDF output
    gp << "set terminal pdfcairo enhanced color font 'CMU Serif,10' linewidth 2 size 6.0in, 5.0in\n";
    gp << "set output '" << file_name << ".pdf'\n";
  
    gp << "set border linewidth 2\n";
    gp << "set xlabel '" << Xlabel << "'\n";
    gp << "set ylabel '" << Ylabel << "'\n";
    gp << "set title '" << Title << "'\n";
    gp << "set ytics scale 2\n";
    gp << "set xtics scale 2\n";

    gp << "plot '-' with points pointtype 7 pointsize 0.75 linecolor rgb" << color << "title '" << Legend << "'\n";
    gp.send1d(boost::make_tuple(x, y));

    std::cout << "Plot saved to " << file_name << ".pdf" << std::endl;
}

std::vector<double> sampling(double &init, double &end, unsigned int &steps){
  
  std::vector<double> array(steps);
  double step_size = (end-init)/(steps);

  for (unsigned int i = 1 ; i < steps+1 ; ++i){
    array[i-1] = end - step_size*i;
  }

  return array;
}

int main(){
  //needed data
  float J = 1.0;
  float H = 0.0;
  unsigned int size = 50;
  unsigned int its = 500000;
  unsigned int steps = 50.0;
  double T_min = 0.01;
  double T_max = 4.0;

  std::vector<double> T(steps , 0.0);
  std::vector<double> mean_energy( steps, 0.0);
  std::vector<double> mean_magnetization(steps, 0.0);
  std::vector<double> meansq_energy( steps, 0.0 );
 //std::vector<double> sqmean_energy( steps, 0.0 );
 //std::vector<double> specificHeat( steps, 0.0 );

  Spin lattice(size);
  T = sampling(T_min,T_max, steps);

  for (unsigned i = 0; i < steps; ++i){
    lattice.configuration_update(1/T[i], J, H, its);
    
    mean_energy[i] = lattice.getEnergy() / (4.0*its);
    mean_magnetization[i] = std::abs(lattice.getMagnetization()) / (2.0*its);
    meansq_energy[i] = std::pow(lattice.getEnergy(), 2) / (4.0*its);


    lattice.configuration_reset();
  }
  
  plot(T, mean_energy, "Energy_plot", "Normalized mean energy per spin", "Temperature (kT)", "Energy", "'red'" ,"mean energy per spin");

  plot(T, mean_magnetization, "Magnetization_plot", "Normalized Magnetization per spin", "Temperature (kT)", "Magnetization", "'blue'" ,"mean magnetization per spin");

  //plot(T, specificHeat, "SpecHeat_plot", "Specific Heat", "Beta", "specific heat", "specific heat");


  return 0;
}
