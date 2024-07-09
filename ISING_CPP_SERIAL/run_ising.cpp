#include "spin.hpp"

void plot(std::vector<double>& x, std::vector<double>& y, std::string file_name, std::string Title, std::string Xlabel, std::string Ylabel, std::string Legend) {
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

    gp << "plot '-' with points pointtype 7 pointsize 0.75 linecolor rgb 'blue' title '" << Legend << "'\n";
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
  unsigned int steps = 100.0;
  double beta_min = 0.02;
  double beta_max = 1.2;

  std::vector<double> beta(steps);
  std::vector<double> mean_energy( steps, 0.0);
  std::vector<double> mean_magnetization(steps, 0.0);
  
  Spin lattice(size);
  beta = sampling(beta_min,beta_max, steps);

  for (unsigned i = 0; i < steps; ++i){
    lattice.configuration_update(beta[i], J, H, its);
    
    mean_magnetization[i] = static_cast<double>(abs(lattice.calcMagnetization()))/ (size*size);

    lattice.configuration_reset();
  }
  

  plot(beta, mean_magnetization, "Magnetization_plot", "Normalized mean Magnetization", "Beta", "mean energy", "magnetization per spin");

  return 0;
}
