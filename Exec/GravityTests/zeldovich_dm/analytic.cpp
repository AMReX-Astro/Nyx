/* Analytic solution for Zel'dovich test (for plots) 
   Units are M_sun, Mpc, and km/s */ 

#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

int main(int argc, char *argv[]) {
   const double PI = M_PI; 
   // Parameters to play with: 
   const int    Np       = 1000;  // Number of points 
   const double box_size = 15.0;  // Box size 
   const double z_i      = 6.0;   // Initial redshift
   const double z_c      = 5.0;   // Redshift of first caustic 
   const double H_0      = 50.0;  // Hubble constant 
   const double lambda   = 5.0;   // Perturbation wavelength 

   // Main stuff: 
   double q, x, vx; 
   const double kk = 2.0*PI/lambda; 

   ofstream Pfile; 
   Pfile.open("solution.dat"); 
   Pfile.precision(8);

   for (int i=0; i<Np; ++i) { 
      q = i*box_size/Np; 
      x = q - (1.0 + z_c) / (1.0 + z_i) * sin(kk * q) / kk; 
      vx = -H_0 * (1.0 + z_c) * sqrt(1.0 + z_i) * sin(kk * q) / kk; 
      vx = vx/(1.0+z_i);  // Nyx convention: a*v
      Pfile << x << " " << vx << endl; 
   }

   Pfile.close();
   return 0;
}
