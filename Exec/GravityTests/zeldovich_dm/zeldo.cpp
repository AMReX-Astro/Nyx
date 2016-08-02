/* Initial conditions for Zel'dovich test 
   Units are M_sun, Mpc, and km/s */ 

#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

int main(int argc, char *argv[]) {
   const double PI = M_PI; 
   // Parameters to play with: 
   const int    Np       = 64;    // Number of particles 
   const double box_size = 15.0;  // Box size 
   const double z_i      = 10.0;  // Initial redshift
   const double z_c      = 5.0;   // Redshift of first caustic 
   const double H_0      = 50.0;  // Hubble constant 
   const double lambda   = 5.0;   // Perturbation wavelength 

   // Main stuff: 
   double q, x, y, z, vx, vy=0, vz=0; 
   const double rho_c = 2.7752e11*(H_0/100.0)*(H_0/100.0); // critical density
   const double Delta = box_size/Np; 
   const double kk = 2.0*PI/lambda; 
   const double p_mass = rho_c * pow((box_size/Np), 3.0); 
   int prnt = Np/10; // for screen prints, no functionality 
   if (Np < 10) prnt=1; 

   ofstream Pfile; 
   Pfile.open("zeldovich_particles.ascii");
   Pfile.precision(8); 
   Pfile << Np*Np*Np << endl; 

   // 1D perturbation (x direction) of positions and velocities 
   cout << endl << "Zeldovich test with " << Np << "^3 particles" << endl;
   for (int i=0; i<Np; ++i) { 
      if (i%prnt == 0) cout << floor(i/prnt)*10 << " percent done." << endl; //rough
      q = (i+0.5)*Delta; 
      x = q - (1.0 + z_c) / (1.0 + z_i) * sin(kk * q) / kk; 
      vx = -H_0 * (1.0 + z_c) * sqrt(1.0 + z_i) * sin(kk * q) / kk; 
      vx = vx/(1.0+z_i);  // Nyx convention: a*v
      for (int j=0; j<Np; ++j) {
         y = (j+0.5)*Delta; 
         for (int k=0; k<Np; ++k) {
            z = (k+0.5)*Delta; 
            Pfile << x << " " << y << " " << z << " " << p_mass << " " 
                  << vx << " " << vy << " " << vz << endl; 
         }
      }
   }

   Pfile.close();
   return 0;
}
