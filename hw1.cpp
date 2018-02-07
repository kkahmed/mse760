/*
 *  hw01.cpp
 *
 *  Created: Feb 2018
 *   Author: Kazi
 *
 *  Usage:
 *      Parameter int N sets the number of cells per side of the NxNxN grid
 *      Parameter bool REDUCED = false: print coordinates in meters vs reduced
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <numeric>
#include <cmath>

// USER PARAMETERS
const int N = 05;           // Single side dimension of lattice i.e. NxNxN = #cells
const bool REDUCED = false; // Print coordinates in reduced units or not

//Set up the non-dim system
double a = 5.26E-10;    //m
double m = 6.60E-26;    //g
double sig = 3.40E-10;  //m
double eps = 0.0104;    //eV
double len = (double)(a/sig)*(N/2.0); //Used for shifting with PBC

/*
 * double pos: [X dir] [Y dir] [Z dir] [4 per cell] [x,y,z coords]
 */
void createLattice(double pos[N][N][N][4][3])
{
  double a2 = a/sig; //Nondimensional
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      for (int k=0; k<N; k++) {
        // l specifies which of the 4 atoms in the cell
        for (int l=0; l<4; l++) {
          // n specifies x y or z coordinate
          for (int n=0; n<3; n++) {
            pos[i][j][k][l][n] = (i+0.5*(l==1 || l==2))*a2*(n==0) +
                                  (j+0.5*(l==1 || l==3))*a2*(n==1) +
                                  (k+0.5*(l==2 || l==3))*a2*(n==2); // :-)
          } //n
        } //l
      } //k
    } //j
  } //i
} //createLattice

/*
 * bool reduced: if reduced units or meters are desired for output
 * double pos: [X dir] [Y dir] [Z dir] [4 per cell] [x,y,z coords]
 */
void printLattice(bool reduced, double pos[N][N][N][4][3])
{
  std::ofstream outfile("lattice.txt");
  outfile << std::scientific << std::setprecision(6);
  double posm[N][N][N][4][3]; //If we want meter output without changing original

  //Traverse atom positions array
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      for (int k=0; k<N; k++) {
        for (int l=0; l<4; l++) {
          if (!reduced) {
            for (int n=0; n<3; n++) {
              posm[i][j][k][l][n] = sig * pos[i][j][k][l][n]; //Scale back to meters
            } //n
            outfile << posm[i][j][k][l][0] << " \t" << posm[i][j][k][l][1] << " \t" << posm[i][j][k][l][2] << " \n";
          } else {
            outfile << pos[i][j][k][l][0] << " \t" << pos[i][j][k][l][1] << " \t" << pos[i][j][k][l][2] << " \n";
          }
        } //l
      } //k
    } //j
  } //i

  if (!reduced) {
    std::cout << "x, y, z coordinates printed to lattice.txt, in meters \n";
  } else {
    std::cout << "x, y, z coordinates printed to lattice.txt, reduced units \n";
  }
} //printLattice

/*
 * bool periodic: if using PBC as well when doing the calculation
 * double r1, r2: xyz coords of the two atoms
 *
 * return: the square of rij
 */
double calcRij2(bool periodic, double r1[3], double r2[3])
{
  if (periodic) {
    double dx = std::abs(r1[0]-r2[0]);
    double dy = std::abs(r1[1]-r2[1]);
    double dz = std::abs(r1[2]-r2[2]);
    if (dx>(len/*+0.0000001*/)) dx -= 2.0*len; //Extra offsets might be needed if
    if (dy>(len/*+0.0000001*/)) dy -= 2.0*len; //there was a machine prec. issue
    if (dz>(len/*+0.0000001*/)) dz -= 2.0*len; //but it seems to work fine now.
    return std::pow(dx,2.0) +  std::pow(dy,2.0) + std::pow(dz,2.0);
  } else {
    return std::pow(r1[0]-r2[0],2.0) +  std::pow(r1[1]-r2[1],2.0) + std::pow(r1[2]-r2[2],2.0);
  }
}

/*
 * double rij2: the square of rij between two atoms
 *
 * return: the LJ potential based on that distance rij
 */
double calcLJ(double rij2)
{
  double temp = std::pow(1.0/rij2,3.0);
  return 4.0*eps*(std::pow(temp,2.0) - temp);
}

/*
 * bool periodic: if using PBC as well when doing the calculation
 * double *pos: pointer to beginning of positions array
 *
 * return: the entire system LJ potential
 */
double calcSystemLJ(bool periodic, double *pos)
{
  int nAtom = (int)std::pow(N,3)*4;
  double r1[3], r2[3];
  double systemLJ = 0.0;

  //Iterate through atom positions
  for (int i=0; i<nAtom-1; i++) {
    for (int l=0; l<3; l++) {
      r1[l] = *(pos + 3*i + l);   //xyz
    }
    for (int j=i+1; j<nAtom; j++) {
      for (int k=0; k<3; k++) {
        r2[k] = *(pos + 3*j + k); //xyz
      }
      systemLJ += calcLJ( calcRij2(periodic,r1,r2) );
    }
  }
  return systemLJ/(double)nAtom;
}

int main(int argc, char const *argv[])
{
  //Set up the arrays
  double pos[N][N][N][4][3]; //(X,Y,Z,atom#,xyz)
  createLattice(pos);

  //Print the coordinates
  std::cout << std::scientific << std::setprecision(6);
  printLattice(REDUCED, pos); //true = reduced units, false = meters

  //Calculate the LJ potential
  double *p;
  p = pos[0][0][0][0];
  std::cout << "Calculating LJ for N = " << N << ", position array in bytes: " << sizeof(pos) << " \n";
  double systemLJ = calcSystemLJ(false,p);
  std::cout << "Non-periodic cohesive energy (eV): \t" << systemLJ << " \n";
  double periodicLJ = calcSystemLJ(true,p);
  std::cout << "Periodic cohesive energy (eV): \t\t" << periodicLJ << " \n";

  return 0;
}
