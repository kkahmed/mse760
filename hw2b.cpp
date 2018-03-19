/*
 *  hw02b.cpp
 *
 *  Created: Mar 2018
 *   Author: Kazi
 *
 *  Usage:
 *      Parameter int N sets the number of cells per side of the NxNxN grid
 *      Parameter int NBINS sets how many histogram bins are used for g(r)
 *      Parameter double T0 is starting temperature in K
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <numeric>
#include <cmath>
//#include <random>

// USER PARAMETERS
const int N = 10;           // Single side dimension of lattice i.e. NxNxN = #cells
const int NBINS = 100;      // Histogram bins for RDF calculation
const double T0 = 500.0;    // In K, will be scaled using kb and eps
const bool REDUCED = true;  // Print coordinates in reduced units or not
const bool ALLFORCES = true;// Print all forces vs 50 forces

//Set up the non-dim system
double a = 5.70E-10;    //m
double m = 6.60E-26;    //kg
double sig = 3.40E-10;  //m
double eps = 0.0104;    //eV
double JeV = 1.60217662E-19;
double tau = sig*std::sqrt(m/(eps*JeV));
double kb = 8.61733248E-5;//eV/K
double fscale = JeV * (eps/sig); // J/m
double cell = (a/sig)*(double)N;
double len = cell/2.0; //Used for shifting with PBC
int nAtom = (int)std::pow(N,3)*4;
double rho = nAtom/std::pow(cell,3.0);

/*
 * double pos: [X dir] [Y dir] [Z dir] [4 per cell] [x,y,z coords]
 * double force: like pos, but just initializes to zero here
 */
void createLattice(double pos[N][N][N][4][3], double force[N][N][N][4][3])
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
            force[i][j][k][l][n] = 0.0;
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
  std::ofstream outfile("lattice10.txt");
  outfile << std::scientific << std::setprecision(6);
  double posm[N][N][N][4][3]; //If we want meter output without changing original

  if (!reduced) {
    std::cout << "x, y, z coordinates printing to lattice10.txt, in meters \n";
    outfile << "x, y, z coordinates, in meters \n";
  } else {
    std::cout << "x, y, z coordinates printing to lattice10.txt, reduced units \n";
    outfile << "x, y, z coordinates, reduced units \n";
  }

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
} //printLattice

/*
 * bool allForces: if printing all forces or just 50
 * bool reduced: if reduced units or meters are desired for output
 * double force: [X dir] [Y dir] [Z dir] [4 per cell] [x,y,z coords]
 */
void printForces(bool allForces, bool reduced, double force[N][N][N][4][3])
{
  std::ofstream outforce("forces10.txt");
  outforce << std::scientific << std::setprecision(6);
  double forceN[N][N][N][4][3]; //If we want Newton output without changing original
  int N1, N2, l1;
  if (!allForces) {
    N1 = 1; N2 = 5; l1 = 2;
    std::cout << "\nPrinting 50 forces \n";
  } else {
    N1 = N; N2 = N; l1=4;
    std::cout << "\nPrinting all forces \n";
  }

  if (!reduced) {
    std::cout << "x, y, z force components printing to forces10.txt, in Newtons \n";
    outforce << "x, y, z force components, in Newtons \n";
  } else {
    std::cout << "x, y, z force components printing to forces10.txt, reduced units \n";
    outforce << "x, y, z force components, reduced units \n";
  }

  //Traverse atom positions array
  for (int i=0; i<N1; i++) {
    for (int j=0; j<N2; j++) {
      for (int k=0; k<N2; k++) {
        for (int l=0; l<l1; l++) {
          if (!reduced) {
            for (int n=0; n<3; n++) {
              forceN[i][j][k][l][n] = fscale * force[i][j][k][l][n]; //Scale back to N
            } //n
            outforce << forceN[i][j][k][l][0] << " \t" << forceN[i][j][k][l][1] << " \t" << forceN[i][j][k][l][2] << " \n";
          } else {
            outforce << force[i][j][k][l][0] << " \t" << force[i][j][k][l][1] << " \t" << force[i][j][k][l][2] << " \n";
          }
        } //l
      } //k
    } //j
  } //i
} //printForces

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
    if (dx>len) dx -= cell;
    if (dy>len) dy -= cell;
    if (dz>len) dz -= cell;
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
 * double rij2: the square of rij between two atoms
 *
 * return: the a coefficient used in force calculation
 */
double calcAccCoeff(double rij2)
{
  double temp = std::pow(1.0/rij2,3.0);
  return 48.0*(std::pow(temp,2.0) - 0.5*temp)/rij2;
}

/*
 * int i: The first atom index
 * int j: The second atom index
 * double r1[3]: The first atom coordinates
 * double r2[3]: The second atom coordinates
 * double ac: a cofficient from calcAccCoeff
 * double *force: pointer to beginning of forces array
 */
void initForce(int i, int j, double r1[3], double r2[3], double ac, double *force)
{
  double r12[3];
  for (int l=0; l<3; l++) {
    //Check each direction and perform required shift before updating force
    r12[l] = r1[l]-r2[l];
    r12[l] = r12[l] - cell*(r12[l] > len) + cell*(r12[l] < -len);
    force[3*i + l] += ac*r12[l]; //Udate force based on atom pair
    force[3*j + l] -= ac*r12[l]; //Newton's law
  }
}

/*
 * double *force: pointer to beginning of forces array
 */
void resetForce(double force[N][N][N][4][3])
{
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      for (int k=0; k<N; k++) {
        // l specifies which of the 4 atoms in the cell
        for (int l=0; l<4; l++) {
          // n specifies x y or z coordinate
          for (int n=0; n<3; n++) {
            force[i][j][k][l][n] = 0.0;
          } //n
        } //l
      } //k
    } //j
  } //i
} //resetForce

/*
 * double *vel: pointer to beginning of velocities array
 *
 * return: kinetic energy per atom
 */
double initVelocity(double *vel)
{
  double rx1, rx2, rx12, rx22;
  bool reject = true;
  double v0 = std::sqrt(3.0*T0*kb/eps);
  double temp = 0.0;
  //Good psuedorandom method, but requires -std=c++11
  //std::uniform_real_distribution<double> unif(-1.0,1.0);
  //std::default_random_engine re;

  //Iterate through atom positions to randomize every velocity
  for (int i=0; i<nAtom; i++) {
    reject = true;
    while (reject) {
      //rx1 = unif(re); rx2 = unif(re); //Good psuedorandom method, but requires -std=c++11
      rx1 = 2.0*std::rand()/RAND_MAX - 1.0; rx2 = 2.0*std::rand()/RAND_MAX - 1.0; //alternative
      rx12 = std::pow(rx1,2.0); rx22 = std::pow(rx2,2.0);
      if (rx12 + rx22 < 1) {
        vel[3*i]     = v0*2.0*rx1*std::sqrt(1 - rx12 - rx22);
        vel[3*i + 1] = v0*2.0*rx2*std::sqrt(1 - rx12 - rx22);
        vel[3*i + 2] = v0*(1.0 - 2.0*(rx12 + rx22));
        reject = false;
        //std::cout <<   vel[3*i] << " \t" <<   vel[3*i + 1] << " \t" <<   vel[3*i + 2] << " \n";
        temp += std::pow(vel[3*i],2.0)+std::pow(vel[3*i+1],2.0)+std::pow(vel[3*i+2],2.0);
      }
    }
  }
  temp = eps*temp/(kb*3.0*(double)nAtom); //Check initial temp
  std::cout << "\nInitialized velocities to \t" << temp << " \t Kelvin \n";
  return (3.0/2.0)*kb*temp;
}

/*
 * bool periodic: if using PBC as well when doing the calculation
 * double *pos: pointer to beginning of positions array
 * double *force: pointer to beginning of forces array
 *
 * return: the entire system LJ potential per atom
 */
double calcSystemLJ(bool periodic, double *pos, double *force)
{
  double r1[3], r2[3];
  double systemLJ = 0.0;
  double rij2, tempLJ;

  //Iterate through atom positions
  for (int i=0; i<nAtom-1; i++) {
    for (int l=0; l<3; l++) {
      r1[l] = *(pos + 3*i + l);   //xyz
    }
    for (int j=i+1; j<nAtom; j++) {
      for (int k=0; k<3; k++) {
        r2[k] = *(pos + 3*j + k); //xyz
      }
      //systemLJ += calcLJ( calcRij2(periodic,r1,r2) );
      //Split the above step so we can utilize the info in updateForce too
      rij2 = calcRij2(periodic,r1,r2);
      systemLJ += calcLJ( rij2 );
      initForce(i,j,r1,r2,calcAccCoeff(rij2),force); //Update force for that pair
    }
  }
  return systemLJ/(double)nAtom;
}

/*
 * double *pos: pointer to beginning of positions array
 * string& name: name of the file to output
 */
void outputRDF(double *pos, std::string name)
{
  int bin;
  double ri = len/(double)NBINS;
  double midpoints[NBINS];  midpoints[0] = ri/2;
  double volumes[NBINS];    volumes[0] = (4.0/3.0)*M_PI*std::pow((double)ri,3.0);
  double gr[NBINS];         gr[0] = 0.0;
  double r1[3], r2[3];
  double rij;

  //Initialize
  for (int a=1; a<NBINS; ++a) {
    midpoints[a] = midpoints[a-1] + ri;
    volumes[a] = (4.0/3.0)*M_PI*(std::pow(ri*((double)a+1.0),3.0) - std::pow(ri*(double)a,3.0));
    gr[a] = 0.0;
  }

  //Iterate through atom positions
  for (int i=0; i<nAtom-1; i++) {
    for (int l=0; l<3; l++) {
      r1[l] = *(pos + 3*i + l);   //xyz
    }
    for (int j=i+1; j<nAtom; j++) {
      for (int k=0; k<3; k++) {
        r2[k] = *(pos + 3*j + k); //xyz
      }
      rij = std::sqrt(calcRij2(true,r1,r2));
      bin = (int)floor(rij/ri);
      if (bin < 100) gr[bin] += 2.0; //Assign to bin (ONLY WITHIN SPHERE!!!!!)
    }
  }

  //Process output
  std::ofstream rdffile(name.c_str());
  rdffile << std::scientific << std::setprecision(6);
  for (int b=0; b<NBINS; ++b) {
    gr[b] = gr[b]/(volumes[b]*(double)nAtom*rho);
    rdffile << midpoints[b] << " \t" << gr[b] << std::endl;//" \n";
  }
  rdffile.flush();
}

/*
 * double *vel: pointer to beginning of velocities array
 * double *force: pointer to beginning of forces array
 * double dt: the time step
 *
 * return: the entire system KE per atom
 */
double velocityStep(double *vel, double *force, double dt)
{
  double temp = 0.0;

  for (int i=0; i<nAtom; i++) {
    for (int l=0; l<3; l++) {
      //Update velocity based on acceleration
      vel[3*i + l] += (dt/2.0) * force[3*i + l];
      temp += std::pow(vel[3*i+l],2.0);
    }//xyz loop
  }//atom loop

  return eps*temp/(2.0*(double)nAtom);
}

/*
 * double *pos: pointer to beginning of positions array
 * double *vel: pointer to beginning of velocities array
 * double dt: the time step
 */
void positionStep(double *pos, double *vel, double dt)
{
  for (int i=0; i<nAtom; i++) {
    for (int l=0; l<3; l++) {
      //Update position based on velocity
      pos[3*i + l] += dt * vel[3*i + l];
      //Account for PBC
      if (pos[3*i + l] > cell) {
        pos[3*i + l] -= cell;
      } else if (pos[3*i + l] < 0.0) {
        pos[3*i + l] += cell;
      }
    }//xyz loop
  }//atom loop
}

int main(int argc, char const *argv[])
{
  //Set up the arrays
  double pos[N][N][N][4][3];    //(X,Y,Z,atom#,xyz)
  double vel[N][N][N][4][3];    //Velocities
  double force[N][N][N][4][3];  //Forces
  double *p;  double *f;  double *v;
  p = pos[0][0][0][0];
  v = vel[0][0][0][0];
  f = force[0][0][0][0];
  std::cout << std::scientific << std::setprecision(6);
  createLattice(pos,force); //Sets xyz coords, and sets forces to zero
  printLattice(!REDUCED, pos); //REDUCED units vs meters
  double KE = initVelocity(v); //Set up the initial random velocity

  //Calculate the LJ potential, calculate forces
  std::cout << "Calculating LJ for N = " << N << ", position array in bytes: " << sizeof(pos) << " \n";
  double PE = calcSystemLJ(true,p,f); //Also performs initForce
  std::cout << "Periodic cohesive energy (eV): \t\t" << PE << " \n";
  printForces(!ALLFORCES, !REDUCED, force); //ALLFORCES vs 50, REDUCED units vs meters

  //Set up for verlet iteration
  double dt[3] = {0.005,0.02,0.001};
  int step[3] = {501,126,2501};
  double out1[4][step[0]]; // At each time step:
  out1[0][0] = 0.0;
  out1[1][0] = KE;
  out1[2][0] = PE;
  out1[3][0]=KE+PE;
  std::cout << "\nKE/atom (eV) start: \t" << KE << " \n";
  std::cout << "PE/atom (eV) start: \t" << PE << " \n";

  std::ofstream out10005("out10005.txt");
  out10005 << std::scientific << std::setprecision(6);
  out10005 << out1[0][0] << "\t" << out1[1][0] << "\t"<< out1[2][0] << "\t"<< out1[3][0] << std::endl;

  outputRDF(p,"RDFinit.txt");

  //Verlet iteration
  std::cout << "\nWorking on dt=0.005 for 2.5LJ time \n";
  for (int k=1; k<step[0]; k++) {
    velocityStep(v,f,dt[0]);              //First velocity halfstep
    positionStep(p,v,dt[0]);              //Set position for next iteration
    resetForce(force);
    out1[2][k] = calcSystemLJ(true,p,f);  //Also performs initForce (update acc)
    out1[1][k] = velocityStep(v,f,dt[0]); //Set velocity for next iteration
    out1[0][k] = out1[0][k-1] + tau*dt[0];//Update time
    out1[3][k] = out1[1][k] + out1[2][k]; //Update total E per atom
    if (k % 50 == 0) {
      out10005 << out1[0][k] << "\t" << out1[1][k] << "\t"<< out1[2][k] << "\t"<< out1[3][k] << std::endl;
      std::cout << k/5 << "\% \t"; std::cout.flush(); //Flush, progress every 10%
    } else {
      out10005 << out1[0][k] << "\t" << out1[1][k] << "\t"<< out1[2][k] << "\t"<< out1[3][k] << " \n";
    }
  }

  outputRDF(p,"RDFfinal.txt");

  std::cout << "\n\nOutput energies printed to outNNttt.txt \n";
  std::cout << "Output format: time (s), Ek/at, Ep/at, Et/at (eV)\n";
  std::cout << "Initial and final g(r) printed to RDFinit.txt and RDFfinal.txt\n";

  return 0;
}
