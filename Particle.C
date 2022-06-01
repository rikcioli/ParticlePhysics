#include <iostream>
#include <vector>
#include <cmath>
#include <TRandom3.h>
#include "Particle.h"
#define E0 180000
using namespace std;

TRandom3* gen = new TRandom3(0);

const double Electron::m0 = 0.511;

void Particle::SetStart(double z0, double thx0, double thy0, double E) {
    posZ = z0;
    theta[0] = gen->Gaus(thx0, 0.00001);
    theta[1] = gen->Gaus(thx0, 0.00001);
    energy = E;
  return;
}

void Particle::Move(double delZ) {
  posZ += delZ;
  return;
}

double Particle::MoveMS(double delZ, double sigma) {
  theta[0]+=gen->Gaus(0, sigma);
  theta[1]+=gen->Gaus(0, sigma);
  double distance = delZ/(cos(theta[0])*cos(theta[1]));
  posZ+=delZ;
  return distance;
}
 
Electron::Electron(int i):Particle() {
    energy = E0;
    posZ = 0;
    theta[0] = gen->Gaus(0, 0.00001);
    theta[1] = gen->Gaus(0, 0.00001);
  return;
} 
