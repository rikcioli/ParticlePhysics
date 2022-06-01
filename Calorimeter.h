#include <vector>
#include <cmath>
#include <TRandom3.h>
#include <TF1.h>
#include <TH1D.h>
#include "Particle.h"

class Calorimeter {
  
 private:
  int N_elec, N_step;
  double lenght, step;
  double rho, Ec, X0, BremsMU, PairMU;
  std::vector<Particle*> vecpart;
  TRandom3* rcal = new TRandom3(0);
  TF1* sigmaPair = new TF1("sigmaPair","1 - 4*x/3 + 4*x*x/3", 0, 1);
  TF1* sigmaBrems = new TF1("sigmaBrems","4/(3*x) - 4.0/3 + x", pow(10,-4), 1);
  TH1D* fracDep;

  void PairProd(Particle* f, Particle* e1, Particle* e2);
  bool Brems(Particle* e, Particle* f);
  void Frac_Dep();
  void Results(double dEdt[], int NdE);
  
 public:
  Calorimeter(double len, double l_step, int N_e);
  void ToyMC();
  ~Calorimeter();

};
