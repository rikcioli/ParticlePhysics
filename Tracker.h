#include <vector>
#include <TH1D.h>
#include <TRandom3.h>

class Tracker {
 private:
  int N_planes = 10, N_elec, detected;
  double lenght;
  double rho, ZoA, X0;
  double thetaMS, d_csi;
  std::vector<double> E_lost;
  std::vector<double> MeanT;
  TH1D* ris;
  TRandom3* r = new TRandom3(0);
  
 public:

  Tracker(int N_p, int N_e, double len);     //Inizializzazione parametri
  double Track(int N_t);                       //ToyMC con calcolo risoluzione
  
};
