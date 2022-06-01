#include <cmath>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include "Elements.h"
#include "Tracker.h"
#include "Particle.h"
#define E0 180000
using namespace std;

Tracker::Tracker(int N_p, int N_e, double len) {
  N_planes = N_p;
  E_lost.reserve(N_p);
  N_elec = N_e;
  MeanT.reserve(N_e);
  lenght = len;
  rho = ANP["Si"].sp_gravity;
  ZoA = ANP["Si"].ZoA;
  X0 = ANP["Si"].rad_lenght;
  thetaMS = (13.6/E0)*sqrt(len*pow(10,-4)/X0)*(1+0.038*log(len*pow(10,-4)/X0));
  d_csi = 0.307075*ZoA*len*pow(10,-4)*rho;
  ris = new TH1D("Energy", "Energy Lost (T. Mean)", floor(sqrt(N_e)), 0.06, 0.25);
  ris->GetXaxis()->SetTitle("E (MeV)");
  ris->GetYaxis()->SetTitle("Number of entries");
  return;
}

double Tracker::Track(int N_t) {
  double Mean = 0;
  double StDev;
  for(int j=0; j<N_elec; j++) {      //ToyMC
    Electron* e = new Electron(0);
    detected=0;
    e->Move(10);
    for(int i=0; i<N_planes; ++i) {         //MS nei piani di silicio                                                    
      double dist = e->MoveMS(0.0003, thetaMS);
      if(r->Uniform() < 0.95) {
	double mp = 3.6*80*dist;
        detected++;
        E_lost[i] = (r->Landau(mp, d_csi));
      }
      else E_lost[i]=0;
      if(i!=N_planes-1) e->Move(0.03);
    }
    for(int k=0; k<N_t; ++k) {          //k=numero di volte in cui si effettua il troncamento                             
      int imax = 0;         //Troncamento                                                                            
      for(int i=1; i<N_planes; ++i) {
        if(E_lost[i]>E_lost[imax]) imax=i;
      }
      E_lost[imax]=0;
      detected--;
    }
    double DeltaE=0;
    for(int i=0; i<N_planes; ++i) DeltaE += E_lost[i];
    MeanT[j] = 1.0*DeltaE/detected;         //Media troncata su M piani utili
    ris->Fill(MeanT[j]);
    Mean += MeanT[j];
  }
  Mean = 1.0*Mean/N_elec;        //Media delle medie troncate
  double SumSq = 0;
  for(int j=0; j<N_elec; ++j) SumSq += pow(MeanT[j]-Mean, 2);
  StDev = sqrt(1.0*SumSq/N_elec);        //Deviazione standard delle medie troncate o Risoluzione assoluta
  TCanvas* c1 = new TCanvas("c1","c1");
  gStyle->SetOptStat("nemrou");
  c1->cd();
  ris->Draw();
  c1->SaveAs("Risoluzione.pdf", "pdf");
  delete c1, ris;
  return StDev/Mean;
}

  


