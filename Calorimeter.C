#include <iostream>
#include <TCanvas.h>
#include <TStyle.h>
#include "Elements.h"
#include "Calorimeter.h"
#define E0 180000
using namespace std;

Calorimeter::Calorimeter(double len, double l_step, int N_e) {
  vecpart.reserve(40000);
  gRandom->SetSeed(0);
  rho = ANP["BGO"].sp_gravity;
  Ec = ANP["BGO"].crit_energy;
  X0 = ANP["BGO"].rad_lenght;
  N_elec = N_e;
  lenght = len;
  step = l_step*X0;
  N_step = len/step;
  BremsMU = -1/X0;
  PairMU = 7*BremsMU/9;
  fracDep = new TH1D("fracDep", "Energy Deposited", 50, 0, 1);
}

void Calorimeter::PairProd(Particle* f, Particle* e1, Particle* e2) {
  double frac = sigmaPair->GetRandom();
  e1->energy = (f->energy)*frac;
  e1->posZ = (f->posZ);
  e1->theta[0] = (f->theta[0]) + (rcal->Gaus(0, Electron::m0/(e1->energy)))/sqrt(2);
  e1->theta[1] = (f->theta[1]) + (rcal->Gaus(0, Electron::m0/(e1->energy)))/sqrt(2);
  e2->energy = (f->energy)*(1-frac);
  e2->posZ = (f->posZ);
  e2->theta[0] = (f->theta[0]) + (rcal->Gaus(0, Electron::m0/(e2->energy)))/sqrt(2);
  e2->theta[1] = (f->theta[1]) + (rcal->Gaus(0, Electron::m0/(e2->energy)))/sqrt(2);
  f->energy = 0;
}

bool Calorimeter::Brems(Particle* e, Particle* f) {
  double frac = sigmaBrems->GetRandom();
  bool hard = true;
  if(frac < 1.0/3) hard = false;
  f->energy = (e->energy)*frac;
  f->posZ = e->posZ;
  e->energy *= (1-frac);
  f->theta[0] = (e->theta[0]) + (rcal->Gaus(0, Electron::m0/(e->energy)))/sqrt(2);
  f->theta[1] = (e->theta[1]) + (rcal->Gaus(0, Electron::m0/(e->energy)))/sqrt(2);
  e->theta[0] += (rcal->Gaus(0, Electron::m0/(e->energy)))/sqrt(2);
  e->theta[1] += (rcal->Gaus(0, Electron::m0/(e->energy)))/sqrt(2);
  if((e->energy < 10.5) || (hard == true)) return true;
  return false;
}

void Calorimeter::Frac_Dep() {
  double E_lost = 0;       //Calcolo dell'energia persa in media
  double frac_E;
  for(int j=0; j<vecpart.size(); ++j) {
    double tmp = vecpart[j]->Enrg();
    if(tmp>Ec) E_lost+=tmp;
  }
  frac_E = E_lost/E0;
  fracDep->Fill(1-frac_E);     //frazione di energia depositata
  return;
}

void Calorimeter::Results(double dEdt[], int NdE) {
  gStyle->SetOptStat("nemr");
  //Risultati sulla frazione di energia depositata in media
  double alpha = 0.5*log(E0/Ec)+0.75;
  double frac_MC, frac_t;
  TF1* dEdt_Func = new TF1("dEdt_Func", "0.5*(((0.5*x)^([0]-1))*exp(-0.5*x))/TMath::Gamma([0])", 0, lenght/X0);
  dEdt_Func->SetParameter(0, alpha);
  TCanvas* c2 = new TCanvas("c2", "c2");
  c2->cd();
  fracDep->GetXaxis()->SetTitle("Fraction deposited");
  fracDep->GetYaxis()->SetTitle("Number of entries");
  fracDep->Draw();
  frac_MC = fracDep->GetMean();
  frac_t = dEdt_Func->Integral(0, lenght/X0);
  c2->SaveAs("Frazione.pdf", "pdf");
  
  //Risultati sull'energia per unità di lunghezza di radiazione
  double xAxis[NdE], yAxis[NdE];
  TCanvas* c3 = new TCanvas("c3", "c3");
  TH1D* dEdt_ToyMC = new TH1D("dEdt_ToyMC", "dE/dt", floor(lenght/X0), 0, floor(lenght/X0));
  TH1D* dEdt_FuncH = new TH1D("dEdt_FuncH", "dE/dt", floor(lenght/X0), 0, floor(lenght/X0));
  for(int i=1; i<=NdE; ++i) {
    xAxis[i-1] = i-1;
    yAxis[i-1] = 0.5*pow(0.5*i, alpha-1)*exp(-0.5*i)/TMath::Gamma(alpha);
  }
  c3->cd();
  dEdt_ToyMC->FillN(NdE, xAxis, dEdt, 1);
  dEdt_ToyMC->GetYaxis()->SetRangeUser(0, 0.1);
  dEdt_ToyMC->GetXaxis()->SetTitle("t (rad lenghts)");
  dEdt_ToyMC->GetYaxis()->SetTitle("Fraction deposited");
  dEdt_ToyMC->Draw("hist");
  dEdt_FuncH->FillN(NdE, xAxis, yAxis, 1);
  dEdt_FuncH->SetLineColor(2);
  dEdt_FuncH->Draw("hist same");
  c3->SaveAs("dEdt.pdf", "pdf");
  delete c2, c3, dEdt_Func, dEdt_ToyMC, dEdt_FuncH;

  cout<<"\nProgramma terminato. \n\nI seguenti grafici evidenziano determinate quantità, calcolate ad ogni iterazione:"<<endl;
  cout<<"Risoluzione.pdf -> Quantità di energia persa nel tracker di silicio;"<<endl;
  cout<<"Frazione.pdf -> Frazione di energia depositata nel calorimetro elettromagnetico;"<<endl;
  cout<<"dEdt.pdf -> Frazione di energia depositata in relazione alla profondità esaminata;"<<endl;
  cout<<"\nLa frazione di energia depositata, valutata con ToyMC (grafico Frazione), vale "<<frac_MC<<";"<<endl;
  cout<<"Rispetto al valore teorico di "<<frac_t<<", c'è uno scarto pari a "<<abs(frac_MC-frac_t)<<"."<<endl;
  return;
}

void Calorimeter::ToyMC() {
  int NdE = floor(lenght/X0);    //N° bin nel confronto finale tra le derivate
  int bincount=1;
  int equivalent=ceil(X0*N_step/lenght);
  double dEcol[3] = {0};     //Array per l'energia persa in ogni step: dE[0] per ionizzazione, dE[1] sotto Ec nello step attuale, dE[2] step precedente
  double dEdt[NdE] = {0};
  
  
  for(int q=0; q<N_elec; ++q) {
    vecpart.push_back(new Electron());
    vecpart[0]->SetStart(10.373, 0, 0, E0);
    for(int i=0; i<N_step; ++i) {     //Avanti di uno step
      int vSize = vecpart.size();
      for(int j=0; j<vSize; ++j) {    //Scorrimento delle particelle
        double tmp = vecpart[j]->Enrg();
        if(i==bincount*equivalent && tmp > 0 && tmp <= Ec) dEcol[1]+=tmp;
        if(tmp > Ec) {
          double dist = vecpart[j]->MoveIon(step, rho);
	  //Differenziazione in base al tipo di particella
	  if(vecpart[j]->Type() == 'e') {
	    double Ion_loss = 1.6*rho*dist*100;
	    if(bincount <= NdE) dEcol[0] += Ion_loss;
	    if(tmp-Ion_loss > Ec) {
	      double prob = 1 - exp(BremsMU*dist);
	      if(rcal->Uniform() < prob) {       //Condizione di avvenuta Bremsstrahlung
		vecpart.push_back(new Photon());
		bool finish;
		do {
		  finish = this->Brems(vecpart[j], vecpart[vecpart.size()-1]);
		  if(finish==false) vecpart.push_back(new Photon());
		} while(finish==false);
	      }
	    }
	  }
	  else {
	    double prob = 1 - exp(PairMU*dist);
	    if(rcal->Uniform() < prob) {        //Condizione di avvenuta produzione di coppie
	      vecpart.push_back(new Electron());
	      vecpart.push_back(new Electron());
	      this->PairProd(vecpart[j], vecpart[vecpart.size()-2], vecpart[vecpart.size()-1]);
	    }
	  }
	}
      }
      if(i>bincount*equivalent) bincount++;       //Fill dell'energia persa allo step i
      if(i==bincount*equivalent && bincount <= NdE) {
	dEdt[bincount-1] += (dEcol[0]+dEcol[1]-dEcol[2])/(E0*N_elec);
	dEcol[2] = dEcol[1];
	dEcol[0] = 0; dEcol[1] = 0;      //Passaggio di staffetta e azzeramento di E degli step attuali
      }
      if(i==N_step-1) dEcol[2]=0;
    }
    this->Frac_Dep();
    if((q+1)%5==0) cout<<"In esecuzione "<<q+1<<" di "<<N_elec<<endl;
    for(int j=0; j<vecpart.size(); ++j) delete vecpart[j];      //Reset e deallocazione per il prossimo elettrone
    vecpart.clear();
    bincount=1;
  }
  this->Results(dEdt, NdE);
  return;
}

Calorimeter::~Calorimeter() {
  delete rcal, sigmaPair, sigmaBrems, fracDep;
}
