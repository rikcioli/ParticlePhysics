#include <iostream>
#include "Elements.h"
#include "Tracker.h"
#include "Calorimeter.h"
using namespace std;

int main() {
  void ANP_Init();
  ANP_Init();
  Tracker tracker(10, 10000, 300);
  Calorimeter calBGO(0.15, 0.01, 100);
  double ris;
  ris = tracker.Track(2);
  calBGO.ToyMC();
  cout<<"La risoluzione in carica del tracker vale "<<ris<<"."<<endl;
  return 0;
}
