#include "Elements.h"

map<string,Element> ANP;

void ANP_Init() {
  Element Silicon = {"Si", 14/28.0855, 2.329, 0.09370, 40.19};
  Element Copper = {"Cu", 29/63.546, 8.960, 0.01436, 19.42};
  Element BGO = {"BGO", 0.42065, 7.130, 0.01118, 10.50};
  ANP.insert( pair<string,Element>("Si", Silicon) );
  ANP.insert( pair<string,Element>("Cu", Copper) );
  ANP.insert( pair<string,Element>("BGO", BGO) );
  return;
}
