extern TRandom3* gen;

class Particle{
  friend class Calorimeter;
  
 protected:
  double posZ;
  double theta[2];
  double energy;
  
 public:
  double Enrg() {
    return energy;
  }
  void SetStart(double z0, double thx0, double thy0, double E);
  void Move(double delZ);
  double MoveMS(double delZ, double sigma);
  virtual double MoveIon(double delZ, double rho) {
    return 0;
  }
  virtual char Type() {
    return 'p';
  }
  
};
  
class Electron:public Particle{

 public:
  static const double m0;

  Electron()=default;
  Electron(int i);
  virtual double MoveIon(double delZ, double rho) {
    double distance = delZ/(cos(theta[0])*cos(theta[1]));
    energy-=1.6*rho*distance*100;
    posZ+=delZ;
    return distance;
  }
  virtual char Type() {
    return 'e';
  }
  
};


class Photon:public Particle{
  
 public:
  virtual double MoveIon(double delZ, double rho) {
    posZ+=delZ;
    double distance = delZ/(cos(theta[0])*cos(theta[1]));
    return distance;
  }
  virtual char Type() {
    return 'f';
  }
  
};
