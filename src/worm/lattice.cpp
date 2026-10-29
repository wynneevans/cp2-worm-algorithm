#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <ranlxd.h>
#include <format>
#include "lattice.h"
#include "config.h"
#include "input_error.h"

using namespace std;

static inline void loadbar(int x, int n, int w = 50)
{
    if ( (x != n) && (x % (n/100+1) != 0) ) return;
    
    float ratio  =  x/(float)n;
    int   c      =  ratio * w;
    
    cout << setw(3) << (int)(ratio*100) << "% [";
    for (int x=0; x<c; x++) cout << "=";
    for (int x=c; x<w; x++) cout << " ";
    cout << "]\r" << flush;
    if (x == n) {
    cout << endl;
    }
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

inline int modulo(int a, int b) {
    const int result = a % b;
    return result >= 0 ? result : result + b;
}


static double rvec[1];
inline double ran() {ranlxd(rvec,1);return rvec[0];} //Generates random number between 0 and 1
inline int ranspin() {ranlxd(rvec,1);return int(2*rvec[0]);} //Generates random number 0 or 1
inline int ranint(int max) {ranlxd(rvec,1);return int(max*rvec[0]);} //Generates integer random number from 0 to max-1

class Lattice {
  ofstream Mfile_, Cfile_;
  string outDir_, BC_;

  int uYopenPlaqPattE02_[9], dYopenPlaqPattE02_[9], sYopenPlaqPattE02_[9], uYopenPlaqPattE13_[9], dYopenPlaqPattE13_[9], sYopenPlaqPattE13_[9], uYopenPlaqPattO02_[9], dYopenPlaqPattO02_[9], sYopenPlaqPattO02_[9], uYopenPlaqPattO13_[9], dYopenPlaqPattO13_[9], sYopenPlaqPattO13_[9];
  double uYopenProbPattE02_[9], dYopenProbPattE02_[9], sYopenProbPattE02_[9], uYopenProbPattE13_[9], dYopenProbPattE13_[9], sYopenProbPattE13_[9], uYopenProbPattO02_[9], dYopenProbPattO02_[9], sYopenProbPattO02_[9], uYopenProbPattO13_[9], dYopenProbPattO13_[9], sYopenProbPattO13_[9];
  int uPlaqPattE0_[30], dPlaqPattE0_[30], sPlaqPattE0_[30], uPlaqPattE1_[30], dPlaqPattE1_[30], sPlaqPattE1_[30], uPlaqPattE2_[30], dPlaqPattE2_[30], sPlaqPattE2_[30], uPlaqPattE3_[30], dPlaqPattE3_[30], sPlaqPattE3_[30];
  int uPlaqPattO0_[30], dPlaqPattO0_[30], sPlaqPattO0_[30], uPlaqPattO1_[30], dPlaqPattO1_[30], sPlaqPattO1_[30], uPlaqPattO2_[30], dPlaqPattO2_[30], sPlaqPattO2_[30], uPlaqPattO3_[30], dPlaqPattO3_[30], sPlaqPattO3_[30];
  double uProbPattE0_[30], dProbPattE0_[30], sProbPattE0_[30], uProbPattE1_[30], dProbPattE1_[30], sProbPattE1_[30], uProbPattE2_[30], dProbPattE2_[30], sProbPattE2_[30], uProbPattE3_[30], dProbPattE3_[30], sProbPattE3_[30];
  double uProbPattO0_[30], dProbPattO0_[30], sProbPattO0_[30], uProbPattO1_[30], dProbPattO1_[30], sProbPattO1_[30], uProbPattO2_[30], dProbPattO2_[30], sProbPattO2_[30], uProbPattO3_[30], dProbPattO3_[30], sProbPattO3_[30];
  int Nx_, Ny_, Nt_, nTherm_, N_meas_, l_, lp_, lm_, plaq_, pos_, Nplaq_, worm_n_, worm_np1_, worm_start_, xStart_, yStart_, tStart_, startP_, startColor_, originalColor_, wormColor_, q3_, q8_;
  int t_dir_, sim_, t_bob_, x_bob_, y_bob_, w3xAdd_, w8xAdd_;
 
  vector<int> xytp1_, xytm1_, plaq_record_, plaq_fwd_, plaq_bwd_, next_fwd_, next_bwd_, diag_fwd_, diag_bwd_, chiUp_, chiUm_, chiTp_, chiTm_, chiVp_, chiVm_;
  double beta_, dt_, mu3_, mu8_, p1_, p2_, p3_, p4_, p5_, p6_, p7_, p8_, p9_, p10_, p11_, p12_, p13_, p14_, p15_, p16_, p17_, p18_, p19_, p20_, p21_, p22_, p23_, p24_, p25_, p26_, p27_, p28_, p29_, p30_, p31_, p32_, p33_, p34_, p35_, p36_, p37_, p38_, p39_, p40_, p41_, p42_, p43_, p44_, p45_, p46_, p47_, p48_;
  vector<double> cXTt0tp1_, cXTtp1t0_, cXTt0tm1_, cXTtm1t0_, cXTu0up1_, cXTup1u0_, cXTu0um1_, cXTum1u0_, cXTv0vp1_, cXTv0vm1_, cXTvp1v0_, cXTvm1v0_, doubleNull_; 

  int *pntPlqPatE0_, *pntPlqPatE1_, *pntPlqPatE2_, *pntPlqPatE3_, *pntPlqPatO0_, *pntPlqPatO1_, *pntPlqPatO2_, *pntPlqPatO3_;
  double *pntPrbPatE0_, *pntPrbPatE1_, *pntPrbPatE2_, *pntPrbPatE3_, *pntPrbPatO0_, *pntPrbPatO1_, *pntPrbPatO2_, *pntPrbPatO3_;
  int *pntPlqPatE02_, *pntPlqPatE13_, *pntPlqPatO02_, *pntPlqPatO13_;
  double *pntPrbPatE02_, *pntPrbPatE13_, *pntPrbPatO02_, *pntPrbPatO13_;
  double *corr1, *corr2;
  string tail_;

    
public:
  Lattice(string runDir, const string& tag, string BC, int N_meas, double beta, int Nx, int Ny, int Nt, double mu3, double mu8, const std::array<double,48>& p, int sim);
  int volume() {return Nx_*Ny_*Nt_;}
  int xyt(int x, int y, int t) { return x + y*Nx_ + t*Nx_*Ny_; }
  int x(int i) { return i%Nx_; }
  int y(int i) { return (i/Nx_)%Ny_; }
  int t(int i) { return i/(Nx_*Ny_); }
  void set_neighbours();
  void set_plaq_origins();
  void setPatts();
  void check();
  void set_AFM_groundstate();
  void wormUpdate();
  void wormUpdateYopen();
  void decisionStructure(int *, double *, int, int, int);
  void decisionStructureYopen(int *, double *, int, int, int);  
  void measure();
  void dirMove();
  void positionTool(int, int, int, int);
  void timelike();
  void spacelike();
  void diag();
  void setWormColorU();
  void setWormColorD();
  void setWormColorS();
  void bounceTool();
  void selectCorrelators();
  void buildCorrelator(double*, double*);
  [[noreturn, gnu::cold]] void impossible(const char* what);
  bool yopenEdge1() { if( t(worm_n_)%4==3 && ( y(worm_n_)==Ny_-1 || y(worm_n_)==0 ) ) { return true; } else { return false; } }
  bool yopenEdge2() { if( t(worm_n_)%4==0 && ( y(worm_n_)==Ny_-1 || y(worm_n_)==0 ) ) { return true; } else { return false; } }
  bool yopenOddEdge1() { if( ( t(worm_n_)%4==1 && y(worm_n_)==Ny_-1 ) || ( t(worm_n_)%4==3 && y(worm_n_)==0 ) ) { return true; } else { return false; } }
  bool yopenOddEdge2() { if( ( t(worm_n_)%4==2 && y(worm_n_)==Ny_-1 ) || ( t(worm_n_)%4==0 && y(worm_n_)==0 ) ) { return true; } else { return false; } }
  

  private:
  typedef void (Lattice::*fptr)();
  fptr runWorm;
  typedef bool (Lattice::*bptr)();
  bptr fwdEdge, bwdEdge;
};

 
// A state no correct simulation reaches: stop the run (its files are then not published).
void Lattice::impossible(const char* what) {
  throw AlgorithmError(std::format("algorithm reached an impossible state in {} (worm at x={}, y={}, t={}; "
                                   "plaquette {}, pos {}, t_dir {}, worm length {})",
                                   what, x(worm_n_), y(worm_n_), t(worm_n_), plaq_, pos_, t_dir_, l_));
}

Lattice::Lattice(string runDir, const string& tag, string BC, int N_meas, double beta, int Nx, int Ny, int Nt, double mu3, double mu8, const std::array<double,48>& pw, int sim) : outDir_(runDir), BC_(BC), beta_(beta), Nx_(Nx), Ny_(Ny), Nt_(Nt), N_meas_(N_meas), l_(0), lp_(0), lm_(0), plaq_(0), pos_(0), Nplaq_(0), worm_n_(0), worm_np1_(0), worm_start_(0), xStart_(0), yStart_(0), tStart_(0), startP_(0), startColor_(0), originalColor_(0), t_dir_(1), xytp1_(Nx_*Ny_*Nt_), xytm1_(Nx_*Ny_*Nt_), plaq_record_(Nx_*(Ny_+1)*Nt_/2), plaq_fwd_(Nx_*Ny_*Nt_), plaq_bwd_(Nx_*Ny_*Nt_), next_fwd_(Nx_*Ny_*Nt_), next_bwd_(Nx_*Ny_*Nt_), diag_fwd_(Nx_*Ny_*Nt_), diag_bwd_(Nx_*Ny_*Nt_), dt_(4*beta_/Nt_), q3_(0), q8_(0), mu3_(mu3), mu8_(mu8), p1_(pw[0]), p2_(pw[1]), p3_(pw[2]), p4_(pw[3]), p5_(pw[4]), p6_(pw[5]), p7_(pw[6]), p8_(pw[7]), p9_(pw[8]), p10_(pw[9]), p11_(pw[10]), p12_(pw[11]), p13_(pw[12]), p14_(pw[13]), p15_(pw[14]), p16_(pw[15]), p17_(pw[16]), p18_(pw[17]), p19_(pw[18]), p20_(pw[19]), p21_(pw[20]), p22_(pw[21]), p23_(pw[22]), p24_(pw[23]), p25_(pw[24]), p26_(pw[25]), p27_(pw[26]), p28_(pw[27]), p29_(pw[28]), p30_(pw[29]), p31_(pw[30]), p32_(pw[31]), p33_(pw[32]), p34_(pw[33]), p35_(pw[34]), p36_(pw[35]), p37_(pw[36]), p38_(pw[37]), p39_(pw[38]), p40_(pw[39]), p41_(pw[40]), p42_(pw[41]), p43_(pw[42]), p44_(pw[43]), p45_(pw[44]), p46_(pw[45]), p47_(pw[46]), p48_(pw[47]), sim_(sim), cXTt0tp1_(Nx_*Nt_/4), cXTtp1t0_(Nx_*Nt_/4), cXTt0tm1_(Nx_*Nt_/4), cXTtm1t0_(Nx_*Nt_/4), cXTu0up1_(Nx_*Nt_/4), cXTup1u0_(Nx_*Nt_/4), cXTu0um1_(Nx_*Nt_/4), cXTum1u0_(Nx_*Nt_/4), cXTv0vp1_(Nx_*Nt_/4), cXTvp1v0_(Nx_*Nt_/4), cXTv0vm1_(Nx_*Nt_/4), cXTvm1v0_(Nx_*Nt_/4), doubleNull_(Nx_*Nt_/4), t_bob_(0), x_bob_(0), y_bob_(0), tail_("empty"), w3xAdd_(0), w8xAdd_(0)
{	
  Mfile_.open(outDir_ + "/M_" + tag + ".dat", ios::out);
  Cfile_.open(outDir_ + "/C_" + tag + ".dat", ios::out);
  if (!Mfile_ || !Cfile_) throw InputError("cannot write the M_ and C_ files in '" + outDir_ + "'");

  if(BC_=="y-open") {
    cout << "y-open" << endl;
    runWorm = &Lattice::wormUpdateYopen;
    if(Ny_%2==1) {
      cout << "odd" << endl;
      fwdEdge = &Lattice::yopenOddEdge1;
      bwdEdge = &Lattice::yopenOddEdge2;
    }
    else if(Ny_%2==0) {
      cout << "even" << endl;
      fwdEdge = &Lattice::yopenEdge1;
      bwdEdge = &Lattice::yopenEdge2;
    }
    else { impossible("Lattice: Ly is neither even nor odd"); }
  }
  else if(BC_=="closed") {
    cout << "closed" << endl;
    runWorm = &Lattice::wormUpdate;
  }
  else { impossible("Lattice: boundary condition is neither closed nor y-open"); }

  
 }

void Lattice::set_neighbours() {
    int n = 0;
    for (int k = 0; k < Nt_; k++ ) {
            for (int j = 0; j < Ny_; j++) {
                for (int i = 0; i < Nx_; i++) {
	
		  xytp1_[n] = xyt(i,j,(k+1)%Nt_);
		  xytm1_[n] = xyt(i,j,(k-1+Nt_)%Nt_);
                  
		  n++;
            }
        }
    }
}

void Lattice::set_AFM_groundstate() {
 
  if(Ny_%2==0) {
    for (int nPlaq=0; nPlaq < volume()/2; nPlaq++) {
      int n,x,y,t,tSand,tSandBott,excess;
      tSand = nPlaq/(2*Ny_*Nx_); // plaquette number / plaquettes in sandwich = sandwich number
      tSandBott = (2*Ny_*Nx_)*tSand; // sandwich number * number of plaquettes in sandwich = start plaquette of sandwich
      excess = nPlaq - tSandBott; // excess is number of plaquettes higher than start plaquette of sandwich
      if(excess < (Nx_/2)*Ny_) {
	if((excess/(Nx_/2))%2==0) {
	  plaq_record_[nPlaq]=13;
	}
	else {
	  plaq_record_[nPlaq]=113;
	}
      }
      else if(excess < Nx_*Ny_) {
	excess = excess - (Nx_/2)*Ny_;
	if((excess/(Ny_/2))%2==0) {
	  plaq_record_[nPlaq]=13;
	}
	else {
	  plaq_record_[nPlaq]=113;
	}
      }
      else if(excess < Nx_*Ny_ + (Ny_*Nx_/2)) {
	excess = excess - Nx_*Ny_;
	if((excess/(Nx_/2))%2==0) {
	  plaq_record_[nPlaq]=113;
	}
	else {
	  plaq_record_[nPlaq]=13;
	}	
      }
      else if(excess < 2*Nx_*Ny_)  {
	excess = excess - Nx_*Ny_ + (Ny_*Nx_/2);
	if((excess/(Ny_/2))%2==0) {
	  plaq_record_[nPlaq]=113;
	}
	else {
	  plaq_record_[nPlaq]=13;
	}	
      }
      else { impossible("set_AFM_groundstate: plaquette outside its time slice"); }
    }


  }
  else {
    for (int nPlaq=0; nPlaq < ((2*Ny_+1)*Nx_)*(Nt_/4); nPlaq++) {
      int n,x,y,t,tSand,tSandBott,excess;
      tSand = nPlaq/((2*Ny_+1)*Nx_); // plaquette number / plaquettes in sandwich = sandwich number
      tSandBott = ((2*Ny_+1)*Nx_)*tSand; // sandwich number * number of plaquettes in sandwich = start plaquette of sandwich
      excess = nPlaq - tSandBott; // excess is number of plaquettes higher than start plaquette of sandwich
      if(excess < (Nx_/2)*Ny_) {
	if((excess/(Nx_/2))%2==0) {
	  plaq_record_[nPlaq]=13;
	}
	else {
	  plaq_record_[nPlaq]=113;
	}
      }
      else if(excess < (Nx_/2)*Ny_ + ((Ny_+1)*Nx_/2)) {
	excess = excess - (Nx_/2)*Ny_;
	if((excess/((Ny_+1)/2))%2==0) {
	  plaq_record_[nPlaq]=13;
	}
	else {
	  plaq_record_[nPlaq]=113;
	}
      }
      else if(excess < Nx_*Ny_ + ((Ny_+1)*Nx_/2)) {
	excess = excess - (Nx_/2)*Ny_ + ((Ny_+1)*Nx_/2);
	if((excess/(Nx_/2))%2==0) {
	  plaq_record_[nPlaq]=113;
	}
	else {
	  plaq_record_[nPlaq]=13;
	}	
      }
      else if(excess < Nx_*Ny_ + (Ny_+1)*Nx_) {
	excess = excess - Nx_*Ny_ + ((Ny_+1)*Nx_/2);
	if((excess/((Ny_+1)/2))%2==0) {
	  plaq_record_[nPlaq]=113;
	}
	else {
	  plaq_record_[nPlaq]=13;
	}	
      }
      else { impossible("set_AFM_groundstate: plaquette outside its time slice"); }
    }
  }
}

void Lattice::check() {
  cout << endl;
  for (int n = 0; n < volume(); n++) {
    printf("n: %d, plaq_fwd_[n]: %d, plaq_bwd_[n]: %d, plaq_record_[plaq_fwd_[n]]: %d, plaq_record_[plaq_bwd_[n]]: %d.\n", n, plaq_fwd_[n], plaq_bwd_[n], plaq_record_[plaq_fwd_[n]], plaq_record_[plaq_bwd_[n]]) ; 
  }
}

void Lattice::diag() {
  dirMove();
  if (t_dir_==1) {
    worm_n_ = diag_fwd_[worm_n_];
    t_bob_++;
    if(wormColor_==1) {
      if(plaq_ == 1 || plaq_ == 4 || plaq_ == 7 || plaq_ == 101 || plaq_ == 104 || plaq_ == 107) {
	setWormColorD();
      }
      else if(plaq_ == 3 || plaq_ == 6 || plaq_ == 12 || plaq_ == 103 || plaq_ == 106 || plaq_ == 112) {
	setWormColorS();
      }
      else { impossible("diag (forward, colour u): no diagonal move across this plaquette"); }
    }
    else if(wormColor_==-1) {
      if(plaq_ == 1 || plaq_ == 4 || plaq_ == 10 || plaq_ == 101 || plaq_ == 104 || plaq_ == 110) {
	setWormColorU();
      }
      else if(plaq_ == 2 || plaq_ == 5 || plaq_ == 11 || plaq_ == 102 || plaq_ == 105 || plaq_ == 111) {
	setWormColorS();
      }
      else { impossible("diag (forward, colour d): no diagonal move across this plaquette"); }
    }
    else if(wormColor_==0) {
      if(plaq_ == 3 || plaq_ == 6 || plaq_ == 9 || plaq_ == 103 || plaq_ == 106 || plaq_ == 109) {
	setWormColorU();
      }
      else if(plaq_ == 2 || plaq_ == 5 || plaq_ == 8 || plaq_ == 102 || plaq_ == 105 || plaq_ == 108) {
	setWormColorD();
      }
      else { impossible("diag (forward, colour s): no diagonal move across this plaquette"); }
    }
    else { impossible("diag (forward): worm colour is not u, d or s"); }
  }
  else if (t_dir_==-1) {
    worm_n_ = diag_bwd_[worm_n_];
    t_bob_--;
    if(wormColor_==1) {
      if(plaq_ == 1 || plaq_ == 4 || plaq_ == 10 || plaq_ == 101 || plaq_ == 104 || plaq_ == 110) {
	setWormColorD();
      }
      else if(plaq_ == 3 || plaq_ == 6 || plaq_ == 9 || plaq_ == 103 || plaq_ == 106 || plaq_ == 109) {
	setWormColorS();
      }
      else { impossible("diag (backward, colour u): no diagonal move across this plaquette"); }
    }
    else if(wormColor_==-1) {
      if(plaq_ == 1 || plaq_ == 4 || plaq_ == 7 || plaq_ == 101 || plaq_ == 104 || plaq_ == 107) {
	setWormColorU();
      }
      else if(plaq_ == 2 || plaq_ == 5 || plaq_ == 8 || plaq_ == 102 || plaq_ == 105 || plaq_ == 108) {
	setWormColorS();
      }
      else { impossible("diag (backward, colour d): no diagonal move across this plaquette"); }
    }
    else if(wormColor_==0) {
      if(plaq_ == 3 || plaq_ == 6 || plaq_ == 12 || plaq_ == 103 || plaq_ == 106 || plaq_ == 112) {
	setWormColorU();
      }
      else if(plaq_ == 2 || plaq_ == 5 || plaq_ == 11 || plaq_ == 102 || plaq_ == 105 || plaq_ == 111) {
	setWormColorD();
      }
      else { impossible("diag (backward, colour s): no diagonal move across this plaquette"); }
    }
    else { impossible("diag (backward): worm colour is not u, d or s"); }
  }
  else {
    impossible("diag: t_dir is not +1 or -1");
  }
}

void Lattice::timelike() {
  if (t_dir_==1) {
     worm_n_ = xytp1_[worm_n_];
     t_bob_++;
  }
  else if (t_dir_==-1) {
    worm_n_ = xytm1_[worm_n_];
    t_bob_--;
  }
  else {
    impossible("timelike: t_dir is not +1 or -1");
  }
}

void Lattice::spacelike() {
  dirMove();
  if (t_dir_==1) {
    worm_n_ = next_fwd_[worm_n_];
  }
  else if (t_dir_==-1) {
    worm_n_ = next_bwd_[worm_n_];
  }
  else {
    impossible("spacelike: t_dir is not +1 or -1");
  }
  t_dir_=-t_dir_;
}

void Lattice::dirMove() {
  
  int xcoord = x(worm_n_);
  int ycoord = y(worm_n_);
  int tcoord = t(worm_n_);

  if (t_dir_==1) {
    if (tcoord%4==3) {
      if (ycoord%2==0) { y_bob_--; }
      else if (ycoord%2==1) { y_bob_++; }
      else { impossible("dirMove (forward): coordinate is neither even nor odd"); }
    }
    else if (tcoord%4==2) {
      if (xcoord%2==0) { x_bob_--; }
      else if (xcoord%2==1) { x_bob_++; }
      else { impossible("dirMove (forward): coordinate is neither even nor odd"); }
    }
    else if (tcoord%4==1) {       
      if (ycoord%2==0) { y_bob_++; }
      else if (ycoord%2==1) { y_bob_--; } 
      else { impossible("dirMove (forward): coordinate is neither even nor odd"); }
    }
    else if (tcoord%4==0) {
      if (xcoord%2==0) { x_bob_++; }
      else if (xcoord%2==1) { x_bob_--; }
      else { impossible("dirMove (forward): coordinate is neither even nor odd"); }
    }
    else {
      impossible("dirMove (forward): t mod 4 is not 0-3");
    }   
  }
  else {
    if (tcoord%4==3) {
      if (xcoord%2==0) { x_bob_--; }
      else if (xcoord%2==1) { x_bob_++; }
      else { impossible("dirMove (backward): coordinate is neither even nor odd"); }
    }
    else if (tcoord%4==2) {
      if (ycoord%2==0) { y_bob_++; }
      else if (ycoord%2==1) { y_bob_--; } 
      else { impossible("dirMove (backward): coordinate is neither even nor odd"); }
    }
    else if (tcoord%4==1) {       
      if (xcoord%2==0) { x_bob_++; }
      else if (xcoord%2==1) { x_bob_--; }
      else { impossible("dirMove (backward): coordinate is neither even nor odd"); }
    }
    else if (tcoord%4==0) {
      if (ycoord%2==0) { y_bob_--; }
      else if (ycoord%2==1) { y_bob_++; }
      else { impossible("dirMove (backward): coordinate is neither even nor odd"); }
    }
    else {
      impossible("dirMove (backward): t mod 4 is not 0-3");
    }   
  }
}



void Lattice::set_plaq_origins() {
  int n = 0;
    
  for (int k = 0; k < Nt_; k++ ) {
    for (int j = 0; j < Ny_; j++) {
      for (int i = 0; i < Nx_; i++) {
                
	if (k%4==3) {
	  if(Ny_%2==1) {
	    plaq_fwd_[n] = (i*(Ny_+1)/2) + (((j-1+Ny_)%Ny_)/2) + ((Ny_+1)*Nx_/2) + (Ny_*Nx_) + ((k/4)*(2*Ny_+1)*Nx_);
	    plaq_bwd_[n] = (((i-1+Nx_)%Nx_)/2) + (j*Nx_/2) + ((Ny_+1)*Nx_/2) + (Ny_*Nx_/2) + (k/4)*(2*Ny_+1)*Nx_;
	  }
	  else {
	    plaq_fwd_[n] = (i*Ny_/2) + (((j-1+Ny_)%Ny_)/2) + (k*Nx_*Ny_/2); 
	    plaq_bwd_[n] = (((i-1+Nx_)%Nx_)/2) + (j*Nx_/2) + (((k-1+Nt_)%Nt_)*Nx_*Ny_/2); 
	  }
		
	  if (j%2==0) { next_fwd_[n] = xyt(i,(j-1+Ny_)%Ny_,k); }
	  if (j%2==1) { next_fwd_[n] = xyt(i,(j+1)%Ny_,k); } //xyp1t_[n]
	  if (i%2==0) { next_bwd_[n] = xyt((i-1+Nx_)%Nx_,j,k); }
	  if (i%2==1) { next_bwd_[n] = xyt((i+1)%Nx_,j,k); } //xp1yt_[n]
	  if (j%2==0) { diag_fwd_[n] = xyt(i,(j-1+Ny_)%Ny_,(k+1)%Nt_); }
	  if (j%2==1) { diag_fwd_[n] = xyt(i,(j+1)%Ny_,(k+1)%Nt_); }
	  if (i%2==0) { diag_bwd_[n] = xyt((i-1+Nx_)%Nx_,j,(k-1+Nt_)%Nt_); }
	  if (i%2==1) { diag_bwd_[n] = xyt((i+1)%Nx_,j,(k-1+Nt_)%Nt_); }		
	}
                
	if (k%4==2) {
	  if(Ny_%2==1) {
	    plaq_fwd_[n] = (((i-1+Nx_)%Nx_)/2) + (j*Nx_/2) + ((Ny_+1)*Nx_/2) + (Ny_*Nx_/2) + ((k/4)*(2*Ny_+1)*Nx_);
	    plaq_bwd_[n] = (i*(Ny_+1)/2) + (j/2) + (Ny_*Nx_/2) + (k/4)*(2*Ny_+1)*Nx_;
	  }
	  else {
	    plaq_fwd_[n] = (((i-1+Nx_)%Nx_)/2) + (j*Nx_/2) + (k*Nx_*Ny_/2);
	    plaq_bwd_[n] = (i*Ny_/2) + (j/2) + (((k-1+Nt_)%Nt_)*Nx_*Ny_/2);
	  }
	  if (i%2==0) { next_fwd_[n] = xyt((i-1+Nx_)%Nx_,j,k); }
	  if (i%2==1) { next_fwd_[n] = xyt((i+1)%Nx_,j,k); } //xp1yt_[n]
	  if (j%2==0) { next_bwd_[n] = xyt(i,(j+1)%Ny_,k); } //xyp1t_[n]
	  if (j%2==1) { next_bwd_[n] = xyt(i,(j-1+Ny_)%Ny_,k); }
	  if (i%2==0) { diag_fwd_[n] = xyt((i-1+Nx_)%Nx_,j,(k+1)%Nt_); }
	  if (i%2==1) { diag_fwd_[n] = xyt((i+1)%Nx_,j,(k+1)%Nt_); }
	  if (j%2==0) { diag_bwd_[n] = xyt(i,(j+1)%Ny_,(k-1+Nt_)%Nt_); }
	  if (j%2==1) { diag_bwd_[n] = xyt(i,(j-1+Ny_)%Ny_,(k-1+Nt_)%Nt_); }
	}
	      
	if (k%4==1) {
	  if(Ny_%2==1) {
	    plaq_fwd_[n] = (i*(Ny_+1)/2) + (j/2) + (Ny_*Nx_/2) + ((k/4)*(2*Ny_+1)*Nx_);
	    plaq_bwd_[n] = (i/2) + (j*Nx_/2) + (k/4)*(2*Ny_+1)*Nx_;
	  }
	  else {
	    plaq_fwd_[n] = (i*Ny_/2) + (j/2) + (k*Nx_*Ny_/2);
	    plaq_bwd_[n] = (i/2) + (j*Nx_/2) + (((k-1+Nt_)%Nt_)*Nx_*Ny_/2);
	  }
	  if (j%2==0) { next_fwd_[n] = xyt(i,(j+1)%Ny_,k); } //xyp1t_[n]
	  if (j%2==1) { next_fwd_[n] = xyt(i,(j-1+Ny_)%Ny_,k); }
	  if (i%2==0) { next_bwd_[n] = xyt((i+1)%Nx_,j,k); } //xp1yt_[n]
	  if (i%2==1) { next_bwd_[n] = xyt((i-1+Nx_)%Nx_,j,k); }
	  if (j%2==0) { diag_fwd_[n] = xyt(i,(j+1)%Ny_,(k+1)%Nt_); }
	  if (j%2==1) { diag_fwd_[n] = xyt(i,(j-1+Ny_)%Ny_,(k+1)%Nt_); }
	  if (i%2==0) { diag_bwd_[n] = xyt((i+1)%Nx_,j,(k-1+Nt_)%Nt_); }
	  if (i%2==1) { diag_bwd_[n] = xyt((i-1+Nx_)%Nx_,j,(k-1+Nt_)%Nt_); }
	}

	if (k%4==0) {
	  if(Ny_%2==1) {
	    plaq_fwd_[n] = (i/2) + (j*Nx_/2) + (k/4)*(2*Ny_+1)*Nx_;
	    plaq_bwd_[n] = (i*(Ny_+1)/2) + (((j-1+Ny_)%Ny_)/2) + ((Ny_+1)*Nx_/2) + (Ny_*Nx_) + ((((k/4)-1+(Nt_/4))%(Nt_/4))*(2*Ny_+1)*Nx_);
	  }
	  else {
	    plaq_fwd_[n] = (i/2) + (j*Nx_/2) + (k*Nx_*Ny_/2);
	    plaq_bwd_[n] = (i*Ny_/2) + (((j-1+Ny_)%Ny_)/2) + (((k-1+Nt_)%Nt_)*Nx_*Ny_/2);
	  }
	  if (i%2==0) { next_fwd_[n] = xyt((i+1)%Nx_,j,k); } //xp1yt_[n]
	  if (i%2==1) { next_fwd_[n] = xyt((i-1+Nx_)%Nx_,j,k); }
	  if (j%2==0) { next_bwd_[n] = xyt(i,(j-1+Ny_)%Ny_,k); }
	  if (j%2==1) { next_bwd_[n] = xyt(i,(j+1)%Ny_,k); } //xyp1t_[n]
	  if (i%2==0) { diag_fwd_[n] = xyt((i+1)%Nx_,j,(k+1)%Nt_); }
	  if (i%2==1) { diag_fwd_[n] = xyt((i-1+Nx_)%Nx_,j,(k+1)%Nt_); }
	  if (j%2==0) { diag_bwd_[n] = xyt(i,(j-1+Ny_)%Ny_,(k-1+Nt_)%Nt_); }
	  if (j%2==1) { diag_bwd_[n] = xyt(i,(j+1)%Ny_,(k-1+Nt_)%Nt_); }
	}
		
	n++;
                
      }
    }
  }
}

void Lattice::decisionStructure(int *pl, double *pr, int s, int t, int d) { 
  double myRan = ran();
  if (plaq_==13 || plaq_==113) {
     if (myRan < *(pr+24)) {
      timelike();
      pos_=t;
      plaq_record_[Nplaq_]=*(pl+24);
     }
    else if (myRan > *(pr+24) && myRan <= (*(pr+24)+*(pr+25))) {
      spacelike();
      pos_=s;
      plaq_record_[Nplaq_]=*(pl+25);
     }
    else {
      bounceTool();
      t_dir_=-t_dir_;
      }
  }
  else if (plaq_==14 || plaq_==114) {
    if (myRan < *(pr+26)) {
      timelike();
      pos_=t;
      plaq_record_[Nplaq_]=*(pl+26);
     }
    else if (myRan > *(pr+26) && myRan <= (*(pr+26)+*(pr+27))) {
      spacelike();
      pos_=s;
      plaq_record_[Nplaq_]=*(pl+27);
     }
    else {
      bounceTool();
      t_dir_=-t_dir_;
      }
  }
  else if (plaq_==15 || plaq_==115) {
     if (myRan < *(pr+28)) {
      timelike();
      pos_=t;
      plaq_record_[Nplaq_]=*(pl+28);
     }
    else if (myRan > *(pr+28) && myRan <= (*(pr+28)+*(pr+29))) {
      spacelike();
      pos_=s;
      plaq_record_[Nplaq_]=*(pl+29);
     }
    else {
      bounceTool();
      t_dir_=-t_dir_;
     }
  }
  else if (plaq_==7 || plaq_==107) {
    if (myRan < *(pr+13)) {
      spacelike();
      pos_=s;
      plaq_record_[Nplaq_]=*(pl+13);
    }
    else if (myRan > *(pr+13) && myRan <= (*(pr+13)+*(pr+12))) {
      diag();
      pos_=d;
      plaq_record_[Nplaq_]=*(pl+12);
    }
    else {
      bounceTool();
      t_dir_=-t_dir_;
    }
  }
  else if (plaq_==8 || plaq_==108) {
    if (myRan < *(pr+15)) {
      spacelike();
      pos_=s; 
      plaq_record_[Nplaq_]=*(pl+15);
    }
    else if (myRan > *(pr+15) && myRan <= (*(pr+15)+*(pr+14))) {
      diag();
      pos_=d;
      plaq_record_[Nplaq_]=*(pl+14);
    }
    else {
      bounceTool();
      t_dir_=-t_dir_;
    }
  }
  else if (plaq_==9 || plaq_==109) {
    if (myRan < *(pr+17)) {
      spacelike();
      pos_=s; 
      plaq_record_[Nplaq_]=*(pl+17);
    }
    else if (myRan > *(pr+17) && myRan <= (*(pr+17)+*(pr+16))) {
      diag();
      pos_=d;
      plaq_record_[Nplaq_]=*(pl+16);
    }
    else {
      bounceTool();
      t_dir_=-t_dir_;
    }
  }
  else if (plaq_==10 || plaq_==110) {
     if (myRan < *(pr+19)) {
      spacelike();
      pos_=s; 
      plaq_record_[Nplaq_]=*(pl+19);
    }
    else if (myRan > *(pr+19) && myRan <= (*(pr+19)+*(pr+18))) {
      diag();
      pos_=d;
      plaq_record_[Nplaq_]=*(pl+18);
    }
    else {
      bounceTool();
      t_dir_=-t_dir_;
    }
  }
  else if (plaq_==11 || plaq_==111) {
    if (myRan < *(pr+21)) {
      spacelike();
      pos_=s; 
      plaq_record_[Nplaq_]=*(pl+21);
     }
    else if (myRan > *(pr+21) && myRan <= (*(pr+21)+*(pr+20))) {
      diag();
      pos_=d;
      plaq_record_[Nplaq_]=*(pl+20);
    }
    else {
      bounceTool();
      t_dir_=-t_dir_;
     }
  }
  else if (plaq_==12 || plaq_==112) {
    if (myRan < *(pr+23)) {
      spacelike();
      pos_=s; 
      plaq_record_[Nplaq_]=*(pl+23);
     }
    else if (myRan > *(pr+23) && myRan <= (*(pr+23)+*(pr+22))) {
      diag();
      pos_=d;
      plaq_record_[Nplaq_]=*(pl+22);
    }
    else {
      bounceTool();
      t_dir_=-t_dir_;
      }
  }
  else if (plaq_==1 || plaq_==101) {
   if (myRan < *(pr+1)) {
      timelike();
      pos_=t;
      plaq_record_[Nplaq_]=*(pl+1);
    }
    else if (myRan > *(pr+1) && myRan <= (*(pr+1) + *pr)) {
      diag();
      pos_=d;
      plaq_record_[Nplaq_]=*pl;
    }
    else {
      bounceTool();
      t_dir_=-t_dir_;
    }
  }
  else if (plaq_==2 || plaq_==102) {
   if (myRan < *(pr+3)) {
      timelike();
      pos_=t;
      plaq_record_[Nplaq_]=*(pl+3);
    }
    else if (myRan > *(pr+3) && myRan <= (*(pr+3) + *(pr+2))) {
      diag();
      pos_=d;
      plaq_record_[Nplaq_]=*(pl+2);
    }
    else {
      bounceTool();
      t_dir_=-t_dir_;
    }
  }
  else if (plaq_==3 || plaq_==103) {
   if (myRan < *(pr+5)) {
      timelike();
      pos_=t;
      plaq_record_[Nplaq_]=*(pl+5);
    }
    else if (myRan > *(pr+5) && myRan <= (*(pr+5) + *(pr+4))) {
      diag();
      pos_=d;
      plaq_record_[Nplaq_]=*(pl+4);
    }
    else {
      bounceTool();
      t_dir_=-t_dir_;
    }
  }
  else if (plaq_==4 || plaq_==104) {
      if (myRan < *(pr+7)) {
      timelike();
      pos_=t;
      plaq_record_[Nplaq_]=*(pl+7);
    }
    else if (myRan > *(pr+7) && myRan <= (*(pr+7)+*(pr+6))) {
      diag();
      pos_=d;
      plaq_record_[Nplaq_]=*(pl+6);
    }
    else {
      bounceTool();
      t_dir_=-t_dir_;
    }
  }
  else if (plaq_==5 || plaq_==105) {
   if (myRan < *(pr+9)) {
      timelike();
      pos_=t;
      plaq_record_[Nplaq_]=*(pl+9);
    }
    else if (myRan > *(pr+9) && myRan <= (*(pr+9)+*(pr+8))) {
      diag();
      pos_=d;
      plaq_record_[Nplaq_]=*(pl+8);
    }
    else {
      bounceTool();
      t_dir_=-t_dir_;
    }
  }
  else if (plaq_==6 || plaq_==106) {
    if (myRan < *(pr+11)) {
      timelike();
      pos_=t;
      plaq_record_[Nplaq_]=*(pl+11);
    }
    else if (myRan > *(pr+11) && myRan <= (*(pr+11)+*(pr+10))) {
      diag();
      pos_=d;
      plaq_record_[Nplaq_]=*(pl+10);
    }
    else {
      bounceTool();
      t_dir_=-t_dir_;
    }
  }
  else {
    impossible("decisionStructure: plaquette is not 1-15 or 101-115");
  }
} 

void Lattice::decisionStructureYopen(int *pl, double *pr, int s, int t, int d) { 
  double myRan = ran();
  if (plaq_==1 || plaq_==101) {
    if (myRan < *pr) {
      timelike();
      pos_=t;
      plaq_record_[Nplaq_]=*pl;
    }
    else {
      bounceTool();
      t_dir_=-t_dir_;
    }
  }
  else if (plaq_==2 || plaq_==102) {
    if (myRan < *(pr+1)) {
      timelike();
      pos_=t;
      plaq_record_[Nplaq_]=*(pl+1);
    }
    else {
      bounceTool();
      t_dir_=-t_dir_;
    }
  }
  else if (plaq_==3 || plaq_==103) {
    if (myRan < *(pr+2)) {
      timelike();
      pos_=t;
      plaq_record_[Nplaq_]=*(pl+2);
    }
    else {
      bounceTool();
      t_dir_=-t_dir_;
    }
  }
  else if (plaq_==4 || plaq_==104) {
    if (myRan < *(pr+3)) {
      timelike();
      pos_=t;
      plaq_record_[Nplaq_]=*(pl+3);
    }
    else {
      bounceTool();
      t_dir_=-t_dir_;
    }
  }
  else if (plaq_==5 || plaq_==105) {
    if (myRan < *(pr+4)) {
      timelike();
      pos_=t;
      plaq_record_[Nplaq_]=*(pl+4);
    }
    else {
      bounceTool();
      t_dir_=-t_dir_;
    } 
  }
  else if (plaq_==6 || plaq_==106) {
    if (myRan < *(pr+5)) {
      timelike();
      pos_=t;
      plaq_record_[Nplaq_]=*(pl+5);
    }
    else {
      bounceTool();
      t_dir_=-t_dir_;
    }
  }
  else if (plaq_==13 || plaq_==113) {
    if (myRan < *(pr+6)) {
      timelike();
      pos_=t;
      plaq_record_[Nplaq_]=*(pl+6);
    }
    else {
      bounceTool();
      t_dir_=-t_dir_;
    }
  }
  else if (plaq_==14 || plaq_==114) {
    if (myRan < *(pr+7)) {
      timelike();
      pos_=t;
      plaq_record_[Nplaq_]=*(pl+7);
    }
    else {
      bounceTool();
      t_dir_=-t_dir_;
    }
  }
  else if (plaq_==15 || plaq_==115) {
    if (myRan < *(pr+8)) {
      timelike();
      pos_=t;
      plaq_record_[Nplaq_]=*(pl+8);
    }
    else {
      bounceTool();
      t_dir_=-t_dir_;
    }
  }
  else {
    impossible("decisionStructureYopen: plaquette is not 1-6, 13-15, 101-106 or 113-115");
  }
}

void Lattice::bounceTool() {
  if(plaq_==13 || plaq_==113) {
    setWormColorU();
  }
  else if(plaq_==14 || plaq_==114) {
    setWormColorD();
  }
  else if(plaq_==15 || plaq_==115) {
    setWormColorS();
  }
  else if(t_dir_==1) {
    if(plaq_==9 || plaq_==109 || plaq_==10 || plaq_==110) {
      setWormColorU();
    }
    else if(plaq_==7 || plaq_==107 || plaq_==8 || plaq_==108) {
      setWormColorD();
    }
    else if(plaq_==11 || plaq_==111 || plaq_==12 || plaq_==112) {
      setWormColorS();
    }
    else if(pos_==0 || pos_==2 || pos_==3 || pos_==4 ) { // 0
      if(plaq_==1 || plaq_==101 || plaq_==6 || plaq_==106) {
	setWormColorU();
      }
      else if(plaq_==2 || plaq_==102 || plaq_==4 || plaq_==104) {
	setWormColorD();
      }
      else if(plaq_==3 || plaq_==103 || plaq_==5 || plaq_==105) {
	setWormColorS();
      }
    }
    else if(pos_==1 || pos_==5 || pos_==6 || pos_==7 || pos_==9 || pos_==11 || pos_==14 || pos_==15 ) { // 1
      if(plaq_==3 || plaq_==103 || plaq_==4 || plaq_==104) {
	setWormColorU();
      }
      else if(plaq_==1 || plaq_==101 || plaq_==5 || plaq_==105) {
	setWormColorD();
      }
      else if(plaq_==2 || plaq_==102 || plaq_==6 || plaq_==106) {
	setWormColorS();
      }
    }
    else { impossible("bounceTool (forward): no colour for this plaquette and position"); }
  }  
  else if(t_dir_==-1) {
     if(plaq_==7 || plaq_==107 || plaq_==12 || plaq_==112) {
      setWormColorU();
    }
    else if(plaq_==10 || plaq_==110 || plaq_==11 || plaq_==111) {
      setWormColorD();
    }
    else if(plaq_==8 || plaq_==108 || plaq_==9 || plaq_==109) {
      setWormColorS();
    }            
    else if(pos_==2 || pos_==4 || pos_==5 || pos_==6 || pos_==8 || pos_==9 || pos_==12 || pos_==14 ) { // 2
      if(plaq_==1 || plaq_==101 || plaq_==6 || plaq_==106) {
	setWormColorU();
      }
      else if(plaq_==2 || plaq_==102 || plaq_==4 || plaq_==104) {
	setWormColorD();
      }
      else if(plaq_==3 || plaq_==103 || plaq_==5 || plaq_==105) {
	setWormColorS();
      }
    }             
    else if(pos_==0 || pos_==1 || pos_==3 || pos_==7 || pos_==10 || pos_==11 || pos_==13 || pos_==15 ) { // 3
      if(plaq_==3 || plaq_==103 || plaq_==4 || plaq_==104) {
	setWormColorU();
      }
      else if(plaq_==1 || plaq_==101 || plaq_==5 || plaq_==105) {
	setWormColorD();
      }
      else if(plaq_==2 || plaq_==102 || plaq_==6 || plaq_==106) {
	setWormColorS();
      }
    }
    else { impossible("bounceTool (backward): no colour for this plaquette and position"); }
  }
}

void Lattice::positionTool(int a, int b, int c, int d) {
  if (y(worm_n_)%2==0) { 
    if (x(worm_n_)%2==0) {
      pos_=a;
    }
    else if (x(worm_n_)%2==1) {
      pos_=b;
    }
    else {
        impossible("positionTool: x or y is neither even nor odd");
    }
  }
  else if (y(worm_n_)%2==1) { 
    if (x(worm_n_)%2==0) {
      pos_=c; 
    }
    else if (x(worm_n_)%2==1) { 
      pos_=d;
    }
    else {
      impossible("positionTool: x or y is neither even nor odd");
    } 
  }
  else { impossible("positionTool: x or y is neither even nor odd"); }
 
  int startPlaq;
  if(t_dir_==1) {
    startPlaq=plaq_record_[plaq_fwd_[worm_n_]];
    if(startPlaq < 100) { // must be on even plaquette
      if(pos_==0 || pos_==2) { // on even and 0
	if(startPlaq==1 || startPlaq==6 || startPlaq==9 || startPlaq==10 || startPlaq==13) {
	  originalColor_=1;
	  if(ran() < 0.5) {
	    setWormColorD();
	    tail_="tm";
	    w3xAdd_=-2;
	    w8xAdd_=0;
	    
	  }
	  else {
	    setWormColorS();
	    tail_="vm";
	    w3xAdd_=-1;
	    w8xAdd_=-1;
	  }
	}
	else if(startPlaq==2 || startPlaq==4 || startPlaq==7 || startPlaq==8 || startPlaq==14) {
	  originalColor_=-1;
	  if(ran() < 0.5) {
	    setWormColorS();
	    tail_="um";
	    w3xAdd_=1;
	    w8xAdd_=-1;
	  }
	  else {
	    setWormColorU();
	    tail_="tp";
	    w3xAdd_=2;
	    w8xAdd_=0;
	  }
	}
	else if(startPlaq==3 || startPlaq==5 || startPlaq==11 || startPlaq==12 || startPlaq==15) {
	  originalColor_=0;
	  if(ran() < 0.5) {
	    setWormColorU();
	    tail_="vp";
	    w3xAdd_=1;
	    w8xAdd_=1;
	  }
	  else {
	    setWormColorD();  
	    tail_="up";
	    w3xAdd_=-1;
	    w8xAdd_=1;
	  }
	}
	else { impossible("positionTool: start plaquette not possible at this position"); }
      }
      else if(pos_==1 || pos_==6) { // on even and 1
	if(startPlaq==3 || startPlaq==4 || startPlaq==9 || startPlaq==10 || startPlaq==13) {
	  originalColor_=11;
	  if(ran() < 0.5) {
	    setWormColorD();
	    tail_="tp";
	      w3xAdd_=2; w8xAdd_=0;
	   
	      }
	  else {
	    setWormColorS(); 
	    tail_="vp";
	    w3xAdd_=1;
	    w8xAdd_=1;
	  }
	}
	else if(startPlaq==1 || startPlaq==5 || startPlaq==7 || startPlaq==8 || startPlaq==14) {
	  originalColor_=-11;
	  if(ran() < 0.5) {
	    setWormColorS();
	    tail_="up";
	    w3xAdd_=-1;
	    w8xAdd_=1;
	  }
	  else {
	    setWormColorU(); 
	    tail_="tm";
	    w3xAdd_=-2;
	    w8xAdd_=0;
	  }
	}
	else if(startPlaq==2 || startPlaq==6 || startPlaq==11 || startPlaq==12 || startPlaq==15) {
	  originalColor_=10;
	  if(ran() < 0.5) {
	    setWormColorU();
	    tail_="vm";
	    w3xAdd_=-1;
	    w8xAdd_=-1;
	  }
	  else {
	    setWormColorD(); 
	    tail_="um";
	    w3xAdd_=1;
	    w8xAdd_=-1;
	  }
	}
	else { impossible("positionTool: start plaquette not possible at this position"); }
      }
      else { impossible("positionTool: position not possible on a plaquette of this parity"); }
    }
    else if(startPlaq > 100) {
      if(pos_==3 || pos_==4) { // on odd and 0
	if(startPlaq==101 || startPlaq==106 || startPlaq==109 || startPlaq==110 || startPlaq==113) {
	  originalColor_=11;
	  if(ran() < 0.5) {
	    setWormColorD();
	    tail_="tp";
	    w3xAdd_=2;
	    w8xAdd_=0;
	  }
	  else {
	    setWormColorS();   
	    tail_="vp";
	    w3xAdd_=1;
	    w8xAdd_=1;
	  }
	}
	else if(startPlaq==102 || startPlaq==104 || startPlaq==107 || startPlaq==108 || startPlaq==114) {
	  originalColor_=-11;
	  if(ran() < 0.5) {
	    setWormColorS();
	    tail_="up";
	    w3xAdd_=-1;
	    w8xAdd_=1;
	  }
	  else {
	    setWormColorU();
	    tail_="tm";
	    w3xAdd_=-2;
	    w8xAdd_=0;
	  }
	}
	else if(startPlaq==103 || startPlaq==105 || startPlaq==111 || startPlaq==112 || startPlaq==115) {
	  originalColor_=10;
	  if(ran() < 0.5) {
	    setWormColorU();
	    tail_="vm";
	    w3xAdd_=-1;
	    w8xAdd_=-1;
	  }
	  else {
	    setWormColorD();  
	    tail_="um";
	    w3xAdd_=1;
	    w8xAdd_=-1;
	  }
	}
	else { impossible("positionTool: start plaquette not possible at this position"); }
      }
      else if(pos_==5 || pos_==7) { // on odd and 1
	if(startPlaq==103 || startPlaq==104 || startPlaq==109 || startPlaq==110 || startPlaq==113) {
	  originalColor_=1;
	  if(ran() < 0.5) {
	    setWormColorD();
	    tail_="tm";
	    w3xAdd_=-2;
	    w8xAdd_=0;
	  }
	  else {
	    setWormColorS();    
	    tail_="vm";
	    w3xAdd_=-1;
	    w8xAdd_=-1;
	  }
	}
	else if(startPlaq==101 || startPlaq==105 || startPlaq==107 || startPlaq==108 || startPlaq==114) {
	  originalColor_=-1;
	  if(ran() < 0.5) {
	    setWormColorS();
	    tail_="um";
	    w3xAdd_=1;
	    w8xAdd_=-1;
	  }
	  else {
	    setWormColorU(); 
	    tail_="tp";
	    w3xAdd_=2;
	    w8xAdd_=0;
	  }
	}
	else if(startPlaq==102 || startPlaq==106 || startPlaq==111 || startPlaq==112 || startPlaq==115) {
	  originalColor_=0;
	  if(ran() < 0.5) {
	    setWormColorU();
	    tail_="vp";
	    w3xAdd_=1;
	    w8xAdd_=1;
	  }
	  else {
	    setWormColorD(); 
	    tail_="up";
	    w3xAdd_=-1;
	    w8xAdd_=1;
	  }
	}
	else { impossible("positionTool: start plaquette not possible at this position"); }
      }
      else { impossible("positionTool: position not possible on a plaquette of this parity"); }
    }
    else { impossible("positionTool: start plaquette is 100"); }    
  }
  else if(t_dir_==-1) {
    startPlaq=plaq_record_[plaq_bwd_[worm_n_]];
    if(startPlaq < 100) { // must be on even plaquette
      if(pos_==2 || pos_==5) { // on even and 2
	if(startPlaq==1 || startPlaq==6 || startPlaq==7 || startPlaq==12 || startPlaq==13) {
	  originalColor_=1;
	  if(ran() < 0.5) {
	    setWormColorD();
	    tail_="tp";
	    w3xAdd_=2;
	    w8xAdd_=0;
	  }
	  else {
	    setWormColorS();  
	    tail_="vp";
	    w3xAdd_=1;
	    w8xAdd_=1;
	  }
	}
	else if(startPlaq==2 || startPlaq==4 || startPlaq==10 || startPlaq==11 || startPlaq==14) {
	  originalColor_=-1;
	  if(ran() < 0.5) {
	    setWormColorS();
	    tail_="up";
	    w3xAdd_=-1;
	    w8xAdd_=1;
	  }
	  else {
	    setWormColorU(); 
	    tail_="tm";
	    w3xAdd_=-2;
	    w8xAdd_=0;
	  }
	}
	else if(startPlaq==3 || startPlaq==5 || startPlaq==8 || startPlaq==9 || startPlaq==15) {
	  originalColor_=0;
	  if(ran() < 0.5) {
	    setWormColorU();
	    tail_="vm";
	    w3xAdd_=-1;
	    w8xAdd_=-1;
	  }
	  else {
	    setWormColorD();  
	    tail_="um";
	    w3xAdd_=1;
	    w8xAdd_=-1;
	  }
	}
	else { impossible("positionTool: start plaquette not possible at this position"); }
      }
      else if(pos_==1 || pos_==3) { // on even and 3
	if(startPlaq==3 || startPlaq==4 || startPlaq==7 || startPlaq==12 || startPlaq==13) {
	  originalColor_=11;
	  if(ran() < 0.5) {
	    setWormColorD();
	    tail_="tm";
	    w3xAdd_=-2;
	    w8xAdd_=0;
	  }
	  else {
	    setWormColorS(); 
	    tail_="vm";
	    w3xAdd_=-1;
	    w8xAdd_=-1;
	  }
	}
	else if(startPlaq==1 || startPlaq==5 || startPlaq==10 || startPlaq==11 || startPlaq==14) {
	  originalColor_=-11;
	  if(ran() < 0.5) {
	    setWormColorS();
	    tail_="um";
	    w3xAdd_=1;
	    w8xAdd_=-1;
	  }
	  else {
	    setWormColorU();   
	    tail_="tp";
	    w3xAdd_=2;
	    w8xAdd_=0;
	  }
	}
	else if(startPlaq==2 || startPlaq==6 || startPlaq==8 || startPlaq==9 || startPlaq==15) {
	  originalColor_=10;
	  if(ran() < 0.5) {
	    setWormColorU();
	    tail_="vp";
	    w3xAdd_=1;
	    w8xAdd_=1;
	  }
	  else {
	    setWormColorD();
	    tail_="up";
	    w3xAdd_=-1;
	    w8xAdd_=1;
	  }
	}
	else { impossible("positionTool: start plaquette not possible at this position"); }
      }
      else { impossible("positionTool: position not possible on a plaquette of this parity"); }
    }
    else if(startPlaq > 100) {
      if(pos_==4 || pos_==6) { // on odd and 2
	if(startPlaq==101 || startPlaq==106 || startPlaq==107 || startPlaq==112 || startPlaq==113) {
	  originalColor_=11;
	  if(ran() < 0.5) {
	    setWormColorD();
	    tail_="tm";
	    w3xAdd_=-2;
	    w8xAdd_=0;
	  }
	  else {
	    setWormColorS(); 
	    tail_="vm";
	    w3xAdd_=-1;
	    w8xAdd_=-1;
	  }
	}
	else if(startPlaq==102 || startPlaq==104 || startPlaq==110 || startPlaq==111 || startPlaq==114) {
	  originalColor_=-11;
	  if(ran() < 0.5) {
	    setWormColorS();
	    tail_="um";
	    w3xAdd_=1;
	    w8xAdd_=-1;
	  }
	  else {
	    setWormColorU();  
	    tail_="tp";
	    w3xAdd_=2;
	    w8xAdd_=0;
	  }
	}
	else if(startPlaq==103 || startPlaq==105 || startPlaq==108 || startPlaq==109 || startPlaq==115) {
	  originalColor_=10;
	  if(ran() < 0.5) {
	    setWormColorU();
	    tail_="vp";
	    w3xAdd_=1;
	    w8xAdd_=1;
	  }
	  else {
	    setWormColorD();  
	    tail_="up";
	    w3xAdd_=-1;
	    w8xAdd_=1;
	  }
	}
	else { impossible("positionTool: start plaquette not possible at this position"); }
      }
      else if(pos_==0 || pos_==7) { // on odd and 3
	if(startPlaq==103 || startPlaq==104 || startPlaq==107 || startPlaq==112 || startPlaq==113) {
	  originalColor_=1;
	  if(ran() < 0.5) {
	    setWormColorD();
	    tail_="tp";
	    w3xAdd_=2;
	    w8xAdd_=0;
	  }
	  else {
	    setWormColorS(); 
	    tail_="vp";
	    w3xAdd_=1;
	    w8xAdd_=1;
	  }
	}
	else if(startPlaq==101 || startPlaq==105 || startPlaq==110 || startPlaq==111 || startPlaq==114) {
	  originalColor_=-1;
	  if(ran() < 0.5) {
	    setWormColorS();
	    tail_="up";
	    w3xAdd_=-1;
	    w8xAdd_=1;
	  }
	  else {
	    setWormColorU();  
	    tail_="tm";
	    w3xAdd_=-2;
	    w8xAdd_=0;
	  }
	}
	else if(startPlaq==102 || startPlaq==106 || startPlaq==108 || startPlaq==109 || startPlaq==115) {
	  originalColor_=0;
	  if(ran() < 0.5) {
	    setWormColorU();
	    tail_="vm";
	    w3xAdd_=-1;
	    w8xAdd_=-1;
	  }
	  else {
	    setWormColorD(); 
	    tail_="um";
	    w3xAdd_=1;
	    w8xAdd_=-1;
	  }
	}
	else { impossible("positionTool: start plaquette not possible at this position"); }
      }
      else { impossible("positionTool: position not possible on a plaquette of this parity"); }
    }
    else { impossible("positionTool: start plaquette is 100"); }   
  }
  else { impossible("positionTool: t_dir is not +1 or -1"); }
 
}

void Lattice::setWormColorU() {
  wormColor_=1;

  pntPlqPatE0_=uPlaqPattE0_;
  pntPlqPatE1_=uPlaqPattE1_;
  pntPlqPatE2_=uPlaqPattE2_;
  pntPlqPatE3_=uPlaqPattE3_;
  pntPlqPatO0_=uPlaqPattO0_;
  pntPlqPatO1_=uPlaqPattO1_;
  pntPlqPatO2_=uPlaqPattO2_;
  pntPlqPatO3_=uPlaqPattO3_;
  
  pntPrbPatE0_=uProbPattE0_;
  pntPrbPatE1_=uProbPattE1_;
  pntPrbPatE2_=uProbPattE2_;
  pntPrbPatE3_=uProbPattE3_;
  pntPrbPatO0_=uProbPattO0_;
  pntPrbPatO1_=uProbPattO1_;
  pntPrbPatO2_=uProbPattO2_;
  pntPrbPatO3_=uProbPattO3_;

  pntPlqPatE02_=uYopenPlaqPattE02_;
  pntPlqPatE13_=uYopenPlaqPattE13_;
  pntPlqPatO02_=uYopenPlaqPattO02_;
  pntPlqPatO13_=uYopenPlaqPattO13_;

  pntPrbPatE02_=uYopenProbPattE02_;
  pntPrbPatE13_=uYopenProbPattE13_;
  pntPrbPatO02_=uYopenProbPattO02_;
  pntPrbPatO13_=uYopenProbPattO13_;

}

void Lattice::setWormColorD() {
  wormColor_ = -1;

  pntPlqPatE0_=dPlaqPattE0_;
  pntPlqPatE1_=dPlaqPattE1_;
  pntPlqPatE2_=dPlaqPattE2_;
  pntPlqPatE3_=dPlaqPattE3_;
  pntPlqPatO0_=dPlaqPattO0_;
  pntPlqPatO1_=dPlaqPattO1_;
  pntPlqPatO2_=dPlaqPattO2_;
  pntPlqPatO3_=dPlaqPattO3_;
  
  pntPrbPatE0_=dProbPattE0_;
  pntPrbPatE1_=dProbPattE1_;
  pntPrbPatE2_=dProbPattE2_;
  pntPrbPatE3_=dProbPattE3_;
  pntPrbPatO0_=dProbPattO0_;
  pntPrbPatO1_=dProbPattO1_;
  pntPrbPatO2_=dProbPattO2_;
  pntPrbPatO3_=dProbPattO3_;

  pntPlqPatE02_=dYopenPlaqPattE02_;
  pntPlqPatE13_=dYopenPlaqPattE13_;
  pntPlqPatO02_=dYopenPlaqPattO02_;
  pntPlqPatO13_=dYopenPlaqPattO13_;

  pntPrbPatE02_=dYopenProbPattE02_;
  pntPrbPatE13_=dYopenProbPattE13_;
  pntPrbPatO02_=dYopenProbPattO02_;
  pntPrbPatO13_=dYopenProbPattO13_;

}

void Lattice::setWormColorS() {
  wormColor_ = 0;

  pntPlqPatE0_=sPlaqPattE0_;
  pntPlqPatE1_=sPlaqPattE1_;
  pntPlqPatE2_=sPlaqPattE2_;
  pntPlqPatE3_=sPlaqPattE3_;
  pntPlqPatO0_=sPlaqPattO0_;
  pntPlqPatO1_=sPlaqPattO1_;
  pntPlqPatO2_=sPlaqPattO2_;
  pntPlqPatO3_=sPlaqPattO3_;
  
  pntPrbPatE0_=sProbPattE0_;
  pntPrbPatE1_=sProbPattE1_;
  pntPrbPatE2_=sProbPattE2_;
  pntPrbPatE3_=sProbPattE3_;
  pntPrbPatO0_=sProbPattO0_;
  pntPrbPatO1_=sProbPattO1_;
  pntPrbPatO2_=sProbPattO2_;
  pntPrbPatO3_=sProbPattO3_;

  pntPlqPatE02_=sYopenPlaqPattE02_;
  pntPlqPatE13_=sYopenPlaqPattE13_;
  pntPlqPatO02_=sYopenPlaqPattO02_;
  pntPlqPatO13_=sYopenPlaqPattO13_;

  pntPrbPatE02_=sYopenProbPattE02_;
  pntPrbPatE13_=sYopenProbPattE13_;
  pntPrbPatO02_=sYopenProbPattO02_;
  pntPrbPatO13_=sYopenProbPattO13_;

}

void Lattice::setPatts() {
  int uPlaqPattE0Temp[30]={0,0,0,6,9,13,10,13,0,1,0,0,1,13,0,9,0,0,0,0,0,10,6,13,0,0,1,10,6,9};
  int dPlaqPattE0Temp[30]={7,14,0,0,0,4,0,0,8,14,0,2,0,0,0,0,0,8,4,14,2,14,0,7,4,7,0,0,2,8};
  int sPlaqPattE0Temp[30]={0,5,11,15,0,0,0,3,0,0,12,15,0,12,5,15,3,15,0,11,0,0,0,0,3,12,5,11,0,0};

  int uPlaqPattE1Temp[30]={10,13,0,4,0,0,0,0,0,3,9,13,4,13,0,9,0,0,0,0,0,10,3,13,0,0,4,10,3,9};
  int dPlaqPattE1Temp[30]={0,0,8,14,0,5,7,14,0,0,0,1,0,0,0,0,0,8,1,14,5,14,0,7,1,7,0,0,5,8};
  int sPlaqPattE1Temp[30]={0,6,0,0,12,15,0,2,11,15,0,0,0,12,2,15,6,15,0,11,0,0,0,0,6,12,2,11,0,0};

  int uPlaqPattE2Temp[30]={0,0,0,6,12,13,7,13,0,1,0,0,0,0,0,7,6,13,1,13,0,12,0,0,0,0,1,7,6,12};
  int dPlaqPattE2Temp[30]={10,14,0,0,0,4,0,0,11,14,0,2,4,14,2,14,0,10,0,0,0,0,0,11,4,10,0,0,2,11};
  int sPlaqPattE2Temp[30]={0,5,8,15,0,0,0,3,0,0,9,15,0,8,0,0,0,0,0,9,5,15,3,15,3,9,5,8,0,0};

  int uPlaqPattE3Temp[30]={7,13,0,4,0,0,0,0,0,3,12,13,0,0,0,7,3,13,4,13,0,12,0,0,0,0,4,7,3,12};
  int dPlaqPattE3Temp[30]={0,0,11,14,0,5,10,14,0,0,0,1,1,14,5,14,0,10,0,0,0,0,0,11,1,10,0,0,5,11};
  int sPlaqPattE3Temp[30]={0,6,0,0,9,15,0,2,8,15,0,0,0,8,0,0,0,0,0,9,2,15,6,15,6,9,2,8,0,0};

  int uPlaqPattO0Temp[30]={100,100,100,106,109,113,110,113,100,101,100,100,101,113,100,109,100,100,100,100,100,110,106,113,100,100,101,110,106,109};
  int dPlaqPattO0Temp[30]={107,114,100,100,100,104,100,100,108,114,100,102,100,100,100,100,100,108,104,114,102,114,100,107,104,107,100,100,102,108};
  int sPlaqPattO0Temp[30]={100,105,111,115,100,100,100,103,100,100,112,115,100,112,105,115,103,115,100,111,100,100,100,100,103,112,105,111,100,100};

  int uPlaqPattO1Temp[30]={110,113,100,104,100,100,100,100,100,103,109,113,104,113,100,109,100,100,100,100,100,110,103,113,100,100,104,110,103,109};
  int dPlaqPattO1Temp[30]={100,100,108,114,100,105,107,114,100,100,100,101,100,100,100,100,100,108,101,114,105,114,100,107,101,107,100,100,105,108};
  int sPlaqPattO1Temp[30]={100,106,100,100,112,115,100,102,111,115,100,100,100,112,102,115,106,115,100,111,100,100,100,100,106,112,102,111,100,100};

  int uPlaqPattO2Temp[30]={100,100,100,106,112,113,107,113,100,101,100,100,100,100,100,107,106,113,101,113,100,112,100,100,100,100,101,107,106,112};
  int dPlaqPattO2Temp[30]={110,114,100,100,100,104,100,100,111,114,100,102,104,114,102,114,100,110,100,100,100,100,100,111,104,110,100,100,102,111};
  int sPlaqPattO2Temp[30]={100,105,108,115,100,100,100,103,100,100,109,115,100,108,100,100,100,100,100,109,105,115,103,115,103,109,105,108,100,100};

  int uPlaqPattO3Temp[30]={107,113,100,104,100,100,100,100,100,103,112,113,100,100,100,107,103,113,104,113,100,112,100,100,100,100,104,107,103,112};
  int dPlaqPattO3Temp[30]={100,100,111,114,100,105,110,114,100,100,100,101,101,114,105,114,100,110,100,100,100,100,100,111,101,110,100,100,105,111};
  int sPlaqPattO3Temp[30]={100,106,100,100,109,115,100,102,108,115,100,100,100,108,100,100,100,100,100,109,102,115,106,115,106,109,102,108,100,100};

  int uYopenPlaqPattTempE02[9]={0,6,13,13,1,0,0,1,6};
  int dYopenPlaqPattTempE02[9]={14,0,4,0,14,2,4,0,2};
  int sYopenPlaqPattTempE02[9]={5,15,0,3,0,15,3,5,0};

  int uYopenPlaqPattTempE13[9]={13,4,0,0,3,13,0,4,3};
  int dYopenPlaqPattTempE13[9]={0,14,5,14,0,1,1,0,5};
  int sYopenPlaqPattTempE13[9]={6,0,15,2,15,0,6,2,0};

  int uYopenPlaqPattTempO02[9]={0,106,113,113,101,0,0,101,106};
  int dYopenPlaqPattTempO02[9]={114,0,104,0,114,102,104,0,102};
  int sYopenPlaqPattTempO02[9]={105,115,0,103,0,115,103,105,0};

  int uYopenPlaqPattTempO13[9]={113,104,0,0,103,113,0,104,103};
  int dYopenPlaqPattTempO13[9]={0,114,105,114,0,101,101,0,105};
  int sYopenPlaqPattTempO13[9]={106,0,115,102,115,0,106,102,0};

  double uYopenProbPattTempE02[9]={0,p44_,p47_,p44_,p47_,0,0,p44_,p47_};
  double dYopenProbPattTempE02[9]={p43_,0,p46_,0,p46_,p43_,p43_,0,p46_};
  double sYopenProbPattTempE02[9]={p48_,p45_,0,p45_,0,p48_,p48_,p45_,0};

  double uYopenProbPattTempE13[9]={p43_,p48_,0,0,p43_,p48_,0,p43_,p48_};
  double dYopenProbPattTempE13[9]={0,p45_,p44_,p44_,0,p45_,p44_,0,p45_};
  double sYopenProbPattTempE13[9]={p46_,0,p47_,p47_,p46_,0,p47_,p46_,0};

  double uYopenProbPattTempO02[9]={0,p43_,p48_,p43_,p48_,0,0,p43_,p48_};
  double dYopenProbPattTempO02[9]={p44_,0,p45_,0,p45_,p44_,p44_,0,p45_};
  double sYopenProbPattTempO02[9]={p47_,p46_,0,p46_,0,p47_,p47_,p46_,0};

  double uYopenProbPattTempO13[9]={p44_,p47_,0,0,p44_,p47_,0,p44_,p47_};
  double dYopenProbPattTempO13[9]={0,p46_,p43_,p43_,0,p46_,p43_,0,p46_};
  double sYopenProbPattTempO13[9]={p45_,0,p48_,p48_,p45_,0,p48_,p45_,0};

  double uProbPattE0Temp[30]={0,0,0,p39_,p29_,p30_,p11_,p12_,0,p38_,0,0,p2_,p1_,0,1,0,0,0,0,0,1,p32_,p31_,0,0,p4_,p3_,p34_,p33_};
  double dProbPattE0Temp[30]={p5_,p6_,0,0,0,p41_,0,0,p23_,p24_,0,p40_,0,0,0,0,0,1,p8_,p7_,p14_,p13_,0,1,p10_,p9_,0,0,p16_,p15_};
  double sProbPattE0Temp[30]={0,p37_,p17_,p18_,0,0,0,p42_,0,0,p35_,p36_,0,1,p20_,p19_,p26_,p25_,0,1,0,0,0,0,p28_,p27_,p22_,p21_,0,0};

  double uProbPattE1Temp[30]={p5_,p6_,0,p37_,0,0,0,0,0,p40_,p35_,p36_,p8_,p7_,0,1,0,0,0,0,0,1,p26_,p25_,0,0,p10_,p9_,p28_,p27_};
  double dProbPattE1Temp[30]={0,0,p17_,p18_,0,p39_,p11_,p12_,0,0,0,p42_,0,0,0,0,0,1,p2_,p1_,p20_,p19_,0,1,p4_,p3_,0,0,p22_,p21_};
  double sProbPattE1Temp[30]={0,p41_,0,0,p29_,p30_,0,p38_,p23_,p24_,0,0,0,1,p14_,p13_,p32_,p31_,0,1,0,0,0,0,p34_,p33_,p16_,p15_,0,0};

  double uProbPattE2Temp[30]={0,0,0,p39_,p29_,p30_,p11_,p12_,0,p38_,0,0,0,0,0,1,p32_,p31_,p2_,p1_,0,1,0,0,0,0,p4_,p3_,p34_,p33_};
  double dProbPattE2Temp[30]={p5_,p6_,0,0,0,p41_,0,0,p23_,p24_,0,p40_,p8_,p7_,p14_,p13_,0,1,0,0,0,0,0,1,p10_,p9_,0,0,p16_,p15_};
  double sProbPattE2Temp[30]={0,p37_,p17_,p18_,0,0,0,p42_,0,0,p35_,p36_,0,1,0,0,0,0,0,1,p20_,p19_,p26_,p25_,p28_,p27_,p22_,p21_,0,0};

  double uProbPattE3Temp[30]={p5_,p6_,0,p37_,0,0,0,0,0,p40_,p35_,p36_,0,0,0,1,p26_,p25_,p8_,p7_,0,1,0,0,0,0,p10_,p9_,p28_,p27_};
  double dProbPattE3Temp[30]={0,0,p17_,p18_,0,p39_,p11_,p12_,0,0,0,p42_,p2_,p1_,p20_,p19_,0,1,0,0,0,0,0,1,p4_,p3_,0,0,p22_,p21_};
  double sProbPattE3Temp[30]={0,p41_,0,0,p29_,p30_,0,p38_,p23_,p24_,0,0,0,1,0,0,0,0,0,1,p14_,p13_,p32_,p31_,p34_,p33_,p16_,p15_,0,0}; 

  double uProbPattO0Temp[30]={0,0,0,p40_,p35_,p36_,p5_,p6_,0,p37_,0,0,p8_,p7_,0,1,0,0,0,0,0,1,p26_,p25_,0,0,p10_,p9_,p28_,p27_};
  double dProbPattO0Temp[30]={p11_,p12_,0,0,0,p42_,0,0,p17_,p18_,0,p39_,0,0,0,0,0,1,p2_,p1_,p20_,p19_,0,1,p4_,p3_,0,0,p22_,p21_};
  double sProbPattO0Temp[30]={0,p38_,p23_,p24_,0,0,0,p41_,0,0,p29_,p30_,0,1,p14_,p13_,p32_,p31_,0,1,0,0,0,0,p34_,p33_,p16_,p15_,0,0};

  double uProbPattO1Temp[30]={p11_,p12_,0,p38_,0,0,0,0,0,p39_,p29_,p30_,p2_,p1_,0,1,0,0,0,0,0,1,p32_,p31_,0,0,p4_,p3_,p34_,p33_};
  double dProbPattO1Temp[30]={0,0,p23_,p24_,0,p40_,p5_,p6_,0,0,0,p41_,0,0,0,0,0,1,p8_,p7_,p14_,p13_,0,1,p10_,p9_,0,0,p16_,p15_};
  double sProbPattO1Temp[30]={0,p42_,0,0,p35_,p36_,0,p37_,p17_,p18_,0,0,0,1,p20_,p19_,p26_,p25_,0,1,0,0,0,0,p28_,p27_,p22_,p21_,0,0};

  double uProbPattO2Temp[30]={0,0,0,p40_,p35_,p36_,p5_,p6_,0,p37_,0,0,0,0,0,1,p26_,p25_,p8_,p7_,0,1,0,0,0,0,p10_,p9_,p28_,p27_};
  double dProbPattO2Temp[30]={p11_,p12_,0,0,0,p42_,0,0,p17_,p18_,0,p39_,p2_,p1_,p20_,p19_,0,1,0,0,0,0,0,1,p4_,p3_,0,0,p22_,p21_};
  double sProbPattO2Temp[30]={0,p38_,p23_,p24_,0,0,0,p41_,0,0,p29_,p30_,0,1,0,0,0,0,0,1,p14_,p13_,p32_,p31_,p34_,p33_,p16_,p15_,0,0};

  double uProbPattO3Temp[30]={p11_,p12_,0,p38_,0,0,0,0,0,p39_,p29_,p30_,0,0,0,1,p32_,p31_,p2_,p1_,0,1,0,0,0,0,p4_,p3_,p34_,p33_};
  double dProbPattO3Temp[30]={0,0,p23_,p24_,0,p40_,p5_,p6_,0,0,0,p41_,p8_,p7_,p14_,p13_,0,1,0,0,0,0,0,1,p10_,p9_,0,0,p16_,p15_};
  double sProbPattO3Temp[30]={0,p42_,0,0,p35_,p36_,0,p37_,p17_,p18_,0,0,0,1,0,0,0,0,0,1,p20_,p19_,p26_,p25_,p28_,p27_,p22_,p21_,0,0}; 


  for(int i = 0; i < 9; i++) {
    uYopenPlaqPattE02_[i]=uYopenPlaqPattTempE02[i];
    dYopenPlaqPattE02_[i]=dYopenPlaqPattTempE02[i];
    sYopenPlaqPattE02_[i]=sYopenPlaqPattTempE02[i];
   
    uYopenPlaqPattE13_[i]=uYopenPlaqPattTempE13[i];
    dYopenPlaqPattE13_[i]=dYopenPlaqPattTempE13[i];
    sYopenPlaqPattE13_[i]=sYopenPlaqPattTempE13[i];

    uYopenPlaqPattO02_[i]=uYopenPlaqPattTempO02[i];
    dYopenPlaqPattO02_[i]=dYopenPlaqPattTempO02[i];
    sYopenPlaqPattO02_[i]=sYopenPlaqPattTempO02[i];
   
    uYopenPlaqPattO13_[i]=uYopenPlaqPattTempO13[i];
    dYopenPlaqPattO13_[i]=dYopenPlaqPattTempO13[i];
    sYopenPlaqPattO13_[i]=sYopenPlaqPattTempO13[i];


    uYopenProbPattE02_[i]=uYopenProbPattTempE02[i];
    dYopenProbPattE02_[i]=dYopenProbPattTempE02[i];
    sYopenProbPattE02_[i]=sYopenProbPattTempE02[i];

    uYopenProbPattE13_[i]=uYopenProbPattTempE13[i];
    dYopenProbPattE13_[i]=dYopenProbPattTempE13[i];
    sYopenProbPattE13_[i]=sYopenProbPattTempE13[i];

    uYopenProbPattO02_[i]=uYopenProbPattTempO02[i];
    dYopenProbPattO02_[i]=dYopenProbPattTempO02[i];
    sYopenProbPattO02_[i]=sYopenProbPattTempO02[i];

    uYopenProbPattO13_[i]=uYopenProbPattTempO13[i];
    dYopenProbPattO13_[i]=dYopenProbPattTempO13[i];
    sYopenProbPattO13_[i]=sYopenProbPattTempO13[i];
  }

  for(int i = 0; i < 30; i++) {
    
    uPlaqPattE0_[i]=uPlaqPattE0Temp[i];
    dPlaqPattE0_[i]=dPlaqPattE0Temp[i];
    sPlaqPattE0_[i]=sPlaqPattE0Temp[i];
    uPlaqPattE1_[i]=uPlaqPattE1Temp[i];
    dPlaqPattE1_[i]=dPlaqPattE1Temp[i];
    sPlaqPattE1_[i]=sPlaqPattE1Temp[i];
    uPlaqPattE2_[i]=uPlaqPattE2Temp[i];
    dPlaqPattE2_[i]=dPlaqPattE2Temp[i];
    sPlaqPattE2_[i]=sPlaqPattE2Temp[i];
    uPlaqPattE3_[i]=uPlaqPattE3Temp[i];
    dPlaqPattE3_[i]=dPlaqPattE3Temp[i];
    sPlaqPattE3_[i]=sPlaqPattE3Temp[i];
    
    uPlaqPattO0_[i]=uPlaqPattO0Temp[i];
    dPlaqPattO0_[i]=dPlaqPattO0Temp[i];
    sPlaqPattO0_[i]=sPlaqPattO0Temp[i];
    uPlaqPattO1_[i]=uPlaqPattO1Temp[i];
    dPlaqPattO1_[i]=dPlaqPattO1Temp[i];
    sPlaqPattO1_[i]=sPlaqPattO1Temp[i];
    uPlaqPattO2_[i]=uPlaqPattO2Temp[i];
    dPlaqPattO2_[i]=dPlaqPattO2Temp[i];
    sPlaqPattO2_[i]=sPlaqPattO2Temp[i];
    uPlaqPattO3_[i]=uPlaqPattO3Temp[i];
    dPlaqPattO3_[i]=dPlaqPattO3Temp[i];
    sPlaqPattO3_[i]=sPlaqPattO3Temp[i];
    
    uProbPattE0_[i]=uProbPattE0Temp[i];
    dProbPattE0_[i]=dProbPattE0Temp[i];
    sProbPattE0_[i]=sProbPattE0Temp[i];
    uProbPattE1_[i]=uProbPattE1Temp[i];
    dProbPattE1_[i]=dProbPattE1Temp[i];
    sProbPattE1_[i]=sProbPattE1Temp[i];
    uProbPattE2_[i]=uProbPattE2Temp[i];
    dProbPattE2_[i]=dProbPattE2Temp[i];
    sProbPattE2_[i]=sProbPattE2Temp[i];
    uProbPattE3_[i]=uProbPattE3Temp[i];
    dProbPattE3_[i]=dProbPattE3Temp[i];
    sProbPattE3_[i]=sProbPattE3Temp[i];
    
    uProbPattO0_[i]=uProbPattO0Temp[i];
    dProbPattO0_[i]=dProbPattO0Temp[i];
    sProbPattO0_[i]=sProbPattO0Temp[i];
    uProbPattO1_[i]=uProbPattO1Temp[i];
    dProbPattO1_[i]=dProbPattO1Temp[i];
    sProbPattO1_[i]=sProbPattO1Temp[i];
    uProbPattO2_[i]=uProbPattO2Temp[i];
    dProbPattO2_[i]=dProbPattO2Temp[i];
    sProbPattO2_[i]=sProbPattO2Temp[i];
    uProbPattO3_[i]=uProbPattO3Temp[i];
    dProbPattO3_[i]=dProbPattO3Temp[i];
    sProbPattO3_[i]=sProbPattO3Temp[i];
  
  } 
    
}


void Lattice::wormUpdate() {
  while(true) {
        
    l_++;
    if(t_bob_>=0) {
      lp_++;
    }
    if(t_bob_<=0) {
      lm_++;
    }
  
    buildCorrelator(corr1,corr2);

    if (t_dir_==1) {
      
      Nplaq_=plaq_fwd_[worm_n_];
      plaq_=plaq_record_[Nplaq_];
    
      if (pos_==0) {      // 0 e
	decisionStructure(pntPlqPatE0_,pntPrbPatE0_,1,2,3);
      }
      else if (pos_==1) { // 1 e
	decisionStructure(pntPlqPatE1_,pntPrbPatE1_,0,3,2);
      }
      else if (pos_==2) { // 0 e
	decisionStructure(pntPlqPatE0_,pntPrbPatE0_,6,5,1);
      }
      else if (pos_==3) { // 0 o
	decisionStructure(pntPlqPatO0_,pntPrbPatO0_,7,4,0);
      } 
      else if (pos_==4) { // 0 o
	decisionStructure(pntPlqPatO0_,pntPrbPatO0_,5,6,7);
      }
      else if (pos_==5) { // 1 o
	decisionStructure(pntPlqPatO1_,pntPrbPatO1_,4,7,6);
      }
      else if (pos_==6) { // 1 e
	decisionStructure(pntPlqPatE1_,pntPrbPatE1_,2,1,5);
      }
      else if (pos_==7) { // 1 o
	decisionStructure(pntPlqPatO1_,pntPrbPatO1_,3,0,4);
      } 
      else { impossible("wormUpdate: pos is not 0-7"); }      
    }
    else if (t_dir_==-1) {
      
      Nplaq_=plaq_bwd_[worm_n_];
      plaq_=plaq_record_[Nplaq_];
      
      if (pos_==0) {      // 3 o
	decisionStructure(pntPlqPatO3_,pntPrbPatO3_,4,7,3);
      }
      else if (pos_==1) { // 3 e
	decisionStructure(pntPlqPatE3_,pntPrbPatE3_,5,6,2);
      }
      else if (pos_==2) { // 2 e
	decisionStructure(pntPlqPatE2_,pntPrbPatE2_,3,0,1);
      }
      else if (pos_==3) { // 3 e
	decisionStructure(pntPlqPatE3_,pntPrbPatE3_,2,1,0);
      }
      else if (pos_==4) { // 2 o
	decisionStructure(pntPlqPatO2_,pntPrbPatO2_,0,3,7);
      }
      else if (pos_==5) { // 2 e
	decisionStructure(pntPlqPatE2_,pntPrbPatE2_,1,2,6);
      }
      else if (pos_==6) { // 2 o
	decisionStructure(pntPlqPatO2_,pntPrbPatO2_,7,4,5);
      }
      else if (pos_==7) { // 3 o
	decisionStructure(pntPlqPatO3_,pntPrbPatO3_,6,5,4);
      }
      else { impossible("wormUpdate: pos is not 0-7"); }      
    }
    else { impossible("wormUpdate: t_dir is not +1 or -1"); }
    
    if(worm_n_ == worm_start_) {
      break;
    }
     
  }
  
}

void Lattice::wormUpdateYopen() {
  while(true) {

    l_++;
    if(t_bob_>=0) {
      lp_++;
    }
    if(t_bob_<=0) {
      lm_++;
    }    
    
    buildCorrelator(corr1,corr2);
    
    if (t_dir_==1) {
      
      Nplaq_=plaq_fwd_[worm_n_];
      plaq_=plaq_record_[Nplaq_];
      
      if( (this->*fwdEdge)() ) { 
	if (pos_==0) {      // 0 e
	  decisionStructureYopen(pntPlqPatE02_,pntPrbPatE02_,1,2,3);
	}
	else if (pos_==1) { // 1 e
	  decisionStructureYopen(pntPlqPatE13_,pntPrbPatE13_,0,3,2);
	}
	else if (pos_==2) { // 0 e
	  decisionStructureYopen(pntPlqPatE02_,pntPrbPatE02_,6,5,1);
	}
	else if (pos_==3) { // 0 o
	  decisionStructureYopen(pntPlqPatO02_,pntPrbPatO02_,7,4,0);
	} 
	else if (pos_==4) { // 0 o
	  decisionStructureYopen(pntPlqPatO02_,pntPrbPatO02_,5,6,7);
	}
	else if (pos_==5) { // 1 o
	  decisionStructureYopen(pntPlqPatO13_,pntPrbPatO13_,4,7,6);
	}
	else if (pos_==6) { // 1 e
	  decisionStructureYopen(pntPlqPatE13_,pntPrbPatE13_,2,1,5);
	}
	else if (pos_==7) { // 1 o
	  decisionStructureYopen(pntPlqPatO13_,pntPrbPatO13_,3,0,4);
	} 
      }
      else if (pos_==0) {      // 0 e
	decisionStructure(pntPlqPatE0_,pntPrbPatE0_,1,2,3);
      }
      else if (pos_==1) { // 1 e
	decisionStructure(pntPlqPatE1_,pntPrbPatE1_,0,3,2);
      }
      else if (pos_==2) { // 0 e
	decisionStructure(pntPlqPatE0_,pntPrbPatE0_,6,5,1);
      }
      else if (pos_==3) { // 0 o
	decisionStructure(pntPlqPatO0_,pntPrbPatO0_,7,4,0);
      } 
      else if (pos_==4) { // 0 o
	decisionStructure(pntPlqPatO0_,pntPrbPatO0_,5,6,7);
      }
      else if (pos_==5) { // 1 o
	decisionStructure(pntPlqPatO1_,pntPrbPatO1_,4,7,6);
      }
      else if (pos_==6) { // 1 e
	decisionStructure(pntPlqPatE1_,pntPrbPatE1_,2,1,5);
      }
      else if (pos_==7) { // 1 o
	decisionStructure(pntPlqPatO1_,pntPrbPatO1_,3,0,4);
      }  
    }
    else if (t_dir_==-1) {
      
      Nplaq_=plaq_bwd_[worm_n_];
      plaq_=plaq_record_[Nplaq_];
    
      if( (this->*bwdEdge)() ) { 
	if (pos_==0) {      // 3 o
	  decisionStructureYopen(pntPlqPatO13_,pntPrbPatO13_,4,7,3);
	}
	else if (pos_==1) { // 3 e
	  decisionStructureYopen(pntPlqPatE13_,pntPrbPatE13_,5,6,2);
	}
	else if (pos_==2) { // 2 e
	  decisionStructureYopen(pntPlqPatE02_,pntPrbPatE02_,3,0,1);
	}
	else if (pos_==3) { // 3 e
	  decisionStructureYopen(pntPlqPatE13_,pntPrbPatE13_,2,1,0);
	}
	else if (pos_==4) { // 2 o
	  decisionStructureYopen(pntPlqPatO02_,pntPrbPatO02_,0,3,7);
	}
	else if (pos_==5) { // 2 e
	  decisionStructureYopen(pntPlqPatE02_,pntPrbPatE02_,1,2,6);
	}
	else if (pos_==6) { // 2 o
	  decisionStructureYopen(pntPlqPatO02_,pntPrbPatO02_,7,4,5);
	}
	else if (pos_==7) { // 3 o
	  decisionStructureYopen(pntPlqPatO13_,pntPrbPatO13_,6,5,4);
	}
      }
      else if (pos_==0) {      // 3 o
	decisionStructure(pntPlqPatO3_,pntPrbPatO3_,4,7,3);
      }
      else if (pos_==1) { // 3 e
	decisionStructure(pntPlqPatE3_,pntPrbPatE3_,5,6,2);
      }
      else if (pos_==2) { // 2 e
	decisionStructure(pntPlqPatE2_,pntPrbPatE2_,3,0,1);
      }
      else if (pos_==3) { // 3 e
	decisionStructure(pntPlqPatE3_,pntPrbPatE3_,2,1,0);
      }
      else if (pos_==4) { // 2 o
	decisionStructure(pntPlqPatO2_,pntPrbPatO2_,0,3,7);
      }
      else if (pos_==5) { // 2 e
	decisionStructure(pntPlqPatE2_,pntPrbPatE2_,1,2,6);
      }
      else if (pos_==6) { // 2 o
	decisionStructure(pntPlqPatO2_,pntPrbPatO2_,7,4,5);
      }
      else if (pos_==7) { // 3 o
	decisionStructure(pntPlqPatO3_,pntPrbPatO3_,6,5,4);
      }
    }
    else { impossible("wormUpdateYopen: t_dir is not +1 or -1"); }
    
    if(worm_n_ == worm_start_) {
      break;
    }
 
  }
  
}



void Lattice::buildCorrelator(double *gram1, double *gram2) {
  if(t_bob_<Nt_ && t_bob_>-Nt_) {
    int nX = x(worm_n_);
    int nT = t(worm_n_);

    int cX = modulo(nX-xStart_,Nx_);
    int cT;  

    if(t_bob_>=0) {
      cT = modulo(nT-tStart_,Nt_)/4;
      *(gram1 + cX + (cT*Nx_)) += 1;
    }
    if(t_bob_<=0) {
      cT = modulo(tStart_-nT,Nt_)/4;
      *(gram2 + cX + (cT*Nx_)) += 1;
    }
  }
}

void Lattice::selectCorrelators() {
  if (q8_==0) {
    if(strcmp(tail_.c_str(),"tp")==0) {
      if(q3_==0) {
	corr1=&cXTt0tp1_[0];
	corr2=&cXTt0tm1_[0];	 
      }
      else if(q3_==2) {
	corr1=&doubleNull_[0];
	corr2=&cXTtp1t0_[0];
      }
      else if(q3_==-2) {
	corr1=&cXTtm1t0_[0];
	corr2=&doubleNull_[0];
      }
    }
    else if(strcmp(tail_.c_str(),"tm")==0) {
      if(q3_==0) {
	corr1=&cXTt0tm1_[0];
	corr2=&cXTt0tp1_[0];
      }
      else if(q3_==2) {
	corr1=&cXTtp1t0_[0];
	corr2=&doubleNull_[0];
      }
      else if(q3_==-2) {
	corr1=&doubleNull_[0];
	corr2=&cXTtm1t0_[0];	 
      }
    }
  }
  if(q8_<2 && q8_>-2) {
    if(strcmp(tail_.c_str(),"up")==0) {
      if(q3_==0) {
	corr1=&cXTu0up1_[0];
	corr2=&cXTu0um1_[0];	
      }
      else if(q3_==-1) {
	corr1=&doubleNull_[0];
	corr2=&cXTup1u0_[0];	 
      }
      else if(q3_==1) {
	corr1=&cXTum1u0_[0];
	corr2=&doubleNull_[0];
      }
    }
    else if(strcmp(tail_.c_str(),"um")==0) {
      if(q3_==0) {
	corr1=&cXTu0um1_[0];
	corr2=&cXTu0up1_[0];	
      }
      else if(q3_==-1) {
	corr1=&cXTup1u0_[0];
	corr2=&doubleNull_[0];
      }
      else if(q3_==1) {
	corr1=&doubleNull_[0];
	corr2=&cXTum1u0_[0];	
      }
    }
    else if(strcmp(tail_.c_str(),"vp")==0) {
      if(q3_==0) {	
	corr1=&cXTv0vp1_[0];	 
	corr2=&cXTv0vm1_[0];	
      }
      else if(q3_==1) {
	corr1=&doubleNull_[0];
	corr2=&cXTvp1v0_[0];	 
      }
      else if(q3_==-1) {
	corr1=&cXTvm1v0_[0];
	corr2=&doubleNull_[0];
      }
    }
    else if(strcmp(tail_.c_str(),"vm")==0) {
      if(q3_==0) {
	corr1=&cXTv0vm1_[0];	 
	corr2=&cXTv0vp1_[0];	
      }
      else if(q3_==1) {
	corr1=&cXTvp1v0_[0];
	corr2=&doubleNull_[0];
      }
      else if(q3_==-1) {
	corr1=&doubleNull_[0];	
	corr2=&cXTvm1v0_[0];	 
      }
    }
  } 
}


void Lattice::measure() {
  cout << "Measuring..." << endl;
  
  int wX = 0, wY = 0, wT = 0, w3x = 0, w8x = 0, w3y = 0, w8y = 0;
  double valt0tp1, valt0tm1, valu0up1, valu0um1, valv0vp1, valv0vm1;
  double valtp1t0, valtm1t0, valup1u0, valum1u0, valvp1v0, valvm1v0; 
  
  for (int i = 1; i <= N_meas_; i++) {
    loadbar(i, N_meas_); 
    
    t_bob_ = 0; x_bob_ = 0; y_bob_ = 0; l_ = 0, lp_ = 0, lm_ = 0;
    worm_start_= ran()*volume();
    xStart_=x(worm_start_); yStart_=y(worm_start_); tStart_=t(worm_start_);
    worm_n_ = worm_start_;
  
    if (t(worm_n_)%4==3) {
      positionTool(7,6,3,2);
    }
    else if (t(worm_n_)%4==2) {
      positionTool(5,4,1,0);
    }
    else if (t(worm_n_)%4==1) {   
      positionTool(2,3,6,7);
    }
    else if (t(worm_n_)%4==0) {
      positionTool(0,1,4,5);
    }
    else {
      impossible("measure: t mod 4 is not 0-3");
    }

    selectCorrelators();

    (this->*runWorm)(); //  worm update;

    wX = x_bob_/Nx_;	
    wY = y_bob_/Ny_;
    wT = t_bob_/Nt_;

    w3x += wX*w3xAdd_;
    w8x += wX*w8xAdd_;

    w3y += wY*w3xAdd_;
    w8y += wY*w8xAdd_;

    q3_ += wT*w3xAdd_; 
    q8_ += wT*w8xAdd_;
    
    Mfile_ << std::format("{:<15},{:<15},{:<15},{:<15},{:<15},{:<15},{:<15},{:<15},{:<15},{:<15}\n", w3x, w8x, w3y, w8y, q3_, q8_, l_, lp_, lm_, tail_);   
    
    if (i%(N_meas_/20) == 0) {
     
      for (unsigned int k = 0; k < Nt_/4; k++) {
      	for (unsigned int j = 0; j < Nx_; j++) {	    
           
      	  valt0tp1 = cXTt0tp1_[j + (k*Nx_)];
      	  valt0tm1 = cXTt0tm1_[j + (k*Nx_)];
      	  valtp1t0 = cXTtp1t0_[j + (k*Nx_)];
      	  valtm1t0 = cXTtm1t0_[j + (k*Nx_)];
      	  
      	  valu0up1 = cXTu0up1_[j + (k*Nx_)];
      	  valu0um1 = cXTu0um1_[j + (k*Nx_)];
      	  valup1u0 = cXTup1u0_[j + (k*Nx_)];
      	  valum1u0 = cXTum1u0_[j + (k*Nx_)];

      	  valv0vp1 = cXTv0vp1_[j + (k*Nx_)];
      	  valv0vm1 = cXTv0vm1_[j + (k*Nx_)];
      	  valvp1v0 = cXTvp1v0_[j + (k*Nx_)];
      	  valvm1v0 = cXTvm1v0_[j + (k*Nx_)];
      	  
      	  Cfile_ << std::format("{:<20.15g},{:<20.15g},{:<20.15g},{:<20.15g},{:<20.15g},{:<20.15g},{:<20.15g},{:<20.15g},{:<20.15g},{:<20.15g},{:<20.15g},{:<20.15g}\n", valt0tp1, valt0tm1, valtp1t0, valtm1t0, valu0up1, valu0um1, valup1u0, valum1u0, valv0vp1, valv0vm1, valvp1v0, valvm1v0); 

      	}
      }

      
      fill(cXTt0tp1_.begin(), cXTt0tp1_.end(), 0);
      fill(cXTt0tm1_.begin(), cXTt0tm1_.end(), 0);
      fill(cXTtp1t0_.begin(), cXTtp1t0_.end(), 0);
      fill(cXTtm1t0_.begin(), cXTtm1t0_.end(), 0);
   
      fill(cXTu0up1_.begin(), cXTu0up1_.end(), 0);
      fill(cXTu0um1_.begin(), cXTu0um1_.end(), 0);
      fill(cXTup1u0_.begin(), cXTup1u0_.end(), 0);
      fill(cXTum1u0_.begin(), cXTum1u0_.end(), 0);

      fill(cXTv0vp1_.begin(), cXTv0vp1_.end(), 0);
      fill(cXTv0vm1_.begin(), cXTv0vm1_.end(), 0);
      fill(cXTvp1v0_.begin(), cXTvp1v0_.end(), 0);
      fill(cXTvm1v0_.begin(), cXTvm1v0_.end(), 0);


    }

  }

    

Mfile_.close();
Cfile_.close();

}

void simulate(const RunConfig& run, const std::array<double, kNumWeights>& p, const std::string& dir,
              const std::string& tag) {
  rlxd_init(2, static_cast<int>(*run.seed));
  Lattice system(dir, tag, run.bc, run.measurements, run.beta, run.Lx, run.Ly, run.Lt, run.mu3, run.mu8, p, 1);
  system.set_neighbours();
  system.set_plaq_origins();
  system.set_AFM_groundstate();
  system.setPatts();
  system.measure();
}
