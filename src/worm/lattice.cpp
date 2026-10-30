#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <ranlxd.h>
#include <format>
#include <initializer_list>
#include <iterator>
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

inline int modulo(int a, int b) {
    const int result = a % b;
    return result >= 0 ? result : result + b;
}


// Random numbers in [0,1) from ranlxd, fetched in blocks (the same sequence as one at a time).
// simulate() empties the block when it seeds ranlxd, so no run uses another run's numbers.
constexpr int kRanBlock = 96;
static double ranBlock[kRanBlock];
static int ranNext = kRanBlock;
inline double ran() {
  if (ranNext == kRanBlock) {
    ranlxd(ranBlock, kRanBlock);
    ranNext = 0;
  }
  return ranBlock[ranNext++];
}

// The moves the worm head can make across a plaquette (see decisionStructure).
enum class Move { time, space, diag };

// The plaquette type 1-15 of a plaquette value (k on one sublattice, k + 100 on the other), or 0 if
// the value is neither.
static inline int plaquetteType(int plaq) {
  if (plaq >= 1 && plaq <= 15) return plaq;
  if (plaq >= 101 && plaq <= 115) return plaq - 100;
  return 0;
}

// A worm's colour, also its index in the per-colour tables (kNone: no such move in kDiagColour).
enum Colour : signed char { kU, kD, kS, kNone = -1 };

// One step of the worm for each direction (0 forward, t_dir_ = +1; 1 backward) and position pos_:
// the patterns to use (parity of the plaquette, 0 even or 1 odd, and its bulk slot; at a y-open edge
// the edge slot is slot % 2) and the positions the head lands on after a spacelike, timelike or
// diagonal move.
struct Step {
  int par, slot, s, t, d;
};
constexpr Step kSteps[2][8] = {
  {{0, 0, 1, 2, 3}, {0, 1, 0, 3, 2}, {0, 0, 6, 5, 1}, {1, 0, 7, 4, 0}, {1, 0, 5, 6, 7}, {1, 1, 4, 7, 6}, {0, 1, 2, 1, 5}, {1, 1, 3, 0, 4}},  // forward
  {{1, 3, 4, 7, 3}, {0, 3, 5, 6, 2}, {0, 2, 3, 0, 1}, {0, 3, 2, 1, 0}, {1, 2, 0, 3, 7}, {0, 2, 1, 2, 6}, {1, 2, 7, 4, 5}, {1, 3, 6, 5, 4}},  // backward
};

// The position of a site (the pos_ a worm starting there gets), by t mod 4 and the parity of x and y:
// {x even and y even, x odd and y even, x even and y odd, x odd and y odd}.
constexpr int kStartPos[4][4] = {{0, 1, 4, 5}, {2, 3, 6, 7}, {5, 4, 1, 0}, {7, 6, 3, 2}};

// The spatial step from a site in time slice t (by t mod 4), forward (t_dir +1) and backward: along
// x or y, by `even` from an even coordinate and by -even from an odd one.
struct SpaceStep {
  bool alongY;
  int even;
};
constexpr SpaceStep kSpaceStep[2][4] = {
  {{false, +1}, {true, +1}, {false, -1}, {true, -1}},  // forward
  {{true, -1}, {false, +1}, {true, +1}, {false, -1}},  // backward
};

// The correlators, in the order of the C_ file's columns; kNull collects the counts of worms that
// add to no correlator and is not written.
enum Corr { kT0Tp1, kT0Tm1, kTp1T0, kTm1T0, kU0Up1, kU0Um1, kUp1U0, kUm1U0, kV0Vp1, kV0Vm1, kVp1V0, kVm1V0,
            kNull, kNumCorr };
constexpr int kNumCorrOut = kNull;  // the correlators written to the C_ file

// The y-open edge flags of a site (Lattice::edge_).
constexpr unsigned char kFwdEdge = 1, kBwdEdge = 2;

// The tail a worm leaves when it starts, which fixes the winding it adds (w3xAdd_, w8xAdd_; see
// positionTool).
struct Tail {
  const char* name;
  int w3, w8;
};

class Lattice {
  ofstream Mfile_, Cfile_;
  string outDir_, BC_;

  int Nx_, Ny_, Nt_, N_meas_;
  int area_, volume_;  // sites per time slice, and in all
  bool yOpen_;         // y-open boundaries (else closed)

  // The worm: its head (worm_n_, pos_, t_dir_, the plaquette plaq_ = plaq_record_[Nplaq_] it is
  // crossing), where it started, its length (l_, lp_, lm_) and how far it has wound (t_bob_ etc.).
  int worm_n_ = 0, pos_ = 0, t_dir_ = 1, plaq_ = 0, Nplaq_ = 0;
  int worm_start_ = 0, xStart_ = 0, yStart_ = 0, tStart_ = 0;
  int l_ = 0, lp_ = 0, lm_ = 0, t_bob_ = 0, x_bob_ = 0, y_bob_ = 0;
  int w3xAdd_ = 0, w8xAdd_ = 0, q3_ = 0, q8_ = 0;
  Colour colour_ = kU;  // the worm's colour (setWormColour)

  vector<int> plaq_fwd_, plaq_bwd_, next_fwd_, next_bwd_;
  vector<unsigned char> plaq_record_;  // the configuration: each plaquette's value, 0-15 or 100-115
  vector<unsigned char> edge_;  // y-open only: kFwdEdge, kBwdEdge for each site (set_edges)
  std::array<double, kNumWeights> p_;  // the run's update weights p01 ... p48
  std::array<vector<double>, kNumCorr> corr_;  // the correlators (Corr), each Lx * Lt/4

  // Probability patterns of this run ([colour][parity][slot][entry], see kBulkWeight), and the
  // current worm's patterns ([parity][slot][entry]), set with its colour by setWormColour.
  double bulkProbAll_[3][2][4][30], edgeProbAll_[3][2][2][9];
  const int (*bulkPlaq_)[4][30];
  const double (*bulkProb_)[4][30];
  const int (*edgePlaq_)[2][9];
  const double (*edgeProb_)[2][9];
  double *corr1 = nullptr, *corr2 = nullptr;
  string tail_ = "empty";


public:
  Lattice(string runDir, const string& tag, string BC, int N_meas, int Nx, int Ny, int Nt, const std::array<double,48>& p);
  int volume() const { return volume_; }
  // The site one time slice later or earlier, with periodic time.
  int tp1(int n) const { const int m = n + area_; return m >= volume_ ? m - volume_ : m; }
  int tm1(int n) const { const int m = n - area_; return m < 0 ? m + volume_ : m; }
  int xyt(int x, int y, int t) { return x + y*Nx_ + t*Nx_*Ny_; }
  int x(int i) { return i%Nx_; }
  int y(int i) { return (i/Nx_)%Ny_; }
  int t(int i) { return i/(Nx_*Ny_); }
  void set_edges();
  void set_plaq_origins();
  int spaceNeighbour(int i, int j, int k, int dir);
  void setPatts();
  void set_AFM_groundstate();
  template <bool yOpen> void wormLoop();
  template <bool yOpen, int dir, int pos> void wormStep();
  void decisionStructure(const int *, const double *, int, int, int);
  void decisionStructureYopen(const int *, const double *, int, int, int);
  // The current worm's patterns for a plaquette of parity par (0 even, 1 odd) in a slot.
  const int* bulkPlaq(int par, int slot) const { return bulkPlaq_[par][slot]; }
  const double* bulkProb(int par, int slot) const { return bulkProb_[par][slot]; }
  const int* edgePlaq(int par, int slot) const { return edgePlaq_[par][slot]; }
  const double* edgeProb(int par, int slot) const { return edgeProb_[par][slot]; }
  void setWormColour(Colour c);
  void makeMove(Move, int, int, int);
  void startWorm(Colour, const Tail&);
  void measure();
  void dirMove();
  int sitePosition(int n) { return kStartPos[t(n)%4][x(n)%2 + 2*(y(n)%2)]; }
  void positionTool();
  void timelike();
  void spacelike();
  void diag();
  void bounceTool();
  void selectCorrelators();
  void buildCorrelator(double*, double*);
  [[noreturn]] void impossible(const char* what);
  // Whether a forward (1) or backward (2) step from a site at (y, t) is at a y-open edge, for even
  // and odd Ly.
  bool yopenEdge1(int y, int t) const { return t%4==3 && (y==Ny_-1 || y==0); }
  bool yopenEdge2(int y, int t) const { return t%4==0 && (y==Ny_-1 || y==0); }
  bool yopenOddEdge1(int y, int t) const { return (t%4==1 && y==Ny_-1) || (t%4==3 && y==0); }
  bool yopenOddEdge2(int y, int t) const { return (t%4==2 && y==Ny_-1) || (t%4==0 && y==0); }
};

// A state no correct simulation reaches: stop the run (its files are then not published).
void Lattice::impossible(const char* what) {
  throw AlgorithmError(std::format("algorithm reached an impossible state in {} (worm at x={}, y={}, t={}; "
                                   "plaquette {}, pos {}, t_dir {}, worm length {})",
                                   what, x(worm_n_), y(worm_n_), t(worm_n_), plaq_, pos_, t_dir_, l_));
}

// -------------------------------------------------------------------------------------------------
// Setup: the lattice geometry, the starting configuration and the run's update patterns.

Lattice::Lattice(string runDir, const string& tag, string BC, int N_meas, int Nx, int Ny, int Nt, const std::array<double,48>& pw)
    : outDir_(runDir), BC_(BC), Nx_(Nx), Ny_(Ny), Nt_(Nt), N_meas_(N_meas), area_(Nx*Ny), volume_(Nx*Ny*Nt),
      yOpen_(BC == "y-open"), plaq_fwd_(volume_), plaq_bwd_(volume_), next_fwd_(volume_), next_bwd_(volume_),
      plaq_record_(Nx*(Ny+1)*Nt/2), p_(pw)
{
  Mfile_.open(outDir_ + "/M_" + tag + ".dat", ios::out);
  Cfile_.open(outDir_ + "/C_" + tag + ".dat", ios::out);
  if (!Mfile_ || !Cfile_) throw InputError("cannot write the M_ and C_ files in '" + outDir_ + "'");
  for (auto& c : corr_) c.assign(Nx_*Nt_/4, 0.0);

  if(BC_=="y-open") {
    cout << "y-open" << endl;
    if(Ny_%2==1) {
      cout << "odd" << endl;
    }
    else if(Ny_%2==0) {
      cout << "even" << endl;
    }
    else { impossible("Lattice: Ly is neither even nor odd"); }
  }
  else if(BC_=="closed") {
    cout << "closed" << endl;
  }
  else { impossible("Lattice: boundary condition is neither closed nor y-open"); }


 }

// The y-open edge flags of each site: whether a step forward (kFwdEdge) or backward (kBwdEdge) from
// it is at a y-open edge.
void Lattice::set_edges() {
  if (!yOpen_) return;
  const bool odd = Ny_%2==1;
  edge_.assign(volume(), 0);
  for (int n = 0; n < volume(); n++) {
    const int yn = y(n), tn = t(n);
    const bool fwd = odd ? yopenOddEdge1(yn, tn) : yopenEdge1(yn, tn);
    const bool bwd = odd ? yopenOddEdge2(yn, tn) : yopenEdge2(yn, tn);
    edge_[n] = (fwd ? kFwdEdge : 0) | (bwd ? kBwdEdge : 0);
  }
}

// The spatial neighbour of site (i, j, k) that a step forward (dir 0) or backward (dir 1) leads to
// (kSpaceStep); y is periodic here too, the y-open edges are handled by edge_.
int Lattice::spaceNeighbour(int i, int j, int k, int dir) {
  const SpaceStep& s = kSpaceStep[dir][k%4];
  if (s.alongY) return xyt(i, (j + (j%2==0 ? s.even : -s.even) + Ny_)%Ny_, k);
  return xyt((i + (i%2==0 ? s.even : -s.even) + Nx_)%Nx_, j, k);
}

// The plaquettes ahead of and behind each site (plaq_fwd_, plaq_bwd_), and its spatial neighbours.
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
        }

        next_fwd_[n] = spaceNeighbour(i, j, k, 0);
        next_bwd_[n] = spaceNeighbour(i, j, k, 1);
        n++;

      }
    }
  }
}

// The starting configuration: plaquette type 13 everywhere, on the sublattice each plaquette belongs
// to (13 even, 113 odd), which is the parity the steps from its sites use (kSteps).
void Lattice::set_AFM_groundstate() {
  for (int n = 0; n < volume(); n++) {
    const int pos = sitePosition(n);
    plaq_record_[plaq_fwd_[n]] = 13 + 100*kSteps[0][pos].par;
    plaq_record_[plaq_bwd_[n]] = 13 + 100*kSteps[1][pos].par;
  }
}

// Update patterns. For each worm colour (u, d, s), sublattice parity of the plaquette (even, odd) and
// slot (the plaquette's place in the time direction: 0-3 in the bulk; 02 and 13 at a y-open edge),
// entry i of a pattern belongs to one move of decisionStructure (offsets in kPlaquetteMoves) or of
// decisionStructureYopen (kYopenEdgeOffset).
//
// kBulkPlaq / kEdgePlaq: the plaquette value the move leaves behind. Written out for the even
// sublattice ([colour][slot][entry]); a plaquette on the odd sublattice has the same value + 100
// (withOddSublattice). Entries of 0 belong to moves of weight 0, which are never made.
template <int S, int N>
struct PlaqPatterns {
  int v[3][2][S][N];  // [colour][parity][slot][entry]
};
template <int S, int N>
constexpr PlaqPatterns<S, N> withOddSublattice(const int (&even)[3][S][N]) {
  PlaqPatterns<S, N> out{};
  for (int c = 0; c < 3; c++)
    for (int s = 0; s < S; s++)
      for (int i = 0; i < N; i++) {
        out.v[c][0][s][i] = even[c][s][i];
        out.v[c][1][s][i] = even[c][s][i] == 0 ? 0 : even[c][s][i] + 100;
      }
  return out;
}

// Bulk (closed boundaries, or away from a y-open edge).
constexpr int kBulkPlaqEven[3][4][30] = {
  {{  0,   0,   0,   6,   9,  13,  10,  13,   0,   1,   0,   0,   1,  13,   0,   9,   0,   0,   0,   0,   0,  10,   6,  13,   0,   0,   1,  10,   6,   9},
   { 10,  13,   0,   4,   0,   0,   0,   0,   0,   3,   9,  13,   4,  13,   0,   9,   0,   0,   0,   0,   0,  10,   3,  13,   0,   0,   4,  10,   3,   9},
   {  0,   0,   0,   6,  12,  13,   7,  13,   0,   1,   0,   0,   0,   0,   0,   7,   6,  13,   1,  13,   0,  12,   0,   0,   0,   0,   1,   7,   6,  12},
   {  7,  13,   0,   4,   0,   0,   0,   0,   0,   3,  12,  13,   0,   0,   0,   7,   3,  13,   4,  13,   0,  12,   0,   0,   0,   0,   4,   7,   3,  12}},  // u
  {{  7,  14,   0,   0,   0,   4,   0,   0,   8,  14,   0,   2,   0,   0,   0,   0,   0,   8,   4,  14,   2,  14,   0,   7,   4,   7,   0,   0,   2,   8},
   {  0,   0,   8,  14,   0,   5,   7,  14,   0,   0,   0,   1,   0,   0,   0,   0,   0,   8,   1,  14,   5,  14,   0,   7,   1,   7,   0,   0,   5,   8},
   { 10,  14,   0,   0,   0,   4,   0,   0,  11,  14,   0,   2,   4,  14,   2,  14,   0,  10,   0,   0,   0,   0,   0,  11,   4,  10,   0,   0,   2,  11},
   {  0,   0,  11,  14,   0,   5,  10,  14,   0,   0,   0,   1,   1,  14,   5,  14,   0,  10,   0,   0,   0,   0,   0,  11,   1,  10,   0,   0,   5,  11}},  // d
  {{  0,   5,  11,  15,   0,   0,   0,   3,   0,   0,  12,  15,   0,  12,   5,  15,   3,  15,   0,  11,   0,   0,   0,   0,   3,  12,   5,  11,   0,   0},
   {  0,   6,   0,   0,  12,  15,   0,   2,  11,  15,   0,   0,   0,  12,   2,  15,   6,  15,   0,  11,   0,   0,   0,   0,   6,  12,   2,  11,   0,   0},
   {  0,   5,   8,  15,   0,   0,   0,   3,   0,   0,   9,  15,   0,   8,   0,   0,   0,   0,   0,   9,   5,  15,   3,  15,   3,   9,   5,   8,   0,   0},
   {  0,   6,   0,   0,   9,  15,   0,   2,   8,  15,   0,   0,   0,   8,   0,   0,   0,   0,   0,   9,   2,  15,   6,  15,   6,   9,   2,   8,   0,   0}},  // s
};
constexpr auto kBulkPlaq = withOddSublattice(kBulkPlaqEven);

// At a y-open edge.
constexpr int kEdgePlaqEven[3][2][9] = {
  {{  0,   6,  13,  13,   1,   0,   0,   1,   6},
   { 13,   4,   0,   0,   3,  13,   0,   4,   3}},  // u
  {{ 14,   0,   4,   0,  14,   2,   4,   0,   2},
   {  0,  14,   5,  14,   0,   1,   1,   0,   5}},  // d
  {{  5,  15,   0,   3,   0,  15,   3,   5,   0},
   {  6,   0,  15,   2,  15,   0,   6,   2,   0}},  // s
};
constexpr auto kEdgePlaq = withOddSublattice(kEdgePlaqEven);

// kBulkWeight / kEdgeWeight: the probability of the move, as an index into the run's weights:
// 0 means probability 0, k = 1-48 means p_k (p01 ... p48 of the table), kOne means probability 1.
constexpr int kOne = 49;

// Bulk.
constexpr unsigned char kBulkWeight[3][2][4][30] = {
  {  // u
    {{ 0,  0,  0, 39, 29, 30, 11, 12,  0, 38,  0,  0,  2,  1,  0, 49,  0,  0,  0,  0,  0, 49, 32, 31,  0,  0,  4,  3, 34, 33},
     { 5,  6,  0, 37,  0,  0,  0,  0,  0, 40, 35, 36,  8,  7,  0, 49,  0,  0,  0,  0,  0, 49, 26, 25,  0,  0, 10,  9, 28, 27},
     { 0,  0,  0, 39, 29, 30, 11, 12,  0, 38,  0,  0,  0,  0,  0, 49, 32, 31,  2,  1,  0, 49,  0,  0,  0,  0,  4,  3, 34, 33},
     { 5,  6,  0, 37,  0,  0,  0,  0,  0, 40, 35, 36,  0,  0,  0, 49, 26, 25,  8,  7,  0, 49,  0,  0,  0,  0, 10,  9, 28, 27}},  // even
    {{ 0,  0,  0, 40, 35, 36,  5,  6,  0, 37,  0,  0,  8,  7,  0, 49,  0,  0,  0,  0,  0, 49, 26, 25,  0,  0, 10,  9, 28, 27},
     {11, 12,  0, 38,  0,  0,  0,  0,  0, 39, 29, 30,  2,  1,  0, 49,  0,  0,  0,  0,  0, 49, 32, 31,  0,  0,  4,  3, 34, 33},
     { 0,  0,  0, 40, 35, 36,  5,  6,  0, 37,  0,  0,  0,  0,  0, 49, 26, 25,  8,  7,  0, 49,  0,  0,  0,  0, 10,  9, 28, 27},
     {11, 12,  0, 38,  0,  0,  0,  0,  0, 39, 29, 30,  0,  0,  0, 49, 32, 31,  2,  1,  0, 49,  0,  0,  0,  0,  4,  3, 34, 33}},  // odd
  },
  {  // d
    {{ 5,  6,  0,  0,  0, 41,  0,  0, 23, 24,  0, 40,  0,  0,  0,  0,  0, 49,  8,  7, 14, 13,  0, 49, 10,  9,  0,  0, 16, 15},
     { 0,  0, 17, 18,  0, 39, 11, 12,  0,  0,  0, 42,  0,  0,  0,  0,  0, 49,  2,  1, 20, 19,  0, 49,  4,  3,  0,  0, 22, 21},
     { 5,  6,  0,  0,  0, 41,  0,  0, 23, 24,  0, 40,  8,  7, 14, 13,  0, 49,  0,  0,  0,  0,  0, 49, 10,  9,  0,  0, 16, 15},
     { 0,  0, 17, 18,  0, 39, 11, 12,  0,  0,  0, 42,  2,  1, 20, 19,  0, 49,  0,  0,  0,  0,  0, 49,  4,  3,  0,  0, 22, 21}},  // even
    {{11, 12,  0,  0,  0, 42,  0,  0, 17, 18,  0, 39,  0,  0,  0,  0,  0, 49,  2,  1, 20, 19,  0, 49,  4,  3,  0,  0, 22, 21},
     { 0,  0, 23, 24,  0, 40,  5,  6,  0,  0,  0, 41,  0,  0,  0,  0,  0, 49,  8,  7, 14, 13,  0, 49, 10,  9,  0,  0, 16, 15},
     {11, 12,  0,  0,  0, 42,  0,  0, 17, 18,  0, 39,  2,  1, 20, 19,  0, 49,  0,  0,  0,  0,  0, 49,  4,  3,  0,  0, 22, 21},
     { 0,  0, 23, 24,  0, 40,  5,  6,  0,  0,  0, 41,  8,  7, 14, 13,  0, 49,  0,  0,  0,  0,  0, 49, 10,  9,  0,  0, 16, 15}},  // odd
  },
  {  // s
    {{ 0, 37, 17, 18,  0,  0,  0, 42,  0,  0, 35, 36,  0, 49, 20, 19, 26, 25,  0, 49,  0,  0,  0,  0, 28, 27, 22, 21,  0,  0},
     { 0, 41,  0,  0, 29, 30,  0, 38, 23, 24,  0,  0,  0, 49, 14, 13, 32, 31,  0, 49,  0,  0,  0,  0, 34, 33, 16, 15,  0,  0},
     { 0, 37, 17, 18,  0,  0,  0, 42,  0,  0, 35, 36,  0, 49,  0,  0,  0,  0,  0, 49, 20, 19, 26, 25, 28, 27, 22, 21,  0,  0},
     { 0, 41,  0,  0, 29, 30,  0, 38, 23, 24,  0,  0,  0, 49,  0,  0,  0,  0,  0, 49, 14, 13, 32, 31, 34, 33, 16, 15,  0,  0}},  // even
    {{ 0, 38, 23, 24,  0,  0,  0, 41,  0,  0, 29, 30,  0, 49, 14, 13, 32, 31,  0, 49,  0,  0,  0,  0, 34, 33, 16, 15,  0,  0},
     { 0, 42,  0,  0, 35, 36,  0, 37, 17, 18,  0,  0,  0, 49, 20, 19, 26, 25,  0, 49,  0,  0,  0,  0, 28, 27, 22, 21,  0,  0},
     { 0, 38, 23, 24,  0,  0,  0, 41,  0,  0, 29, 30,  0, 49,  0,  0,  0,  0,  0, 49, 14, 13, 32, 31, 34, 33, 16, 15,  0,  0},
     { 0, 42,  0,  0, 35, 36,  0, 37, 17, 18,  0,  0,  0, 49,  0,  0,  0,  0,  0, 49, 20, 19, 26, 25, 28, 27, 22, 21,  0,  0}},  // odd
  },
};

// At a y-open edge.
constexpr unsigned char kEdgeWeight[3][2][2][9] = {
  {  // u
    {{ 0, 44, 47, 44, 47,  0,  0, 44, 47},
     {43, 48,  0,  0, 43, 48,  0, 43, 48}},  // even
    {{ 0, 43, 48, 43, 48,  0,  0, 43, 48},
     {44, 47,  0,  0, 44, 47,  0, 44, 47}},  // odd
  },
  {  // d
    {{43,  0, 46,  0, 46, 43, 43,  0, 46},
     { 0, 45, 44, 44,  0, 45, 44,  0, 45}},  // even
    {{44,  0, 45,  0, 45, 44, 44,  0, 45},
     { 0, 46, 43, 43,  0, 46, 43,  0, 46}},  // odd
  },
  {  // s
    {{48, 45,  0, 45,  0, 48, 48, 45,  0},
     {46,  0, 47, 47, 46,  0, 47, 46,  0}},  // even
    {{47, 46,  0, 46,  0, 47, 47, 46,  0},
     {45,  0, 48, 48, 45,  0, 48, 45,  0}},  // odd
  },
};

// Fills the probability patterns from the run's weights (kBulkWeight, kEdgeWeight).
void Lattice::setPatts() {
  auto weight = [&](int k) { return k == 0 ? 0.0 : k == kOne ? 1.0 : p_[k - 1]; };
  for (int c = 0; c < 3; c++)
    for (int par = 0; par < 2; par++) {
      for (int slot = 0; slot < 4; slot++)
        for (int i = 0; i < 30; i++) bulkProbAll_[c][par][slot][i] = weight(kBulkWeight[c][par][slot][i]);
      for (int slot = 0; slot < 2; slot++)
        for (int i = 0; i < 9; i++) edgeProbAll_[c][par][slot][i] = weight(kEdgeWeight[c][par][slot][i]);
    }
}

// Sets the worm's colour and points its patterns at those of that colour.
void Lattice::setWormColour(Colour c) {
  colour_ = c;
  bulkPlaq_ = kBulkPlaq.v[c];
  bulkProb_ = bulkProbAll_[c];
  edgePlaq_ = kEdgePlaq.v[c];
  edgeProb_ = edgeProbAll_[c];
}

// -------------------------------------------------------------------------------------------------
// Starting a worm.

// The tails a worm can start with, and the winding each adds.
constexpr Tail kTp{"tp", 2, 0}, kTm{"tm", -2, 0}, kVp{"vp", 1, 1}, kVm{"vm", -1, -1}, kUp{"up", -1, 1},
    kUm{"um", 1, -1};

// The bit mask of a set of plaquette types.
constexpr unsigned types(std::initializer_list<int> list) {
  unsigned mask = 0;
  for (int type : list) mask |= 1u << type;
  return mask;
}

// Start plaquette types that start a worm the same way: the two equally likely starts (colour and
// tail).
struct StartGroup {
  unsigned types;
  Colour colour1;
  Tail tail1;
  Colour colour2;
  Tail tail2;
};

// The start plaquettes possible for one direction, sublattice parity (start plaquette value above
// 100) and pair of positions, in three groups.
struct StartCase {
  int tDir;
  bool odd;
  int pos1, pos2;
  StartGroup groups[3];
};

constexpr StartCase kStartCases[8] = {
  {1, false, 0, 2, {{types({1, 6, 9, 10, 13}),  kD, kTm, kS, kVm},
                    {types({2, 4, 7, 8, 14}),   kS, kUm, kU, kTp},
                    {types({3, 5, 11, 12, 15}), kU, kVp, kD, kUp}}},
  {1, false, 1, 6, {{types({3, 4, 9, 10, 13}),  kD, kTp, kS, kVp},
                    {types({1, 5, 7, 8, 14}),   kS, kUp, kU, kTm},
                    {types({2, 6, 11, 12, 15}), kU, kVm, kD, kUm}}},
  {1, true,  3, 4, {{types({1, 6, 9, 10, 13}),  kD, kTp, kS, kVp},
                    {types({2, 4, 7, 8, 14}),   kS, kUp, kU, kTm},
                    {types({3, 5, 11, 12, 15}), kU, kVm, kD, kUm}}},
  {1, true,  5, 7, {{types({3, 4, 9, 10, 13}),  kD, kTm, kS, kVm},
                    {types({1, 5, 7, 8, 14}),   kS, kUm, kU, kTp},
                    {types({2, 6, 11, 12, 15}), kU, kVp, kD, kUp}}},
  {-1, false, 2, 5, {{types({1, 6, 7, 12, 13}),  kD, kTp, kS, kVp},
                     {types({2, 4, 10, 11, 14}), kS, kUp, kU, kTm},
                     {types({3, 5, 8, 9, 15}),   kU, kVm, kD, kUm}}},
  {-1, false, 1, 3, {{types({3, 4, 7, 12, 13}),  kD, kTm, kS, kVm},
                     {types({1, 5, 10, 11, 14}), kS, kUm, kU, kTp},
                     {types({2, 6, 8, 9, 15}),   kU, kVp, kD, kUp}}},
  {-1, true,  4, 6, {{types({1, 6, 7, 12, 13}),  kD, kTm, kS, kVm},
                     {types({2, 4, 10, 11, 14}), kS, kUm, kU, kTp},
                     {types({3, 5, 8, 9, 15}),   kU, kVp, kD, kUp}}},
  {-1, true,  0, 7, {{types({3, 4, 7, 12, 13}),  kD, kTp, kS, kVp},
                     {types({1, 5, 10, 11, 14}), kS, kUp, kU, kTm},
                     {types({2, 6, 8, 9, 15}),   kU, kVm, kD, kUm}}},
};

// Sets pos_ from the worm's start site (kStartPos), then chooses the worm's colour and tail from
// the plaquette it starts on (kStartCases).
void Lattice::positionTool() {
  pos_ = sitePosition(worm_n_);

  if (t_dir_ != 1 && t_dir_ != -1) impossible("positionTool: t_dir is not +1 or -1");
  const int startPlaq = plaq_record_[t_dir_==1 ? plaq_fwd_[worm_n_] : plaq_bwd_[worm_n_]];
  if (startPlaq == 100) impossible("positionTool: start plaquette is 100");
  const bool odd = startPlaq > 100;
  const StartCase* start = nullptr;
  for (const StartCase& c : kStartCases) {
    if (c.tDir == t_dir_ && c.odd == odd && (pos_ == c.pos1 || pos_ == c.pos2)) start = &c;
  }
  if (!start) impossible("positionTool: position not possible on a plaquette of this parity");
  const int type = plaquetteType(startPlaq);
  for (const StartGroup& g : start->groups) {
    if (!((g.types >> type) & 1u)) continue;
    if (ran() < 0.5) startWorm(g.colour1, g.tail1);
    else startWorm(g.colour2, g.tail2);
    return;
  }
  impossible("positionTool: start plaquette not possible at this position");
}

// Gives the new worm its colour and tail.
void Lattice::startWorm(Colour colour, const Tail& tail) {
  setWormColour(colour);
  tail_ = tail.name;
  w3xAdd_ = tail.w3;
  w8xAdd_ = tail.w8;
}

// -------------------------------------------------------------------------------------------------
// The worm update: one step after another until the head is back at the tail.

// One step of the worm in direction dir from position pos (kSteps). With yOpen, a step at a y-open
// edge (edge_) uses the edge patterns.
template <bool yOpen, int dir, int pos>
void Lattice::wormStep() {
  constexpr Step step = kSteps[dir][pos];
  bool edge = false;
  if constexpr (yOpen) edge = edge_[worm_n_] & (dir == 0 ? kFwdEdge : kBwdEdge);
  if (edge) {
    decisionStructureYopen(edgePlaq(step.par,step.slot%2),edgeProb(step.par,step.slot%2),step.s,step.t,step.d);
  }
  else {
    decisionStructure(bulkPlaq(step.par,step.slot),bulkProb(step.par,step.slot),step.s,step.t,step.d);
  }
}

// Moves the worm until its head is back at the tail. With yOpen, steps at a y-open edge (edge_) use
// the edge patterns.
template <bool yOpen>
void Lattice::wormLoop() {
  while(true) {

    l_++;
    if(t_bob_>=0) {
      lp_++;
    }
    if(t_bob_<=0) {
      lm_++;
    }

    buildCorrelator(corr1,corr2);

    if (t_dir_ != 1 && t_dir_ != -1) {
      impossible(yOpen ? "wormLoop (y-open): t_dir is not +1 or -1" : "wormLoop (closed): t_dir is not +1 or -1");
    }
    const bool forward = t_dir_==1;
    Nplaq_ = forward ? plaq_fwd_[worm_n_] : plaq_bwd_[worm_n_];
    plaq_=plaq_record_[Nplaq_];
    // A switch, so that each step's patterns and positions are compile-time constants (a run-time
    // lookup in kSteps here made the simulator 2-5% slower).
    if (forward) {
      switch (pos_) {
        case 0: wormStep<yOpen, 0, 0>(); break;
        case 1: wormStep<yOpen, 0, 1>(); break;
        case 2: wormStep<yOpen, 0, 2>(); break;
        case 3: wormStep<yOpen, 0, 3>(); break;
        case 4: wormStep<yOpen, 0, 4>(); break;
        case 5: wormStep<yOpen, 0, 5>(); break;
        case 6: wormStep<yOpen, 0, 6>(); break;
        case 7: wormStep<yOpen, 0, 7>(); break;
        default: impossible(yOpen ? "wormLoop (y-open): pos is not 0-7" : "wormLoop (closed): pos is not 0-7");
      }
    }
    else {
      switch (pos_) {
        case 0: wormStep<yOpen, 1, 0>(); break;
        case 1: wormStep<yOpen, 1, 1>(); break;
        case 2: wormStep<yOpen, 1, 2>(); break;
        case 3: wormStep<yOpen, 1, 3>(); break;
        case 4: wormStep<yOpen, 1, 4>(); break;
        case 5: wormStep<yOpen, 1, 5>(); break;
        case 6: wormStep<yOpen, 1, 6>(); break;
        case 7: wormStep<yOpen, 1, 7>(); break;
        default: impossible(yOpen ? "wormLoop (y-open): pos is not 0-7" : "wormLoop (closed): pos is not 0-7");
      }
    }

    if(worm_n_ == worm_start_) {
      break;
    }
  }
}


// For each plaquette type (closed boundaries): the two moves the worm can make across it, with the
// offsets in the pattern arrays (pr: probability, pl: the plaquette value it leaves behind) of the
// first and second move. Any other outcome is a bounce. Entry 0 is unused.
struct PlaquetteMoves {
  Move first, second;
  int a, b;
};
constexpr PlaquetteMoves kPlaquetteMoves[16] = {
  {},
  {Move::time,  Move::diag,   1,  0}, {Move::time,  Move::diag,   3,  2}, {Move::time,  Move::diag,   5,  4},
  {Move::time,  Move::diag,   7,  6}, {Move::time,  Move::diag,   9,  8}, {Move::time,  Move::diag,  11, 10},
  {Move::space, Move::diag,  13, 12}, {Move::space, Move::diag,  15, 14}, {Move::space, Move::diag,  17, 16},
  {Move::space, Move::diag,  19, 18}, {Move::space, Move::diag,  21, 20}, {Move::space, Move::diag,  23, 22},
  {Move::time,  Move::space, 24, 25}, {Move::time,  Move::space, 26, 27}, {Move::time,  Move::space, 28, 29},
};

// At a y-open edge the only move is timelike: the offset of its entries for each plaquette type, or
// -1 for the types that cannot occur there.
constexpr int kYopenEdgeOffset[16] = {-1, 0, 1, 2, 3, 4, 5, -1, -1, -1, -1, -1, -1, 6, 7, 8};

// Moves the worm head across the current plaquette and sets pos_ to the position it lands on
// (s, t, d for a spacelike, timelike or diagonal move).
inline void Lattice::makeMove(Move move, int s, int t, int d) {
  switch (move) {
    case Move::time:  timelike();  pos_=t; break;
    case Move::space: spacelike(); pos_=s; break;
    case Move::diag:  diag();      pos_=d; break;
  }
}

// One step of the worm across the plaquette plaq_: the first move, the second move or a bounce,
// with the probabilities pr and new plaquette values pl of the pattern for this plaquette type.
void Lattice::decisionStructure(const int *pl, const double *pr, int s, int t, int d) {
  double myRan = ran();
  const int type = plaquetteType(plaq_);
  if (type == 0) impossible("decisionStructure: plaquette is not 1-15 or 101-115");
  const PlaquetteMoves& m = kPlaquetteMoves[type];
  if (myRan < pr[m.a]) {
    makeMove(m.first, s, t, d);
    plaq_record_[Nplaq_]=pl[m.a];
  }
  else if (myRan > pr[m.a] && myRan <= (pr[m.a]+pr[m.b])) {
    makeMove(m.second, s, t, d);
    plaq_record_[Nplaq_]=pl[m.b];
  }
  else {
    bounceTool();
    t_dir_=-t_dir_;
  }
}

// The same at a y-open edge, where the only move is timelike.
void Lattice::decisionStructureYopen(const int *pl, const double *pr, int s, int t, int d) {
  double myRan = ran();
  const int a = kYopenEdgeOffset[plaquetteType(plaq_)];
  if (a < 0) impossible("decisionStructureYopen: plaquette is not 1-6, 13-15, 101-106 or 113-115");
  if (myRan < pr[a]) {
    timelike();
    pos_=t;
    plaq_record_[Nplaq_]=pl[a];
  }
  else {
    bounceTool();
    t_dir_=-t_dir_;
  }
}

// -------------------------------------------------------------------------------------------------
// The moves of the worm head, and the colour it takes.

void Lattice::timelike() {
  if (t_dir_==1) {
    worm_n_ = tp1(worm_n_);
     t_bob_++;
  }
  else if (t_dir_==-1) {
    worm_n_ = tm1(worm_n_);
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

// The colour a worm of colour u, d or s takes after a diagonal move forward (t_dir +1) or backward
// across a plaquette of type 1-15 (entry 0 unused).
constexpr Colour kDiagColour[2][3][16] = {
  {  // forward    1      2      3      4      5      6      7      8      9     10     11     12  13-15
    {kNone,    kD, kNone,    kS,    kD, kNone,    kS,    kD, kNone, kNone, kNone, kNone,    kS, kNone, kNone, kNone},  // u
    {kNone,    kU,    kS, kNone,    kU,    kS, kNone, kNone, kNone, kNone,    kU,    kS, kNone, kNone, kNone, kNone},  // d
    {kNone, kNone,    kD,    kU, kNone,    kD,    kU, kNone,    kD,    kU, kNone, kNone, kNone, kNone, kNone, kNone},  // s
  },
  {  // backward
    {kNone,    kD, kNone,    kS,    kD, kNone,    kS, kNone, kNone,    kS,    kD, kNone, kNone, kNone, kNone, kNone},  // u
    {kNone,    kU,    kS, kNone,    kU,    kS, kNone,    kU,    kS, kNone, kNone, kNone, kNone, kNone, kNone, kNone},  // d
    {kNone, kNone,    kD,    kU, kNone,    kD,    kU, kNone, kNone, kNone, kNone,    kD,    kU, kNone, kNone, kNone},  // s
  },
};

// Moves the worm head diagonally across the plaquette plaq_, changing its colour (kDiagColour).
void Lattice::diag() {
  dirMove();
  int dir;
  if (t_dir_==1) {
    worm_n_ = tp1(next_fwd_[worm_n_]);  // the spatial neighbour, one time slice later
    t_bob_++;
    dir = 0;
  }
  else if (t_dir_==-1) {
    worm_n_ = tm1(next_bwd_[worm_n_]);  // the spatial neighbour, one time slice earlier
    t_bob_--;
    dir = 1;
  }
  else {
    impossible("diag: t_dir is not +1 or -1");
  }
  const Colour c = kDiagColour[dir][colour_][plaquetteType(plaq_)];
  if (c == kNone) impossible("diag: no diagonal move across this plaquette for this colour");
  setWormColour(c);
}

// Records the spatial step the worm head is about to take (x_bob_, y_bob_; kSpaceStep).
// Records the spatial step the worm head is about to take (x_bob_, y_bob_). These are the steps of
// kSpaceStep, written out as branches: looking them up in kSpaceStep here made the simulator 1-2%
// slower (2026-10-30).
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

// The colour a worm takes when it bounces back from a plaquette of type 1-15, forward (t_dir +1) or
// backward, in one of two groups of positions (kBounceGroup; only types 1-6 depend on the group).
constexpr Colour kBounceColour[2][2][16] = {
  {  // forward    1   2   3   4   5   6   7   8   9  10  11  12  13  14  15
    {kNone,       kU, kD, kS, kD, kS, kU, kD, kD, kU, kU, kS, kS, kU, kD, kS},
    {kNone,       kD, kS, kU, kU, kD, kS, kD, kD, kU, kU, kS, kS, kU, kD, kS},
  },
  {  // backward
    {kNone,       kU, kD, kS, kD, kS, kU, kU, kS, kS, kD, kD, kU, kU, kD, kS},
    {kNone,       kD, kS, kU, kU, kD, kS, kU, kS, kS, kD, kD, kU, kU, kD, kS},
  },
};

// The group in kBounceColour of each position 0-7, forward and backward.
constexpr int kBounceGroup[2][8] = {
  {0, 1, 0, 0, 0, 1, 1, 1},
  {1, 1, 0, 1, 0, 0, 0, 1},
};

// The worm bounces back from the plaquette plaq_ and takes the colour kBounceColour gives it.
// plaq_ was checked by decisionStructure(Yopen), t_dir_ (+1 or -1) and pos_ (0-7) by wormLoop.
void Lattice::bounceTool() {
  const int dir = t_dir_==1 ? 0 : 1;
  setWormColour(kBounceColour[dir][kBounceGroup[dir][pos_]][plaquetteType(plaq_)]);
}

// -------------------------------------------------------------------------------------------------
// Correlators and measurements.

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

// The correlators a worm adds to (corr1, corr2), by its tail and the winding numbers q3, q8 so far.
// A worm whose tail and windings match no row leaves corr1 and corr2 as they were.
struct CorrChoice {
  const char* tail;
  int maxQ8;  // the row applies only if |q8| <= maxQ8
  int q3;
  Corr corr1, corr2;
};
constexpr CorrChoice kCorrChoices[] = {
  {"tp", 0,  0, kT0Tp1, kT0Tm1}, {"tp", 0,  2, kNull,  kTp1T0}, {"tp", 0, -2, kTm1T0, kNull},
  {"tm", 0,  0, kT0Tm1, kT0Tp1}, {"tm", 0,  2, kTp1T0, kNull},  {"tm", 0, -2, kNull,  kTm1T0},
  {"up", 1,  0, kU0Up1, kU0Um1}, {"up", 1, -1, kNull,  kUp1U0}, {"up", 1,  1, kUm1U0, kNull},
  {"um", 1,  0, kU0Um1, kU0Up1}, {"um", 1, -1, kUp1U0, kNull},  {"um", 1,  1, kNull,  kUm1U0},
  {"vp", 1,  0, kV0Vp1, kV0Vm1}, {"vp", 1,  1, kNull,  kVp1V0}, {"vp", 1, -1, kVm1V0, kNull},
  {"vm", 1,  0, kV0Vm1, kV0Vp1}, {"vm", 1,  1, kVp1V0, kNull},  {"vm", 1, -1, kNull,  kVm1V0},
};

static const CorrChoice* findCorrChoice(const string& tail, int q3, int q8) {
  for (const CorrChoice& c : kCorrChoices) {
    if (tail == c.tail && q3 == c.q3 && q8 >= -c.maxQ8 && q8 <= c.maxQ8) return &c;
  }
  return nullptr;
}

void Lattice::selectCorrelators() {
  if (const CorrChoice* c = findCorrChoice(tail_, q3_, q8_)) {
    corr1 = corr_[c->corr1].data();
    corr2 = corr_[c->corr2].data();
  }
}


void Lattice::measure() {
  cout << "Measuring..." << endl;

  int wX = 0, wY = 0, wT = 0, w3x = 0, w8x = 0, w3y = 0, w8y = 0;
  string line;

  for (int i = 1; i <= N_meas_; i++) {
    loadbar(i, N_meas_);

    t_bob_ = 0; x_bob_ = 0; y_bob_ = 0; l_ = 0, lp_ = 0, lm_ = 0;
    worm_start_= ran()*volume();
    xStart_=x(worm_start_); yStart_=y(worm_start_); tStart_=t(worm_start_);
    worm_n_ = worm_start_;

    positionTool();

    selectCorrelators();

    if (yOpen_) wormLoop<true>(); else wormLoop<false>();

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
      // One line per x and t/4: the correlators, comma separated.
      for (int k = 0; k < Nt_/4; k++) {
        for (int j = 0; j < Nx_; j++) {
          const int m = j + k*Nx_;
          Cfile_ << std::format("{:<20.15g},{:<20.15g},{:<20.15g},{:<20.15g},{:<20.15g},{:<20.15g},{:<20.15g},{:<20.15g},{:<20.15g},{:<20.15g},{:<20.15g},{:<20.15g}\n",
                                corr_[0][m], corr_[1][m], corr_[2][m], corr_[3][m], corr_[4][m], corr_[5][m], corr_[6][m], corr_[7][m], corr_[8][m], corr_[9][m], corr_[10][m], corr_[11][m]);
        }
      }
      for (int c = 0; c < kNumCorrOut; c++) fill(corr_[c].begin(), corr_[c].end(), 0);


    }

  }



Mfile_.close();
Cfile_.close();

}

void simulate(const RunConfig& run, const std::array<double, kNumWeights>& p, const std::string& dir,
              const std::string& tag) {
  rlxd_init(2, static_cast<int>(*run.seed));
  ranNext = kRanBlock;
  Lattice system(dir, tag, run.bc, run.measurements, run.Lx, run.Ly, run.Lt, p);
  system.set_edges();
  system.set_plaq_origins();
  system.set_AFM_groundstate();
  system.setPatts();
  system.measure();
}
