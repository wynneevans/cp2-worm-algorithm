#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <cstring>
#include <cstdio>
#include <ctime>
#include <sstream>
#include <boost/format.hpp>
#include <cp2worm/jackknife.h>   

using namespace std;
using namespace boost;

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, elems);
  return elems;
}

void printJackknifedCorrelator(int C, double* slice, int nL, double* cLC, double* cL, double* cLerr) {
  
  for (int v = 0; v < nL; v++) {
    JK jacker(C, 1, 0);
    for (int i = 0; i < C; i++) {
      *(slice+i) = *(cLC+((i*nL)+v));
    }
    double norm;
    jacker.jackknife(slice);
    if(v==0) {
    norm = jacker.jackMean();
    }
    *(cL+v) = jacker.jackMean()/norm;
    *(cLerr+v) = jacker.jackSE()/norm;
    
  }
}

int main(int argc, char** argv) {
  string dir = argv[1];
  string file = argv[2];
  vector<string> splitFile = split(argv[2], '_');
  string fileType = splitFile[0];
  cout << "fileType: " << fileType << endl;
  string model = splitFile[1];
  cout << "model: " << model << endl;
  string BC = splitFile[3];
  cout << "BC: " << BC << endl;
  int nMeas = stoi(splitFile[5]);
  cout << "nMeas: " << nMeas << endl;
  double beta = stod(splitFile[7]);
  cout << "beta: " << beta << endl;
  int nX = stoi(splitFile[9]);
  cout << "nX: " << nX << endl;
  int nY = stoi(splitFile[10]);
  cout << "nY: " << nY << endl;
  int nT = stoi(splitFile[11])/4;
  cout << "nT: " << nT << endl;
  double mu3 = stof(splitFile[13]);
  cout << "mu3: " << mu3 << endl;
  double mu8 = stof(splitFile[15]);
  cout << "mu8: " << mu8 << endl;
  
  string line;
  double pI=3.14159, dubb1, dubb2, dubb3, dubb4, dubb5, dubb6, dubb7, dubb8, dubb9, dubb10, dubb11, dubb12;
  int nCtotal = atoi(argv[3]);
  int nCskip = atoi(argv[4]);
  int nC = nCtotal - nCskip;
  int vol = nX*nT*nC;
  int alwaysSkipOne = nCskip*nX*nT;
  
  vector<double> corrMeasSlice(nC);
  vector<double> cXTt0tp1(vol), cXTt0tm1(vol), cXTtp1t0(vol), cXTtm1t0(vol), cXTu0up1(vol), cXTu0um1(vol), cXTup1u0(vol), cXTum1u0(vol), cXTv0vp1(vol), cXTv0vm1(vol), cXTvp1v0(vol), cXTvm1v0(vol);  
  vector<double> cXt0tp1(nX*nC), cXt0tm1(nX*nC), cXtp1t0(nX*nC), cXtm1t0(nX*nC), cXu0up1(nX*nC), cXu0um1(nX*nC), cXup1u0(nX*nC), cXum1u0(nX*nC), cXv0vp1(nX*nC), cXv0vm1(nX*nC), cXvp1v0(nX*nC), cXvm1v0(nX*nC); 
  vector<double> cTt0tp1(nT*nC), cTt0tm1(nT*nC), cTtp1t0(nT*nC), cTtm1t0(nT*nC), cTu0up1(nT*nC), cTu0um1(nT*nC), cTup1u0(nT*nC), cTum1u0(nT*nC), cTv0vp1(nT*nC), cTv0vm1(nT*nC), cTvp1v0(nT*nC), cTvm1v0(nT*nC); 

  vector<double> cXt0tp1_(nX), cXt0tm1_(nX), cXtp1t0_(nX), cXtm1t0_(nX), cXu0up1_(nX), cXu0um1_(nX), cXup1u0_(nX), cXum1u0_(nX), cXv0vp1_(nX), cXv0vm1_(nX), cXvp1v0_(nX), cXvm1v0_(nX);
  vector<double> cTt0tp1_(nT), cTt0tm1_(nT), cTtp1t0_(nT), cTtm1t0_(nT), cTu0up1_(nT), cTu0um1_(nT), cTup1u0_(nT), cTum1u0_(nT), cTv0vp1_(nT), cTv0vm1_(nT), cTvp1v0_(nT), cTvm1v0_(nT);

  vector<double> cXt0tp1err_(nX), cXt0tm1err_(nX), cXtp1t0err_(nX), cXtm1t0err_(nX), cXu0up1err_(nX), cXu0um1err_(nX), cXup1u0err_(nX), cXum1u0err_(nX), cXv0vp1err_(nX), cXv0vm1err_(nX), cXvp1v0err_(nX), cXvm1v0err_(nX); 
  vector<double> cTt0tp1err_(nT), cTt0tm1err_(nT), cTtp1t0err_(nT), cTtm1t0err_(nT), cTu0up1err_(nT), cTu0um1err_(nT), cTup1u0err_(nT), cTum1u0err_(nT), cTv0vp1err_(nT), cTv0vm1err_(nT), cTvp1v0err_(nT), cTvm1v0err_(nT);   

  stringstream oOpen;
  oOpen << dir << "/" << file;// << ".dat";
  
  ifstream data(oOpen.str());
  cout << oOpen.str() << endl;

  for(int m = 0; m < alwaysSkipOne; m++) {
    getline(data, line);
  }
  
  for(int m = 0; m < vol; m++) {
    getline(data, line);
    stringstream ss( line );
    vector<string> record;
    
    while ( ss.good() )
      {
        string substr;
        getline( ss, substr, ',' );
        record.push_back( substr );
      }
      
    // data >> dubb1 >> dubb2 >> dubb3 >> dubb4 >> dubb5 >> dubb6 >> dubb7 >> dubb8 >> dubb9 >> dubb10 >> dubb11 >> dubb12;
    
    dubb1 = std::stod(record[0]);
    dubb2 = std::stod(record[1]);
    dubb3 = std::stod(record[2]);
    dubb4 = std::stod(record[3]);
    dubb5 = std::stod(record[4]);
    dubb6 = std::stod(record[5]);
    dubb7 = std::stod(record[6]);
    dubb8 = std::stod(record[7]);
    dubb9 = std::stod(record[8]);
    dubb10 = std::stod(record[9]);
    dubb11 = std::stod(record[10]);
    dubb12 = std::stod(record[11]);

    cout << ' ' << dubb1 << ' ' << dubb2 << ' ' << dubb3 << ' ' << dubb4 << ' ' << dubb5 << ' ' << dubb6 << ' ' << dubb7 << ' ' << dubb8 << ' ' << dubb9 << ' ' << endl;

    cXTt0tp1[m] = dubb1;
    cXTt0tm1[m] = dubb2;
    cXTtp1t0[m] = dubb3;
    cXTtm1t0[m] = dubb4;
    cXTu0up1[m] = dubb5; 
    cXTu0um1[m] = dubb6;
    cXTup1u0[m] = dubb7;
    cXTum1u0[m] = dubb8; 
    cXTv0vp1[m] = dubb9;
    cXTv0vm1[m] = dubb10; 
    cXTvp1v0[m] = dubb11; 
    cXTvm1v0[m] = dubb12;

  } 
  
  data.close();


  for (int c=0; c<nC; c++) {
    for (int x=0; x<nX; x++) {
      for (int t=0; t<nT; t++) {

	cXt0tp1[x+c*nX] += cXTt0tp1[x+t*nX+c*nX*nT];
	cXt0tm1[x+c*nX] += cXTt0tm1[x+t*nX+c*nX*nT];
	cXtp1t0[x+c*nX] += cXTtp1t0[x+t*nX+c*nX*nT];
	cXtm1t0[x+c*nX] += cXTtm1t0[x+t*nX+c*nX*nT];
	cXu0up1[x+c*nX] += cXTu0up1[x+t*nX+c*nX*nT]; 
	cXu0um1[x+c*nX] += cXTu0um1[x+t*nX+c*nX*nT];
	cXup1u0[x+c*nX] += cXTup1u0[x+t*nX+c*nX*nT];
	cXum1u0[x+c*nX] += cXTum1u0[x+t*nX+c*nX*nT]; 
	cXv0vp1[x+c*nX] += cXTv0vp1[x+t*nX+c*nX*nT];
	cXv0vm1[x+c*nX] += cXTv0vm1[x+t*nX+c*nX*nT]; 
	cXvp1v0[x+c*nX] += cXTvp1v0[x+t*nX+c*nX*nT]; 
	cXvm1v0[x+c*nX] += cXTvm1v0[x+t*nX+c*nX*nT];
        
      }
    }
  }

  printJackknifedCorrelator(nC, &corrMeasSlice[0], nX, &cXt0tp1[0], &cXt0tp1_[0], &cXt0tp1err_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nX, &cXt0tm1[0], &cXt0tm1_[0], &cXt0tm1err_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nX, &cXtp1t0[0], &cXtp1t0_[0], &cXtp1t0err_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nX, &cXtm1t0[0], &cXtm1t0_[0], &cXtm1t0err_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nX, &cXu0up1[0], &cXu0up1_[0], &cXu0up1err_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nX, &cXu0um1[0], &cXu0um1_[0], &cXu0um1err_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nX, &cXup1u0[0], &cXup1u0_[0], &cXup1u0err_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nX, &cXum1u0[0], &cXum1u0_[0], &cXum1u0err_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nX, &cXv0vp1[0], &cXv0vp1_[0], &cXv0vp1err_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nX, &cXv0vm1[0], &cXv0vm1_[0], &cXv0vm1err_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nX, &cXvp1v0[0], &cXvp1v0_[0], &cXvp1v0err_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nX, &cXvm1v0[0], &cXvm1v0_[0], &cXvm1v0err_[0]);

  for (int c=0; c<nC; c++) {
    for (int t=0; t<nT; t++) {
      for (int x=0; x<nX; x++) {

	cTt0tp1[t+c*nT] += cXTt0tp1[x+t*nX+c*nX*nT];
	cTt0tm1[t+c*nT] += cXTt0tm1[x+t*nX+c*nX*nT];
	cTtp1t0[t+c*nT] += cXTtp1t0[x+t*nX+c*nX*nT];
	cTtm1t0[t+c*nT] += cXTtm1t0[x+t*nX+c*nX*nT];
	cTu0up1[t+c*nT] += cXTu0up1[x+t*nX+c*nX*nT]; 
	cTu0um1[t+c*nT] += cXTu0um1[x+t*nX+c*nX*nT];
	cTup1u0[t+c*nT] += cXTup1u0[x+t*nX+c*nX*nT];
	cTum1u0[t+c*nT] += cXTum1u0[x+t*nX+c*nX*nT]; 
	cTv0vp1[t+c*nT] += cXTv0vp1[x+t*nX+c*nX*nT];
	cTv0vm1[t+c*nT] += cXTv0vm1[x+t*nX+c*nX*nT]; 
	cTvp1v0[t+c*nT] += cXTvp1v0[x+t*nX+c*nX*nT]; 
	cTvm1v0[t+c*nT] += cXTvm1v0[x+t*nX+c*nX*nT];

      }
    }
  }

  printJackknifedCorrelator(nC, &corrMeasSlice[0], nT, &cTt0tp1[0], &cTt0tp1_[0], &cTt0tp1err_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nT, &cTt0tm1[0], &cTt0tm1_[0], &cTt0tm1err_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nT, &cTtp1t0[0], &cTtp1t0_[0], &cTtp1t0err_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nT, &cTtm1t0[0], &cTtm1t0_[0], &cTtm1t0err_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nT, &cTu0up1[0], &cTu0up1_[0], &cTu0up1err_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nT, &cTu0um1[0], &cTu0um1_[0], &cTu0um1err_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nT, &cTup1u0[0], &cTup1u0_[0], &cTup1u0err_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nT, &cTum1u0[0], &cTum1u0_[0], &cTum1u0err_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nT, &cTv0vp1[0], &cTv0vp1_[0], &cTv0vp1err_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nT, &cTv0vm1[0], &cTv0vm1_[0], &cTv0vm1err_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nT, &cTvp1v0[0], &cTvp1v0_[0], &cTvp1v0err_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nT, &cTvm1v0[0], &cTvm1v0_[0], &cTvm1v0err_[0]);
  
  stringstream wWrite;
  wWrite << format("%s/xiWorm_nCskip_%f_%s") % dir % nCskip % file;

  cout << wWrite.str() << endl;

  ofstream outData(wWrite.str());

  outData << "cX" << endl;
  for(int x=0; x<nX; x++) {
    outData << format("%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d\n") % cXt0tp1_[x] % cXt0tp1err_[x] % cXt0tm1_[x] % cXt0tm1err_[x] % cXtp1t0_[x] % cXtp1t0err_[x] % cXtm1t0_[x] % cXtm1t0err_[x] % cXu0up1_[x] % cXu0up1err_[x] % cXu0um1_[x] % cXu0um1err_[x] % cXup1u0_[x] % cXup1u0err_[x] % cXum1u0_[x] % cXum1u0err_[x] % cXv0vp1_[x] % cXv0vp1err_[x] % cXv0vm1_[x] % cXv0vm1err_[x] % cXvp1v0_[x] % cXvp1v0err_[x] % cXvm1v0_[x] % cXvm1v0err_[x];
  }
 
  outData << "cT" << endl;
  for(int t=0; t<nT; t++) {
    outData << format("%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d%-20.12d\n") % cTt0tp1_[t] % cTt0tp1err_[t] % cTt0tm1_[t] % cTt0tm1err_[t] % cTtp1t0_[t] % cTtp1t0err_[t] % cTtm1t0_[t] % cTtm1t0err_[t] % cTu0up1_[t] % cTu0up1err_[t] % cTu0um1_[t] % cTu0um1err_[t] % cTup1u0_[t] % cTup1u0err_[t] % cTum1u0_[t] % cTum1u0err_[t] %  cTv0vp1_[t] % cTv0vp1err_[t] % cTv0vm1_[t] % cTv0vm1err_[t] %  cTvp1v0_[t] % cTvp1v0err_[t] % cTvm1v0_[t] % cTvm1v0err_[t];
  }
 
  outData.close();

}

 
