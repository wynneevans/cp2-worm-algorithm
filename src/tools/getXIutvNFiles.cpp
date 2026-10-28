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
    
    //cout  << *(cLC+v) << " " << *(cLC+v+nL) << " " << *(cLC+v+(2*nL)) << " "  << *(cLC+v+(3*nL)) << " " << *(cLC+v+(4*nL)) << " " << " " << *(cL+v) << " " << *(cLerr+v) << endl; //
  
  }
  
}


int main(int argc, char** argv) {
  //Wsqrd_CP2AFM_BC_closed_nTherm_10000_nMeas_100000_beta_30.60_XYT_50_50_12240_nDef_0_nSim_1.dat
  string dir = argv[1];
  string file = argv[2];
  vector<string> splitFile = split(argv[2], '_');
  string fileType = splitFile[0];
  cout << "fileType: " << fileType << endl;
  string model = splitFile[1];
  cout << "model: " << model << endl;
  string BC = splitFile[3];
  cout << "BC: " << BC << endl;
  int nTherm = stoi(splitFile[5]);
  cout << "nTherm: " << nTherm << endl;
  int nMeas = stoi(splitFile[7]);
  cout << "nMeas: " << nMeas << endl;
  double beta = stod(splitFile[9]);
  cout << "beta: " << beta << endl;
  int nX = stoi(splitFile[11]);
  cout << "nX: " << nX << endl;
  int nY = stoi(splitFile[12]);
  cout << "nY: " << nY << endl;
  int nT = stoi(splitFile[13])/4;
  cout << "nT: " << nT << endl;
  if(strcmp(model.c_str(),"wormCP2AFMmu2dplaq")==0) {
    double mu3 = stof(splitFile[15]);
    cout << "mu3: " << mu3 << endl;
    double mu8 = stof(splitFile[17]);
    cout << "mu8: " << mu8 << endl;
    //    int nSim = stoi(splitFile[19]);
    //    cout << "nSim :" << nSim << endl;
  }
  else {
    int nDef = stoi(splitFile[15]);
    cout << "nDef :" << nDef << endl;
    //    int nSim = stoi(splitFile[17]);
    //    cout << "nSim :" << nSim << endl;
  }
  
  string line;
  double pI=3.14159, dubb1, dubb2, dubb3, dubb4, dubb5, dubb6;
  int corrPF = atoi(argv[3]);
  int nSims = atoi(argv[4]);
  int nC = corrPF*nSims;
  int volPF = nX*nT*corrPF;
  int vol = nX*nT*nC;
  vector<double> corrMeasSlice(nC);
  vector<double> cXTupun(vol), cXTtptn(vol), cXTvpvn(vol);
  vector<double> cXupun(nX*nC), cXtptn(nX*nC), cXvpvn(nX*nC);
  vector<double> cTupun(nT*nC), cTtptn(nT*nC), cTvpvn(nT*nC);

  vector<double> cXTup(vol), cXTun(vol), cXTtp(vol), cXTtn(vol), cXTvp(vol), cXTvn(vol);
  vector<double> cXup(nX*nC), cXun(nX*nC), cXtp(nX*nC), cXtn(nX*nC), cXvp(nX*nC), cXvn(nX*nC);
  vector<double> cTup(nT*nC), cTun(nT*nC), cTtp(nT*nC), cTtn(nT*nC), cTvp(nT*nC), cTvn(nT*nC);

  vector<double> c0xupun(nC), c0xtptn(nC), c0xvpvn(nC);
  vector<double> cLxupun(nC), cLxtptn(nC), cLxvpvn(nC);
  vector<double> c0tupun(nC), c0ttptn(nC), c0tvpvn(nC);
  vector<double> cLtupun(nC), cLttptn(nC), cLtvpvn(nC);
  vector<double> xiXupun(nC), xiXtptn(nC), xiXvpvn(nC);
  vector<double> xiTupun(nC), xiTtptn(nC), xiTvpvn(nC);
  
  vector<double> c0xup(nC), c0xun(nC), c0xtp(nC), c0xtn(nC), c0xvp(nC), c0xvn(nC);
  vector<double> cLxup(nC), cLxun(nC), cLxtp(nC), cLxtn(nC), cLxvp(nC), cLxvn(nC);
  vector<double> c0tup(nC), c0tun(nC), c0ttp(nC), c0ttn(nC), c0tvp(nC), c0tvn(nC);
  vector<double> cLtup(nC), cLtun(nC), cLttp(nC), cLttn(nC), cLtvp(nC), cLtvn(nC);
  vector<double> xiXup(nC), xiXun(nC), xiXtp(nC), xiXtn(nC), xiXvp(nC), xiXvn(nC);
  vector<double> xiTup(nC), xiTun(nC), xiTtp(nC), xiTtn(nC), xiTvp(nC), xiTvn(nC);
  
  fill(cXupun.begin(),cXupun.end(),0); fill(cXtptn.begin(),cXtptn.end(),0); fill(cXvpvn.begin(),cXvpvn.end(),0);
  fill(cTupun.begin(),cTupun.end(),0);  fill(cTtptn.begin(),cTtptn.end(),0);  fill(cTvpvn.begin(),cTvpvn.end(),0);
  fill(c0xupun.begin(),c0xupun.end(),0); fill(c0xtptn.begin(),c0xtptn.end(),0);  fill(c0xvpvn.begin(),c0xvpvn.end(),0);
  fill(cLxupun.begin(),cLxupun.end(),0);  fill(cLxtptn.begin(),cLxtptn.end(),0);  fill(cLxvpvn.begin(),cLxvpvn.end(),0);
  fill(c0tupun.begin(),c0tupun.end(),0);  fill(c0ttptn.begin(),c0ttptn.end(),0);  fill(c0tvpvn.begin(),c0tvpvn.end(),0);
  fill(cLtupun.begin(),cLtupun.end(),0); fill(cLttptn.begin(),cLttptn.end(),0);  fill(cLtvpvn.begin(),cLtvpvn.end(),0);
  fill(xiXupun.begin(),xiXupun.end(),0);  fill(xiXtptn.begin(),xiXtptn.end(),0);  fill(xiXvpvn.begin(),xiXvpvn.end(),0);
  fill(xiTupun.begin(),xiTupun.end(),0); fill(xiTtptn.begin(),xiTtptn.end(),0);  fill(xiTvpvn.begin(),xiTvpvn.end(),0);

  fill(cXup.begin(),cXup.end(),0); fill(cXtp.begin(),cXtp.end(),0); fill(cXvp.begin(),cXvp.end(),0);
  fill(cXun.begin(),cXun.end(),0); fill(cXtn.begin(),cXtn.end(),0); fill(cXvn.begin(),cXvn.end(),0);

  fill(cTup.begin(),cTup.end(),0);  fill(cTtp.begin(),cTtp.end(),0);  fill(cTvp.begin(),cTvp.end(),0);
  fill(cTun.begin(),cTun.end(),0);  fill(cTtn.begin(),cTtn.end(),0);  fill(cTvn.begin(),cTvn.end(),0);

  fill(c0xup.begin(),c0xup.end(),0); fill(c0xtp.begin(),c0xtp.end(),0);  fill(c0xvp.begin(),c0xvp.end(),0);
  fill(c0xun.begin(),c0xun.end(),0); fill(c0xtn.begin(),c0xtn.end(),0);  fill(c0xvn.begin(),c0xvn.end(),0);

  fill(cLxup.begin(),cLxup.end(),0);  fill(cLxtp.begin(),cLxtp.end(),0);  fill(cLxvp.begin(),cLxvp.end(),0);
  fill(cLxun.begin(),cLxun.end(),0);  fill(cLxtn.begin(),cLxtn.end(),0);  fill(cLxvn.begin(),cLxvn.end(),0);

  fill(c0tup.begin(),c0tup.end(),0);  fill(c0ttp.begin(),c0ttp.end(),0);  fill(c0tvp.begin(),c0tvp.end(),0);
  fill(c0tun.begin(),c0tun.end(),0);  fill(c0ttn.begin(),c0ttn.end(),0);  fill(c0tvn.begin(),c0tvn.end(),0);

  fill(cLtup.begin(),cLtup.end(),0); fill(cLttp.begin(),cLttp.end(),0);  fill(cLtvp.begin(),cLtvp.end(),0);
  fill(cLtun.begin(),cLtun.end(),0); fill(cLttn.begin(),cLttn.end(),0);  fill(cLtvn.begin(),cLtvn.end(),0);

  fill(xiXup.begin(),xiXup.end(),0);  fill(xiXtp.begin(),xiXtp.end(),0);  fill(xiXvp.begin(),xiXvp.end(),0);
  fill(xiXun.begin(),xiXun.end(),0);  fill(xiXtn.begin(),xiXtn.end(),0);  fill(xiXvn.begin(),xiXvn.end(),0);

  fill(xiTup.begin(),xiTup.end(),0); fill(xiTtp.begin(),xiTtp.end(),0);  fill(xiTvp.begin(),xiTvp.end(),0);
  fill(xiTun.begin(),xiTun.end(),0); fill(xiTtn.begin(),xiTtn.end(),0);  fill(xiTvn.begin(),xiTvn.end(),0);

  vector<double> cXupun_(nX), cXtptn_(nX), cXvpvn_(nX);
  vector<double> cXupunerr_(nX), cXtptnerr_(nX), cXvpvnerr_(nX);
  vector<double> cTupun_(nT), cTtptn_(nT), cTvpvn_(nT);
  vector<double> cTupunerr_(nT), cTtptnerr_(nT), cTvpvnerr_(nT);

  vector<double> sQupun_(nX*nT), sQtptn_(nX*nT), sQvpvn_(nX*nT);
  vector<double> sQupunerr_(nX*nT), sQtptnerr_(nX*nT), sQvpvnerr_(nX*nT);

  vector<double> cXup_(nX), cXtp_(nX), cXvp_(nX);
  vector<double> cXuperr_(nX), cXtperr_(nX), cXvperr_(nX);
  vector<double> cTup_(nT), cTtp_(nT), cTvp_(nT);
  vector<double> cTuperr_(nT), cTtperr_(nT), cTvperr_(nT);

  vector<double> cXun_(nX), cXtn_(nX), cXvn_(nX);
  vector<double> cXunerr_(nX), cXtnerr_(nX), cXvnerr_(nX);
  vector<double> cTun_(nT), cTtn_(nT), cTvn_(nT);
  vector<double> cTunerr_(nT), cTtnerr_(nT), cTvnerr_(nT);

  double dubb;
  
  for(int ns = 0; ns < nSims; ns++) {
 
    stringstream oOpen;
    oOpen << dir << "/" << file << "_" << ns+1 << ".dat";
    
    ifstream data(oOpen.str());
    cout << oOpen.str() << endl;
    
    for(int m = 0; m < volPF; m++) {
      
      //data >> dubb1 >> dubb2 >> dubb3 >> dubb4 >> dubb5 >> dubb6;
      data >> dubb >> dubb >> dubb >> dubb >> dubb >> dubb >> dubb >> dubb >> dubb >> dubb >> dubb >> dubb >> dubb1 >> dubb2 >> dubb3 >> dubb4 >> dubb5 >> dubb6;
      
      cXTupun[m+ns*volPF]= dubb1 + dubb2;
      cXTtptn[m+ns*volPF]= dubb3 + dubb4;
      cXTvpvn[m+ns*volPF]= dubb5 + dubb6;

      cXTup[m+ns*volPF] = dubb1;
      cXTtp[m+ns*volPF] = dubb3;
      cXTvp[m+ns*volPF] = dubb5;

      cXTun[m+ns*volPF] = dubb2;
      cXTtn[m+ns*volPF] = dubb4;
      cXTvn[m+ns*volPF] = dubb6;
      
    } 
  
    data.close();
    
  }

  for (int c=0; c<nC; c++) {
    for (int x=0; x<nX; x++) {
      for (int t=0; t<nT; t++) {

        cXupun[x+c*nX]+=cXTupun[x+t*nX+c*nX*nT];
	cXtptn[x+c*nX]+=cXTtptn[x+t*nX+c*nX*nT];
	cXvpvn[x+c*nX]+=cXTvpvn[x+t*nX+c*nX*nT];

        cXup[x+c*nX]+=cXTup[x+t*nX+c*nX*nT];
	cXtp[x+c*nX]+=cXTtp[x+t*nX+c*nX*nT];
	cXvp[x+c*nX]+=cXTvp[x+t*nX+c*nX*nT];

        cXun[x+c*nX]+=cXTun[x+t*nX+c*nX*nT];
	cXtn[x+c*nX]+=cXTtn[x+t*nX+c*nX*nT];
	cXvn[x+c*nX]+=cXTvn[x+t*nX+c*nX*nT];
        
      }
    }
  }

  printJackknifedCorrelator(nC, &corrMeasSlice[0], nX, &cXupun[0], &cXupun_[0], &cXupunerr_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nX, &cXtptn[0], &cXtptn_[0], &cXtptnerr_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nX, &cXvpvn[0], &cXvpvn_[0], &cXvpvnerr_[0]);

  printJackknifedCorrelator(nC, &corrMeasSlice[0], nX, &cXup[0], &cXup_[0], &cXuperr_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nX, &cXtp[0], &cXtp_[0], &cXtperr_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nX, &cXvp[0], &cXvp_[0], &cXvperr_[0]);

  printJackknifedCorrelator(nC, &corrMeasSlice[0], nX, &cXun[0], &cXun_[0], &cXunerr_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nX, &cXtn[0], &cXtn_[0], &cXtnerr_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nX, &cXvn[0], &cXvn_[0], &cXvnerr_[0]);

  for (int c=0; c<nC; c++) {
    for (int t=0; t<nT; t++) {
      for (int x=0; x<nX; x++) {
       cTupun[t+c*nT]+=cXTupun[x+t*nX+c*nX*nT];
       cTtptn[t+c*nT]+=cXTtptn[x+t*nX+c*nX*nT];
       cTvpvn[t+c*nT]+=cXTvpvn[x+t*nX+c*nX*nT];

       cTup[t+c*nT]+=cXTup[x+t*nX+c*nX*nT];
       cTtp[t+c*nT]+=cXTtp[x+t*nX+c*nX*nT];
       cTvp[t+c*nT]+=cXTvp[x+t*nX+c*nX*nT];

       cTun[t+c*nT]+=cXTun[x+t*nX+c*nX*nT];
       cTtn[t+c*nT]+=cXTtn[x+t*nX+c*nX*nT];
       cTvn[t+c*nT]+=cXTvn[x+t*nX+c*nX*nT];

       }
    }
  }

  printJackknifedCorrelator(nC, &corrMeasSlice[0], nT, &cTupun[0], &cTupun_[0], &cTupunerr_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nT, &cTtptn[0], &cTtptn_[0], &cTtptnerr_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nT, &cTvpvn[0], &cTvpvn_[0], &cTvpvnerr_[0]);

  printJackknifedCorrelator(nC, &corrMeasSlice[0], nT, &cTup[0], &cTup_[0], &cTuperr_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nT, &cTtp[0], &cTtp_[0], &cTtperr_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nT, &cTvp[0], &cTvp_[0], &cTvperr_[0]);

  printJackknifedCorrelator(nC, &corrMeasSlice[0], nT, &cTun[0], &cTun_[0], &cTunerr_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nT, &cTtn[0], &cTtn_[0], &cTtnerr_[0]);
  printJackknifedCorrelator(nC, &corrMeasSlice[0], nT, &cTvn[0], &cTvn_[0], &cTvnerr_[0]);
 
  
  for (int c=0; c<nC; c++) {
    for (int x=0; x<nX; x++) {
      c0xupun[c] += cXupun[x+c*nX];
      cLxupun[c] += cos(2*pI*x/nX)*cXupun[x+c*nX];
      c0xtptn[c] += cXtptn[x+c*nX];
      cLxtptn[c] += cos(2*pI*x/nX)*cXtptn[x+c*nX];
      c0xvpvn[c] += cXvpvn[x+c*nX];
      cLxvpvn[c] += cos(2*pI*x/nX)*cXvpvn[x+c*nX];

      c0xup[c] += cXup[x+c*nX];
      cLxup[c] += cos(2*pI*x/nX)*cXup[x+c*nX];
      c0xtp[c] += cXtp[x+c*nX];
      cLxtp[c] += cos(2*pI*x/nX)*cXtp[x+c*nX];
      c0xvp[c] += cXvp[x+c*nX];
      cLxvp[c] += cos(2*pI*x/nX)*cXvp[x+c*nX];

      c0xun[c] += cXun[x+c*nX];
      cLxun[c] += cos(2*pI*x/nX)*cXun[x+c*nX];
      c0xtn[c] += cXtn[x+c*nX];
      cLxtn[c] += cos(2*pI*x/nX)*cXtn[x+c*nX];
      c0xvn[c] += cXvn[x+c*nX];
      cLxvn[c] += cos(2*pI*x/nX)*cXvn[x+c*nX];
    }
  }


  for (int c=0; c<nC; c++) {
    for (int t=0; t<nT; t++) {
      c0tupun[c] += cTupun[t+c*nT];
      cLtupun[c] += cos(2*pI*t/nT)*cTupun[t+c*nT];
      c0ttptn[c] += cTtptn[t+c*nT];
      cLttptn[c] += cos(2*pI*t/nT)*cTtptn[t+c*nT];
      c0tvpvn[c] += cTvpvn[t+c*nT];
      cLtvpvn[c] += cos(2*pI*t/nT)*cTvpvn[t+c*nT];

      c0tup[c] += cTup[t+c*nT];
      cLtup[c] += cos(2*pI*t/nT)*cTup[t+c*nT];
      c0ttp[c] += cTtp[t+c*nT];
      cLttp[c] += cos(2*pI*t/nT)*cTtp[t+c*nT];
      c0tvp[c] += cTvp[t+c*nT];
      cLtvp[c] += cos(2*pI*t/nT)*cTvp[t+c*nT];

      c0tun[c] += cTun[t+c*nT];
      cLtun[c] += cos(2*pI*t/nT)*cTun[t+c*nT];
      c0ttn[c] += cTtn[t+c*nT];
      cLttn[c] += cos(2*pI*t/nT)*cTtn[t+c*nT];
      c0tvn[c] += cTvn[t+c*nT];
      cLtvn[c] += cos(2*pI*t/nT)*cTvn[t+c*nT];
    }
  }
  
  for (int c=0; c<nC; c++) {
    xiXupun[c]=(nX/(2*pI))*sqrt((c0xupun[c]/cLxupun[c])-1);
    xiXtptn[c]=(nX/(2*pI))*sqrt((c0xtptn[c]/cLxtptn[c])-1);
    xiXvpvn[c]=(nX/(2*pI))*sqrt((c0xvpvn[c]/cLxvpvn[c])-1);

    xiXup[c]=(nX/(2*pI))*sqrt((c0xup[c]/cLxup[c])-1);
    xiXtp[c]=(nX/(2*pI))*sqrt((c0xtp[c]/cLxtp[c])-1);
    xiXvp[c]=(nX/(2*pI))*sqrt((c0xvp[c]/cLxvp[c])-1);

    xiXun[c]=(nX/(2*pI))*sqrt((c0xun[c]/cLxun[c])-1);
    xiXtn[c]=(nX/(2*pI))*sqrt((c0xtn[c]/cLxtn[c])-1);
    xiXvn[c]=(nX/(2*pI))*sqrt((c0xvn[c]/cLxvn[c])-1);
  }

  
  //cLx/=cX[0];
  for (int c=0; c<nC; c++) {	
    xiTupun[c]=(nT/(2*pI))*sqrt((c0tupun[c]/cLtupun[c])-1); 
    xiTtptn[c]=(nT/(2*pI))*sqrt((c0ttptn[c]/cLttptn[c])-1); 
    xiTvpvn[c]=(nT/(2*pI))*sqrt((c0tvpvn[c]/cLtvpvn[c])-1); 

    xiTup[c]=(nT/(2*pI))*sqrt((c0tup[c]/cLtup[c])-1); 
    xiTtp[c]=(nT/(2*pI))*sqrt((c0ttp[c]/cLttp[c])-1); 
    xiTvp[c]=(nT/(2*pI))*sqrt((c0tvp[c]/cLtvp[c])-1); 

    xiTun[c]=(nT/(2*pI))*sqrt((c0tun[c]/cLtun[c])-1); 
    xiTtn[c]=(nT/(2*pI))*sqrt((c0ttn[c]/cLttn[c])-1); 
    xiTvn[c]=(nT/(2*pI))*sqrt((c0tvn[c]/cLtvn[c])-1); 

   }

  JK xJackupun(nC,1,0), tJackupun(nC,1,0);
  JK xJacktptn(nC,1,0), tJacktptn(nC,1,0);
  JK xJackvpvn(nC,1,0), tJackvpvn(nC,1,0);

  JK xJackup(nC,1,0), tJackup(nC,1,0);
  JK xJacktp(nC,1,0), tJacktp(nC,1,0);
  JK xJackvp(nC,1,0), tJackvp(nC,1,0);

  JK xJackun(nC,1,0), tJackun(nC,1,0);
  JK xJacktn(nC,1,0), tJacktn(nC,1,0);
  JK xJackvn(nC,1,0), tJackvn(nC,1,0);

  xJackupun.jackknife(&xiXupun[0]);
  double xiXupun_ = xJackupun.jackMean();
  double xiXupunerr_ = xJackupun.jackSE();
  xJacktptn.jackknife(&xiXtptn[0]);
  double xiXtptn_ = xJacktptn.jackMean();
  double xiXtptnerr_ = xJacktptn.jackSE();
  xJackvpvn.jackknife(&xiXvpvn[0]);
  double xiXvpvn_ = xJackvpvn.jackMean();
  double xiXvpvnerr_ = xJackvpvn.jackSE();

  xJackup.jackknife(&xiXup[0]);
  double xiXup_ = xJackup.jackMean();
  double xiXuperr_ = xJackup.jackSE();
  xJacktp.jackknife(&xiXtp[0]);
  double xiXtp_ = xJacktp.jackMean();
  double xiXtperr_ = xJacktp.jackSE();
  xJackvp.jackknife(&xiXvp[0]);
  double xiXvp_ = xJackvp.jackMean();
  double xiXvperr_ = xJackvp.jackSE();

  xJackun.jackknife(&xiXun[0]);
  double xiXun_ = xJackun.jackMean();
  double xiXunerr_ = xJackun.jackSE();
  xJacktn.jackknife(&xiXtn[0]);
  double xiXtn_ = xJacktn.jackMean();
  double xiXtnerr_ = xJacktn.jackSE();
  xJackvn.jackknife(&xiXvn[0]);
  double xiXvn_ = xJackvn.jackMean();
  double xiXvnerr_ = xJackvn.jackSE();


  tJackupun.jackknife(&xiTupun[0]);
  double xiTupun_ = tJackupun.jackMean();
  double xiTupunerr_ = tJackupun.jackSE();
  tJacktptn.jackknife(&xiTtptn[0]);
  double xiTtptn_ = tJacktptn.jackMean();
  double xiTtptnerr_ = tJacktptn.jackSE();
  tJackvpvn.jackknife(&xiTvpvn[0]);
  double xiTvpvn_ = tJackvpvn.jackMean();
  double xiTvpvnerr_ = tJackvpvn.jackSE();

  tJackup.jackknife(&xiTup[0]);
  double xiTup_ = tJackup.jackMean();
  double xiTuperr_ = tJackup.jackSE();
  tJacktp.jackknife(&xiTtp[0]);
  double xiTtp_ = tJacktp.jackMean();
  double xiTtperr_ = tJacktp.jackSE();
  tJackvp.jackknife(&xiTvp[0]);
  double xiTvp_ = tJackvp.jackMean();
  double xiTvperr_ = tJackvp.jackSE();

  tJackun.jackknife(&xiTun[0]);
  double xiTun_ = tJackun.jackMean();
  double xiTunerr_ = tJackun.jackSE();
  tJacktn.jackknife(&xiTtn[0]);
  double xiTtn_ = tJacktn.jackMean();
  double xiTtnerr_ = tJacktn.jackSE();
  tJackvn.jackknife(&xiTvn[0]);
  double xiTvn_ = tJackvn.jackMean();
  double xiTvnerr_ = tJackvn.jackSE();

  stringstream wWrite;
  wWrite << boost::format("%s/xi_All_%ss_%d.dat") % dir % file % nSims;

  cout << wWrite.str() << endl;

  ofstream outData(wWrite.str());
  outData << boost::format("%f   %f   %f\n upun:     %.4f   %.4f   %.4f   %.4f\n tptn:     %.4f   %.4f   %.4f   %.4f\n vpvn:     %.4f   %.4f   %.4f   %.4f\n up:     %.4f   %.4f   %.4f   %.4f\n tp:     %.4f   %.4f   %.4f   %.4f\n vp:     %.4f   %.4f   %.4f   %.4f\n un:     %.4f   %.4f   %.4f   %.4f\n tn:     %.4f   %.4f   %.4f   %.4f\n vn:     %.4f   %.4f   %.4f   %.4f\n") % beta % nX % nT % xiXupun_ % xiXupunerr_ % xiTupun_ % xiTupunerr_ % xiXtptn_ % xiXtptnerr_ % xiTtptn_ % xiTtptnerr_ % xiXvpvn_ % xiXvpvnerr_ % xiTvpvn_ % xiTvpvnerr_ % xiXup_ % xiXuperr_ % xiTup_ % xiTuperr_ % xiXtp_ % xiXtperr_ % xiTtp_ % xiTtperr_ % xiXvp_ % xiXvperr_ % xiTvp_ % xiTvperr_ % xiXun_ % xiXunerr_ % xiTun_ % xiTunerr_ % xiXtn_ % xiXtnerr_ % xiTtn_ % xiTtnerr_ % xiXvn_ % xiXvnerr_ % xiTvn_ % xiTvnerr_;


 cout << boost::format("%f   %f   %f\n upun:     %.4f   %.4f   %.4f   %.4f\n tptn:     %.4f   %.4f   %.4f   %.4f\n vpvn:     %.4f   %.4f   %.4f   %.4f\n up:     %.4f   %.4f   %.4f   %.4f\n tp:     %.4f   %.4f   %.4f   %.4f\n vp:     %.4f   %.4f   %.4f   %.4f\n un:     %.4f   %.4f   %.4f   %.4f\n tn:     %.4f   %.4f   %.4f   %.4f\n vn:     %.4f   %.4f   %.4f   %.4f\n") % beta % nX % nT % xiXupun_ % xiXupunerr_ % xiTupun_ % xiTupunerr_ % xiXtptn_ % xiXtptnerr_ % xiTtptn_ % xiTtptnerr_ % xiXvpvn_ % xiXvpvnerr_ % xiTvpvn_ % xiTvpvnerr_ % xiXup_ % xiXuperr_ % xiTup_ % xiTuperr_ % xiXtp_ % xiXtperr_ % xiTtp_ % xiTtperr_ % xiXvp_ % xiXvperr_ % xiTvp_ % xiTvperr_ % xiXun_ % xiXunerr_ % xiTun_ % xiTunerr_ % xiXtn_ % xiXtnerr_ % xiTtn_ % xiTtnerr_ % xiXvn_ % xiXvnerr_ % xiTvn_ % xiTvnerr_;

 outData << "cX" << endl;
 for(int x=0; x<nX; x++) {
   outData << boost::format("%.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f\n") % cXupun_[x] % cXupunerr_[x] % cXtptn_[x] % cXtptnerr_[x]  % cXvpvn_[x] % cXvpvnerr_[x] % cXup_[x] % cXuperr_[x] % cXtp_[x] % cXtperr_[x]  % cXvp_[x] % cXvperr_[x] % cXun_[x] % cXunerr_[x] % cXtn_[x] % cXtnerr_[x]  % cXvn_[x] % cXvnerr_[x];
 }
 
 outData << "cT" << endl;
 for(int t=0; t<nT; t++) {
   outData << boost::format("%.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f\n") % cTupun_[t] % cTupunerr_[t] % cTtptn_[t] % cTtptnerr_[t] % cTvpvn_[t] % cTvpvnerr_[t] % cTup_[t] % cTuperr_[t] % cTtp_[t] % cTtperr_[t] % cTvp_[t] % cTvperr_[t]  % cTun_[t] % cTunerr_[t] % cTtn_[t] % cTtnerr_[t] % cTvn_[t] % cTvnerr_[t];
 }

 
 outData.close();
  

}
