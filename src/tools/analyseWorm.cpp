#include <algorithm> 
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <array>
#include <string>
#include <cstdio>
#include <ctime>
#include <sstream>
#include <format>
#include <cp2worm/jackknife.h>   

using namespace std;

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

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void printHistPoints(string name, int L, int* hist, ofstream* out) {

  int val;
  *out << endl << name << endl;
  
  for(int j = 0; j < L; j++) {
    for(int i = 0; i < L; i++) {

      val = *(hist + i + (j*L));
      
      if(val>0) {
	*out << std::format("{:<20}{:<20}{:<20}\n", val, i, j);
      }
    }
  }
}

void printHistSum(string name, int L, int* hist1, int* hist2, ofstream* out) {

  int val;
  *out << endl << name << endl;
  
  for(int j = 0; j < L; j++) {
    for(int i = 0; i < L; i++) {

      val = *(hist1 + i + (j*L)) + *(hist2 + i + (j*L));
      
      if(val>0) {
	*out << std::format("{:<20}{:<20}{:<20}\n", val, i, j);
      }
    }
  }
}

int main(int argc, char** argv) {

  if(argc==1) {
    cout << "   string dir = argv[1];\n  string file = argv[2];\n  int skip = atoi(argv[3]);" << endl;
    return 0;
  }

  string dir = argv[1];
  string file = argv[2];
  int skip = atoi(argv[3]);
  int binJK = atoi(argv[4]);
  int nMeasMax = atoi(argv[5]);
  int histX = atoi(argv[6]);
  int histHalf = (histX - 1)/2;

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
  int nT = stoi(splitFile[11]);
  cout << "nT: " << nT << endl;

  int nn = nMeasMax-skip;
  string line, shifter11;
  double dubb1, dubb2, dubb3, dubb4, dubb5, dubb6, dubb7, dubb8, dubb9, dubb10;
  
  vector<double> vec1(nn), vec2(nn), vec3(nn), vec4(nn), vec5(nn), vec6(nn), vec7(nn), vec8(nn), vec9(nn), vec10(nn), vec11(nn), vec12(nn), vec13(nn);
  int int12, int13, q3prev = 0, q8prev = 0, x3prev = 0, x8prev = 0, n = 0, p = 0;
  stringstream oOpen;

  int xs = 0, ys = 0, xps = 0, yps = 0, xpps = 0, ypps = 0, Xxpps = 0, Xxs = 0, Xys = 0, Xypps = 0, flag = 0, token = 0;
  vector<int> chiTp_(histX*histX), chiTm_(histX*histX), chiVp_(histX*histX), chiVm_(histX*histX), chiUp_(histX*histX), chiUm_(histX*histX);
  vector<int> stateHist(histX*histX), stateHistx3(histX*histX), stateHistx8(histX*histX);
      
  oOpen << dir << "/" << argv[2];

  ifstream data(oOpen.str());
   
  cout << oOpen.str() << endl;

  for(int m = 0; m < skip; m++) {
    getline(data, line);
  }
  for(int m = 0; m < nn; m++) {
    // Wx3 Wx8 Wy3 Wy8 Q3 Q8 l lp lm shift_operator
    getline(data, line);
    stringstream ss( line );
    vector<string> record;
    
    while ( ss.good() )
      {
        string substr;
        getline( ss, substr, ',' );
        record.push_back( substr );
      }

    // data >> ',' >> dubb1 >> ',' >> dubb2 >> ',' >> dubb3 >> ',' >> dubb4 >> ',' >> dubb5 >> ',' >> dubb6 >> ',' >> dubb7 >> ',' >> dubb8 >> ',' >> dubb9 >> ',' >> shifter11;
    
    dubb1 = std::stod(record[0]);
    dubb2 = std::stod(record[1]);
    dubb3 = std::stod(record[2]);
    dubb4 = std::stod(record[3]);
    dubb6 = std::stod(record[5]);
    dubb7 = std::stod(record[6]);
    dubb8 = std::stod(record[7]);
    dubb9 = std::stod(record[8]);
    shifter11 = record[9];

    cout << ' ' << dubb1 << ' ' << dubb2 << ' ' << dubb3 << ' ' << dubb4 << ' ' << dubb5 << ' ' << dubb6 << ' ' << dubb7 << ' ' << dubb8 << ' ' << dubb9 << ' ' << shifter11 << endl;

    vec1[m]=dubb1*dubb1;
    vec2[m]=dubb2*dubb2;
    vec3[m]=dubb3*dubb3;
    vec4[m]=dubb4*dubb4;
    vec5[m]=dubb5*dubb5;
    vec6[m]=dubb6*dubb6;
    vec7[m]=dubb7;
    vec8[m]=dubb5;
    vec9[m]=dubb6;   

    if(q3prev!=int(dubb5) || q8prev!=int(dubb6)) {
      	vec12[n]=dubb1*dubb1;
	vec13[n]=dubb2*dubb2;
	p++;
      if(q3prev==1) {
	vec10[n]=dubb1*dubb1;
	vec11[n]=dubb2*dubb2;
	n++;
      }
    }

    xs = q3prev;
    ys = q8prev;
 
    q3prev = int(dubb5);
    q8prev = int(dubb6);
  
    xps = xs + histHalf;
    yps = ys + histHalf;

    if(xps > histX || xps < 0 || yps > histX || yps < 0) {
      flag = 1;
    }

    if(shifter11=="tp"){
      chiTp_[xps + (yps*histX)] += dubb8;
      chiTm_[xps + (yps*histX)] += dubb9;
    }
    else if(shifter11=="tm"){
      chiTm_[xps + (yps*histX)] += dubb8;
      chiTp_[xps + (yps*histX)] += dubb9;		  
    }
    else if(shifter11=="vp"){
      chiVp_[xps + (yps*histX)] += dubb8;
      chiVm_[xps + (yps*histX)] += dubb9;
    }
    else if(shifter11=="vm"){
      chiVm_[xps + (yps*histX)] += dubb8;
      chiVp_[xps + (yps*histX)] += dubb9;
    }
    else if(shifter11=="up"){
      chiUp_[xps + (yps*histX)] += dubb8;
      chiUm_[xps + (yps*histX)] += dubb9;
    }
    else if(shifter11=="um"){
      chiUm_[xps + (yps*histX)] += dubb8;
      chiUp_[xps + (yps*histX)] += dubb9;
    }
  
    xpps = q3prev + histHalf;
    ypps = q8prev + histHalf;

    if(xpps > histX || xpps < 0 || ypps > histX || ypps < 0) {
      token = 1;
    }
        
    stateHist[xpps + (ypps*histX)]++;
    stateHistx3[xpps + (ypps*histX)] += dubb1*dubb1;
    stateHistx8[xpps + (ypps*histX)] += dubb2*dubb2;
  }
 
  JK jack1(nn, binJK, 0), jack2(nn, binJK, 0), jack3(nn, binJK, 0), jack4(nn, binJK, 0), jack5(nn, binJK, 0), jack6(nn, binJK, 0), jack7(nn, binJK, 0), jack8(nn, binJK, 0),  jack9(nn, binJK, 0), jack10(n, 1, 0), jack11(n, 1, 0), jack12(p, 1, 0), jack13(p, 1, 0);

  jack1.jackknife(&vec1[0]);
  jack2.jackknife(&vec2[0]);
  jack3.jackknife(&vec3[0]);
  jack4.jackknife(&vec4[0]);
  jack5.jackknife(&vec5[0]);
  jack6.jackknife(&vec6[0]);
  jack7.jackknife(&vec7[0]);
  jack8.jackknife(&vec8[0]);
  jack9.jackknife(&vec9[0]);
  jack10.jackknife(&vec10[0]);
  jack11.jackknife(&vec11[0]);
  jack12.jackknife(&vec12[0]);
  jack13.jackknife(&vec13[0]);

  
  double W1 = jack1.jackMean();
  double W1err = jack1.jackSE();
  
  double W2 = 3*jack2.jackMean();
  double W2err = 3*jack2.jackSE();
  
  double w3 = jack3.jackMean();
  double w3err = jack3.jackSE();

  double w4 = 3*jack4.jackMean();
  double w4err = 3*jack4.jackSE();

  double w5 = jack5.jackMean();
  double w5err = jack5.jackSE();

  double w6 =  3*jack6.jackMean();
  double w6err = 3*jack6.jackSE();
  
  double w7 = jack7.jackMean();
  double w7err = jack7.jackSE();
  
  double w8 = jack8.jackMean();
  double w8err = jack8.jackSE();
  
  double w9 = sqrt(3)*jack9.jackMean();
  double w9err = sqrt(3)*jack9.jackSE();

  double w10 = jack10.jackMean();
  double w10err = jack10.jackSE();
  
  double w11 = 3*jack11.jackMean();
  double w11err = 3*jack11.jackSE();

  double w12 = jack12.jackMean();
  double w12err = jack12.jackSE();
  
  double w13 = 3*jack13.jackMean();
  double w13err = 3*jack13.jackSE();

  data.close();

  double total10 = 0.0, total11 = 0.0;
  for(int i = 0; i<n; i++) {
    total10 += vec10[i];
    total11 += vec11[i];
  }

  stringstream wWrite;
  wWrite << std::format("{}/analysis_skip_{}_nMeasMax_{}_{}", dir, skip, nMeasMax, file);

  cout << wWrite.str() << endl;

  ofstream outData(wWrite.str());
  outData << std::format("{:<13.6f}{:<13.6f}{:<13.6f}{:<13.6f}{:<13.6f}{:<13.6f}{:<13.6f}{:<13.6f}{:<13.6f}{:<13.6f}{:<13.6f}{:<13.6f}{:<13.6f}{:<13.6f}{:<13.6f}{:<13.6f}{:<13.0f}{:<13.0f}{:<13.6f}{:<13.6f}{:<13.6f}{:<13.6f}{:<13.6f}{:<13.6f}{:<13.6f}{:<13.6f}{:<13}{:<13.1f}{:<13.1f}\n", w5, w5err, w6, w6err, w8, w8err, w9, w9err, W1, W1err, W2, W2err, w3, w3err, w4, w4err, w7, w7err, w10, w10err, w11, w11err, w12, w12err, w13, w13err, n, total10, total11);

  outData << endl;

  auto minmax5 = minmax_element(vec8.begin(),vec8.end());
  auto minmax6 = minmax_element(vec9.begin(),vec9.end());

  outData << "q3 min: " << vec8[(minmax5.first - vec8.begin())] << ", q3 max: " << vec8[(minmax5.second - vec8.begin())] << endl;
  outData << "q8 min: " << vec9[(minmax6.first - vec9.begin())] << ", q8 max: " << vec9[(minmax6.second - vec9.begin())] << endl;

  outData << endl;
 
  if(token==1) { outData << "token histogram too small" << endl; }
  if(flag==1) { outData << "flag histogram too small" << endl; }
 
  outData << endl << "Histogram Points" << endl;

  printHistPoints("States", histX, &stateHist[0], &outData);
  printHistPoints("Tp", histX, &chiTp_[0], &outData);
  printHistPoints("Tm", histX, &chiTm_[0], &outData);
  printHistPoints("Up", histX, &chiUp_[0], &outData);
  printHistPoints("Um", histX, &chiUm_[0], &outData);
  printHistPoints("Vp", histX, &chiVp_[0], &outData);
  printHistPoints("Vm", histX, &chiVm_[0], &outData);
  printHistPoints("statesX3", histX, &stateHistx3[0], &outData);
  printHistPoints("statesX8", histX, &stateHistx8[0], &outData);
  printHistSum("TpAddTm", histX, &chiTp_[0], &chiTm_[0], &outData);

  outData.close();

}
