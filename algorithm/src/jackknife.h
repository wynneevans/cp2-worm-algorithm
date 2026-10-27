#ifndef JACKKNIFE
#define JACKKNIFE

using namespace std;

class JK {
    int nMeas_, binSize_, skip_, nBin_, nBinSkip_;
    double jackMean_, jackSE_, total_, binTotal_;
    vector<double> means_;
    public:
    JK();
    JK(int nMeas, int binSize, int skip);
    ~JK();
    void jackknife(double *);
    double jackMean();
    double jackSE();
};

JK::JK() : nMeas_(0), binSize_(0), skip_(0), nBin_(0), nBinSkip_(0), jackMean_(0), jackSE_(0), total_(0), means_(0) { }

JK::JK(int nMeas, int binSize, int skip) : nMeas_(nMeas), binSize_(binSize), skip_(skip), nBin_(nMeas_/binSize_), nBinSkip_(skip_/binSize_), jackMean_(0), jackSE_(0), total_(0), means_(nBin_)
{
    fill(means_.begin(), means_.end(), 0);
}

JK::~JK(void) {
    //cout << "Destructor called" << endl;
}

void JK::jackknife(double *O) {
    O = O + skip_;
    //calculate sum of all measurements
    for (int i = 0; i < nMeas_-skip_; i++) {
        total_ = total_ + (*O);
        O++;
    }
    
    //reset the pointer
    O = O - nMeas_ + skip_;
    
    //calculate the jackknife means by subtracting the ith measurement from total
    //and calculatethe overall jackknife mean.
    for (int i = 0; i < nBin_-nBinSkip_; i++) {
        binTotal_ = 0;
        for (int j = 0; j < binSize_; j++) {
            binTotal_ = binTotal_ + (*O);
            O++;
        }
        
        
        means_[i] = (total_ - binTotal_)/((nBin_-1-nBinSkip_)*binSize_);
	//	cout << means_[i] << endl;
        jackMean_ += means_[i];
       
    
    }
    jackMean_ = jackMean_/(nBin_-nBinSkip_);
    //    cout << "b: " << jackMean_ << endl;
    
    //calculate the jackknife standard error
    for (int i = 0; i < nBin_-nBinSkip_; i++) {
        jackSE_ += (means_[i]-jackMean_)*(means_[i]-jackMean_);
    }
    jackSE_ = sqrt(((nBin_-1-nBinSkip_)*jackSE_)/(nBin_-nBinSkip_));
    
}

double JK::jackMean() {
    return jackMean_;
}

double JK::jackSE() {
  return jackSE_;
}


#endif
