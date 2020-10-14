#include<assert.h>

void calcYields(vector<TString> filename, vector<double> weights, TString hist)
{
  assert (filename.size()==weights.size());
  for (size_t i=0; i<filename.size(); ++i){
    TFile *f = new TFile(filename[i]);
    TH1F* h = (TH1F*)f->Get(hist);
    
    //double sum = h->Integral(6,40);
    std::cout << filename[i] << " yield: 1 vtx: " << (h->GetBinContent(h->FindBin(1)))*weights[i] << " 2vtx: " << (h->GetBinContent(h->FindBin(2)))*weights[i] << std::endl;
    
  }
  
}




void getYeilds()
{
  const vector<TString> filename {
                               "ggToNN_800M_1mm.root",
                               "QCD_HT700to1000.root",
                               "QCD_HT1000to1500.root",
                               "QCD_HT1500to2000.root",
                               "QCD_HT2000toInf.root",
                               "TTJets_HT600To800.root",
                               "TTJets_HT800To1200.root",
                               "TTJets_HT1200To2500.root",
                               "TTJets_HT2500ToInf.root"
  };
  //const vector<double> weights {4.15e-04, 5.49, 2.70, 0.353, 0.141, 9.25e-04, 7.76e-04, 4.13e-04, 1.14e-05}; //2017
  const vector<double> weights {5.97e-04, 8.74, 4.33, 0.535, 0.218, 7.6e-03, 4.3e-03, 2.7e-03, 5.8e-05}; //2018
  for (size_t i=0; i<filename.size(); ++i){
    TFile *f = new TFile(filename[i]);
    TH1F* h = (TH1F*)f->Get("vtx_tkSize");
    double sum = h->Integral(6,40);
    std::cout << filename[i] << " # vertices: unweightaed: " << sum << " weighted: " << weights[i]*sum << std::endl;
  }
}
