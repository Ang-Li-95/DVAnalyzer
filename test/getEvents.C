#include<assert.h>

void calcYields(TString path17, TString path18, vector<TString> filename, vector<double> weights17, vector<double> weights18, TString hist)
{
  assert (filename.size()==weights.size());
  for (size_t i=0; i<filename.size(); ++i){
    TFile *f17 = new TFile(path17 + filename[i]);
    TH1F* h17 = (TH1F*)f17->Get(hist);
    TFile *f18 = new TFile(path18 + filename[i]);
    TH1F* h18 = (TH1F*)f18->Get(hist);

    
    //double sum = h->Integral(6,40);
    std::cout << filename[i] 
      << " yield: 1 vtx: " << (h17->GetBinContent(h17->FindBin(1)))*weights17[i] + (h18->GetBinContent(h18->FindBin(1)))*weights18[i]
      << " 2vtx: " << (h17->GetBinContent(h17->FindBin(2)))*weights17[i] + (h18->GetBinContent(h18->FindBin(2)))*weights18[i]
      << std::endl;
    
  }
  
}




void getEvents()
{
  const TString path17 = "../../../../../CMSSW_9_4_15/src/LLPAnalyzer/DVAnalyzer/test/";
  const TString path18 = "./";
  const vector<TString> filename {
                               "ggToNN_800M_1mm.root",
                               "QCD_HT700to1000.root",
                               "QCD_HT1000to1500.root",
                               "QCD_HT1500to2000.root",
                               "QCD_HT2000toInf.root",
                               "TTJets600To800.root",
                               "TTJets800To1200.root",
                               "TTJets1200To2500.root",
                               "TTJets2500ToInf.root"
  };
  const vector<double> weights17 {4.15e-04, 5.49, 2.70, 0.353, 0.141, 9.25e-04, 7.76e-04, 4.13e-04, 1.14e-05}; //2017
  const vector<double> weights18 {5.97e-04, 8.74, 4.33, 0.535, 0.218, 7.6e-03, 4.3e-03, 2.7e-03, 5.8e-05}; //2018
  const vector<TString> hists {"nvtx_per_event_3tks", "nvtx_per_event_4tks", "nvtx_per_event_5tks"};
  for (auto &hist:hists){
    std::cout << hist << std::endl;
    calcYields(path17, path18, filename, weights17, weights18, hist);
    std::cout << std::endl;
  }
}
