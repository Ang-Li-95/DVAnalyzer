void TwoPlotsOnCanvas(TFile* f_1, TFile* f_2, TString name, bool is1Signal=true, double dlimit=-999, double ulimit=-999)
{
  TH1F *h_1 = (TH1F*)f_1->Get(name);
  TH1F *h_2 = (TH1F*)f_2->Get(name);
  TCanvas *c = new TCanvas("c"+name, "c"+name, 800,600);
  c->cd();
  h_1->Scale(1./(h_1->Integral()));
  h_2->Scale(1./(h_2->Integral()));
  if(is1Signal){
    h_2->SetLineColor(kBlue);
    h_1->SetLineColor(kRed);
  }
  else{
    h_1->SetLineColor(kBlue);
    h_2->SetLineColor(kRed);
  }
  if(ulimit!=-999 && dlimit!=-999){
    h_1->GetXaxis()->SetRangeUser(dlimit, ulimit);
    h_2->GetXaxis()->SetRangeUser(dlimit, ulimit);
  }
  h_1->Draw("hist");
  h_2->Draw("hist same");
  //c->BuildLegend();
  auto l = new TLegend(0.8,0.7,0.98,0.90);
  if(is1Signal){
    l->AddEntry(h_1,"signal");
    l->AddEntry(h_2,"background");
  }
  else{
    l->AddEntry(h_1, "background");
    l->AddEntry(h_2, "signal");
  }
  l->Draw();
}

void TwoPlots2D(TFile* f_1, TFile* f_2, TString name)
{
  TH2F *h_1 = (TH2F*)f_1->Get(name);
  TH2F *h_2 = (TH2F*)f_2->Get(name);
  TCanvas *c = new TCanvas("c"+name, "c"+name, 800,800);
  c->cd();
  TH2F *ts = new TH2F("ts","vertex y vs. x", 80, -4, 4, 80, -4, 4);
  ts->Add(h_1);
  ts->Add(h_2);
  ts->Draw("colz");
}



void ComparePlot_vtx()
{
  gStyle->SetOptStat(0);
  TString sig = "ggToNN_800M_1mm.root";
  TString bkg = "QCD_HT1000to1500.root";
  TFile *f_sig = new TFile(sig);
  TFile *f_bkg = new TFile(bkg);
  TwoPlotsOnCanvas(f_bkg, f_sig, "vtx_tkSize", false);
  TwoPlotsOnCanvas(f_bkg, f_sig, "vtx_dBV", false);
  TwoPlotsOnCanvas(f_sig, f_bkg, "vtx_sigma_dBV", true);
  TwoPlots2D(f_bkg, f_sig, "vtx_xy");
}
