#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TCanvas.h>
#include <TChain.h>

void IllustrateStubs()
{
  TChain* t_L1TrackTrigger = new TChain("L1TrackTrigger");
  t_L1TrackTrigger->Add("/Users/souvik/CMSPhase2Upgrades/Samples/SLHC/GEN/612_SLHC6_MU/MU_612_SLHC6.root");
  
  int STUB_n;
  std::vector<float> *STUB_x=0;
  std::vector<float> *STUB_y=0;
  std::vector<float> *STUB_z=0;
  std::vector<int> *STUB_layer=0;
  
  t_L1TrackTrigger->SetBranchAddress("STUB_n", &STUB_n);
  t_L1TrackTrigger->SetBranchAddress("STUB_x", &STUB_x);
  t_L1TrackTrigger->SetBranchAddress("STUB_z", &STUB_z);
  t_L1TrackTrigger->SetBranchAddress("STUB_y", &STUB_y);
  t_L1TrackTrigger->SetBranchAddress("STUB_layer", &STUB_layer);
  
  t_L1TrackTrigger->SetBranchAddress("STUB_layer", &STUB_layer);
  
  TH2F *h_Barrel_STUB_x_y=new TH2F("h_STUB_x_y", "STUB x vs y", 1000, -120, 120, 1000, -120, 120);
  TH2F *h_STUB_r_z=new TH2F("h_STUB_r_z", "h_STUB_r_z", 1000, -300, 300, 1000, 0, 120);
  // TH3F *h_STUB_3D=new TH3F("h_STUB_3D", "h_STUB_3D", 1000, -300, 300, 1000, -120, 120, 1000, -120, 120);
  
  TH2F *h_STUB_r_z_hiDeg=(TH2F*)h_STUB_r_z->Clone("h_STUB_r_z_hiDeg");
  
  TH1F *h_STUB_n=new TH1F("h_STUB_n", "Number of Stubs in an Event; Number of Stubs", 100, 0., 100.);
  
  int nEvents=t_L1TrackTrigger->GetEntries();
  for (int i=0; i<1e5; ++i)
  {
    t_L1TrackTrigger->GetEvent(i);
    
    h_STUB_n->Fill(STUB_n);
    
    for (unsigned int j=0; j<STUB_x->size(); ++j)
    {
      if (
          STUB_layer->at(j)==5 ||
          STUB_layer->at(j)==6 ||
          STUB_layer->at(j)==7 ||
          STUB_layer->at(j)==8 ||
          STUB_layer->at(j)==9 ||
          STUB_layer->at(j)==10
         ) h_Barrel_STUB_x_y->Fill(STUB_x->at(j), STUB_y->at(j));
    
      double r=pow(STUB_x->at(j)*STUB_x->at(j)+STUB_y->at(j)*STUB_y->at(j), 0.5);
      h_STUB_r_z->Fill(STUB_z->at(j),r);
      // h_STUB_3D->Fill(STUB_x->at(j), STUB_y->at(j), STUB_z->at(j));
      
      if (STUB_x->size()>10) h_STUB_r_z_hiDeg->Fill(STUB_z->at(j),r);
      
    }
  }
  
  TCanvas *c_STUB_n=new TCanvas("c_STUB_n", "c_STUB_n", 700, 700);
  h_STUB_n->Draw();
  
  TCanvas *c_Barrel_STUB_x_y=new TCanvas("c_Barrel_STUB_x_y", "c_Barrel_STUB_x_y", 700, 700);
  h_Barrel_STUB_x_y->Draw();
  
  TCanvas *c_STUB_r_z=new TCanvas("c_STUB_r_z", "c_STUB_r_z", 1500, 700);
  h_STUB_r_z->Draw();
  h_STUB_r_z_hiDeg->SetMarkerColor(kRed); h_STUB_r_z_hiDeg->SetMarkerStyle(3);
  h_STUB_r_z_hiDeg->Draw("same");
  
  // TCanvas *c_STUB_3D=new TCanvas("c_STUB_3D", "c_STUB_3D", 700, 700);
  // h_STUB_3D->Draw();
  
  delete t_L1TrackTrigger;
}
  
  
