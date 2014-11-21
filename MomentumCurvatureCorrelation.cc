#include <TROOT.h>
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TF1.h"
#include "TAxis.h"

#include <iostream>
#include <vector>
#include <assert.h>

double pi=3.14159265358979;
double deltaPhi(double phi1, double phi2) { 
  double result = phi1 - phi2;
  while (result > pi) result -= 2*pi;
  while (result <= -pi) result += 2*pi;
  return result;
}

void MomentumCurvatureCorrelation()
{
  gROOT->ProcessLine(".L ../loader.C+");
  
  TFile *file=new TFile("tracks.root");
  TTree *tree=(TTree*)file->Get("ntupler/tree");
  
  std::vector<float> *tracks_rho=0;
  std::vector<float> *tracks_phi=0;
  std::vector<float> *genParts_pt=0;
  std::vector<float> *genParts_phi=0;
  
  tree->SetBranchAddress("AMTTTracks_rho", &(tracks_rho));
  tree->SetBranchAddress("AMTTTracks_phi0", &(tracks_phi));
  tree->SetBranchAddress("genParts_pt", &(genParts_pt));
  tree->SetBranchAddress("genParts_phi", &(genParts_phi));
  
  std::vector<float> track_rho;
  std::vector<float> track_pT;
  
  unsigned int nEvents=tree->GetEntries();
  std::cout<<"Using nEvents = "<<nEvents<<std::endl;
  for (unsigned int i=0; i<nEvents; ++i)
  {
    tree->GetEntry(i);
    
    unsigned int ntracks=tracks_rho->size();
    unsigned int nparts=genParts_pt->size();
    if (nparts!=1) {std::cout<<"Event "<<i<<" had no gen particles!"<<std::endl; continue;}
    
    double genphi=genParts_phi->at(0);
    double minDeltaPhi=0.001;
    int i_track=-1;
    for (unsigned int j=0; j<ntracks; ++j)
    {
      double trackphi=tracks_phi->at(j);
      double trackDeltaPhi=fabs(deltaPhi(trackphi,genphi));
      if (trackDeltaPhi<minDeltaPhi)
      {
        minDeltaPhi=trackDeltaPhi;
        i_track=j;
      }
    }
    if (i_track!=-1)
    {
      track_rho.push_back(fabs(tracks_rho->at(i_track)));
      track_pT.push_back(genParts_pt->at(0));
    }
    else
    {
      // std::cout<<"No reliable L1 track for event "<<i<<std::endl;
    }
  }
  
  TGraph *g_pT_rho=new TGraph(track_rho.size(), &(track_rho.at(0)), &(track_pT.at(0)));
  g_pT_rho->SetTitle("Track p_{T} vs #rho; #rho (cm); p_{T} (GeV)"); g_pT_rho->GetYaxis()->SetTitleOffset(1.3);
  
  TF1 *f_line=new TF1("f_line", "pol2");
  g_pT_rho->Fit(f_line);
  
  TCanvas *c_pT_rho=new TCanvas("c_pT_rho", "c_pT_rho", 700, 700);
  g_pT_rho->Draw("A*");
  c_pT_rho->SetLogy(); c_pT_rho->SetLogx();
  c_pT_rho->SaveAs("c_pT_rho.png");
}
      
      
      
      
      
