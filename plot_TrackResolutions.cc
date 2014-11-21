#include <TROOT.h>
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TEllipse.h"
#include "TGAxis.h"
#include "TF1.h"
#include "TPaveText.h"
#include "TLine.h"

#include <iostream>
#include <vector>
#include <assert.h>
#include <math.h>

std::string sign;

std::string itoa(int i) 
{
  char res[10];
  sprintf(res, "%d", i);
  std::string ret(res);
  return ret;
}

std::string ftoa(double i) 
{
  char res[10];
  sprintf(res, "%.2f", i);
  std::string ret(res);
  return ret;
}

double pi=3.14159265358979;

double deltaPhi(double phi1, double phi2) 
{ 
  double result = phi1 - phi2;
  while (result > pi) result -= 2*pi;
  while (result <= -pi) result += 2*pi;
  return result;
}

double convert_rhoTopT(double rho)
{
  double pT=-99;
  if (sign=="Minus")     pT = (0.025778) + (0.011428)*rho + (1.96861e-08)*(rho*rho);
  else if (sign=="Plus") pT = (0.0242589) + (0.0114315)*rho + (-2.17911e-08)*(rho*rho);
  return pT;
}

double convert_cotthetaToeta(double cottheta)
{
  double theta=fabs(atan(1./cottheta));
  double eta=-log(tan(theta/2.))*cottheta/fabs(cottheta);
  return eta;
}

Double_t yzSinusoid(Double_t *z, Double_t *par)
{
  // par[0] = rho, par[1] = phi0, par[2] = cottheta, par[3] = z0
  Double_t r=par[0]*sin((z[0]-par[3])/(par[0]*par[2]));
  Double_t phi=par[1]-asin(r/(2.*par[0]));
  Double_t y=r*sin(phi);
  return y;
}

Double_t xzSinusoid(Double_t *z, Double_t *par)
{
  // par[0] = rho, par[1] = phi0, par[2] = cottheta, par[3] = z0
  Double_t r=par[0]*sin((z[0]-par[3])/(par[0]*par[2]));
  Double_t phi=par[1]-asin(r/(2.*par[0]));
  Double_t x=r*cos(phi);
  return x;
}

double bestTrackPhi(double phi0, double rho, double ri)
{
  double phi=phi0-asin(ri/(2.*rho));
  return phi;
}

double bestTrackz(double rho, double cottheta, double z0, double ri)
{
  double z=rho*cottheta*asin(ri/rho)+z0;
  return z;
}

void plot_TrackResolutions(std::string signString)
{
  gROOT->ProcessLine(".L loader.C+");
  sign=signString;
  
  TFile *file=new TFile(("tracks_Mu"+sign+".root").c_str());
  TTree *tree=(TTree*)file->Get("ntupler/tree");
  
  std::vector<float>               *tracks_rho=0;
  std::vector<float>               *tracks_phi0=0;
  std::vector<float>               *tracks_cottheta0=0;
  std::vector<float>               *tracks_z0=0;
  std::vector<float>               *tracks_chi2=0;
  std::vector<float>               *tracks_chi2Red=0;
  std::vector<unsigned int>        *tracks_roadIndex=0;
  std::vector<float>               *genParts_pt=0;
  std::vector<float>               *genParts_phi=0;
  std::vector<float>               *genParts_eta=0;
  std::vector<float>               *genParts_z0=0;
  // std::vector<std::vector<float> > *hits_x = 0;
  // std::vector<std::vector<float> > *hits_y = 0;
  std::vector<std::vector<float> > *hits_z = 0;
  std::vector<std::vector<float> > *hits_r = 0;
  std::vector<std::vector<float> > *hits_phi = 0;
  std::vector<std::vector<float> > *hits_stubWidths = 0;
  
  tree->SetBranchAddress("AMTTTracks_rho", &(tracks_rho));
  tree->SetBranchAddress("AMTTTracks_phi0", &(tracks_phi0));
  tree->SetBranchAddress("AMTTTracks_cottheta0", &(tracks_cottheta0));
  tree->SetBranchAddress("AMTTTracks_z0", &(tracks_z0));
  tree->SetBranchAddress("AMTTTracks_chi2", &(tracks_chi2));
  tree->SetBranchAddress("AMTTTracks_chi2Red", &(tracks_chi2Red));
  tree->SetBranchAddress("AMTTTracks_roadIndex", &(tracks_roadIndex));
  tree->SetBranchAddress("genParts_pt", &(genParts_pt));
  tree->SetBranchAddress("genParts_phi", &(genParts_phi));
  tree->SetBranchAddress("genParts_eta", &(genParts_eta));
  tree->SetBranchAddress("genParts_vz", &(genParts_z0));
  // tree->SetBranchAddress("AMTTRoads_hitXs", &(hits_x));
  // tree->SetBranchAddress("AMTTRoads_hitYs", &(hits_y));
  tree->SetBranchAddress("AMTTRoads_hitZs", &(hits_z));
  tree->SetBranchAddress("AMTTRoads_hitRs", &(hits_r));
  tree->SetBranchAddress("AMTTRoads_hitPhis", &(hits_phi));
  tree->SetBranchAddress("AMTTRoads_hitStubWidths", &(hits_stubWidths));
  
  // chi^2 plots
  TH1F *h_chi2=new TH1F("h_chi2", "Fit #chi^{2}; #chi^{2}", 100, 0, 20);
  TH1F *h_chi2Red=new TH1F("h_chi2Red", "Fit #chi^{2}/nDOF; #chi^{2}/nDOF", 100, 0, 10);
  
  // Residual plots
  TH1F *h_residual_phi_1=new TH1F("h_residual_phi_1", "#phi Measurement Residual Layer 1; #phi_{1}^{m} - #phi_{1}", 100, -0.005, 0.005);
  TH1F *h_residual_phi_2=new TH1F("h_residual_phi_2", "#phi Measurement Residual Layer 2; #phi_{2}^{m} - #phi_{2}", 100, -0.005, 0.005);
  TH1F *h_residual_phi_3=new TH1F("h_residual_phi_3", "#phi Measurement Residual Layer 3; #phi_{3}^{m} - #phi_{3}", 100, -0.005, 0.005);
  TH1F *h_residual_phi_4=new TH1F("h_residual_phi_4", "#phi Measurement Residual Layer 4; #phi_{4}^{m} - #phi_{4}", 100, -0.005, 0.005);
  TH1F *h_residual_phi_5=new TH1F("h_residual_phi_5", "#phi Measurement Residual Layer 5; #phi_{5}^{m} - #phi_{5}", 100, -0.005, 0.005);
  TH1F *h_residual_phi_6=new TH1F("h_residual_phi_6", "#phi Measurement Residual Layer 6; #phi_{6}^{m} - #phi_{6}", 100, -0.01, 0.01);
  TH1F *h_residual_z_1=new TH1F("h_residual_z_1", "z Measurement Residual Layer 1; z_{1}^{m} - z_{1} (cm)", 100, -1, 1);
  TH1F *h_residual_z_2=new TH1F("h_residual_z_2", "z Measurement Residual Layer 2; z_{2}^{m} - z_{2} (cm)", 100, -1, 1);
  TH1F *h_residual_z_3=new TH1F("h_residual_z_3", "z Measurement Residual Layer 3; z_{3}^{m} - z_{3} (cm)", 100, -1, 1);
  TH1F *h_residual_z_4=new TH1F("h_residual_z_4", "z Measurement Residual Layer 4; z_{4}^{m} - z_{4} (cm)", 100, -5, 5);
  TH1F *h_residual_z_5=new TH1F("h_residual_z_5", "z Measurement Residual Layer 5; z_{5}^{m} - z_{5} (cm)", 100, -7, 7);
  TH1F *h_residual_z_6=new TH1F("h_residual_z_6", "z Measurement Residual Layer 6; z_{6}^{m} - z_{6} (cm)", 100, -7, 7);
  
  // Stub widths in various layers
  TH1F *h_stubWidth_1=new TH1F("h_stubWidth_1", "Stub Width Layer 1; #Deltas", 100, -10, 10);
  TH1F *h_stubWidth_2=new TH1F("h_stubWidth_2", "Stub Width Layer 2; #Deltas", 100, -10, 10);
  TH1F *h_stubWidth_3=new TH1F("h_stubWidth_3", "Stub Width Layer 3; #Deltas", 100, -10, 10);
  TH1F *h_stubWidth_4=new TH1F("h_stubWidth_4", "Stub Width Layer 4; #Deltas", 100, -10, 10);
  TH1F *h_stubWidth_5=new TH1F("h_stubWidth_5", "Stub Width Layer 5; #Deltas", 100, -10, 10);
  TH1F *h_stubWidth_6=new TH1F("h_stubWidth_6", "Stub Width Layer 6; #Deltas", 100, -10, 10);
  
  // Resolution plots
  TH1F *h_res_pT=new TH1F("h_res_pT", "; p_{T}^{track} - p_{T}^{gen} (GeV)", 400, -0.5, 0.5); h_res_pT->Sumw2();
  TH2F *h_res_gen_pT=new TH2F("h_res_gen_pT", "; gen p_{T} (GeV); p_{T}^{track} - p_{T}^{gen} (GeV)", 200, 0, 200, 100, -15, 15); h_res_gen_pT->Sumw2();
  TH2F *h_res_pT_chi2=new TH2F("h_res_pT_chi2", "; #chi^{2}; p_{T}^{track} - p_{T}^{gen} (GeV)", 200, 0, 20, 100, -10, 10); h_res_pT_chi2->Sumw2();
  
  TH1F *h_resFrac_pT=new TH1F("h_resFrac_pT", "; (p_{T}^{track} - p_{T}^{gen})/p_{T}", 400, -0.1, 0.1); h_resFrac_pT->Sumw2();
  TH2F *h_resFrac_gen_pT=new TH2F("h_resFrac_gen_pT", "; gen p_{T} (GeV); (p_{T}^{track} - p_{T}^{gen})/p_{T}^{gen} (GeV)", 200, 0, 200, 400, -0.1, 0.1); h_resFrac_gen_pT->Sumw2();
  
  TH1F *h_res_phi=new TH1F("h_res_phi", "; #phi_{0}^{track} - #phi_{0}^{gen}", 400, -0.01, 0.01); h_res_phi->Sumw2();
  TH2F *h_res_phi_pT=new TH2F("h_res_phi_pT", "; gen p_{T} (GeV); #phi^{track} - #phi^{gen} (GeV)", 200, 0, 200, 400, -0.005, 0.005); h_res_phi_pT->Sumw2();
  
  TH1F *h_res_eta=new TH1F("h_res_eta", "; #eta^{track} - #eta^{gen}", 400, -0.02, 0.02); h_res_eta->Sumw2();
  TH2F *h_res_eta_pT=new TH2F("h_res_eta_pT", "; gen p_{T} (GeV); #eta^{track} - #eta^{gen}", 200, 0, 200, 400, -0.02, 0.02);
  
  TH1F *h_res_z=new TH1F("h_res_z", "; z_{0}^{track} - z_{0}^{gen} (cm)", 400, -0.5, 0.5); h_res_z->Sumw2();
  TH2F *h_res_z_pT=new TH2F("h_res_z_pT", "; gen p_{T} (GeV); z_{0}^{track} - z_{0}^{gen} (cm)", 200, 0, 200, 400, -0.5, 0.5);
  
  TGraph *g_y_x=0, *g_y_z=0, *g_x_z=0;
  
  unsigned int nEvents=tree->GetEntries();
  std::cout<<"nEvents = "<<nEvents<<std::endl;
  unsigned int nplots=0;
  for (unsigned int i=0; i<nEvents; ++i)
  {
    tree->GetEntry(i);
    
    unsigned int ntracks=tracks_rho->size();
    unsigned int nparts=genParts_pt->size();
    // assert(nparts==1); // There has to be exactly one gen track for this
    if (nparts!=1) {std::cout<<"Event "<<i<<" had no gen particles!"<<std::endl; continue;}
    
    // Take the gen track, find the closest L1 track in \phi and use it to measure track resolutions
    double genphi=genParts_phi->at(0);
    double genpT=genParts_pt->at(0);
    double geneta=genParts_eta->at(0);
    double genz0=genParts_z0->at(0);
    double minDeltaPhi=999;
    
    int i_track=-1;
    for (unsigned int j=0; j<ntracks; ++j)
    {
      double trackphi=tracks_phi0->at(j);
      double trackDeltaPhi=fabs(deltaPhi(trackphi, genphi));
      if (trackDeltaPhi<minDeltaPhi)
      {
        minDeltaPhi=trackDeltaPhi;
        i_track=j;
      }
    }
    if (i_track!=-1)
    {
    
      double track_rho=tracks_rho->at(i_track);
      double track_phi=tracks_phi0->at(i_track);
      double track_cottheta=tracks_cottheta0->at(i_track);
      double track_z0=tracks_z0->at(i_track);
      double track_chi2=tracks_chi2->at(i_track);
      double track_chi2Red=tracks_chi2Red->at(i_track);
      
      double trackpT=convert_rhoTopT(fabs(tracks_rho->at(i_track)));
      double trackphi=tracks_phi0->at(i_track);
      double tracketa=convert_cotthetaToeta(tracks_cottheta0->at(i_track));
      
      unsigned int roadIndex=tracks_roadIndex->at(i_track);
      
      /*std::cout<<"roadIndex = "<<roadIndex<<std::endl;
      std::cout<<"hits_r->at(roadIndex).at(0) = "<<hits_r->at(roadIndex).at(0)<<std::endl;
      std::cout<<"hits_r->at(roadIndex).at(1) = "<<hits_r->at(roadIndex).at(1)<<std::endl;
      std::cout<<"hits_r->at(roadIndex).at(2) = "<<hits_r->at(roadIndex).at(2)<<std::endl;
      std::cout<<"hits_r->at(roadIndex).at(3) = "<<hits_r->at(roadIndex).at(3)<<std::endl;
      std::cout<<"hits_r->at(roadIndex).at(4) = "<<hits_r->at(roadIndex).at(4)<<std::endl;
      std::cout<<"hits_r->at(roadIndex).at(5) = "<<hits_r->at(roadIndex).at(5)<<std::endl;*/
      
      // chi^2 plots
      h_chi2->Fill(track_chi2);
      h_chi2Red->Fill(track_chi2Red);
      
      // Residual plots
      if (hits_phi->at(roadIndex).size()==6) // && genpT>100)
      {
        double residual_phi_1 = hits_phi->at(roadIndex).at(0) - bestTrackPhi(track_phi, track_rho, hits_r->at(roadIndex).at(0));
        double residual_phi_2 = hits_phi->at(roadIndex).at(1) - bestTrackPhi(track_phi, track_rho, hits_r->at(roadIndex).at(1));
        double residual_phi_3 = hits_phi->at(roadIndex).at(2) - bestTrackPhi(track_phi, track_rho, hits_r->at(roadIndex).at(2));
        double residual_phi_4 = hits_phi->at(roadIndex).at(3) - bestTrackPhi(track_phi, track_rho, hits_r->at(roadIndex).at(3));
        double residual_phi_5 = hits_phi->at(roadIndex).at(4) - bestTrackPhi(track_phi, track_rho, hits_r->at(roadIndex).at(4));
        double residual_phi_6 = hits_phi->at(roadIndex).at(5) - bestTrackPhi(track_phi, track_rho, hits_r->at(roadIndex).at(5));
        double residual_z_1 = hits_z->at(roadIndex).at(0) - bestTrackz(track_rho, track_cottheta, track_z0, hits_r->at(roadIndex).at(0));
        double residual_z_2 = hits_z->at(roadIndex).at(1) - bestTrackz(track_rho, track_cottheta, track_z0, hits_r->at(roadIndex).at(1));
        double residual_z_3 = hits_z->at(roadIndex).at(2) - bestTrackz(track_rho, track_cottheta, track_z0, hits_r->at(roadIndex).at(2));
        double residual_z_4 = hits_z->at(roadIndex).at(3) - bestTrackz(track_rho, track_cottheta, track_z0, hits_r->at(roadIndex).at(3));
        double residual_z_5 = hits_z->at(roadIndex).at(4) - bestTrackz(track_rho, track_cottheta, track_z0, hits_r->at(roadIndex).at(4));
        double residual_z_6 = hits_z->at(roadIndex).at(5) - bestTrackz(track_rho, track_cottheta, track_z0, hits_r->at(roadIndex).at(5));
      
        h_residual_phi_1->Fill(residual_phi_1);
        h_residual_phi_2->Fill(residual_phi_2);
        h_residual_phi_3->Fill(residual_phi_3);
        h_residual_phi_4->Fill(residual_phi_4);
        h_residual_phi_5->Fill(residual_phi_5);
        h_residual_phi_6->Fill(residual_phi_6);
        h_residual_z_1->Fill(residual_z_1);
        h_residual_z_2->Fill(residual_z_2);
        h_residual_z_3->Fill(residual_z_3);
        h_residual_z_4->Fill(residual_z_4);
        h_residual_z_5->Fill(residual_z_5);
        h_residual_z_6->Fill(residual_z_6);
        
        h_stubWidth_1->Fill(hits_stubWidths->at(roadIndex).at(0));
        h_stubWidth_2->Fill(hits_stubWidths->at(roadIndex).at(1));
        h_stubWidth_3->Fill(hits_stubWidths->at(roadIndex).at(2));
        h_stubWidth_4->Fill(hits_stubWidths->at(roadIndex).at(3));
        h_stubWidth_5->Fill(hits_stubWidths->at(roadIndex).at(4));
        h_stubWidth_6->Fill(hits_stubWidths->at(roadIndex).at(5));
        // std::cout<<"stub width = "<<hits_stubWidths->at(roadIndex).at(5)<<std::endl;
      }
      
      // Resolution plots
      h_res_pT->Fill(trackpT-genpT);
      h_resFrac_pT->Fill((trackpT-genpT)/genpT);
      h_res_gen_pT->Fill(genpT, trackpT-genpT);
      h_res_pT_chi2->Fill(track_chi2, trackpT-genpT);
      h_resFrac_gen_pT->Fill(genpT, (trackpT-genpT)/genpT);
      
      h_res_phi->Fill(deltaPhi(trackphi, genphi));
      h_res_phi_pT->Fill(genpT, deltaPhi(trackphi, genphi));
      
      h_res_z->Fill(track_z0-genz0);
      h_res_z_pT->Fill(genpT, track_z0-genz0);
      
      h_res_eta->Fill(tracketa-geneta);
      h_res_eta_pT->Fill(genpT, tracketa-geneta);
      
      // Fill the stub (x, y) corresponding to this track
      if (genpT>2.0 && genpT<2.1 && nplots<10)
      {
        unsigned int nStubs=hits_z->at(roadIndex).size();
        std::vector<float> hits_x, hits_y;
        for (unsigned int j=0; j<nStubs; ++j)
        {
          hits_x.push_back(hits_r->at(roadIndex).at(j)*cos(hits_phi->at(roadIndex).at(j)));
          hits_y.push_back(hits_r->at(roadIndex).at(j)*sin(hits_phi->at(roadIndex).at(j)));
        }
        g_y_x=new TGraph(nStubs, &(hits_x.at(0)), &(hits_y.at(0)));
        g_y_z=new TGraph(nStubs, &(hits_z->at(roadIndex).at(0)), &(hits_y.at(0)));
        g_x_z=new TGraph(nStubs, &(hits_z->at(roadIndex).at(0)), &(hits_x.at(0)));
        
        // Draw it right away
        TCanvas *c_TrackFit=new TCanvas("c_TrackFit", "c_TrackFit", 1400, 700);
        c_TrackFit->Divide(2,1);
        c_TrackFit->cd(1);
        g_y_x->SetTitle("Transverse View of Track Fit; x (cm); y (cm)");
        g_y_x->GetXaxis()->SetLimits(-120, 120); 
        g_y_x->GetYaxis()->SetRangeUser(-120, 120); g_y_x->GetYaxis()->SetTitleOffset(1.4);
        g_y_x->Draw("A*");
        TPad *p_y_x=(TPad*)c_TrackFit->GetPad(1);
        p_y_x->Range(-120, -120, 120, 120);
        double track_center_x=(track_rho)*sin(track_phi);
        double track_center_y=-(track_rho)*cos(track_phi);
        TEllipse *layer1=new TEllipse(track_center_x, track_center_y, fabs(track_rho)); 
        layer1->SetFillStyle(4000); layer1->SetLineColor(kRed);
        layer1->Draw("same");
        TGaxis *new_xaxis=new TGaxis(-120, 0, 120, 0, -120, 120, 10, ""); new_xaxis->Draw();
        TGaxis *new_yaxis=new TGaxis(0, -120, 0, 120, -120, 120, 10, ""); new_yaxis->Draw();
        c_TrackFit->cd();
        // TPaveText *pave_trackFit=new TPaveText(20, -20, 100, -100);
        TPaveText *pave_event=new TPaveText(0.40, 0.92, 0.50, 0.98);
        pave_event->AddText(("Event # "+itoa(i)).c_str());
        pave_event->Draw();
        TPaveText *pave_trackFit=new TPaveText(0.45, 0.51, 0.53, 0.9);
        pave_trackFit->AddText("Track Parameters");
        pave_trackFit->SetTextAlign(12);
        pave_trackFit->AddText(("|#rho| = "+ftoa(fabs(track_rho))+" cm").c_str());
        pave_trackFit->AddText(("#phi_{0} = "+ftoa(track_phi)).c_str());
        pave_trackFit->AddText(("cot(#theta_{0}) = "+ftoa(track_cottheta)).c_str());
        pave_trackFit->AddText(("z_{0} = "+ftoa(track_z0)+" cm").c_str());
        pave_trackFit->AddText(("#chi^{2}/nDOF = "+ftoa(track_chi2Red)).c_str());
        pave_trackFit->Draw();
        TPaveText *pave_trackKin=new TPaveText(0.45, 0.1, 0.53, 0.49);
        pave_trackKin->AddText("Track Kinematics");
        pave_trackKin->SetTextAlign(12);
        pave_trackKin->AddText(("p_{T} = "+ftoa(convert_rhoTopT(fabs(track_rho)))+" GeV").c_str());
        pave_trackKin->AddText(("gen p_{T} = "+ftoa(genpT)+" GeV").c_str());
        pave_trackKin->AddText(("#eta = "+ftoa(convert_cotthetaToeta(track_cottheta))).c_str());
        pave_trackKin->AddText(("gen #eta = "+ftoa(geneta)).c_str());
        pave_trackKin->AddText(("#phi = "+ftoa(track_phi)).c_str());
        pave_trackKin->AddText(("gen #phi = "+ftoa(genphi)).c_str());
        pave_trackKin->Draw();
        TPad *p_Long=(TPad*)c_TrackFit->GetPad(2);
        p_Long->Divide(1, 2);
        TPad *p_y_z=(TPad*)p_Long->GetPad(1);
        p_y_z->cd();
        g_y_z->SetTitle("Longitudinal View of Track Fit; z (cm); y (cm)");
        g_y_z->GetXaxis()->SetLimits(-120, 120);
        g_y_z->GetYaxis()->SetRangeUser(-120, 120); g_y_z->GetYaxis()->SetTitleOffset(0.8);
        g_y_z->Draw("A*");
        TF1 *f_y_z=new TF1("f_y_z", yzSinusoid, -2.*fabs(track_rho*track_cottheta), 2.*fabs(track_rho*track_cottheta), 4);
        f_y_z->SetParameters(track_rho, track_phi, track_cottheta, track_z0);
        f_y_z->SetLineWidth(0);
        f_y_z->Draw("same");
        p_y_z->Range(-120, -120, 120, 120);
        TGaxis *new_zaxis=new TGaxis(-120, 0, 120, 0, -120, 120, 10, ""); new_zaxis->Draw();
        new_yaxis->Draw();
        TPad *p_x_z=(TPad*)p_Long->GetPad(2);
        p_x_z->cd();
        g_x_z->SetTitle("Longitudinal View of Track Fit; z (cm); x (cm)");
        g_x_z->GetXaxis()->SetLimits(-120, 120);
        g_x_z->GetYaxis()->SetRangeUser(-120, 120); g_x_z->GetYaxis()->SetTitleOffset(0.8);
        g_x_z->Draw("A*");
        TF1 *f_x_z=new TF1("f_x_z", xzSinusoid, -2.*fabs(track_rho*track_cottheta), 2.*fabs(track_rho*track_cottheta), 4);
        f_x_z->SetParameters(track_rho, track_phi, track_cottheta, track_z0);
        f_x_z->SetLineWidth(0);
        f_x_z->Draw("same");
        p_x_z->Range(-120, -120, 120, 120);
        new_xaxis->Draw();
        new_yaxis->Draw();
        
        c_TrackFit->Update();
        c_TrackFit->SaveAs(("Visuals_Mu"+sign+"/c_TrackFit_"+itoa(i)+".png").c_str());
        
        delete c_TrackFit;
        nplots+=1;
      }
    }
  }
  
  // Make some nice plots
  /*
  gROOT->SetStyle("Plain");
  // gStyle->SetOptStat(000000000);
  
  TCanvas *c_res_pT=new TCanvas("c_res_pT", "c_res_pT", 700, 700);
  h_res_pT->GetYaxis()->SetTitleOffset(1.2);
  h_res_pT->GetXaxis()->SetRangeUser(-2, 3);
  h_res_pT->Draw("hist");
  TLine *line=new TLine(0, 0, 0, h_res_pT->GetMaximum()); line->Draw();
  c_res_pT->SaveAs("c_res_pT.png");
  
  TCanvas *c_resFrac_pT=new TCanvas("c_resFrac_pT", "c_resFrac_pT", 700, 700);
  h_resFrac_pT->GetYaxis()->SetTitleOffset(1.2);
  h_resFrac_pT->Draw("hist");
  line=new TLine(0, 0, 0, h_resFrac_pT->GetMaximum()); line->Draw();
  c_resFrac_pT->SaveAs("c_resFrac_pT.png");
  
  TCanvas *c_res_phi=new TCanvas("c_res_phi", "c_res_phi", 700, 700);
  h_res_phi->GetYaxis()->SetTitleOffset(1.2);
  h_res_phi->Draw("hist");
  line=new TLine(0, 0, 0, h_res_phi->GetMaximum()); line->Draw();
  c_res_phi->SaveAs("c_res_phi.png");
  
  TCanvas *c_res_eta=new TCanvas("c_res_eta", "c_res_eta", 700, 700);
  h_res_eta->GetYaxis()->SetTitleOffset(1.2);
  h_res_eta->Draw("hist");
  line=new TLine(0, 0, 0, h_res_eta->GetMaximum()); line->Draw();
  c_res_eta->SaveAs("c_res_eta.png");
  
  TCanvas *c_res_z=new TCanvas("c_res_z", "c_res_z", 700, 700);
  h_res_z->GetYaxis()->SetTitleOffset(1.2);
  h_res_z->Draw("hist");
  line=new TLine(0, 0, 0, h_res_z->GetMaximum()); line->Draw();
  c_res_z->SaveAs("c_res_z.png");
  
  TCanvas *c_res_pT_chi2=new TCanvas("c_res_pT_chi2", "c_res_pT_chi2", 700, 700);
  h_res_pT_chi2->GetYaxis()->SetTitleOffset(1.2);
  h_res_pT_chi2->Draw("colz");
  h_res_pT_chi2->ProfileX()->Draw("same e1");
  line=new TLine(0, 0, 200, 0); line->Draw();
  c_res_pT_chi2->SaveAs("c_res_pT_chi2.png");
  
  TCanvas *c_res_gen_pT=new TCanvas("c_res_gen_pT", "c_res_gen_pT", 700, 700);
  h_res_gen_pT->GetYaxis()->SetTitleOffset(1.2);
  h_res_gen_pT->Draw("colz");
  h_res_gen_pT->ProfileX()->Draw("same e1");
  line=new TLine(0, 0, 200, 0); line->Draw();
  c_res_gen_pT->SaveAs("c_res_gen_pT.png");
  
  TCanvas *c_res_gen_phi=new TCanvas("c_res_gen_phi", "c_res_gen_phi", 700, 700);
  h_res_phi_pT->GetYaxis()->SetTitleOffset(1.2);
  h_res_phi_pT->Draw("colz");
  h_res_phi_pT->ProfileX()->Draw("same e1");
  line=new TLine(0, 0, 200, 0); line->Draw();
  c_res_gen_phi->SaveAs("c_res_gen_phi.png");
  
  TCanvas *c_chi2=new TCanvas("c_chi2", "c_chi2", 700, 700);
  h_chi2->Draw();
  c_chi2->SetLogy();
  c_chi2->SaveAs("c_chi2.png");
  
  TCanvas *c_chi2Red=new TCanvas("c_chi2Red", "c_chi2Red", 700, 700);
  h_chi2Red->Draw();
  c_chi2Red->SetLogy();
  c_chi2Red->SaveAs("c_chi2Red.png");
  
  // Hit residuals
  TCanvas *c_residual_phi_1=new TCanvas("c_residual_phi_1", "c_residual_phi_1", 700, 700);
  h_residual_phi_1->Draw();
  c_residual_phi_1->SaveAs("c_residual_phi_1.png");
  
  TCanvas *c_residual_phi_2=new TCanvas("c_residual_phi_2", "c_residual_phi_2", 700, 700);
  h_residual_phi_2->Draw();
  c_residual_phi_2->SaveAs("c_residual_phi_2.png");
  
  TCanvas *c_residual_phi_3=new TCanvas("c_residual_phi_3", "c_residual_phi_3", 700, 700);
  h_residual_phi_3->Draw();
  c_residual_phi_3->SaveAs("c_residual_phi_3.png");
  
  TCanvas *c_residual_phi_4=new TCanvas("c_residual_phi_4", "c_residual_phi_4", 700, 700);
  h_residual_phi_4->Draw();
  c_residual_phi_4->SaveAs("c_residual_phi_4.png");
  
  TCanvas *c_residual_phi_5=new TCanvas("c_residual_phi_5", "c_residual_phi_5", 700, 700);
  h_residual_phi_5->Draw();
  c_residual_phi_5->SaveAs("c_residual_phi_5.png");
  
  TCanvas *c_residual_phi_6=new TCanvas("c_residual_phi_6", "c_residual_phi_6", 700, 700);
  h_residual_phi_6->Draw();
  c_residual_phi_6->SaveAs("c_residual_phi_6.png");
  
  TCanvas *c_residual_z_1=new TCanvas("c_residual_z_1", "c_residual_z_1", 700, 700);
  h_residual_z_1->Draw();
  c_residual_z_1->SaveAs("c_residual_z_1.png");
  
  TCanvas *c_residual_z_2=new TCanvas("c_residual_z_2", "c_residual_z_2", 700, 700);
  h_residual_z_2->Draw();
  c_residual_z_2->SaveAs("c_residual_z_2.png");
  
  TCanvas *c_residual_z_3=new TCanvas("c_residual_z_3", "c_residual_z_3", 700, 700);
  h_residual_z_3->Draw();
  c_residual_z_3->SaveAs("c_residual_z_3.png");
  
  TCanvas *c_residual_z_4=new TCanvas("c_residual_z_4", "c_residual_z_4", 700, 700);
  h_residual_z_4->Draw();
  c_residual_z_4->SaveAs("c_residual_z_4.png");
  
  TCanvas *c_residual_z_5=new TCanvas("c_residual_z_5", "c_residual_z_5", 700, 700);
  h_residual_z_5->Draw();
  c_residual_z_5->SaveAs("c_residual_z_5.png");
  
  TCanvas *c_residual_z_6=new TCanvas("c_residual_z_6", "c_residual_z_6", 700, 700);
  h_residual_z_6->Draw();
  c_residual_z_6->SaveAs("c_residual_z_6.png");
  
  // Stub widths
  TCanvas *c_stubWidth_1=new TCanvas("c_stubWidth_1", "c_stubWidth_1", 700, 700);
  h_stubWidth_1->Draw();
  c_stubWidth_1->SaveAs("c_stubWidth_1.png");
  */
  
  // Save the plots in histograms for a dedicated displayer
  TFile *outfile=new TFile(("TrackHistograms_Mu"+sign+".root").c_str(), "recreate");
  h_chi2->Write();
  h_chi2Red->Write();
  h_residual_phi_1->Write();
  h_residual_phi_2->Write();
  h_residual_phi_3->Write();
  h_residual_phi_4->Write();
  h_residual_phi_5->Write();
  h_residual_phi_6->Write();
  h_residual_z_1->Write();
  h_residual_z_2->Write();
  h_residual_z_3->Write();
  h_residual_z_4->Write();
  h_residual_z_5->Write();
  h_residual_z_6->Write();
  h_stubWidth_1->Write();
  h_stubWidth_2->Write();
  h_stubWidth_3->Write();
  h_stubWidth_4->Write();
  h_stubWidth_5->Write();
  h_stubWidth_6->Write();
  h_res_pT->Write();
  h_resFrac_pT->Write();
  h_res_phi->Write();
  h_res_eta->Write();
  h_res_eta_pT->Write();
  h_res_z->Write();
  h_res_z_pT->Write();
  h_res_pT_chi2->Write();
  h_resFrac_gen_pT->Write();
  h_res_gen_pT->Write();
  h_res_phi_pT->Write();
  outfile->Close();
}
  
