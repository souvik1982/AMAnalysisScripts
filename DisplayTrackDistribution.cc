#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLine.h"
#include "TProfile.h"
#include "TGAxis.h"
#include <iostream>

#include "TDRStyle.h"

TH1D* RMSY(TH2F *h)
{
  TH1D* h_RMS=h->ProjectionX();
  h_RMS->SetName(("res_"+std::string(h->GetName())).c_str());
  for (int i=1; i<h_RMS->GetNbinsX(); i+=1)
  {
    h_RMS->SetBinContent(i, h->ProjectionY("_px", i, i)->GetRMS());
    h_RMS->SetBinError(i, h->ProjectionY("_px", i, i)->GetRMSError());
  }
  return h_RMS;
}

TH1D* RMSY_for_pT(TH2F *h)
{
  TH1D* h_RMS=h->ProjectionX();
  h_RMS->SetName(("res_"+std::string(h->GetName())).c_str());
  for (int i=1; i<h_RMS->GetNbinsX(); i+=1)
  {
    h_RMS->SetBinContent(i, (h->ProjectionY("_px", i, i)->GetRMS())/h->GetBinCenter(i));
    h_RMS->SetBinError(i, h->ProjectionY("_px", i, i)->GetRMSError()/h->GetBinCenter(i));
  }
  return h_RMS;
}

void DisplayTrackDistribution()
{
  TFile *f_muMinus=new TFile("TrackHistograms_MuMinus.root");
  TFile *f_muPlus=new TFile("TrackHistograms_MuPlus.root");
  
  // Mu minus plots
  TH1F *h_residual_phi_1_muMinus=(TH1F*)f_muMinus->Get("h_residual_phi_1");
  TH1F *h_residual_phi_2_muMinus=(TH1F*)f_muMinus->Get("h_residual_phi_2");
  TH1F *h_residual_phi_3_muMinus=(TH1F*)f_muMinus->Get("h_residual_phi_3");
  TH1F *h_residual_phi_4_muMinus=(TH1F*)f_muMinus->Get("h_residual_phi_4");
  TH1F *h_residual_phi_5_muMinus=(TH1F*)f_muMinus->Get("h_residual_phi_5");
  TH1F *h_residual_phi_6_muMinus=(TH1F*)f_muMinus->Get("h_residual_phi_6");
  
  TH1F *h_residual_z_1_muMinus=(TH1F*)f_muMinus->Get("h_residual_z_1");
  TH1F *h_residual_z_2_muMinus=(TH1F*)f_muMinus->Get("h_residual_z_2");
  TH1F *h_residual_z_3_muMinus=(TH1F*)f_muMinus->Get("h_residual_z_3");
  TH1F *h_residual_z_4_muMinus=(TH1F*)f_muMinus->Get("h_residual_z_4");
  TH1F *h_residual_z_5_muMinus=(TH1F*)f_muMinus->Get("h_residual_z_5");
  TH1F *h_residual_z_6_muMinus=(TH1F*)f_muMinus->Get("h_residual_z_6");
  
  TH1F *h_res_pT_MuMinus=(TH1F*)f_muMinus->Get("h_res_pT");
  TH2F *h_res_gen_pT_MuMinus=(TH2F*)f_muMinus->Get("h_res_gen_pT");
  TH1F *h_resFrac_pT_MuMinus=(TH1F*)f_muMinus->Get("h_resFrac_pT");
  TH2F *h_resFrac_gen_pT_MuMinus=(TH2F*)f_muMinus->Get("h_resFrac_gen_pT");
  TH1F *h_res_phi_MuMinus=(TH1F*)f_muMinus->Get("h_res_phi");
  TH2F *h_res_phi_pT_MuMinus=(TH2F*)f_muMinus->Get("h_res_phi_pT");
  TH1F *h_res_z_MuMinus=(TH1F*)f_muMinus->Get("h_res_z");
  TH2F *h_res_z_pT_MuMinus=(TH2F*)f_muMinus->Get("h_res_z_pT");
  TH1F *h_res_eta_MuMinus=(TH1F*)f_muMinus->Get("h_res_eta");
  TH2F *h_res_eta_pT_MuMinus=(TH2F*)f_muMinus->Get("h_res_eta_pT");
  
  // Mu plus plots
  TH1F *h_residual_phi_1_muPlus =(TH1F*)f_muPlus->Get("h_residual_phi_1");
  TH1F *h_residual_phi_2_muPlus =(TH1F*)f_muPlus->Get("h_residual_phi_2");
  TH1F *h_residual_phi_3_muPlus =(TH1F*)f_muPlus->Get("h_residual_phi_3");
  TH1F *h_residual_phi_4_muPlus =(TH1F*)f_muPlus->Get("h_residual_phi_4");
  TH1F *h_residual_phi_5_muPlus =(TH1F*)f_muPlus->Get("h_residual_phi_5");
  TH1F *h_residual_phi_6_muPlus =(TH1F*)f_muPlus->Get("h_residual_phi_6");
  
  TH1F *h_residual_z_1_muPlus =(TH1F*)f_muPlus->Get("h_residual_z_1");
  TH1F *h_residual_z_2_muPlus =(TH1F*)f_muPlus->Get("h_residual_z_2");
  TH1F *h_residual_z_3_muPlus =(TH1F*)f_muPlus->Get("h_residual_z_3");
  TH1F *h_residual_z_4_muPlus =(TH1F*)f_muPlus->Get("h_residual_z_4");
  TH1F *h_residual_z_5_muPlus =(TH1F*)f_muPlus->Get("h_residual_z_5");
  TH1F *h_residual_z_6_muPlus =(TH1F*)f_muPlus->Get("h_residual_z_6");
  
  TH1F *h_res_pT_MuPlus=(TH1F*)f_muPlus->Get("h_res_pT");
  TH2F *h_res_gen_pT_MuPlus=(TH2F*)f_muPlus->Get("h_res_gen_pT");
  TH1F *h_resFrac_pT_MuPlus=(TH1F*)f_muPlus->Get("h_resFrac_pT");
  TH2F *h_resFrac_gen_pT_MuPlus=(TH2F*)f_muPlus->Get("h_resFrac_gen_pT");
  TH1F *h_res_phi_MuPlus=(TH1F*)f_muPlus->Get("h_res_phi");
  TH2F *h_res_phi_pT_MuPlus=(TH2F*)f_muPlus->Get("h_res_phi_pT");
  TH1F *h_res_z_MuPlus=(TH1F*)f_muPlus->Get("h_res_z");
  TH2F *h_res_z_pT_MuPlus=(TH2F*)f_muPlus->Get("h_res_z_pT");
  TH1F *h_res_eta_MuPlus=(TH1F*)f_muPlus->Get("h_res_eta");
  TH2F *h_res_eta_pT_MuPlus=(TH2F*)f_muPlus->Get("h_res_eta_pT");
  
  h_residual_phi_1_muMinus->SetLineColor(kBlue);
  h_residual_phi_2_muMinus->SetLineColor(kBlue);
  h_residual_phi_3_muMinus->SetLineColor(kBlue);
  h_residual_phi_4_muMinus->SetLineColor(kBlue);
  h_residual_phi_5_muMinus->SetLineColor(kBlue);
  h_residual_phi_6_muMinus->SetLineColor(kBlue);
  
  h_residual_z_1_muMinus->SetLineColor(kBlue);
  h_residual_z_2_muMinus->SetLineColor(kBlue);
  h_residual_z_3_muMinus->SetLineColor(kBlue);
  h_residual_z_4_muMinus->SetLineColor(kBlue);
  h_residual_z_5_muMinus->SetLineColor(kBlue);
  h_residual_z_6_muMinus->SetLineColor(kBlue);
  
  h_res_pT_MuMinus->SetLineColor(kBlue);
  h_resFrac_pT_MuMinus->SetLineColor(kBlue);
  h_res_phi_MuMinus->SetLineColor(kBlue);
  h_res_phi_pT_MuMinus->SetLineColor(kBlue);
  h_res_z_MuMinus->SetLineColor(kBlue);
  h_res_z_pT_MuMinus->SetLineColor(kBlue);
  h_res_eta_MuMinus->SetLineColor(kBlue);
  h_res_eta_pT_MuMinus->SetLineColor(kBlue);
  
  h_residual_phi_1_muPlus->SetLineColor(kRed);
  h_residual_phi_2_muPlus->SetLineColor(kRed);
  h_residual_phi_3_muPlus->SetLineColor(kRed);
  h_residual_phi_4_muPlus->SetLineColor(kRed);
  h_residual_phi_5_muPlus->SetLineColor(kRed);
  h_residual_phi_6_muPlus->SetLineColor(kRed);
  
  h_residual_z_1_muPlus->SetLineColor(kRed);
  h_residual_z_2_muPlus->SetLineColor(kRed);
  h_residual_z_3_muPlus->SetLineColor(kRed);
  h_residual_z_4_muPlus->SetLineColor(kRed);
  h_residual_z_5_muPlus->SetLineColor(kRed);
  h_residual_z_6_muPlus->SetLineColor(kRed);
  
  h_res_pT_MuPlus->SetLineColor(kRed);
  h_resFrac_pT_MuPlus->SetLineColor(kRed);
  h_res_phi_MuPlus->SetLineColor(kRed);
  h_res_phi_pT_MuPlus->SetLineColor(kRed);
  h_res_z_MuPlus->SetLineColor(kRed);
  h_res_z_pT_MuPlus->SetLineColor(kRed);
  h_res_eta_MuMinus->SetLineColor(kRed);
  h_res_eta_pT_MuMinus->SetLineColor(kRed);
  
  // Some plots may be added up
  TH2F *h_res_z_pT=(TH2F*)h_res_z_pT_MuMinus->Clone("h_res_z_pT"); h_res_z_pT->Add(h_res_z_pT_MuPlus);
  TH2F *h_res_phi_pT=(TH2F*)h_res_phi_pT_MuMinus->Clone("h_res_phi_pT"); h_res_phi_pT->Add(h_res_phi_pT_MuPlus);
  
  // gROOT->SetStyle("Plain");
  // gStyle->SetOptStat(000000000);
  TStyle *myStyle=setTDRStyle();
  myStyle->cd();
  
  TCanvas *c_residual_phi_1=new TCanvas("c_residual_phi_1", "c_residual_phi_1", 700, 700);
  h_residual_phi_1_muMinus->Draw();
  h_residual_phi_1_muPlus->Draw("same");
  c_residual_phi_1->SaveAs("c_residual_phi_1.png");
  
  TCanvas *c_residual_phi_2=new TCanvas("c_residual_phi_2", "c_residual_phi_2", 700, 700);
  h_residual_phi_2_muMinus->Draw();
  h_residual_phi_2_muPlus->Draw("same");
  c_residual_phi_2->SaveAs("c_residual_phi_2.png");
  
  TCanvas *c_residual_phi_3=new TCanvas("c_residual_phi_3", "c_residual_phi_3", 700, 700);
  h_residual_phi_3_muMinus->Draw();
  h_residual_phi_3_muPlus->Draw("same");
  c_residual_phi_3->SaveAs("c_residual_phi_3.png");
  
  TCanvas *c_residual_phi_4=new TCanvas("c_residual_phi_4", "c_residual_phi_4", 700, 700);
  h_residual_phi_4_muMinus->Draw();
  h_residual_phi_4_muPlus->Draw("same");
  c_residual_phi_4->SaveAs("c_residual_phi_4.png");
  
  TCanvas *c_residual_phi_5=new TCanvas("c_residual_phi_5", "c_residual_phi_5", 700, 700);
  h_residual_phi_5_muMinus->Draw();
  h_residual_phi_5_muPlus->Draw("same");
  c_residual_phi_5->SaveAs("c_residual_phi_5.png");
  
  TCanvas *c_residual_phi_6=new TCanvas("c_residual_phi_6", "c_residual_phi_6", 700, 700);
  h_residual_phi_6_muMinus->Draw();
  h_residual_phi_6_muPlus->Draw("same");
  c_residual_phi_6->SaveAs("c_residual_phi_6.png");
  
  TCanvas *c_residual_z_1=new TCanvas("c_residual_z_1", "c_residual_z_1", 700, 700);
  h_residual_z_1_muMinus->Draw();
  h_residual_z_1_muPlus->Draw("same");
  c_residual_z_1->SaveAs("c_residual_z_1.png");
  
  TCanvas *c_residual_z_2=new TCanvas("c_residual_z_2", "c_residual_z_2", 700, 700);
  h_residual_z_2_muMinus->Draw();
  h_residual_z_2_muPlus->Draw("same");
  c_residual_z_2->SaveAs("c_residual_z_2.png");
  
  TCanvas *c_residual_z_3=new TCanvas("c_residual_z_3", "c_residual_z_3", 700, 700);
  h_residual_z_3_muMinus->Draw();
  h_residual_z_3_muPlus->Draw("same");
  c_residual_z_3->SaveAs("c_residual_z_3.png");
  
  TCanvas *c_residual_z_4=new TCanvas("c_residual_z_4", "c_residual_z_4", 700, 700);
  h_residual_z_4_muMinus->Draw();
  h_residual_z_4_muPlus->Draw("same");
  c_residual_z_4->SaveAs("c_residual_z_4.png");
  
  TCanvas *c_residual_z_5=new TCanvas("c_residual_z_5", "c_residual_z_5", 700, 700);
  h_residual_z_5_muMinus->Draw();
  h_residual_z_5_muPlus->Draw("same");
  c_residual_z_5->SaveAs("c_residual_z_5.png");
  
  TCanvas *c_residual_z_6=new TCanvas("c_residual_z_6", "c_residual_z_6", 700, 700);
  h_residual_z_6_muMinus->Draw();
  h_residual_z_6_muPlus->Draw("same");
  c_residual_z_6->SaveAs("c_residual_z_6.png");
  
  TCanvas *c_res_pT=new TCanvas("h_res_pT", "h_res_pT", 700, 700);
  h_res_pT_MuMinus->Draw("");
  h_res_pT_MuPlus->Draw("same");
  TLine *line=new TLine(0, 0, 0, h_res_pT_MuMinus->GetMaximum()); line->Draw();
  c_res_pT->SaveAs("c_res_pT.png");
  
  TCanvas *c_resFrac_pT=new TCanvas("c_resFrac_pT", "c_resFrac_pT", 700, 700);
  h_resFrac_pT_MuMinus->Draw("");
  h_resFrac_pT_MuPlus->Draw("same");
  line=new TLine(0, 0, 0, h_resFrac_pT_MuMinus->GetMaximum()); line->Draw();
  c_resFrac_pT->SaveAs("c_resFrac_pT.png");
  
  TCanvas *c_resFrac_gen_pT=new TCanvas("c_resFrac_gen_pT", "c_resFrac_gen_pT", 700, 700);
  h_resFrac_gen_pT_MuMinus->RebinX(5);
  h_resFrac_gen_pT_MuPlus->RebinX(5);
  TH1D *h_rmsResFrac_gen_pT_MuMinus=RMSY(h_resFrac_gen_pT_MuMinus);
  TH1D *h_rmsResFrac_gen_pT_MuPlus=RMSY(h_resFrac_gen_pT_MuPlus);
  h_rmsResFrac_gen_pT_MuMinus->SetYTitle("#sigma((p_{T}^{track} - p_{T}^{gen})/p_{T}^{gen})");
  h_rmsResFrac_gen_pT_MuMinus->GetYaxis()->SetRangeUser(0, 0.05);
  ((TGaxis*)h_rmsResFrac_gen_pT_MuMinus->GetYaxis())->SetMaxDigits(3);
  h_rmsResFrac_gen_pT_MuMinus->SetLineColor(kBlue);
  h_rmsResFrac_gen_pT_MuPlus->SetLineColor(kRed);
  h_rmsResFrac_gen_pT_MuMinus->Draw("E0");
  h_rmsResFrac_gen_pT_MuPlus->Draw("E0 same");
  c_resFrac_gen_pT->SaveAs("c_resFrac_gen_pT.png");
  
  TCanvas *c_res_gen_pT=new TCanvas("c_res_gen_pT", "c_res_gen_pT", 700, 700);
  h_res_gen_pT_MuMinus->RebinX(5);
  h_res_gen_pT_MuPlus->RebinX(5);
  TH1D *h_rmsRes_gen_pT_MuMinus=RMSY_for_pT(h_res_gen_pT_MuMinus);
  TH1D *h_rmsRes_gen_pT_MuPlus=RMSY_for_pT(h_res_gen_pT_MuPlus);
  h_rmsRes_gen_pT_MuMinus->SetYTitle("#sigma(p_{T}^{track} - p_{T}^{gen})/p_{T}^{gen}");
  h_rmsRes_gen_pT_MuMinus->GetYaxis()->SetRangeUser(0, 0.05);
  // ((TGaxis*)h_rmsRes_gen_pT_MuMinus->GetYaxis())->SetMaxDigits(3);
  h_rmsRes_gen_pT_MuMinus->SetLineColor(kBlue);
  h_rmsRes_gen_pT_MuPlus->SetLineColor(kRed);
  h_rmsRes_gen_pT_MuMinus->Draw("E0");
  h_rmsRes_gen_pT_MuPlus->Draw("E0 same");
  c_res_gen_pT->SaveAs("c_res_gen_pT.png");
  
  TCanvas *c_res_phi=new TCanvas("c_res_phi", "c_res_phi", 700, 700);
  h_res_phi_MuMinus->Draw();
  h_res_phi_MuPlus->Draw("same");
  c_res_phi->SaveAs("c_res_phi.png");
  
  TCanvas *c_res_phi_pT=new TCanvas("c_res_phi_pT", "c_res_phi_pT", 700, 700);
  h_res_phi_pT_MuMinus->RebinX(5);
  h_res_phi_pT_MuPlus->RebinX(5);
  h_res_phi_pT->RebinX(5);
  TH1D *h_rmsRes_phi_pT_MuMinus=RMSY(h_res_phi_pT_MuMinus);
  TH1D *h_rmsRes_phi_pT_MuPlus=RMSY(h_res_phi_pT_MuPlus);
  TH1D *h_rmsRes_phi_pT=RMSY(h_res_phi_pT);
  h_rmsRes_phi_pT_MuMinus->SetYTitle("#sigma(#phi_{0}^{track} - #phi_{0}^{gen})");
  h_rmsRes_phi_pT_MuMinus->GetYaxis()->SetRangeUser(0, 0.002);
  ((TGaxis*)h_rmsRes_phi_pT_MuMinus->GetYaxis())->SetMaxDigits(4);
  h_rmsRes_phi_pT_MuMinus->SetLineColor(kBlue);
  h_rmsRes_phi_pT_MuPlus->SetLineColor(kRed);
  h_rmsRes_phi_pT->SetMarkerStyle(20);
  h_rmsRes_phi_pT_MuMinus->Draw("E0");
  h_rmsRes_phi_pT_MuPlus->Draw("E0 same");
  h_rmsRes_phi_pT->Draw("E0 same");
  c_res_phi_pT->SaveAs("c_res_phi_pT.png");
  
  TCanvas *c_res_z=new TCanvas("c_res_z", "c_res_z", 700, 700);
  h_res_z_MuMinus->Draw();
  h_res_z_MuPlus->Draw("same");
  c_res_z->SaveAs("c_res_z.png");
  
  TCanvas *c_res_z_pT=new TCanvas("c_res_z_pT", "c_res_z_pT", 700, 700);
  h_res_z_pT_MuMinus->RebinX(5);
  h_res_z_pT_MuPlus->RebinX(5);
  h_res_z_pT->RebinX(5);
  TH1D *h_rmsRes_z_pT_MuMinus=RMSY(h_res_z_pT_MuMinus);
  TH1D *h_rmsRes_z_pT_MuPlus=RMSY(h_res_z_pT_MuPlus);
  TH1D *h_rmsRes_z_pT=RMSY(h_res_z_pT);
  h_rmsRes_z_pT_MuMinus->SetYTitle("#sigma(z_{0}^{track} - z_{0}^{gen}) (cm)");
  h_rmsRes_z_pT_MuMinus->GetYaxis()->SetRangeUser(0, 0.2);
  h_rmsRes_z_pT_MuMinus->SetLineColor(kBlue);
  h_rmsRes_z_pT_MuPlus->SetLineColor(kRed);
  h_rmsRes_z_pT->SetMarkerStyle(20);
  h_rmsRes_z_pT_MuMinus->Draw("E0");
  h_rmsRes_z_pT_MuPlus->Draw("E0 same");
  // h_rmsRes_z_pT->Draw("E0 same");
  c_res_z_pT->SaveAs("c_res_z_pT.png");
  
  TCanvas *c_res_eta=new TCanvas("c_res_eta", "c_res_eta", 700, 700);
  h_res_eta_MuMinus->Draw();
  h_res_eta_MuPlus->Draw("same");
  c_res_eta->SaveAs("c_res_eta.png");
  
  TCanvas *c_res_eta_pT=new TCanvas("c_res_eta_pT", "c_res_eta_pT", 700, 700);
  h_res_eta_pT_MuMinus->RebinX(5);
  h_res_eta_pT_MuPlus->RebinX(5);
  TH1D *h_rmsRes_eta_pT_MuMinus=RMSY(h_res_eta_pT_MuMinus);
  TH1D *h_rmsRes_eta_pT_MuPlus=RMSY(h_res_eta_pT_MuPlus);
  h_rmsRes_eta_pT_MuMinus->SetYTitle("#sigma(#eta_{0}^{track} - #eta_{0}^{gen})");
  h_rmsRes_eta_pT_MuMinus->GetYaxis()->SetRangeUser(0, 0.003);
  h_rmsRes_eta_pT_MuMinus->SetLineColor(kBlue);
  h_rmsRes_eta_pT_MuPlus->SetLineColor(kRed);
  h_rmsRes_eta_pT_MuMinus->Draw("E0");
  h_rmsRes_eta_pT_MuPlus->Draw("E0 same");
  c_res_eta_pT->SaveAs("c_res_eta_pT.png");
  
}

