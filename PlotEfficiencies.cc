#include <TFile.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TGraphAsymmErrors.h>
#include <TStyle.h>

int rebin=1;

void PlotEfficiencies()
{
  TFile *file=new TFile("EfficiencyPlots.root");
  TH1F *h_eff_eta_num=(TH1F*)file->Get("h_eff_eta_num"); h_eff_eta_num->Rebin(4);
  TH1F *h_eff_eta_den=(TH1F*)file->Get("h_eff_eta_den"); h_eff_eta_den->Rebin(4);
  TH1F *h_eff_phi_num=(TH1F*)file->Get("h_eff_phi_num"); h_eff_phi_num->Rebin(1);
  TH1F *h_eff_phi_den=(TH1F*)file->Get("h_eff_phi_den"); h_eff_phi_den->Rebin(1);
  TH1F *h_eff_pT_num=(TH1F*)file->Get("h_eff_pT_num"); h_eff_pT_num->Rebin(2);
  TH1F *h_eff_pT_den=(TH1F*)file->Get("h_eff_pT_den"); h_eff_pT_den->Rebin(2);
  TH2F *h_eff_eta_phi_num=(TH2F*)file->Get("h_eff_eta_phi_num"); h_eff_eta_phi_num->Rebin2D(4,4);
  TH2F *h_eff_eta_phi_den=(TH2F*)file->Get("h_eff_eta_phi_den"); h_eff_eta_phi_den->Rebin2D(4,4);
  TH2F *h_eff_px_py_num=(TH2F*)file->Get("h_eff_px_py_num"); h_eff_px_py_num->Rebin2D(4,4);
  TH2F *h_eff_px_py_den=(TH2F*)file->Get("h_eff_px_py_den"); h_eff_px_py_den->Rebin2D(4,4);
  
  TGraphAsymmErrors *g_eff_eta=new TGraphAsymmErrors(h_eff_eta_num, h_eff_eta_den);
  g_eff_eta->SetTitle("Pattern Recognition Efficiency; Generated track #eta; Efficiency");
  
  TGraphAsymmErrors *g_eff_phi=new TGraphAsymmErrors(h_eff_phi_num, h_eff_phi_den, "cp");
  g_eff_phi->SetTitle("Pattern Recognition Efficiency; Generated track #phi; Efficiency");
  
  TGraphAsymmErrors *g_eff_pT=new TGraphAsymmErrors(h_eff_pT_num, h_eff_pT_den);
  g_eff_pT->SetTitle("Pattern Recognition Efficiency; Generated track p_{T} [GeV]; Efficiency");
  
  TH2F *h_eff_px_py=(TH2F*)h_eff_px_py_num->Clone("h_eff_px_py");
  h_eff_px_py->Divide(h_eff_px_py_den);
  h_eff_px_py->SetTitle("Pattern Recognition Efficiency; Generated Track p_{x} [GeV]; Generated Track p_y [GeV]");
  
  TH2F *h_eff_eta_phi=(TH2F*)h_eff_eta_phi_num->Clone("h_eff_eta_phi");
  h_eff_eta_phi->Divide(h_eff_eta_phi_den);
  h_eff_eta_phi->SetTitle("Pattern Recognition Efficiency; Generated Track #eta; Generated Track #phi");
  
  gStyle->SetErrorX(0);
  gStyle->SetOptStat(00000);
  
  TCanvas *c_eff_eta=new TCanvas("c_eff_eta", "c_eff_eta", 700, 700);
  g_eff_eta->SetLineColor(kBlack);
  g_eff_eta->SetMarkerStyle(20);
  g_eff_eta->Draw("ACP");
  g_eff_eta->SetMinimum(0);
  c_eff_eta->SaveAs("c_eff_eta.png");
  
  TCanvas *c_eff_phi=new TCanvas("c_eff_phi", "c_eff_phi", 700, 700);
  g_eff_phi->SetLineColor(kBlack);
  g_eff_phi->SetMarkerStyle(20);
  g_eff_phi->Draw("ACP");
  g_eff_phi->SetMinimum(0);
  c_eff_phi->SaveAs("c_eff_phi.png");
  
  TCanvas *c_eff_pT=new TCanvas("c_eff_pT", "c_eff_pT", 700, 700);
  g_eff_pT->SetLineColor(kBlack);
  g_eff_pT->SetMarkerStyle(20);
  g_eff_pT->Draw("ACP");
  c_eff_pT->SetLogx();
  c_eff_pT->SaveAs("c_eff_pT.png");
  
  TCanvas *c_eff_eta_phi=new TCanvas("c_eff_eta_phi", "c_eff_eta_phi", 1000, 700);
  h_eff_eta_phi->Draw("colz");
  c_eff_eta_phi->SaveAs("c_eff_eta_phi.png");
  
  TCanvas *c_eff_px_py=new TCanvas("c_eff_px_py", "c_eff_px_py", 700, 700);
  h_eff_px_py->Draw("colz");
  c_eff_px_py->SaveAs("c_eff_px_py.png");
  
}
