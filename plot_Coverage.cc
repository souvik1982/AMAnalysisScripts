#include <TROOT.h>
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"

void plot_Coverage()
{
  TFile *file_ss32=new TFile("plots_BankGeneration_ss32.root");
  TGraph *g_coverage_bankSize_ss32=(TGraph*)file_ss32->Get("g_coverage_bankSize");
  
  TFile *file_ss64=new TFile("plots_BankGeneration_ss64.root");
  TGraph *g_coverage_bankSize_ss64=(TGraph*)file_ss64->Get("g_coverage_bankSize");
  
  TFile *file_ss128=new TFile("plots_BankGeneration_ss128.root");
  TGraph *g_coverage_bankSize_ss128=(TGraph*)file_ss128->Get("g_coverage_bankSize");
  
  TFile *file_ss256=new TFile("plots_BankGeneration_ss256.root");
  TGraph *g_coverage_bankSize_ss256=(TGraph*)file_ss256->Get("g_coverage_bankSize");
  
  TFile *file_ss512=new TFile("plots_BankGeneration_ss512.root");
  TGraph *g_coverage_bankSize_ss512=(TGraph*)file_ss512->Get("g_coverage_bankSize");
  
  TFile *file_ss1024=new TFile("plots_BankGeneration_ss1024.root");
  TGraph *g_coverage_bankSize_ss1024=(TGraph*)file_ss1024->Get("g_coverage_bankSize");
  
  TCanvas *c_coverage_bankSize=new TCanvas("c_coverage_bankSize", "c_coverage_bankSize", 700, 700);
  g_coverage_bankSize_ss32->SetTitle("CMS Preliminary Phase II Simulation; # patterns; running estimate for coverage");
  g_coverage_bankSize_ss32->GetYaxis()->SetTitleOffset(1.4);
  g_coverage_bankSize_ss32->GetYaxis()->SetRangeUser(0, 1.1);
  g_coverage_bankSize_ss32->SetLineColor(kRed); g_coverage_bankSize_ss32->SetLineColor(kRed);     g_coverage_bankSize_ss32->Draw("AC");    
  g_coverage_bankSize_ss64->SetLineWidth(2);    g_coverage_bankSize_ss64->SetLineColor(kViolet);  g_coverage_bankSize_ss64->Draw("C");
  g_coverage_bankSize_ss128->SetLineWidth(2);   g_coverage_bankSize_ss128->SetLineColor(kGreen);  g_coverage_bankSize_ss128->Draw("C");
  g_coverage_bankSize_ss256->SetLineWidth(2);   g_coverage_bankSize_ss256->SetLineColor(kBlue);   g_coverage_bankSize_ss256->Draw("C");
  g_coverage_bankSize_ss512->SetLineWidth(2);   g_coverage_bankSize_ss512->SetLineColor(kCyan);   g_coverage_bankSize_ss512->Draw("C");
  g_coverage_bankSize_ss1024->SetLineWidth(2);  g_coverage_bankSize_ss1024->SetLineColor(kRed+2); g_coverage_bankSize_ss1024->Draw("C");
  TLegend *leg=new TLegend(0.7, 0.5, 0.9, 0.2); leg->SetBorderSize(0); leg->SetFillColor(kWhite);
  leg->AddEntry(g_coverage_bankSize_ss32, "ss32", "l");
  leg->AddEntry(g_coverage_bankSize_ss64, "ss64", "l");
  leg->AddEntry(g_coverage_bankSize_ss128, "ss128", "l");
  leg->AddEntry(g_coverage_bankSize_ss256, "ss256", "l");
  leg->AddEntry(g_coverage_bankSize_ss512, "ss512", "l");
  leg->AddEntry(g_coverage_bankSize_ss1024, "ss1024", "l");
  leg->Draw();
  c_coverage_bankSize->SaveAs("c_coverage_bankSize.png");
  
}
  
