#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <set>

#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TLine.h>

std::string itoa(int i) 
{
  char res[10];
  sprintf(res, "%d", i);
  std::string ret(res);
  return ret;
}

void returnSectorEtaPhi(int eta, int phi, double &eta_min, double &eta_max, double &phi_min, double &phi_max)
{
  eta_min=10;
  eta_max=-10;
  if (eta==1)      {eta_min=-2; eta_max=-0.88;}
  else if (eta==2) {eta_min=-1.85; eta_max=0.17;}
  else if (eta==3) {eta_min=-1.38; eta_max=0.82;}
  else if (eta==4) {eta_min=-0.85; eta_max=1.37;}
  else if (eta==5) {eta_min=-0.19; eta_max=1.85;}
  else if (eta==6) {eta_min=0.97; eta_max=2;}
  
  phi_min=10;
  phi_max=-10;
  if (phi==1) {phi_min=-1.78; phi_max=-0.54;}
  else if (phi==2) {phi_min=-0.96; phi_max=0.25;}
  else if (phi==3) {phi_min=-0.22; phi_max=1.05;}
  else if (phi==4) {phi_min=0.59; phi_max=1.86;}
  else if (phi==5) {phi_min=1.36; phi_max=2.61;}
  else if (phi==6) {phi_min=2.16; phi_max=-2.85;}
  else if (phi==7) {phi_min=2.93; phi_max=-2.10;}
  else if (phi==8) {phi_min=-2.54; phi_max=-1.29;}
}

void SectorRangePlotter()
{ 
  // Read the Sector CSV file into data structures
  std::ifstream sectorFile("/Users/souvik/CMSPhase2Upgrades/trigger_sector_map.csv");
  std::string line;
  std::vector<std::pair<int , int> > sectorIDs;
  std::vector<std::set<int> > moduleIDs;
  char buff[8];
  getline(sectorFile, line);
  while (!sectorFile.eof())
  {
    getline(sectorFile, line);
    std::stringstream lineStream(line);
    int eta, phi;
    lineStream.getline(buff, 7, ','); eta=atoi(buff);
    lineStream.getline(buff, 7, ','); phi=atoi(buff);
    sectorIDs.push_back(make_pair(eta, phi));
    std::set<int> moduleInSector;
    // std::cout<<"Read: ";
    while (lineStream.getline(buff, 7, ','))
    {
      int moduleID=atoi(buff);
      moduleInSector.insert(moduleID);
    }
    moduleIDs.push_back(moduleInSector);
    // std::cout<<" --- "<<std::endl;
  }
  
  // std::cout<<"sectorIDs.size() = "<<sectorIDs.size()<< std::endl;
  // std::cout<<"moduleIDs.size() = "<<moduleIDs.size()<<std::endl;
  // std::cout<<sectorIDs.at(0).first<<std::endl;
  
  // Open a TChain to hold the contents of the ROOT file containing stubs
  TChain *tree = new TChain("L1TrackTrigger");
  tree->Add("/Users/souvik/CMSPhase2Upgrades/Samples/SLHC/GEN/612_SLHC6_MU/MU_612_SLHC6.root");
  
  // Declare and assign variables for the stub information
  int STUB_n;
  std::vector<int> *STUB_layer=0;
  std::vector<int> *STUB_ladder=0;
  std::vector<int> *STUB_module=0;
  std::vector<float> *STUB_pxGEN=0;
  std::vector<float> *STUB_pyGEN=0;
  std::vector<float> *STUB_etaGEN=0;
  tree->SetBranchAddress("STUB_n", &STUB_n);
  tree->SetBranchAddress("STUB_layer", &STUB_layer);
  tree->SetBranchAddress("STUB_ladder", &STUB_ladder);
  tree->SetBranchAddress("STUB_module", &STUB_module);
  tree->SetBranchAddress("STUB_pxGEN", &STUB_pxGEN);
  tree->SetBranchAddress("STUB_pyGEN", &STUB_pyGEN);
  tree->SetBranchAddress("STUB_etaGEN", &STUB_etaGEN);
  
  TH2F *h_eta_phi_Sector[7][9];
  for (unsigned int eta=1; eta<=6; ++eta)
  {
    for (unsigned int phi=1; phi<=8; ++phi)
    {
      std::string name="h_eta_phi_Sector_"+itoa(eta)+"_"+itoa(phi);
      h_eta_phi_Sector[eta][phi]=new TH2F(name.c_str(), name.c_str(), 100, -2.5, 2.5, 100, -3.14, 3.14);
    }
  }
  int nEvents=tree->GetEntries();
  // nEvents=1000;
  for (int i=0; i<nEvents; ++i)
  {
    tree->GetEvent(i);
    
    int moduleID=0;
    if (STUB_layer->size()>0)
    {
      // Make a number from the stub layer, ladder, module
      // layer*10000 + ladder*100 + module (From Stefano's email)
      moduleID=STUB_layer->at(0)*10000+STUB_ladder->at(0)*100+STUB_module->at(0);
    }
    // std::cout<<"moduleID = "<<moduleID<<std::endl;
    for (unsigned int j=0; j<48; ++j)
    {
      if (moduleIDs.at(j).find(moduleID)!=moduleIDs.at(j).end())
      {
        int sectorEta=sectorIDs.at(j).first;
        int sectorPhi=sectorIDs.at(j).second;
        // std::cout<<"Event "<<i<<", module found in Sector (eta, phi) = "<<sectorEta<<", "<<sectorPhi<<", track (eta, phi) = "<<STUB_etaGEN->at(0)<<", "<<atan2(STUB_pyGEN->at(0), STUB_pxGEN->at(0))<<std::endl;
        h_eta_phi_Sector[sectorEta][sectorPhi]->Fill(STUB_etaGEN->at(0), atan2(STUB_pyGEN->at(0), STUB_pxGEN->at(0)));
      }
    }
  }
  
  ofstream outfile;
  outfile.open("Sectoring/Sectoring.html");
  outfile<<"<html>"<<std::endl;
  outfile<<"<head>"<<std::endl;
  outfile<<"</head>"<<std::endl;
  outfile<<"<body>"<<std::endl;
  outfile<<"<h1 align='center'> Sectoring of the Tracker in the BE5D 6x8 Geometry </h1>"<<std::endl;
  outfile<<"<ul>"<<std::endl;
  outfile<<" <li> We iterate over a 1 million event muon gun sample. </li>"<<std::endl;
  outfile<<" <li> For each event, we consider the module ID of the stubs hit and look up the official sector they belong in. </li>"<<std::endl;
  outfile<<" <li> We fill the histogram of that sector with the generator level eta, phi of the track that produced the stubs. </li>"<<std::endl;
  outfile<<"</ul>"<<std::endl;
  outfile<<"<table border='1'>"<<std::endl;
  TFile *outFile=new TFile("Sectoring/Sectoring.root", "recreate");
  outFile->cd();
  for (unsigned int phi=1; phi<=8; ++phi)
  {
    outfile<<"<tr>"<<std::endl;
    for (unsigned int eta=1; eta<=6; ++eta)
    {
      double eta_min, eta_max, phi_min, phi_max;
      returnSectorEtaPhi(eta, phi, eta_min, eta_max, phi_min, phi_max);
      TCanvas *c=new TCanvas("c", "c", 700, 700);
      h_eta_phi_Sector[eta][phi]->Draw("colz");
      TLine eta_min_line(eta_min, -3.14, eta_min, 3.14); eta_min_line.Draw();
      TLine eta_max_line(eta_max, -3.14, eta_max, 3.14); eta_max_line.Draw();
      TLine phi_min_line(-2.5, phi_min, 2.5, phi_min); phi_min_line.Draw();
      TLine phi_max_line(-2.5, phi_max, 2.5, phi_max); phi_max_line.Draw();
      c->SaveAs(("Sectoring/c_eta_phi_"+itoa(eta)+"_"+itoa(phi)+".png").c_str());
      c->SaveAs(("Sectoring/c_eta_phi_"+itoa(eta)+"_"+itoa(phi)+".root").c_str());
      h_eta_phi_Sector[eta][phi]->Write();
      delete c;
      outfile<<" <td> "<<std::endl;
      outfile<<"   Eta = "<<eta<<", Phi = "<<phi<<"</br>"<<std::endl;
      outfile<<"   EtaRange={"<<eta_min<<", "<<eta_max<<"}, PhiRange={"<<phi_min<<", "<<phi_max<<"} <br/>"<<std::endl;
      outfile<<"   <img src='"<<("c_eta_phi_"+itoa(eta)+"_"+itoa(phi)+".png").c_str()<<"'/>"<<std::endl;
      outfile<<" </td>"<<std::endl;
    }
    outfile<<"</tr>"<<std::endl;
  }
  outfile<<"</body>"<<std::endl;
  outfile<<"</html>"<<std::endl;
  outFile->Close();
  
}
  
