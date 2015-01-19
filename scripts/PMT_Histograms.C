// 
// PMT_Histograms.C
// Martin McHugh
// 2014-06-04
//
// For plotting E and dE Energies and PMT responses. 

// Compile with:
//	> g++ -o PMT_Histograms PMT_Histograms.C `root-config --cflags --glibs`
// Run with
//	> ./PMT_Histograms


#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <cstdio>
#include <vector>
#include <map>
#include <utility>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector.h"
#include "TBrowser.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMath.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TString.h"
#include "TApplication.h"
#include "TChain.h"
#include "TSystem.h"
#include "TChainElement.h"
#include "TGraph.h"
#include "TF1.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TMultiGraph.h"
#include "TLeaf.h"

Double_t SkewNormal(Double_t* x, Double_t* par) {
  
  Double_t sqrt2 = TMath::Sqrt(2.0);
  Double_t argErfc = -1.0*par[3]*(x[0]-par[1])/(sqrt2*par[2]);
  Double_t signal = par[0]*TMath::Gaus(x[0], par[1], par[2])*TMath::Erfc(argErfc);
  return signal;

}

int main(Int_t argc, Char_t *argv[]) {

  Double_t MicroAmp = 6.241e15;		// # electrons/second
  Double_t MeVtoJoule = 1.602e-13;	// J/MeV 
  
  const char* FileDir = "/home/mjmchugh/Mott/MottG4/output/July_2014";
  const char* FileStem = "AsymRun";
  
  TChain* pChain = new TChain("Mott");
  for(Int_t i=1; i<101; i++) pChain->Add(Form("%s/%s_%i.root", FileDir, FileStem, i));

  TH1F* hLeft_E = new TH1F("hLeft_E","hLeft_E",1000,0,5);
  TH1F* hLeft_dE = new TH1F("hLeft_dE","hLeft_dE",1000,0,1.5);
  TH1F* hLeft_E_PMT = new TH1F("hLeft_E_PMT","hLeft_E_PMT",1000,0,5000);
  TH1F* hLeft_dE_PMT = new TH1F("hLeft_dE_PMT","hLeft_dE_PMT",1009,0,1500);

  Double_t Left_E, Left_dE;
  Int_t Left_E_PMT, Left_dE_PMT;  
  
  Int_t nEntries = pChain->GetEntries();

  std::cout << "Total number of hits: " << nEntries << std::endl;

  pChain->SetBranchAddress("Left_E",&Left_E);
  pChain->SetBranchAddress("Left_dE",&Left_dE);
  pChain->SetBranchAddress("Left_E_PMT",&Left_E_PMT);
  pChain->SetBranchAddress("Left_dE_PMT",&Left_dE_PMT);
 
  for(Int_t i=0; i<nEntries; i++) {
    pChain->GetEntry(i);
    if(Left_dE_PMT > 0 && Left_E_PMT > 0) {
      hLeft_E->Fill(Left_E);
      hLeft_dE->Fill(Left_dE);
      hLeft_E_PMT->Fill(Left_E_PMT);
      hLeft_dE_PMT->Fill(Left_dE_PMT);
    }
  }

  gStyle->SetOptFit(1); 	//Make sure fit stats appear

  hLeft_E->Fit("gaus","V","E1",4.0,5.0);
  hLeft_E_PMT->Fit("gaus","V","E1",4000,5000);

  //TF1* Fit1 = new TF1("SkewNormal1",SkewNormal,3.0,5.0,4);
  //Fit1->SetParNames("Constant","Mean","Sigma","Alpha");
  //Fit1->SetParameters(160.0,4.7,.03,-1.0);
  //Fit1->SetLineColor(kBlack);
  //Fit1->SetLineWidth(3);

  //TF1* Fit2 = new TF1("SkewNormal2",SkewNormal,3000.0,5000.0,4);
  //Fit2->SetParNames("Constant","Mean","Sigma","Alpha");
  //Fit2->SetParameters(160.0,4700.,90.0,-3.0);  
  //Fit2->SetLineColor(kBlack);
  //Fit2->SetLineWidth(3);

  TCanvas* c1 = new TCanvas("c1", "Left Detector Histograms", 1000, 1000);
  c1->Divide(2,2);
  c1->cd(1);
  //hLeft_E->GetXaxis()->SetTitle("LEFT_E [MeV]");
  hLeft_E->Draw();
  //hLeft_E->Fit("SkewNormal1","r");
  c1->cd(2);
  //hLeft_dE->GetXaxis()->SetTitle("LEFT_dE [MeV]");
  hLeft_dE->Draw();
  c1->cd(3);
  //hLeft_E_PMT->GetXaxis()->SetTitle("LEFT_E PEs");
  hLeft_E_PMT->Draw();
  //hLeft_E_PMT->Fit("SkewNormal2","r");
  c1->cd(4);
  //hLeft_dE_PMT->GetXaxis()->SetTitle("LEFT_dE PEs");
  hLeft_dE_PMT->Draw();
  c1->cd();

  // Save the plots as a ROOTfile to make things better!
  TFile* plotFile = new TFile("LeftDetectorPlots.root", "RECREATE");
  c1->Write();
  plotFile->Close();

  return 0;
}

