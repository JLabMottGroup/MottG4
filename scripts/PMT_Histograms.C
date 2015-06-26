// 
// PMT_Histograms.C
// Martin McHugh
// 2014-06-04
//
// For plotting E and dE Energies and PMT responses. 

// Compile with:
//	> g++ -o PMT_Histograms PMT_Histograms.C `root-config --cflags --glibs`
// Run with
//	> ./PMT_Histograms FileStem


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
  
  //const char* FileDir = "/lustre/expphy/volatile/hallc/qweak/mjmchugh/Mott";
  const char* FileDir = "/u/home/mjmchugh/MottG4/Mott_Polarimeter/build";
  const char* FileStem = argv[1];
  
  TChain* pChain = new TChain("Mott");
  //for(Int_t i=2; i<3; i++) 
  pChain->Add(Form("%s/%s_%i.root", FileDir, FileStem, 2));

  TH1F* hLeft_E = new TH1F("hLeft_E","hLeft_E",1000,0,8);
  TH1F* hRight_E = new TH1F("hRight_E_PMT","hRight_E",1000,0,8);

  Double_t Left_E, Right_E;
  Double_t PrimaryVertexX, PrimaryVertexY, PrimaryVertexZ; 
  Int_t nEntries = pChain->GetEntries();

  std::cout << "Total number of hits: " << nEntries << std::endl;

  pChain->SetBranchAddress("Left_E",&Left_E);
  pChain->SetBranchAddress("Right_E",&Right_E);
 
  for(Int_t i=0; i<nEntries; i++) {
    pChain->GetEntry(i);
    if(Left_E > 4.0) hLeft_E->Fill(Left_E);
    if(Right_E > 4.0) hRight_E->Fill(Right_E);
  }

  //gStyle->SetOptFit(1); 	//Make sure fit stats appear

  //hLeft_E->Fit("gaus","V","E1",4.0,5.0);
  //hLeft_E_PMT->Fit("gaus","V","E1",4000,5000);

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

  TCanvas* c1 = new TCanvas("c1", "E_PMT Detector Histograms", 1000, 500);
  c1->Divide(2,1);
  c1->cd(1);
  //hLeft_E_PMT->Scale(1.0/hLeft_E_PMT->GetEntries());
  hLeft_E->GetXaxis()->SetTitle("LEFT_E");
  hLeft_E->Draw();
  c1->cd(2);
  //hRight_E_PMT->Scale(1.0/hRight_E_PMT->GetEntries());
  hRight_E->GetXaxis()->SetTitle("RIGHT_E");
  hRight_E->Draw();
  c1->cd(3);
  c1->cd();

  //Asymmetry calculation
  Double_t L = hLeft_E->GetEntries();
  Double_t R = hRight_E->GetEntries();
  Double_t A = (L - R)/(L + R);
  Double_t dA = TMath::Abs(A)*TMath::Sqrt((L+R)*(1.0/((L-R)*(L-R)) + 1.0/((L+R)*(L+R))));

  std::cout << L << "\t" << R << "\t" << A << "\t" << dA << std::endl;
  
  // Save the plots as a ROOTfile to make things better!
  TFile* plotFile = new TFile(Form("%s_EDetectorPlots.root", FileStem), "RECREATE");
  c1->Write();
  plotFile->Close();

  return 0;
}

