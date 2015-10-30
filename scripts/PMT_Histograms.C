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
  
  const char* FileDir = "/lustre/expphy/volatile/hallc/qweak/mjmchugh/Mott";
  //const char* FileDir = "/u/home/mjmchugh/MottG4/Mott_Polarimeter/build";
  const char* FileStem = argv[1];
  
  TChain* pChain1 = new TChain("Mott");
  TChain* pChain2 = new TChain("Mott");
  for(Int_t i=1; i<11; i++) {
    pChain1->Add(Form("%s/Au_1um_%i.root", FileDir, i));
    pChain2->Add(Form("%s/Au2_1um_%i.root", FileDir, i));
  }

  TH1F* hLeft_E1 = new TH1F("hLeft_E1","hLeft_E1",1000,0,8);
  TH1F* hRight_E1 = new TH1F("hRight_E1","hRight_E1",1000,0,8);
  TH1F* hLeft_E2 = new TH1F("hLeft_E2","hLeft_E2",1000,0,8);
  TH1F* hRight_E2 = new TH1F("hRight_E2","hRight_E2",1000,0,8);
  TH1F* hCS1L = new TH1F("hCS1L","hCS1L",1000,0,30e-24);
  TH1F* hCS1R = new TH1F("hCS1R","hCS1R",1000,0,30e-24);
  TH1F* hCS2L = new TH1F("hCS2L","hCS2L",1000,0,30e-24);
  TH1F* hCS2R = new TH1F("hCS2R","hCS2R",1000,0,30e-24);
  TH1F* hTheta1L = new TH1F("hTheta1L","hTheta1L",180,0,180);
  TH1F* hTheta1R = new TH1F("hTheta1R","hTheta1R",180,0,180);
  TH1F* hTheta2L = new TH1F("hTheta2L","hTheta2L",180,0,180);
  TH1F* hTheta2R = new TH1F("hTheta2R","hTheta2R",180,0,180);
  TH1F* hPhi1L = new TH1F("hPhi1L","hPhi1L",360,-180,180);
  TH1F* hPhi1R = new TH1F("hPhi1R","hPhi1R",360,-180,180);
  TH1F* hPhi2L = new TH1F("hPhi2L","hPhi2L",360,-180,180);
  TH1F* hPhi2R = new TH1F("hPhi2R","hPhi2R",360,-180,180);

  Double_t Left_E1, Right_E1, Left_E2, Right_E2;
  Double_t PrimaryCrossSection, SecondaryCrossSection;
  Double_t PrimaryVertexTheta, PrimaryVertexPhi, SecondaryVertexTheta, SecondaryVertexPhi;
  Int_t nEntries1 = pChain1->GetEntries();
  Int_t nEntries2 = pChain2->GetEntries();

  std::cout << "Total number of hits: " << nEntries1 << " " << nEntries2 << std::endl;

  pChain1->SetBranchAddress("Left_E",&Left_E1);
  pChain1->SetBranchAddress("Right_E",&Right_E1);
  pChain2->SetBranchAddress("Left_E",&Left_E2);
  pChain2->SetBranchAddress("Right_E",&Right_E2);
  pChain2->SetBranchAddress("PrimaryCrossSection",&PrimaryCrossSection);
  pChain2->SetBranchAddress("SecondaryCrossSection",&SecondaryCrossSection);
  pChain2->SetBranchAddress("PrimaryVertexTheta",&PrimaryVertexTheta);
  pChain2->SetBranchAddress("SecondaryVertexTheta",&SecondaryVertexTheta); 
  pChain2->SetBranchAddress("PrimaryVertexPhi",&PrimaryVertexPhi);
  pChain2->SetBranchAddress("SecondaryVertexPhi",&SecondaryVertexPhi);

  for(Int_t i=0; i<nEntries1; i++) {
    pChain1->GetEntry(i);
    if(Left_E1 > 0) {
      hLeft_E1->Fill(Left_E1);
      //hCSL->Fill(PrimaryCrossSection*SecondaryCrossSection);
    }
    if(Right_E1 > 0) {
      hRight_E1->Fill(Right_E1);
      //hCSR->Fill(PrimaryCrossSection*SecondaryCrossSection);
    }
    if (i%100000 == 0) std::cout << i << std::endl;
  }
  for(Int_t i=0; i<nEntries2; i++) {
    pChain2->GetEntry(i);
    if(Left_E2 > 0) {
      hLeft_E2->Fill(Left_E2);
      hCS1L->Fill(PrimaryCrossSection);
      hCS2L->Fill(SecondaryCrossSection);
      hTheta1L->Fill(PrimaryVertexTheta);
      hTheta2L->Fill(SecondaryVertexTheta);
      hPhi1L->Fill(PrimaryVertexPhi);
      hPhi2L->Fill(SecondaryVertexPhi);
    }
    if(Right_E2 > 0) {
      hRight_E2->Fill(Right_E2);
      hCS1R->Fill(PrimaryCrossSection);
      hCS2R->Fill(SecondaryCrossSection);
      hTheta1R->Fill(PrimaryVertexTheta);
      hTheta2R->Fill(SecondaryVertexTheta);
      hPhi1R->Fill(PrimaryVertexPhi);
      hPhi2R->Fill(SecondaryVertexPhi);
    }
    if (i%100000 == 0) std::cout << i << std::endl;
  }

  TCanvas* c1 = new TCanvas("c1", "Left Detector Histograms", 1200, 800);
  c1->Divide(3,2);
  c1->cd(1);
  hCS1L->Draw("");
  c1->cd(2);
  hTheta1L->Draw("");
  c1->cd(3);
  /*hLeft_E_PMT->Scale(1.0/hLeft_E_PMT->GetEntries());
  hLeft_E1->GetXaxis()->SetTitle("LEFT_E");
  hLeft_E1->Scale(1.0/hLeft_E1->GetEntries());
  hLeft_E1->Draw("");
  hLeft_E2->Scale(1.0/hLeft_E2->GetEntries());
  hLeft_E2->SetLineColor(kRed);
  hLeft_E2->Draw("same");*/
  hPhi1L->Draw("");
  c1->cd(4);
  hCS2L->Draw("");
  c1->cd(5);
  hTheta2L->Draw("");
  c1->cd(6);
  /*hRight_E_PMT->Scale(1.0/hRight_E_PMT->GetEntries());
  hRight_E1->GetXaxis()->SetTitle("RIGHT_E");
  hRight_E1->Scale(1.0/hRight_E1->GetEntries());
  hRight_E1->Draw("");
  hRight_E2->Scale(1.0/hRight_E2->GetEntries());
  hRight_E2->SetLineColor(kRed);
  hRight_E2->Draw("same");*/
  //hPhi2L->Draw("");
  c1->cd();
  c1->Print(Form("Left_Plots.root"));

  TCanvas* c2 = new TCanvas("c2", "Right Detector Histograms", 1200, 800);
  c2->Divide(3,2);
  c2->cd(1);
  hCS1R->Draw("");
  c2->cd(2);
  hTheta1R->Draw("");
  c2->cd(3);
  hPhi1R->Draw("");
  c2->cd(4);
  hCS2R->Draw("");
  c2->cd(5);
  hTheta2R->Draw("");
  c2->cd(6);
  //hPhi2R->Draw("");
  c2->cd();
  c2->Print(Form("Right_Plots.root"));

  /*Asymmetry calculation
  Double_t L = hLeft_E2->GetEntries();
  Double_t R = hRight_E2->GetEntries();
  Double_t A = (L - R)/(L + R);
  Double_t dA = TMath::Abs(A)*TMath::Sqrt((L+R)*(1.0/((L-R)*(L-R)) + 1.0/((L+R)*(L+R))));

  std::cout << L << "\t" << R << "\t" << A << "\t" << dA << std::endl;
  
  /*Rate calculation
  Double_t Lum = 3.789e35;			// [Hz/(cm^3 uA)]
  Double_t thickness = 1.0e-4;			// 1.0 um [cm]
  Int_t hitsL = hCSL->GetEntries();
  Int_t hitsR = hCSR->GetEntries();
  Double_t avg_cs_L = hCSL->GetMean();
  Double_t avg_cs_R = hCSR->GetMean();
  Double_t phase_space = 2*0.00397479;

  std::cout << "Rate L = " << Lum*thickness*avg_cs_L*hitsL*phase_space/10000000.0 << " Hz/uA" <<  std::endl;
  std::cout << "Rate R = " << Lum*thickness*avg_cs_R*hitsR*phase_space/10000000.0 << " Hz/uA" <<  std::endl;
  */
  return 0;
}

