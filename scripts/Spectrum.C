// 
// Spectrum.C
// Martin McHugh
// 2015-08-10
//
// Makes combined spectrum, weighted by rates, to compare to data

// Compile with:
//	> g++ -o Spectrum Spectrum.C `root-config --cflags --glibs`
// Run with
//	> ./Spectrum [target thickness]


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
#include "TF2.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TProfile.h"
#include "TMultiGraph.h"
#include "TLeaf.h"
#include "TProfile2D.h"

int main(Int_t argc, Char_t *argv[]) {

  // Assumes 944 um foil Rates are in Hz/uA
  Double_t Rate_L1 = 91.0;
  Double_t Rate_R1 = 282.0;
  Double_t Rate_L2 = 77.5;
  Double_t Rate_R2 = 38.0;

  //Get one rootfile each for Single and Double Scattering
  const char* FileDir = "/lustre/expphy/volatile/hallc/qweak/mjmchugh/Mott/Round2";
  TChain* pChain1 = new TChain("Mott");
  TChain* pChain2 = new TChain("Mott");
  pChain1->Add(Form("%s/Single_0.944um_1.root",FileDir));
  pChain2->Add(Form("%s/Double_0.944um_1.root",FileDir));

  TH1F* hLeft_E1 = new TH1F("hLeft_E1","hLeft_E1",600,0,6);
  TH1F* hRight_E1 = new TH1F("hRight_E1","hRight_E1",600,0,6);
  TH1F* hLeft_E2 = new TH1F("hLeft_E2","hLeft_E2",600,0,6);
  TH1F* hRight_E2 = new TH1F("hRight_E2","hRight_E2",600,0,6);

  /*
  Int_t nEntries1 = pChain1->GetEntries();
  Int_t nEntries2 = pChain2->GetEntries();
  Double_t Left_E1, Left_dE1;
  Double_t Right_E1, Right_dE1;
  Double_t Left_E2, Left_dE2;
  Double_t Right_E2, Right_dE2;
  pChain1->SetBranchAddress("Left_E",&Left_E1);
  pChain1->SetBranchAddress("Left_dE",&Left_dE1);
  pChain1->SetBranchAddress("Right_E",&Right_E1);
  pChain1->SetBranchAddress("Right_dE",&Right_dE1);
  pChain2->SetBranchAddress("Left_E",&Left_E2);
  pChain2->SetBranchAddress("Left_dE",&Left_dE2);
  pChain2->SetBranchAddress("Right_E",&Right_E2);
  pChain2->SetBranchAddress("Right_dE",&Right_dE2);

  for(Int_t i=0; i<5000000; i++) {
    pChain1->GetEntry(i);
    if(Left_E1>0 && Left_dE1>1) hLeft_E1->Fill(Left_E1);
    if(Right_E1>0 && Right_dE1>1) hRight_E1->Fill(Right_E1);
    if(i%100000==0) std::cout << 1 << "\t" << i << std::endl;
  }

  for(Int_t i=0; i<5000000; i++) {
    pChain2->GetEntry(i);
    if(Left_E2>0 && Left_dE2>1) hLeft_E2->Fill(Left_E2);
    if(Right_E2>0 && Right_dE2>1) hRight_E2->Fill(Right_E2);
    if(i%100000==0) std::cout << 2 << "\t" << i << std::endl;
  }
  */

  TCanvas* c1 = new TCanvas("c1","c1",1000,1000);
  c1->Divide(2,2);
  c1->cd(1);
  pChain1->Draw("Left_E>>hLeft_E1","Left_E!=0&&Left_dE!=0","");
  hLeft_E1->Scale(Rate_L1);
  c1->cd(2);
  pChain1->Draw("Right_E>>hRight_E1","Right_E!=0&&Right_dE!=0","");
  hRight_E1->Scale(Rate_R1);
  c1->cd(3);
  pChain2->Draw("Left_E>>hLeft_E2","Left_E!=0&&Left_dE!=0","");
  hLeft_E2->Scale(Rate_L2);
  c1->cd(4);
  pChain2->Draw("Right_E>>hRight_E2","Right_E!=0&&Right_dE!=0","");
  hRight_E2->Scale(Rate_R2);
  c1->cd();
  c1->Print("ScaledSpectra.root");
  c1->Print("ScaledSpectra.png");
  delete c1;

  TH1F* hLeft_ESum = new TH1F("hLeft_ESum","hLeft_ESum",600,0,6);
  TH1F* hRight_ESum = new TH1F("hRight_ESum","hRight_ESum",600,0,6);
  TH1F* hSum = new TH1F("hSum","hSum",600,0,6);
  TH1F* hDiff = new TH1F("hDiff","hDiff",600,0,6);  
  TH1F* hAsym = new TH1F("hAsym","hAsym",600,0,6);

  hLeft_ESum->Add(hLeft_E1,hLeft_E2);
  hRight_ESum->Add(hRight_E1,hRight_E2);
  hSum->Add(hLeft_ESum,hRight_ESum);
  hDiff->Add(hLeft_ESum,hRight_ESum,1.0,-1.0);
  hAsym->Divide(hDiff,hSum);
  
  TCanvas* c2 = new TCanvas("c2","c2",1000,500);
  c2->Divide(2,1);
  c2->cd(1);
  hLeft_ESum->Draw();
  c2->cd(2);
  hRight_ESum->Draw();
  c2->cd();
  c2->Print("SummedSpectra.root");
  c2->Print("SummedSpectra.png");
  delete c2;

  TCanvas* c3 = new TCanvas("c3","c3",1000,500);
  c3->Divide(2,1);
  c3->cd(1);
  hDiff->Draw();
  c3->cd(2);
  hAsym->Draw();
  c3->cd();
  c3->Print("DiffSpectra.root");
  c3->Print("DiffSpectra.png");
  delete c3;   

  TCanvas* c4 = new TCanvas("c4","c4",1000,1000);
  c4->cd(1);
  hSum->Draw();  
  c4->Draw();
  c4->Print("CombinedSpectra.root");
  c4->Print("CombinedSpectra.png");
  delete c4;

  return 0;
}
