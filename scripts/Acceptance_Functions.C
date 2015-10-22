// 
// Acceptance_Functions.C
// Martin McHugh
// 2015-10-22
//
// Calculates and plots the acceptance function (N_hit/N_thrown)[var] for several
// variables

// Compile with:
//	> g++ -o Acceptance_Functions Acceptance_Functions.C `root-config --cflags --glibs`
// Run with
//	> ./Acceptance_Functions

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

  const Double_t pi = 3.1415926;			// [J/MeV] 
  const Double_t N_Beam = 6.241e12;			// electrons per second per micrioamp [Hz/uA]  
  const Double_t N_A = 6.022e23;			// molcules per mol [1/mol]
  const Double_t rho_Au = 19.30;			// Gold density [g/cm^3]
  const Double_t A_Au = 196.97;				// Gold atomic weight [g/mol]
  const Double_t targetCenterZ = -(3987.8);		// Target center position [um]
  Double_t d = 1.0;					// target_thickness [um]

  const char* FileDir = "/home/mjmchugh/Mott/MottG4/Mott_Polarimeter/build";
  TChain* pChain = new TChain("Mott");
  for(Int_t i=0; i<10; i++) pChain->Add(Form("%s/Single_%i.root", FileDir, i));
  
  // Hits
  TH1F* hHitX = new TH1F("hHitX","hHitX",100,-3.0,3.0);
  TH1F* hHitY = new TH1F("hHitY","hHitY",100,-3.0,3.0);
  TH1F* hHitZ = new TH1F("hHitZ","hHitZ",100,0,1.0);
  TH1F* hHitE = new TH1F("hHitE","hHitE",100,4.0,6.0);
  TH1F* hHitChi = new TH1F("hHitChi","hHitChi",100,170.0,175.0);
  TH1F* hHitPsi = new TH1F("hHitPsi","hHitPsi",100,-10.0,10.0);
  TH2F* hHitXY = new TH2F("hHitXY","hHitXY",100,-3.0,3.0,100,-3.0,3.0);
  TH2F* hHitZE = new TH2F("hHitZE","hHitZE",100,0,1.0,100,4.0,6.0);
  TH2F* hHitChiPsi = new TH2F("hHitChiPsi","hHitChiPsi",100,170.0,175.0,100,-10.0,10.0);
  
  // Thrown
  TH1F* hThrownX = new TH1F("hThrownX","hThrownX",100,-3.0,3.0);
  TH1F* hThrownY = new TH1F("hThrownY","hThrownY",100,-3.0,3.0);
  TH1F* hThrownZ = new TH1F("hThrownZ","hThrownZ",100,0,1.0);
  TH1F* hThrownE = new TH1F("hThrownE","hThrownE",100,4.0,6.0);
  TH1F* hThrownChi = new TH1F("hThrownChi","hThrownChi",100,170.0,175.0);
  TH1F* hThrownPsi = new TH1F("hThrownPsi","hThrownPsi",100,-10.0,10.0);
  TH2F* hThrownXY = new TH2F("hThrownXY","hThrownXY",100,-3.0,3.0,100,-3.0,3.0);
  TH2F* hThrownZE = new TH2F("hThrownZE","hThrownZE",100,0,1.0,100,4.0,6.0);
  TH2F* hThrownChiPsi = new TH2F("hThrownChiPsi","hThrownChiPsi",100,170.0,175.0,100,-10.0,10.0);
  
  // epsilon
  TH1F* hAccX = new TH1F("hAccX","hAccX",100,-3.0,3.0);
  TH1F* hAccY = new TH1F("hAccY","hAccY",100,-3.0,3.0);
  TH1F* hAccZ = new TH1F("hAccZ","hAccZ",100,0,1.0);
  TH1F* hAccE = new TH1F("hAccE","hAccE",100,4.0,6.0);
  TH1F* hAccChi = new TH1F("hAccChi","hAccChi",100,170.0,175.0);
  TH1F* hAccPsi = new TH1F("hAccPsi","hAccPsi",100,-10.0,10.0);
  TH2F* hAccXY = new TH2F("hAccXY","hAccXY",100,-3.0,3.0,100,-3.0,3.0);
  TH2F* hAccZE = new TH2F("hAccZE","hAccZE",100,0,1.0,100,4.0,6.0);
  TH2F* hAccChiPsi = new TH2F("hAccChiPsi","hAccChiPsi",100,170.0,175.0,100,-10.0,10.0);
  
  Int_t nEntries = pChain->GetEntries();
  Double_t X, Y, Z, E;
  Double_t Theta;
  Double_t Phi;
  Double_t Left_E, Left_dE;
  
  std::cout << nEntries << std::endl;
  
  pChain->SetBranchAddress("PrimaryVertexX",&X);
  pChain->SetBranchAddress("PrimaryVertexY",&Y);
  pChain->SetBranchAddress("PrimaryVertexZ",&Z);
  pChain->SetBranchAddress("PrimaryVertexKEprime",&E);
  pChain->SetBranchAddress("PrimaryVertexTheta",&Theta);
  pChain->SetBranchAddress("PrimaryVertexPhi",&Phi);
  pChain->SetBranchAddress("Left_E",&Left_E);
  pChain->SetBranchAddress("Left_dE",&Left_dE);
  
  // Now to fill All of those nice histograms
  Int_t nLeft = 0;
  Int_t nLeftHits = 0;
  for(Int_t i=0; i<nEntries; i++) {
    pChain->GetEntry(i);
    
    Double_t z = 1000.0*Z - (targetCenterZ - d/2.0); 		    
    
    if(-10<=Phi&&Phi<=10) {
      nLeft++;
      hThrownX->Fill(X);
      hThrownY->Fill(Y);
      hThrownZ->Fill(z);            
      hThrownE->Fill(E);
      hThrownChi->Fill(Theta);
      hThrownPsi->Fill(Phi);
      hThrownXY->Fill(X,Y);
      hThrownZE->Fill(z,E);
      hThrownChiPsi->Fill(Theta,Phi);
      if(Left_E>0 && Left_dE>0) {
        nLeftHits++;
        hHitX->Fill(X);
        hHitY->Fill(Y);
        hHitZ->Fill(z);            
        hHitE->Fill(E);
        hHitChi->Fill(Theta);
        hHitPsi->Fill(Phi);
        hHitXY->Fill(X,Y);
        hHitZE->Fill(z,E);
        hHitChiPsi->Fill(Theta,Phi);
      }
    }
    
  }
  // Epsilon is the ratio of the two
  hAccX->Divide(hHitX,hThrownX);
  hAccY->Divide(hHitY,hThrownY);
  hAccZ->Divide(hHitZ,hThrownZ);
  hAccE->Divide(hHitE,hThrownE);      
  hAccChi->Divide(hHitChi,hThrownChi);
  hAccPsi->Divide(hHitPsi,hThrownPsi);
  hAccXY->Divide(hHitXY,hThrownXY);
  hAccZE->Divide(hHitZE,hThrownZE);
  hAccChiPsi->Divide(hHitChiPsi,hThrownChiPsi);          
  
  TCanvas* cHitFunc = new TCanvas("cHitFunc", "Hits", 1200, 1200);
  cHitFunc->Divide(3,3);
  cHitFunc->cd(1);
  hHitX->Draw();
  cHitFunc->cd(4);
  hHitY->Draw();
  cHitFunc->cd(7);
  hHitXY->Draw("CONT");  
  cHitFunc->cd(2);
  hHitZ->Draw();
  cHitFunc->cd(5);
  hHitE->Draw();
  cHitFunc->cd(8);
  hHitZE->Draw("CONT");
  cHitFunc->cd(3);
  hHitChi->Draw();
  cHitFunc->cd(6);
  hHitPsi->Draw();
  cHitFunc->cd(9);
  hHitChiPsi->Draw();
  cHitFunc->cd();
  cHitFunc->Print("Hits.root");

  TCanvas* cThrownFunc = new TCanvas("cThrownFunc", "Thrown", 1200, 1200);
  cThrownFunc->Divide(3,3);
  cThrownFunc->cd(1);
  hThrownX->Draw();
  cThrownFunc->cd(4);
  hThrownY->Draw();
  cThrownFunc->cd(7);
  hThrownXY->Draw("CONT");  
  cThrownFunc->cd(2);
  hThrownZ->Draw();
  cThrownFunc->cd(5);
  hThrownE->Draw();
  cThrownFunc->cd(8);
  hThrownZE->Draw("CONT");
  cThrownFunc->cd(3);
  hThrownChi->Draw();
  cThrownFunc->cd(6);
  hThrownPsi->Draw();
  cThrownFunc->cd(9);
  hThrownChiPsi->Draw();
  cThrownFunc->cd();
  cThrownFunc->Print("Thrown.root");
  
  TCanvas* cAccFunc = new TCanvas("cAccFunc", "Acceptance_Functions", 1200, 800);
  cAccFunc->Divide(3,3);
  cAccFunc->cd(1);
  hAccX->Draw();
  cAccFunc->cd(4);
  hAccY->Draw();
  cAccFunc->cd(7);
  hAccXY->Draw("CONT");  
  cAccFunc->cd(2);
  hAccZ->Draw();
  cAccFunc->cd(5);
  hAccE->Draw();
  cAccFunc->cd(8);
  hAccZE->Draw("CONT");
  cAccFunc->cd(3);
  hAccChi->Draw();
  cAccFunc->cd(6);
  hAccPsi->Draw();
  cAccFunc->cd(9);
  hAccChiPsi->Draw();
  cAccFunc->cd();
  cAccFunc->Print("AcceptanceFunction.root");
  
  return 0;
}
  
