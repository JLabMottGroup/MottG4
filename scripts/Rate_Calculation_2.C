// 
// Rate_Calculation_2.C
// Martin McHugh
// 2015-10-06
//
// Performing Monte Carlo Integration of the detector rates using simulated 
// double-scattering events.

// Compile with:
//	> g++ -o Rate_Calculation_2 Rate_Calculation_2.C `root-config --cflags --glibs`
// Run with
//	> ./Rate_Calculation_2 0.05


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
#include "TProfile2D.h"

Double_t SkewNormal(Double_t* x, Double_t* par) {
  
  Double_t sqrt2 = TMath::Sqrt(2.0);
  Double_t argErfc = -1.0*par[3]*(x[0]-par[1])/(sqrt2*par[2]);
  Double_t signal = par[0]*TMath::Gaus(x[0], par[1], par[2])*TMath::Erfc(argErfc);
  return signal;

}

int main(Int_t argc, Char_t *argv[]) {

  const Double_t pi = 3.1415926;			// [J/MeV] 
  const Double_t N_Beam = 6.241e12;			// electrons per second per micrioamp [Hz/uA]  
  const Double_t N_A = 6.022e23;			// molcules per mol [1/mol]
  const Double_t rho_Au = 19.30;			// Gold density [g/cm^3]
  const Double_t A_Au = 196.97;				// Gold atomic weight [g/mol]
  const Double_t targetCenterZ = -(3.9878);		// Target center position [mm]
  Double_t d = 1.0e-4*atof(argv[1]);			// target_thickness [cm]

  const char* FileDir = "/home/mjmchugh/Mott/MottG4/Mott_Polarimeter/build";
  TChain* pChain = new TChain("Mott");
  pChain->Add(Form("%s/Double_%sum.root", FileDir,argv[1]));

  TH2F* hNhit_Left = new TH2F("hNhit_Left","hNhit_Left",100,170.0,175.0,100,-10.0,10.0);
  TH2F* hNhit_Right = new TH2F("hNhit_Right","hNhit_Right",100,170.0,175.0,100,170.0,190.0);
  TH2F* hNthrown_Left = new TH2F("hNthrown_Left","hNthrown_Left",100,170.0,175.0,100,-10.0,10.0);
  TH2F* hNthrown_Right = new TH2F("hNthrown_Right","hNthrown_Right",100,170.0,175.0,100,170.0,190.0);
  TH2F* hAccFunc_Left = new TH2F("hAccFunc_Left","hAccFunc_Left",100,170.0,175.0,100,-10.0,10.0);
  TH2F* hAccFunc_Right = new TH2F("hAccFunc_Right","hAccFunc_Right",100,170.0,175.0,100,170.0,190.0);

  Double_t Left_E, Left_dE;
  Double_t Right_E, Right_dE;
  Double_t PrimaryCrossSection;
  Double_t PrimaryVertexTheta;
  Double_t PrimaryVertexPhi;
  Double_t PrimaryVertexX, PrimaryVertexY, PrimaryVertexZ;
  Double_t SecondaryCrossSection;
  Double_t SecondaryVertexTheta;
  Double_t SecondaryVertexPhi;
  Double_t SecondaryVertexX, SecondaryVertexY, SecondaryVertexZ;
  Double_t Theta;
  Double_t Phi;
  Int_t nEntries = pChain->GetEntries();

  pChain->SetBranchAddress("Left_E",&Left_E);
  pChain->SetBranchAddress("Left_dE",&Left_dE);
  pChain->SetBranchAddress("Right_E",&Right_E);
  pChain->SetBranchAddress("Right_dE",&Right_dE);
  pChain->SetBranchAddress("PrimaryCrossSection",&PrimaryCrossSection);
  pChain->SetBranchAddress("PrimaryVertexTheta",&PrimaryVertexTheta);
  pChain->SetBranchAddress("PrimaryVertexPhi",&PrimaryVertexPhi);
  pChain->SetBranchAddress("PrimaryVertexX",&PrimaryVertexX);
  pChain->SetBranchAddress("PrimaryVertexY",&PrimaryVertexY);
  pChain->SetBranchAddress("PrimaryVertexZ",&PrimaryVertexZ);
  pChain->SetBranchAddress("SecondaryCrossSection",&SecondaryCrossSection);
  pChain->SetBranchAddress("SecondaryVertexTheta",&SecondaryVertexTheta);
  pChain->SetBranchAddress("SecondaryVertexPhi",&SecondaryVertexPhi);
  pChain->SetBranchAddress("SecondaryVertexX",&SecondaryVertexX);
  pChain->SetBranchAddress("SecondaryVertexY",&SecondaryVertexY);
  pChain->SetBranchAddress("SecondaryVertexZ",&SecondaryVertexZ);
  pChain->SetBranchAddress("Theta",&Theta);
  pChain->SetBranchAddress("Phi",&Phi);

  // First we need to simply get the acceptance function epsilon(theta,phi)
  for(Int_t i=0; i<nEntries; i++) {
    pChain->GetEntry(i);
    
    if(170<=Phi && Phi<180) {
      hNthrown_Right->Fill(Theta,Phi);
      if(Right_E>0 && Right_dE>0) hNhit_Right->Fill(Theta,Phi);
    } else if(-180<=Phi && Phi<=-170) {
      hNthrown_Right->Fill(Theta,Phi+360.0);
      if(Right_E>0 && Right_dE>0) hNhit_Right->Fill(Theta,Phi+360.0);
    } else {
      hNthrown_Left->Fill(Theta,Phi);
      if(Left_E>0 && Left_dE>0) hNhit_Left->Fill(Theta,Phi);
    }

  }
  // Epsilon is the ratio of the two
  hAccFunc_Left->Divide(hNhit_Left,hNthrown_Left);
  hAccFunc_Right->Divide(hNhit_Right,hNthrown_Right);  
  
  // Now we must actually calculate the integral using a Monte Carlo estimator
  Double_t nLeft = 0;
  Double_t nRight = 0;
  Double_t sumLeft = 0;
  Double_t sumRight = 0;
  for(Int_t i=0; i<nEntries; i++) {
    pChain->GetEntry(i);
        
    Double_t z = 0.1*PrimaryVertexZ - (0.1*targetCenterZ - d/2.0); 		// cm (from mm)
    Double_t xi_max = (d - z)/TMath::Abs(TMath::Cos(PrimaryVertexTheta));
    if (xi_max > 0.0157) xi_max = 0.0157;
    Double_t CS1 = PrimaryCrossSection;
    Double_t CS2 = SecondaryCrossSection;
    Double_t epsilon;
     
    // Find if it's an event to the left or the right
    if (-10<=Phi && Phi<=10) {
      nLeft++;
      epsilon = hAccFunc_Left->GetBinContent(hAccFunc_Left->FindBin(Theta,Phi));
      sumLeft += xi_max*CS1*CS2*epsilon;
    } else {
      nRight++;
      epsilon = hAccFunc_Right->GetBinContent(hAccFunc_Right->FindBin(Theta,Phi));
      sumRight += xi_max*CS1*CS2*epsilon;
    }
  }    
  Double_t estLeft = sumLeft/nLeft;
  Double_t estRight = sumRight/nRight;
  
  // Now we must calculate the error using the variance
  Double_t varSumLeft = 0;
  Double_t varSumRight = 0;  
  for(Int_t i=0; i<nEntries; i++) {
    pChain->GetEntry(i);
        
    Double_t z = 0.1*PrimaryVertexZ - (0.1*targetCenterZ - d/2.0); 		// cm (from mm)
    Double_t xi_max = (d - z)/TMath::Abs(TMath::Cos(PrimaryVertexTheta));
    if (xi_max > 0.0157) xi_max = 0.0157;
    Double_t CS1 = PrimaryCrossSection;
    Double_t CS2 = SecondaryCrossSection;
    Double_t epsilon;
     
    // Find if it's an event to the left or the right
    if (-10<=Phi && Phi<=10) {
      epsilon = hAccFunc_Left->GetBinContent(hAccFunc_Left->FindBin(Theta,Phi));
      Double_t varLeft = xi_max*CS1*CS2*epsilon - estLeft;
      varSumLeft += varLeft*varLeft;
    } else {
      epsilon = hAccFunc_Right->GetBinContent(hAccFunc_Right->FindBin(Theta,Phi));
      Double_t varRight = xi_max*CS1*CS2*epsilon - estRight;
      varSumRight += varRight*varRight;
    }
  }
  Double_t errLeft = TMath::Sqrt(varSumLeft)/nLeft;
  Double_t errRight = TMath::Sqrt(varSumRight)/nRight;    

///
////
/////////// *** this is currently Wrong as of 2015-10-06 I'm just fucking around w/ it.
///
  //Now to give the proper units and print 
  Double_t constant = (2.0/9.0)*pi*pi*(TMath::Cos(pi/36.0)-TMath::Cos(pi/18.0))
                      *N_Beam*(N_A*rho_Au/A_Au)*(N_A*rho_Au/A_Au);
  Double_t rateLeft = constant*estLeft;
  Double_t rateRight = constant*estRight;
  errLeft *= constant;
  errRight *= constant;                      
  
  std::cout << "# Left = " << nLeft << " Rate = " << rateLeft << " +- " << errLeft << std::endl;
  std::cout << "# Right = " << nRight << " Rate = " << rateRight << " +- " << errRight << std::endl; 
  
  TCanvas* c1 = new TCanvas("c1", "Histograms", 1200, 800);
  c1->Divide(3,2);
  c1->cd(1);
  hNhit_Left->Draw("CONT4");
  c1->cd(2);
  hNthrown_Left->Draw("CONT4");
  c1->cd(3);
  hAccFunc_Left->Draw("CONT4");
  c1->cd(4);
  hNhit_Right->Draw("CONT4");
  c1->cd(5);
  hNthrown_Right->Draw("CONT4");
  c1->cd(6);
  hAccFunc_Right->Draw("CONT4");
  c1->cd();
  c1->Print("AcceptanceFunction.root");

  return 0;
}

