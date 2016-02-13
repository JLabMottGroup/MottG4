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
//	> ./Rate_Calculation_2 [target thickness] nFiles


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

Double_t gauss2D(Double_t *x, Double_t *par) {
  Double_t r1 = Double_t((x[0]-par[1])/par[2]);
  Double_t r2 = Double_t((x[1]-par[3])/par[4]);
  return par[0]*TMath::Exp(-0.5*(r1*r1+r2*r2));
} 

int main(Int_t argc, Char_t *argv[]) {

  const Double_t pi = 3.1415926;			// [J/MeV] 
  const Double_t N_Beam = 6.241e12;			// electrons per second per micrioamp [Hz/uA]  
  const Double_t N_A = 6.022e23;			// molcules per mol [1/mol]
  const Double_t rho_Au = 19.30;			// Gold density [g/cm^3]
  const Double_t A_Au = 196.97;				// Gold atomic weight [g/mol]
  const Double_t targetCenterZ = -(3.9878);		// Target center position [mm]
  Double_t d = 1.0e-4*atof(argv[1]);			// target_thickness [cm]
  Int_t nFiles = atoi(argv[2]);

  const char* FileDir = "/lustre/expphy/volatile/hallc/qweak/mjmchugh/Mott/Round2";
  TChain* pChain = new TChain("Mott");
  for(Int_t i=1; i<nFiles+1; i++) {
    pChain->Add(Form("%s/Double_%sum_%i.root",FileDir,argv[1],i));
    std::cout << "Added File " << i << std::endl;
  }

  Int_t nBinsTheta = 200;
  Int_t nBinsPhi = 200;
  TH2F* hNhit_Left = new TH2F("hNhit_Left","hNhit_Left",nBinsTheta,170.0,175.0,nBinsPhi,-10.0,10.0);
  TH2F* hNhit_Right = new TH2F("hNhit_Right","hNhit_Right",nBinsTheta,170.0,175.0,nBinsPhi,170.0,190.0);
  TH2F* hNthrown_Left = new TH2F("hNthrown_Left","hNthrown_Left",nBinsTheta,170.0,175.0,nBinsPhi,-10.0,10.0);
  TH2F* hNthrown_Right = new TH2F("hNthrown_Right","hNthrown_Right",nBinsTheta,170.0,175.0,nBinsPhi,170.0,190.0);
  TH2F* hAccFunc_Left = new TH2F("hAccFunc_Left","hAccFunc_Left",nBinsTheta,170.0,175.0,nBinsPhi,-10.0,10.0);
  TH2F* hAccFunc_Right = new TH2F("hAccFunc_Right","hAccFunc_Right",nBinsTheta,170.0,175.0,nBinsPhi,170.0,190.0);

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

  // First we need to simply get the acceptance functions
  for(Int_t i=0; i<nEntries; i++) {
    pChain->GetEntry(i);
    
    if(-10<=Phi&&Phi<=10) {
      hNthrown_Left->Fill(Theta,Phi);
      if( (4.6<Left_E&&Left_E<5.1) && Left_dE>0) hNhit_Left->Fill(Theta,Phi);
    } else {
      hNthrown_Right->Fill(Theta,Phi);
      if( (4.6<Right_E&&Right_E<5.1) && Right_dE>0) hNhit_Right->Fill(Theta,Phi);      
    }
    
  }
  // Epsilon is the ratio of the two
  hAccFunc_Left->Divide(hNhit_Left,hNthrown_Left);
  hAccFunc_Right->Divide(hNhit_Right,hNthrown_Right); 
  
  const Int_t npar = 5;
  Double_t fitParamsLeft[npar] = {1.0,172.5,0.4,0.0,3.0};
  Double_t fitParamsRight[npar] = {1.0,172.5,0.4,180.0,3.0};
  TF2 *fepsilonLeft = new TF2("fepsilonLeft",gauss2D,170.0,175.0,-10.0,10.0,npar);
  TF2 *fepsilonRight = new TF2("fepsilonRight",gauss2D,170.0,175.0,170.0,190.0,npar);
  fepsilonLeft->SetParameters(fitParamsLeft);
  fepsilonRight->SetParameters(fitParamsRight);
  hAccFunc_Left->Fit("fepsilonLeft");
  hAccFunc_Right->Fit("fepsilonRight");
  
  ///////////////////////////
  // Second Method: calculate the integral using a Monte Carlo estimator
  Double_t nLeft = 0;
  Double_t nRight = 0;
  Double_t sumLeft = 0;
  Double_t sumRight = 0;
  Double_t sumSquareLeft = 0;
  Double_t sumSquareRight = 0;
  for(Int_t i=0; i<nEntries; i++) {
    pChain->GetEntry(i);
            
    Double_t sigma1 = PrimaryCrossSection;
    Double_t sigma2 = SecondaryCrossSection;
    Double_t epsilon;
    
    Double_t lambda = 0.182769;					// mm
    Double_t xi = TMath::Sqrt( (SecondaryVertexX-PrimaryVertexX)*(SecondaryVertexX-PrimaryVertexX)
			      +(SecondaryVertexY-PrimaryVertexY)*(SecondaryVertexY-PrimaryVertexY)
			      +(SecondaryVertexZ-PrimaryVertexZ)*(SecondaryVertexZ-PrimaryVertexZ));

    Double_t I_xi = TMath::Exp(-xi/lambda);

    // Find if it's an event to the left or the right
    if (-10<=Phi && Phi<=10) {
      nLeft++;
      epsilon = fepsilonLeft->Eval(Theta,Phi);
      sumLeft += sigma1*sigma2*epsilon*I_xi;
      sumSquareLeft += sigma1*sigma2*epsilon*I_xi*sigma1*sigma2*epsilon*I_xi;
    } else {
      nRight++;
      epsilon = fepsilonRight->Eval(Theta,Phi);
      sumRight += sigma1*sigma2*epsilon*I_xi;
      sumSquareRight += sigma1*sigma2*epsilon*I_xi*sigma1*sigma2*epsilon*I_xi;
    }
  }    
  Double_t estLeft = sumLeft/nLeft;
  Double_t estSquareLeft = sumSquareLeft/nLeft;
  Double_t varLeft = TMath::Sqrt( (estSquareLeft - estLeft*estLeft) / nLeft);
  Double_t estRight = sumRight/nRight;
  Double_t estSquareRight = sumSquareRight/nRight;
  Double_t varRight = TMath::Sqrt( (estSquareRight - estRight*estRight) / nRight);  

  //Now to give the proper units and print 
  Double_t constant = (2.0/9.0)*pi*pi*(TMath::Cos(pi/36.0)-TMath::Cos(pi/18.0))
                      *N_Beam*(N_A*rho_Au*d/A_Au)*(N_A*rho_Au*d/A_Au);
  Double_t rateLeft = constant*estLeft;
  Double_t rateRight = constant*estRight;
  varLeft *= constant;
  varRight *= constant;                      
  
  printf(" \\hline %g & %i & %g $\\pm$ %g & %g $\\pm$ %g \\\\ \n", 
            d*10000000, nEntries, rateLeft/100.0, varLeft/100.0, rateRight/100.0, varRight/100.0);
  return 0;
}

