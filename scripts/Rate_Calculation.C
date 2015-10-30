// 
// Rate_Calculation.C
// Martin McHugh
// 2015-08-10
//
// Performing Monte Carlo Integration of the detector rates using simulated 
// single-scattering events.

// Compile with:
//	> g++ -o Rate_Calculation Rate_Calculation.C `root-config --cflags --glibs`
// Run with
//	> ./Rate_Calculation [target thickness]


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

  const char* FileDir = "/home/mjmchugh/Mott/MottG4/Mott_Polarimeter/build";
  TChain* pChain = new TChain("Mott");
  //pChain->Add(Form("%s/Single_%sum.root", FileDir,argv[1]));
  //for(Int_t i=0; i<atoi(argv[2]); i++) pChain->Add(Form("%s/Single_%i.root", FileDir,i));
  pChain->Add(Form("%s/Single_%s.root", FileDir,argv[1]));
  
  Int_t nBinsTheta = 100;
  Int_t nBinsPhi = 100;
  TProfile2D* pCS_Left = new TProfile2D("pCS_Left","pCS_Left",nBinsTheta,170.0,175.0,nBinsPhi,-10.0,10.0);
  TProfile2D* pCS_Right = new TProfile2D("pCS_Right","pCS_Right",nBinsTheta,170.0,175.0,nBinsPhi,170.0,190.0);
  TH2F* hNhit_Left = new TH2F("hNhit_Left","hNhit_Left",nBinsTheta,170.0,175.0,nBinsPhi,-10.0,10.0);
  TH2F* hNhit_Right = new TH2F("hNhit_Right","hNhit_Right",nBinsTheta,170.0,175.0,nBinsPhi,170.0,190.0);
  TH2F* hNthrown_Left = new TH2F("hNthrown_Left","hNthrown_Left",nBinsTheta,170.0,175.0,nBinsPhi,-10.0,10.0);
  TH2F* hNthrown_Right = new TH2F("hNthrown_Right","hNthrown_Right",nBinsTheta,170.0,175.0,nBinsPhi,170.0,190.0);
  TH2F* hAccFunc_Left = new TH2F("hAccFunc_Left","hAccFunc_Left",nBinsTheta,170.0,175.0,nBinsPhi,-10.0,10.0);
  TH2F* hAccFunc_Right = new TH2F("hAccFunc_Right","hAccFunc_Right",nBinsTheta,170.0,175.0,nBinsPhi,170.0,190.0);

  Double_t Left_E, Left_dE;
  Double_t Right_E, Right_dE;
  Double_t PrimaryCrossSection;
  Double_t Theta;
  Double_t Phi;
  Int_t nEntries = pChain->GetEntries();

  pChain->SetBranchAddress("Left_E",&Left_E);
  pChain->SetBranchAddress("Left_dE",&Left_dE);
  pChain->SetBranchAddress("Right_E",&Right_E);
  pChain->SetBranchAddress("Right_dE",&Right_dE);
  pChain->SetBranchAddress("PrimaryCrossSection",&PrimaryCrossSection);
  pChain->SetBranchAddress("PrimaryVertexTheta",&Theta);
  pChain->SetBranchAddress("PrimaryVertexPhi",&Phi);

  // First we need to simply get the acceptance functions
  // and fill the cross-section diagram
  for(Int_t i=0; i<nEntries; i++) {
    pChain->GetEntry(i);
    
    if(-10<=Phi&&Phi<=10) {
      pCS_Left->Fill(Theta,Phi,PrimaryCrossSection);
      hNthrown_Left->Fill(Theta,Phi);
      if(Left_E>0 && Left_dE>0) hNhit_Left->Fill(Theta,Phi);
    } else {
      pCS_Right->Fill(Theta,Phi,PrimaryCrossSection);
      hNthrown_Right->Fill(Theta,Phi);
      if(Right_E>0 && Right_dE>0) hNhit_Right->Fill(Theta,Phi);      
    }
    
  }
  // Epsilon is the ratio of the two
  hAccFunc_Left->Divide(hNhit_Left,hNthrown_Left);
  hAccFunc_Right->Divide(hNhit_Right,hNthrown_Right);

  TCanvas* c1 = new TCanvas("c1","c1",800,400);
  c1->Divide(2,1);
  c1->cd(1);
  pCS_Left->GetXaxis()->SetTitle("#chi [deg]");
  pCS_Left->GetYaxis()->SetTitle("#psi [deg]");
  pCS_Left->Draw("CONT4");
  c1->cd(2);
  hAccFunc_Left->GetXaxis()->SetTitle("#chi [deg]");
  hAccFunc_Left->GetYaxis()->SetTitle("#psi [deg]");
  hAccFunc_Left->Draw("CONT4");
  c1->cd();
  c1->Print("CSandEpsilon.root");
  delete c1;
  
  const Int_t npar = 5;
  Double_t fitParamsLeft[npar] = {1.0,172.5,0.4,0.0,3.0};
  Double_t fitParamsRight[npar] = {1.0,172.5,0.4,180.0,3.0};
  TF2 *fepsilonLeft = new TF2("fepsilonLeft",gauss2D,170.0,175.0,-10.0,10.0,npar);
  TF2 *fepsilonRight = new TF2("fepsilonRight",gauss2D,170.0,175.0,170.0,190.0,npar);
  fepsilonLeft->SetParameters(fitParamsLeft);
  fepsilonRight->SetParameters(fitParamsRight);
  hAccFunc_Left->Fit("fepsilonLeft");
  hAccFunc_Right->Fit("fepsilonRight");

  ////////////////////////////////
  // First Method ( Qweak Method ) integrate over theta and phi with averaged cross sections
  Double_t Sum_L = 0, Sum_R = 0;
  Double_t dSum_L = 0, dSum_R = 0;
  Double_t X_L = 0, X_R = 0;
  Double_t dX_L = 0, dX_R = 0;
  Double_t Theta_min = 170.0*3.1515926/180.0;
  Double_t Delta_theta = (5.0/nBinsTheta)*(3.1515926/180.0);
  Double_t Delta_phi = 4.0*Delta_theta;
  for(Int_t theta_i=1; theta_i<nBinsTheta+1; theta_i++) {
    for(Int_t phi_j=1; phi_j<nBinsTheta+1; phi_j++) {
      Double_t s_i = TMath::Sin(Theta_min+theta_i*Delta_theta);
      Double_t sigma_L = pCS_Left->GetBinContent(theta_i,phi_j);
      Double_t sigma_R = pCS_Right->GetBinContent(theta_i,phi_j);
      Double_t epsilon_L = hAccFunc_Left->GetBinContent(theta_i,phi_j);
      Double_t epsilon_R = hAccFunc_Right->GetBinContent(theta_i,phi_j);
      Double_t d_sigma_L = pCS_Left->GetBinError(theta_i,phi_j);
      Double_t d_sigma_R = pCS_Right->GetBinError(theta_i,phi_j);
      Double_t d_epsilon_L = hAccFunc_Left->GetBinError(theta_i,phi_j);
      Double_t d_epsilon_R = hAccFunc_Right->GetBinError(theta_i,phi_j);

      X_L = sigma_L*epsilon_L*s_i;
      X_R = sigma_R*epsilon_R*s_i;      
      
      if(sigma_L==0 || epsilon_L==0) {
        dX_L = 0;
      } else {
        dX_L = (d_sigma_L/sigma_L)*(d_sigma_L/sigma_L)
                + (d_epsilon_L/epsilon_L)*(d_epsilon_L/epsilon_L);	// (dX/X)^2 
        dX_L = X_L*X_L*dX_L;	// dX^2
      }
      
      if(sigma_R==0 || epsilon_R==0) { 
        dX_R = 0;
      } else {
        dX_R = (d_sigma_R/sigma_R)*(d_sigma_R/sigma_R) 
                + (d_epsilon_R/epsilon_R)*(d_epsilon_R/epsilon_R);	// (dX/X)^2 
        dX_R = X_R*X_R*dX_R;	// dX^2
      } 
       
      //std::cout << pCS->GetBinContent(theta_i,phi_j) << " " << pCS->GetBinError(theta_i,phi_j) << std::endl;
      Sum_L += X_L;
      Sum_R += X_R;
      dSum_L += dX_L;
      dSum_R += dX_R;
    }
  }
  Sum_L *= (Delta_theta*Delta_phi);
  Sum_R *= (Delta_theta*Delta_phi);
  dSum_L = TMath::Sqrt(dSum_L);
  dSum_L *= (Delta_theta*Delta_phi);
  dSum_R = TMath::Sqrt(dSum_R);
  dSum_R *= (Delta_theta*Delta_phi);

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
            
    Double_t sigma = PrimaryCrossSection;
    Double_t epsilon;
     
    // Find if it's an event to the left or the right
    if (-10<=Phi && Phi<=10) {
      nLeft++;
      epsilon = fepsilonLeft->Eval(Theta,Phi);
      sumLeft += sigma*epsilon;
      sumSquareLeft += sigma*epsilon*sigma*epsilon;
    } else {
      nRight++;
      epsilon = fepsilonRight->Eval(Theta,Phi);
      sumRight += sigma*epsilon;
      sumSquareRight += sigma*epsilon*sigma*epsilon;
    }
  }    
  Double_t estLeft = sumLeft/nLeft;
  Double_t estSquareLeft = sumSquareLeft/nLeft;
  Double_t varLeft = TMath::Sqrt( (estSquareLeft - estLeft*estLeft) / nLeft);
  Double_t estRight = sumRight/nRight;
  Double_t estSquareRight = sumSquareRight/nRight;
  Double_t varRight = TMath::Sqrt( (estSquareRight - estRight*estRight) / nRight);
  
  ////////////////////
  // Results        //
  ////////////////////
  Double_t luminosity = N_Beam*N_A*rho_Au*d/A_Au;
  Sum_L *= luminosity;	dSum_L *= luminosity;
  Sum_R *= luminosity;	dSum_R *= luminosity;
  Double_t dR = TMath::Sqrt(dSum_L*dSum_L + dSum_R+dSum_R)/2.0;

  Double_t constant = luminosity*(pi/9.0)*(TMath::Cos(pi/36.0)-TMath::Cos(pi/18.0));
  estLeft *= constant;	varLeft *= constant;
  estRight *= constant;	varRight *= constant;  
 
  std::cout << "Left Rate1 = " << Sum_L << " +- " << dSum_L << std::endl;  
  std::cout << "Right Rate1 = " << Sum_R << " +- " << dSum_R << std::endl;
  std::cout << "Avg1 = " << (Sum_L + Sum_R)/2.0 << " +- " << dR << std::endl;
  std::cout << "Left Rate2 = " << estLeft << " +- " << varLeft << std::endl;  
  std::cout << "Right Rate2 = " << estRight << " +- " << varRight << std::endl;
  
  printf("  \\\\hline %g & %g $\\pm$ %g & %g $\\pm$ %g & %g $\\pm$ %g &  %g $\\pm$ %g &\\\\ \n", 
            d, Sum_L, dSum_L, Sum_R, dSum_R, estLeft, varLeft, estRight, varRight);
  //std::cout << Sum_L << " & " << Sum_R << " & " << estLeft << " & " << estRight << std::endl;

  return 0;
}

