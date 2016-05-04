//////////////////////////////////////
// DEPRICATED DEPRICATED DEPRICATED //
// DEPRICATED DEPRICATED DEPRICATED //
// DEPRICATED DEPRICATED DEPRICATED //
//////////////////////////////////////




// CrossSectionPlot.C
// Martin McHugh
// 2014-06-18
//
// For plotting CrossSections and spin functions. 

// Compile with:
//	> g++ -o CrossSectionPlot CrossSectionPlot.C `root-config --cflags --glibs`
// Run with
//	> ./CrossSectionPlot Ag 5

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

int main(Int_t argc, Char_t *argv[]) {

  Double_t PI = 3.14159265359;
  
  std::ifstream inFile;
  inFile.open("/home/mjmchugh/Mott/MottG4/CrossSections/Mott-DWBA-Au-5MeV.out");

  // Read in the first line
  std::string headers; 
  for(Int_t i=0; i<6; i++)
    inFile >> headers;
  
  // Read in all the data and determine how long the file is
  std::vector <Double_t> Theta, I, T, U, S;
  Double_t theta, i, t, u, s;
  Int_t nLines = 0;
  while(inFile.good()) {
    inFile >> theta >> i >> t >> u >> s;
    Theta.push_back(theta);
    I.push_back(i);
    T.push_back(t);
    U.push_back(u);
    S.push_back(s);
    nLines++;
  }

  inFile.close();

  Double_t ScAngle[nLines], CS[nLines], SpinT[nLines], SpinU[nLines], Sherman[nLines];
  Double_t CrossSection = 0;
  for(Int_t ii=1; ii<nLines-1; ii++) {
    ScAngle[ii] = Theta[ii];
    CS[ii] = I[ii];
    SpinT[ii] = T[ii];
    SpinU[ii] = U[ii];
    Sherman[ii] = S[ii];
    Double_t dTheta = Theta[ii+1]-Theta[ii];
    CrossSection += 2.0*PI*I[ii]*TMath::Sin(Theta[ii]*PI/180.0)*dTheta*PI/180.0;
  }
  ScAngle[nLines-1] = Theta[nLines-1];
  CS[nLines-1] = I[nLines-1];
  SpinT[nLines-1] = T[nLines-1];
  SpinU[nLines-1] = U[nLines-1];
  Sherman[nLines-1] = S[nLines-1];

  std::cout << "Total Cross Section = " << CrossSection*1.0e24 << " barn" << std::endl;  
  
  // define style here 
  // general parameters
  gStyle->SetOptDate(0);     gStyle->SetOptTitle(0);
  gStyle->SetStatColor(10);  gStyle->SetStatH(0.2);
  gStyle->SetStatW(0.3);     gStyle->SetOptStat(0); 
  
  // canvas parameters
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasColor(0);

  // pads parameters
  //  gStyle->SetPadColor(39); 
  gStyle->SetPadColor(0); 
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadBorderSize(0);
  gStyle->SetPadBottomMargin(0.18);
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadLeftMargin(0.14);
  gStyle->SetLabelSize(0.04,"x");
  gStyle->SetLabelSize(0.04,"y");
  gStyle->SetTitleXSize(0.08);
  gStyle->SetPaperSize(10,12);

  gStyle->SetTitleYOffset(0.8);
  gStyle->SetTitleYSize(0.10);
  gROOT->ForceStyle();  
  TGraph* gI = new TGraph(nLines,ScAngle,CS);
  TGraph* gT = new TGraph(nLines,ScAngle,SpinT);
  TGraph* gU = new TGraph(nLines,ScAngle,SpinU);
  TGraph* gS = new TGraph(nLines,ScAngle,Sherman);
  
  TCanvas* c1 = new TCanvas("c1","c1",800,800);
  c1->Divide(2,2);
  
  c1->cd(1);
  gI->SetLineWidth(2.0);
  gI->SetLineStyle(20);
  gI->Draw("L");
  
  c1->cd(2);
  gT->Draw("L");
  c1->cd(3);
  gU->Draw("L");
  c1->cd(4);
  gS->Draw("L");

  c1->cd();
  c1->Print("Plots.gif");
  //TFile* plotFile = new TFile("XaviDataPlots.root", "RECREATE");
  //c1->Write();
  //plotFile->Close();

  return 0;
}

Double_t CrossSection(Double_t theta) {

}


