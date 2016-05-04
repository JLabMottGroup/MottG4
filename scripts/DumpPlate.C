//////////////////////////////////////
// DEPRICATED DEPRICATED DEPRICATED //
// DEPRICATED DEPRICATED DEPRICATED //
// DEPRICATED DEPRICATED DEPRICATED //
//////////////////////////////////////


// DumpPlate.C
// Martin McHugh
// Winter 2013-14
//
// For analyzing GEANT4 rootfiles on the farm machines

// Compile with:
//	> g++ -o DumpPlate DumpPlate.C `root-config --cflags --glibs`
// Run with
//	> ./DumpPlate [runnumber]


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

using namespace std;

int main(Int_t argc, Char_t *argv[]) {

  Double_t PI = 3.141592653;
  Double_t RadialStep = 5.08; 		// mm
  Double_t LayerStep = 127.0/120.0; 	// mm
  Double_t MicroAmp = 6.241e15;		// # electrons/second
  Double_t MeVtoJoule = 1.602e-13;	// J/MeV 
  
  const char* FileDir = "/lustre/expphy/volatile/hallc/qweak/mjmchugh/Jan_2014";
  const char* FileStem = "MottSim";
  
  TChain* pChain = new TChain("Mott");

  Int_t run = 0;

  if ( argc > 1 ) {
    run = atoi( argv[1] );
  } 
    
  pChain->Add(Form("%s/%s_%i.root", FileDir, FileStem, run));

  Int_t nEntries = pChain->GetEntries();
  
  Int_t nBeSegments = 120;
  Int_t nCuSegments = 360;
  
  vector < Double_t > BeDumpSegments;
  vector < Double_t > CuDumpSegments;  
  
  Double_t eDep = 0, TotalE = 0, GrandTotalE = 0;
  //Int_t nEvent = 0;
  Int_t nSegment = 0;
  
  cout << "Reading Data for " << nEntries << " Entries" << endl;
  
  for(nSegment=0; nSegment<nBeSegments; nSegment++) {
    //pChain->SetBranchAddress("Event_ID", &nEvent);
    pChain->SetBranchAddress(Form("BeSegment%i", nSegment), &eDep);
    for(Int_t i=0; i<nEntries; i++) {
      pChain->GetEntry(i);
      TotalE += eDep;
    }
    GrandTotalE += TotalE;
    Int_t Circle = nSegment % 20;
    Int_t Layer = (nSegment-Circle)/20;
    Double_t Volume = PI*RadialStep*RadialStep*(2*Circle+1)*LayerStep;
    TotalE = (TotalE)/(Volume*nEntries);
    cout << "Read Beryllium Layer " << Layer << " Circle " << Circle << ": E_seg = " << TotalE << " E_tot = " << GrandTotalE << endl;
    BeDumpSegments.push_back(TotalE);
    pChain->ResetBranchAddresses();
    TotalE = 0;
  }
  
  for(nSegment=0; nSegment<nCuSegments; nSegment++) {
    //pChain->SetBranchAddress("Event_ID", &nEvent);
    pChain->SetBranchAddress(Form("CuSegment%i", nSegment), &eDep);
    for(Int_t i=0; i<nEntries; i++) {
      pChain->GetEntry(i);
      TotalE += eDep;
    }
    GrandTotalE += TotalE;
    Int_t Circle = nSegment % 20;
    Int_t Layer = (nSegment-Circle)/20;
    Double_t Volume = PI*RadialStep*RadialStep*(2*Circle+1)*LayerStep;
    TotalE = (TotalE)/(Volume*nEntries);
    cout << "Read Copper Layer " << Layer << " Circle " << Circle << ": E_seg = " << TotalE << " E_tot = " << GrandTotalE << endl;
    CuDumpSegments.push_back(TotalE);
    pChain->ResetBranchAddresses();
    TotalE = 0;
  }
  
  cout << "Data Read" << endl;
  
  cout << endl << "======== Beryllium Data ===========" << endl << endl; 
  
  for(Int_t i=0; i<20; i++) {
    for(Int_t j=0; j<6; j++) {
      cout << BeDumpSegments[j*20+i] << ",";
    }
    cout << endl;
  }
  
  cout << endl << "======== Copper Data ===========" << endl << endl; 
  
  for(Int_t i=0; i<20; i++) {
    for(Int_t j=0; j<18; j++) {
      cout << CuDumpSegments[j*20+i] << ",";
    }
    cout << endl;
  }
  
  return 0;
}
