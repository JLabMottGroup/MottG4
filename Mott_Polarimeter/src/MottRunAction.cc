//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: MottRunAction.cc,v 3.9 2013/12/05 19:56:26 mjmchugh Exp $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "MottRunAction.hh"
#include "MottRunActionMessenger.hh"
#include "MottAnalysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include <sstream>
#include <stdlib.h>     /* getenv */

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottRunAction::MottRunAction()
{
  G4cout << "\tEntering MottRunAction::MottRunAction()" <<G4endl; 

  //Default File names and locations.
  G4String defaultPath = getenv("MOTTG4DIR");
  rootFileName = defaultPath + "/Mott_Polarimeter/MottSim.root";
  myMessenger = new MottRunActionMessenger(this);
  
  G4cout << "\tLeaving MottRunAction::MottRunAction()" <<G4endl; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottRunAction::~MottRunAction()
{
  G4cout << "\tEntering MottRunAction::~MottRunAction()" <<G4endl;
  
  if(myMessenger) delete myMessenger; 
  
  G4cout << "\tEntering MottRunAction::~MottRunAction()" <<G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MottRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "\tEntering MottRunAction::BeginOfRunAction()" << G4endl;

  G4int runID = aRun->GetRunID();
  std::stringstream runNo;
  runNo << runID;
  
  G4cout << "### Run " << runID << " start." << G4endl;
  
  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in MottAnalysis.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() 
         << " analysis manager" << G4endl;

  // Open an output ROOTfile
  // G4String fileName = "MottSim_" + runNo.str() + ".root";
  // For the farm
  G4String fileName = rootFileName;
  analysisManager->OpenFile(fileName);
  analysisManager->SetFirstHistoId(1);
  
  // Creating Ntuple
  analysisManager->CreateNtuple("Mott", "Edep and TrackL");
  
  // UP Detector (0-1)
  analysisManager->CreateNtupleDColumn("Up_E");
  analysisManager->CreateNtupleDColumn("Up_dE");
  //analysisManager->CreateNtupleDColumn("Up_E_dl");
  //analysisManager->CreateNtupleDColumn("Up_dE_dl");
  
  // DOWN Detector (2-3)
  analysisManager->CreateNtupleDColumn("Down_E");
  analysisManager->CreateNtupleDColumn("Down_dE");
  //analysisManager->CreateNtupleDColumn("Down_E_dl");
  //analysisManager->CreateNtupleDColumn("Down_dE_dl");
  
  // LEFT Detector (4-5)
  analysisManager->CreateNtupleDColumn("Left_E");
  analysisManager->CreateNtupleDColumn("Left_dE");
  //analysisManager->CreateNtupleDColumn("Left_E_dl");
  //analysisManager->CreateNtupleDColumn("Left_dE_dl");
  
  // RIGHT Detector (6-7)
  analysisManager->CreateNtupleDColumn("Right_E");
  analysisManager->CreateNtupleDColumn("Right_dE");
  //analysisManager->CreateNtupleDColumn("Right_E_dl");
  //analysisManager->CreateNtupleDColumn("Right_dE_dl");
  
  // "PMT" response.
  analysisManager->CreateNtupleIColumn("Up_E_PMT");	// 8
  analysisManager->CreateNtupleIColumn("Up_dE_PMT");	// 9
  analysisManager->CreateNtupleIColumn("Down_E_PMT");	// 10
  analysisManager->CreateNtupleIColumn("Down_dE_PMT");  // 11
  analysisManager->CreateNtupleIColumn("Left_E_PMT");   // 12
  analysisManager->CreateNtupleIColumn("Left_dE_PMT");  // 13
  analysisManager->CreateNtupleIColumn("Right_E_PMT");  // 14
  analysisManager->CreateNtupleIColumn("Right_dE_PMT"); // 15
  
  // Event #
  analysisManager->CreateNtupleIColumn("Event_ID");     // 16  
  
  // Primary Vertex Info
  analysisManager->CreateNtupleDColumn("PrimaryVertexKEprime");	// 17
  analysisManager->CreateNtupleDColumn("PrimaryVertexX");	// 18	// Location
  analysisManager->CreateNtupleDColumn("PrimaryVertexY");	// 19
  analysisManager->CreateNtupleDColumn("PrimaryVertexZ");	// 20
  analysisManager->CreateNtupleDColumn("PrimaryVertexTheta");	// 21
  analysisManager->CreateNtupleDColumn("PrimaryVertexPhi");	// 22
  analysisManager->CreateNtupleDColumn("PrimaryVertexPX");	// 23	// Polarization
  analysisManager->CreateNtupleDColumn("PrimaryVertexPY");	// 24
  analysisManager->CreateNtupleDColumn("PrimaryVertexPZ");	// 25

  // Primary Vertex Dynamics
  analysisManager->CreateNtupleDColumn("PrimaryCrossSection");	// 26
  analysisManager->CreateNtupleDColumn("PrimarySherman");	// 27
  analysisManager->CreateNtupleDColumn("PrimarySpinT");		// 28
  analysisManager->CreateNtupleDColumn("PrimarySpinU");		// 29

  // Secondary Vertex Info
  analysisManager->CreateNtupleDColumn("SecondaryVertexKEprime");	// 30
  analysisManager->CreateNtupleDColumn("SecondaryVertexX");		// 31	// Location
  analysisManager->CreateNtupleDColumn("SecondaryVertexY");		// 32
  analysisManager->CreateNtupleDColumn("SecondaryVertexZ");		// 33
  analysisManager->CreateNtupleDColumn("SecondaryVertexTheta");		// 34
  analysisManager->CreateNtupleDColumn("SecondaryVertexPhi");		// 35
  analysisManager->CreateNtupleDColumn("SecondaryVertexPX");		// 36	// Polarization
  analysisManager->CreateNtupleDColumn("SecondaryVertexPY");		// 37
  analysisManager->CreateNtupleDColumn("SecondaryVertexPZ");		// 38

  // Secondary Vertex Dynamics
  analysisManager->CreateNtupleDColumn("SecondaryCrossSection");	// 39
  analysisManager->CreateNtupleDColumn("SecondarySherman");		// 40
  analysisManager->CreateNtupleDColumn("SecondarySpinT");		// 41
  analysisManager->CreateNtupleDColumn("SecondarySpinU");		// 42

  // Scattering Direction
  analysisManager->CreateNtupleDColumn("Theta");		// 43
  analysisManager->CreateNtupleDColumn("Phi");			// 44

  /* Dump Plate segments.
  G4int nBeSegments = 20*6;
  G4int nCuSegments = 20*18;
  for(G4int i=0; i<nBeSegments; i++) {
    G4String name = "BeSegment";
    std::ostringstream number;
    number << i;
    name += number.str();
    analysisManager->CreateNtupleDColumn(name);
  }
  for(G4int i=0; i<nCuSegments; i++) {
    G4String name = "CuSegment";
    std::ostringstream number;
    number << i;
    name += number.str();
    analysisManager->CreateNtupleDColumn(name);
  }
  */
  
  analysisManager->FinishNtuple();
  
  G4cout << "\tLeaving MottRunAction::BeginOfRunAction()" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MottRunAction::EndOfRunAction(const G4Run* aRun)
{
  G4cout << "\tEntering MottRunAction::EndOfRunAction()" << G4endl;

  G4int runID = aRun->GetRunID();
  G4cout << "### Run " << runID << "end." << G4endl;  

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Write and Save ROOTfile.
  analysisManager->Write();
  analysisManager->CloseFile(); 

  // Complete cleanup
  delete G4AnalysisManager::Instance();
  
  G4cout << "\tLeaving MottRunAction::EndOfRunAction()" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
