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
#include "MottAnalysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottRunAction::MottRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottRunAction::~MottRunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MottRunAction::BeginOfRunAction(const G4Run* aRun)
{

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
  G4String fileName = "/volatile/hallc/qweak/mjmchugh/Jan_2014/MottSim_" + runNo.str() + ".root";
  analysisManager->OpenFile(fileName);
  analysisManager->SetFirstHistoId(1);
  
  // Creating Ntuple
  analysisManager->CreateNtuple("Mott", "Edep and TrackL");
  
  // UP Detector (0-3)
  analysisManager->CreateNtupleDColumn("up_e");
  analysisManager->CreateNtupleDColumn("up_de");
  analysisManager->CreateNtupleDColumn("up_e_dl");
  analysisManager->CreateNtupleDColumn("up_de_dl");
  
  // DOWN Detector (4-7)
  analysisManager->CreateNtupleDColumn("down_e");
  analysisManager->CreateNtupleDColumn("down_de");
  analysisManager->CreateNtupleDColumn("down_e_dl");
  analysisManager->CreateNtupleDColumn("down_de_dl");
  
  // LEFT Detector (8-11)
  analysisManager->CreateNtupleDColumn("left_e");
  analysisManager->CreateNtupleDColumn("left_de");
  analysisManager->CreateNtupleDColumn("left_e_dl");
  analysisManager->CreateNtupleDColumn("left_de_dl");
  
  // RIGHT Detector (12-15)
  analysisManager->CreateNtupleDColumn("right_e");
  analysisManager->CreateNtupleDColumn("right_de");
  analysisManager->CreateNtupleDColumn("right_e_dl");
  analysisManager->CreateNtupleDColumn("right_de_dl");
  
  // "PMT" response.
  analysisManager->CreateNtupleIColumn("Up_E_PMT");	// 16
  analysisManager->CreateNtupleIColumn("Up_dE_PMT");	// 17
  analysisManager->CreateNtupleIColumn("Down_E_PMT");	// 18
  analysisManager->CreateNtupleIColumn("Down_dE_PMT");  // 19
  analysisManager->CreateNtupleIColumn("Left_E_PMT");   // 20
  analysisManager->CreateNtupleIColumn("Left_dE_PMT");  // 21
  analysisManager->CreateNtupleIColumn("Right_E_PMT");  // 22
  analysisManager->CreateNtupleIColumn("Right_dE_PMT"); // 23
  
  // Event #
  analysisManager->CreateNtupleIColumn("Event_ID");     // 24  
  
  // Dump segments.
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
  
  
  analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MottRunAction::EndOfRunAction(const G4Run*)
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Write and Save ROOTfile.
  analysisManager->Write();
  analysisManager->CloseFile(); 

  // Complete cleanup
  delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



