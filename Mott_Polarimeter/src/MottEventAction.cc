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
// $Id: MottEventAction.cc,v 3.9 2013/12/02 21:53:35 mjmchugh Exp $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#include "MottEventAction.hh"
#include "MottAnalysis.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
MottEventAction::MottEventAction()
 : G4UserEventAction(),
   Edep(),
   dEdep(),
   EtrackL(),
   dEtrackL(),
   NumEPhotons(),
   NumdEPhotons()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
MottEventAction::~MottEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void MottEventAction::BeginOfEventAction(const G4Event*)
{
  for(G4int i=0; i<4; i++) {
    Edep[i] = 0.0;
    dEdep[i] = 0.0;
    EtrackL[i] = 0.0;
    dEtrackL[i] = 0.0;
    NumEPhotons[i] = 0; 
    NumdEPhotons[i] = 0; 
  }
  
  /*
  for(G4int i=0; i<6*20; i++) {
    BeEnergyDeposited[i] = 0;
  }
  
  for(G4int i=0; i<18*20; i++) {
    CuEnergyDeposited[i] = 0;
  }
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void MottEventAction::EndOfEventAction(const G4Event* evt)
{
  // std::cout << "\tEntering MottEventAction::EndOfEventAction()" << std::endl;
  
  G4int event_id = evt->GetEventID();

  // get number of stored trajectories
  //G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  //G4int n_trajectories = 0;
  //if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  
  /* periodic printing
  if (event_id < 100 || event_id%100 == 0) {
    G4cout << ">>> Event " << evt->GetEventID() << G4endl;
    G4cout << "    " << n_trajectories 
           << " trajectories stored in this event." << G4endl;
  }
  */
  
  // Get Analysis Manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  G4int HasBeenHit = 0;
  
  // Check Detector packages for a hit
  for(G4int i=0; i<4; i++) {
    if(Edep[i]!=0) HasBeenHit++;
    if(dEdep[i]!=0) HasBeenHit++;
  }

  /* Check Dump for a hit
  for(G4int i=0; i<6*20; i++) {
    if(BeEnergyDeposited[i]!=0) HasBeenHit++;
  }
  */

  // Fill ntuples
  if (HasBeenHit>0) {
    
    // UP
    analysisManager->FillNtupleDColumn(0, Edep[0]);
    analysisManager->FillNtupleDColumn(1, dEdep[0]);
    analysisManager->FillNtupleDColumn(2, EtrackL[0]);
    analysisManager->FillNtupleDColumn(3, dEtrackL[0]);
    
    // DOWN
    analysisManager->FillNtupleDColumn(4, Edep[1]);
    analysisManager->FillNtupleDColumn(5, dEdep[1]);
    analysisManager->FillNtupleDColumn(6, EtrackL[1]);
    analysisManager->FillNtupleDColumn(7, dEtrackL[1]);
    
    // LEFT
    analysisManager->FillNtupleDColumn(8, Edep[2]);
    analysisManager->FillNtupleDColumn(9, dEdep[2]);
    analysisManager->FillNtupleDColumn(10, EtrackL[2]);
    analysisManager->FillNtupleDColumn(11, dEtrackL[2]);
    
    // RIGHT
    analysisManager->FillNtupleDColumn(12, Edep[3]);
    analysisManager->FillNtupleDColumn(13, dEdep[3]);
    analysisManager->FillNtupleDColumn(14, EtrackL[3]);
    analysisManager->FillNtupleDColumn(15, dEtrackL[3]);
    
    // "PMT" Response
    analysisManager->FillNtupleIColumn(16, NumEPhotons[0]);
    analysisManager->FillNtupleIColumn(17, NumdEPhotons[0]);
    analysisManager->FillNtupleIColumn(18, NumEPhotons[1]);
    analysisManager->FillNtupleIColumn(19, NumdEPhotons[1]);
    analysisManager->FillNtupleIColumn(20, NumEPhotons[2]);
    analysisManager->FillNtupleIColumn(21, NumdEPhotons[2]);
    analysisManager->FillNtupleIColumn(22, NumEPhotons[3]);
    analysisManager->FillNtupleIColumn(23, NumdEPhotons[3]);
    
    // Event_ID    
    analysisManager->FillNtupleIColumn(24, event_id);
    
    /* Be Dump Plate
    for(G4int i=0; i<6*20; i++) {
      analysisManager->FillNtupleDColumn(25+i, BeEnergyDeposited[i]);
    }
    
    // Be Dump Plate
    for(G4int i=0; i<18*20; i++) {
      analysisManager->FillNtupleDColumn(145+i, CuEnergyDeposited[i]);
    }
    */
      
    analysisManager->AddNtupleRow();
  }

  // std::cout << "\tLeaving MottEventAction::EndOfEventAction()" << std::endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
