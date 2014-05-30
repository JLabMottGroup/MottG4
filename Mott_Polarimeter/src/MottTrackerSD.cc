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
// $Id: MottTrackerSD.cc,v 3.8 2013/12/02 21:53:35 mjmchugh Exp $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "MottTrackerSD.hh"
#include "MottEventAction.hh"

#include "G4Material.hh"
#include "G4RunManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

#include "G4OpticalPhoton.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottTrackerSD::MottTrackerSD(G4String name)
:G4VSensitiveDetector(name)
{
  // This method is run before the beam turns on.
  G4String HCname;
  HCname = SensitiveDetectorName;
  HCname += "HitsCollection";
  //G4cout << HCname << " MottTrackerSD::MottTrackerSD(G4String name)" << G4endl;  
  collectionName.insert(HCname);
  // Somewhere in here I have to get the event action pointer created in Mott.cc 
  // and set it equal to pEventAction
  // Solved: See Below!
  pEventAction = (MottEventAction*) G4RunManager::GetRunManager()->GetUserEventAction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottTrackerSD::~MottTrackerSD(){ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MottTrackerSD::Initialize(G4HCofThisEvent* HCE)
{
  // This method runs when a run begins.    
  trackerCollection = new MottTrackerHitsCollection
                          (SensitiveDetectorName,collectionName[0]); 
  // G4cout << "MottTrackerSD::Initialize(G4HCofThisEvent* HCE)" << G4endl;                 
  // G4cout << SensitiveDetectorName << G4endl;                     
  // G4cout << collectionName[0] << G4endl;
  static G4int HCID = -1;
  if(HCID<0)
  { HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]); }
  HCE->AddHitsCollection( HCID, trackerCollection );   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool MottTrackerSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  // std::cout << "\tEntering  MottTrackerSD::ProcessHits()" << std::endl;

  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double dl   = aStep->GetStepLength();
 
  G4StepPoint* point = aStep->GetPreStepPoint();
  G4int copy = point->GetTouchableHandle()->GetCopyNumber();
  
  if(edep==0.) return false;
  
  MottTrackerHit* newHit = new MottTrackerHit();
  newHit->SetCopyNo(copy);
  newHit->SetTrackID(aStep->GetTrack()->GetTrackID());
  newHit->SetEdep(edep);
  newHit->SetTrackLength(dl);
  newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());
  trackerCollection->insert(newHit);
  
  //newHit->Print();
  //newHit->Draw();

  // std::cout << "\tLeaving MottTrackerSD::ProcessHits()" << std::endl;

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MottTrackerSD::EndOfEvent(G4HCofThisEvent*)
{
  // std::cout << "\tEntering MottTrackerSD::EndOfEvent()" << std::endl;

  G4int copy = 0;
  G4int NbHits = trackerCollection->entries();
  G4double totalEdep = 0.0;
  G4double totalLength = 0.0;

  G4int ngamma = GetPhotonCount();
  
  for (G4int i=0; i<NbHits; i++){
    //(*trackerCollection)[i]->Print();
    totalEdep += (*trackerCollection)[i]->GetEdep();
    totalLength += (*trackerCollection)[i]->GetTrackLength();
  }
    
  if (verboseLevel>0) { 
    G4cout << "\n<<<<< Detector " << trackerCollection->GetSDname() << " >>>>>"
           << "\nCopy Number: " << copy
           << "\nNumber of Hits: " << NbHits
           << "\nTotal Energy Deposited: " << G4BestUnit(totalEdep, "Energy")
           << "\nTrack Length: " << G4BestUnit(totalLength, "Length") << "\n" << G4endl;
  }

  // UP
  if ( trackerCollection->GetSDname() == "dE_1" ) { 
    pEventAction->SetdEdep(totalEdep, 0);
    pEventAction->SetdEtrackL(totalLength, 0);
    pEventAction->SetNumdEPhotons(ngamma, 0);
  }
  if ( trackerCollection->GetSDname() == "E_1" ) { 
    pEventAction->SetEdep(totalEdep, 0);
    pEventAction->SetEtrackL(totalLength, 0);
    pEventAction->SetNumEPhotons(ngamma, 0);
  }  
  
  // DOWN
  if ( trackerCollection->GetSDname() == "dE_2" ) { 
    pEventAction->SetdEdep(totalEdep, 1);
    pEventAction->SetdEtrackL(totalLength, 1);
    pEventAction->SetNumdEPhotons(ngamma, 1);
  }
  if ( trackerCollection->GetSDname() == "E_2" ) { 
    pEventAction->SetEdep(totalEdep, 1);
    pEventAction->SetEtrackL(totalLength, 1);
    pEventAction->SetNumEPhotons(ngamma, 1);
  }
  
  // LEFT
  if ( trackerCollection->GetSDname() == "dE_3" ) { 
    pEventAction->SetdEdep(totalEdep, 2);
    pEventAction->SetdEtrackL(totalLength, 2);
    pEventAction->SetNumdEPhotons(ngamma, 2);
  }
  if ( trackerCollection->GetSDname() == "E_3" ) { 
    pEventAction->SetEdep(totalEdep, 2);
    pEventAction->SetEtrackL(totalLength, 2);    
    pEventAction->SetNumEPhotons(ngamma, 2);
  }  
  
  // RIGHT
  if ( trackerCollection->GetSDname() == "dE_4" ) { 
    pEventAction->SetdEdep(totalEdep, 3);
    pEventAction->SetdEtrackL(totalLength, 3);
    pEventAction->SetNumdEPhotons(ngamma, 3);
  }
  if ( trackerCollection->GetSDname() == "E_4" ) { 
    pEventAction->SetEdep(totalEdep, 3);
    pEventAction->SetEtrackL(totalLength, 3);
    pEventAction->SetNumEPhotons(ngamma, 3);
  } 
  
  if ( trackerCollection->GetSDname() == "BeDump" ) {
    for (G4int i=0; i<NbHits; i++) {
      copy = (*trackerCollection)[i]->GetCopyNo();
      G4double edep = (*trackerCollection)[i]->GetEdep();
      pEventAction->AddBeEnergyDeposited(edep, copy);
    }  
  }
  
  if ( trackerCollection->GetSDname() == "CuDump" ) {
    for (G4int i=0; i<NbHits; i++) {
      copy = (*trackerCollection)[i]->GetCopyNo();
      G4double edep = (*trackerCollection)[i]->GetEdep();
      pEventAction->AddCuEnergyDeposited(edep, copy);
    }  
  }
  
  SetPhotonCount(0);
  
  // std::cout << "\tLeaving MottTrackerSD::EndOfEvent()" << std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

