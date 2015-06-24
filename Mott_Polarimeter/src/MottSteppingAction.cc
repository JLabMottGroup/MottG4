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
// $Id: MottSteppingAction.cc,v 3.2 2013/03/12 16:34:43 mjmchugh Exp $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "MottSteppingAction.hh"
#include "MottTrackerSD.hh"

#include "G4SteppingManager.hh"
#include "G4TrackVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4OpticalPhoton.hh"
#include "G4PionMinus.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottSteppingAction::MottSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MottSteppingAction::UserSteppingAction(const G4Step* theStep)
{
  // G4cout << "\tEntering MottSteppingAction::UserSteppingAction()" << G4endl;
  
  // ** all the track info is postStep **
  G4Track*              theTrack     = theStep->GetTrack();
  
  if(!theTrack->GetNextVolume()) return;	// Stop us from segfaulting when OutOfWorld
  
  G4ParticleDefinition* particleType = theTrack->GetDefinition();
  G4StepPoint*          thePostPoint = theStep->GetPostStepPoint();
  G4VPhysicalVolume*    thePostPV    = thePostPoint->GetPhysicalVolume();
  G4String              particleName = theTrack->GetDefinition()->GetParticleName();
  G4ProcessManager*     pm           = particleType->GetProcessManager();
  
  G4OpBoundaryProcessStatus boundaryStatus=Undefined;
  G4OpBoundaryProcess* boundary=NULL;
  
  if(!boundary){
    G4int nprocesses = pm->GetProcessListLength();
    G4ProcessVector* pv = pm->GetProcessList();
    for(G4int i=0; i<nprocesses; i++){
      if((*pv)[i]->GetProcessName()=="OpBoundary"){
        boundary = (G4OpBoundaryProcess*)(*pv)[i];
        break;
      }
    }
  }
  
  if(particleName=="opticalphoton"){
    boundaryStatus=boundary->GetStatus(); 
    if(thePostPoint->GetStepStatus()==fGeomBoundary) {
      switch(boundaryStatus) {
        case Absorption: {
          // G4cout<<"Absorption"<<G4endl;
          break;
        } case Detection: {
          G4SDManager* SDMan = G4SDManager::GetSDMpointer();
          //G4cout << "Detected a photon in "<< thePostPV->GetName() << G4endl;
          if (thePostPV->GetName() == "physiEPhotoUp") {				// E Up
            MottTrackerSD* mottSD = (MottTrackerSD*)SDMan->FindSensitiveDetector("/Mott/E_1");
            mottSD->IncrementPhotonCount();
            theTrack->SetTrackStatus(fStopAndKill);
          }
          if (thePostPV->GetName() == "physidEPhotoUp") {				// dE Up
            MottTrackerSD* mottSD = (MottTrackerSD*)SDMan->FindSensitiveDetector("/Mott/dE_1");
            mottSD->IncrementPhotonCount();
            theTrack->SetTrackStatus(fStopAndKill);
          }
          if (thePostPV->GetName() == "physiEPhotoDown") {				// E Down
            MottTrackerSD* mottSD = (MottTrackerSD*)SDMan->FindSensitiveDetector("/Mott/E_2");
            mottSD->IncrementPhotonCount();
            theTrack->SetTrackStatus(fStopAndKill);
          }
          if (thePostPV->GetName() == "physidEPhotoDown") {				// dE Down
            MottTrackerSD* mottSD = (MottTrackerSD*)SDMan->FindSensitiveDetector("/Mott/dE_2");
            mottSD->IncrementPhotonCount();
            theTrack->SetTrackStatus(fStopAndKill);
          }          
          if (thePostPV->GetName() == "physiEPhotoLeft") {				// E Left
            MottTrackerSD* mottSD = (MottTrackerSD*)SDMan->FindSensitiveDetector("/Mott/E_3");
            mottSD->IncrementPhotonCount();
            theTrack->SetTrackStatus(fStopAndKill);
          }
          if (thePostPV->GetName() == "physidEPhotoLeft") {				// dE Left
            MottTrackerSD* mottSD = (MottTrackerSD*)SDMan->FindSensitiveDetector("/Mott/dE_3");
            mottSD->IncrementPhotonCount();
            theTrack->SetTrackStatus(fStopAndKill);
          }
          if (thePostPV->GetName() == "physiEPhotoRight") {				// E Right
            MottTrackerSD* mottSD = (MottTrackerSD*)SDMan->FindSensitiveDetector("/Mott/E_4");
            mottSD->IncrementPhotonCount();
            theTrack->SetTrackStatus(fStopAndKill);
          }
          if (thePostPV->GetName() == "physidEPhotoRight") {				// dE Right
            MottTrackerSD* mottSD = (MottTrackerSD*)SDMan->FindSensitiveDetector("/Mott/dE_4");
            mottSD->IncrementPhotonCount();
            theTrack->SetTrackStatus(fStopAndKill);
          }          
          break; 
        } case FresnelReflection: {
          // G4cout<<"FresnelReflection"<<G4endl;
          break;
        } case TotalInternalReflection: {
          //G4cout<<"TotalInternalReflection"<<G4endl;
          break;
        } case SpikeReflection: {
          //G4cout<<"SpikeReflection"<<G4endl;
          break;
        } default: {
          //G4cout<<"Undefined"<<G4endl;
          break;
        }
      } // end switch(boundaryStatus)
    } // end (if fGeomBoundary)
  } // end of optical photon process

  // G4cout << "\tLeaving MottSteppingAction::UserSteppingAction()" << G4endl;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

