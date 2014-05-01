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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottSteppingAction::MottSteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MottSteppingAction::UserSteppingAction(const G4Step* theStep)
{
  // std::cout << "\tEntering MottSteppingAction::UserSteppingAction()" << std::endl;
  
  /*
  // ** all the track info is postStep **
  G4Track*              theTrack     = theStep->GetTrack();
  G4ParticleDefinition* particleType = theTrack->GetDefinition();
  G4StepPoint*          thePrePoint  = theStep->GetPreStepPoint();
  G4VPhysicalVolume*    thePrePV     = thePrePoint->GetPhysicalVolume();
  G4StepPoint*          thePostPoint = theStep->GetPostStepPoint();
  G4VPhysicalVolume*    thePostPV    = thePostPoint->GetPhysicalVolume(); 
  G4TouchableHistory*   theTouchable = (G4TouchableHistory*)(thePrePoint->GetTouchable());
  G4int                 ReplicaNo    = thePostPV->GetCopyNo();
  G4String              particleName = theTrack->GetDefinition()->GetParticleName();
  G4ProcessManager*     pm           = particleType->GetProcessManager();

  G4OpBoundaryProcessStatus boundaryStatus=Undefined;
  static G4OpBoundaryProcess* boundary=NULL;

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

  if(particleType==G4OpticalPhoton::OpticalPhotonDefinition()){
    boundaryStatus=boundary->GetStatus(); 
    if(thePostPoint->GetStepStatus()==fGeomBoundary){
      switch(boundaryStatus){
      case Absorption:
        {
          std::cout<<"Absorption"<<std::endl;
          break;
        }
      case Detection:
        {
          std::cout<<"Detected a photon."<<std::endl;
          break; 
        }
      case FresnelReflection:
        {
          std::cout<<"FresnelReflection"<<std::endl;
          break;
        }
      case TotalInternalReflection:
        {
          std::cout<<"TotalInternalReflection"<<std::endl;
          break;
        }
      case SpikeReflection:
        {
          std::cout<<"SpikeReflection"<<std::endl;
          break;
        }
      default:
        {
          std::cout<<"Undefined"<<std::endl;
          break;
        }
      }
    }
  }  // end of optical photon process

  */
  // std::cout << "\tLeaving MottSteppingAction::UserSteppingAction()" << std::endl;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

