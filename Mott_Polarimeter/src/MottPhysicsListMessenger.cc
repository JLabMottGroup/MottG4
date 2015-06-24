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
// $Id: MottPhysicsListMessenger.cc,v 1.1 2014/01/16 23:40:03 mjmchugh Exp mjmchugh $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "MottPhysicsList.hh"
#include "MottPhysicsListMessenger.hh"

#include "globals.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottPhysicsListMessenger::MottPhysicsListMessenger(MottPhysicsList* myMPGA)
{

  myPhysicsList = myMPGA;

  physicsDir = new G4UIdirectory("/PhysicsList/");
  physicsDir->SetGuidance("Adjust what physics is portrayed");

  switchOpticalPhotonsCmd = new G4UIcmdWithAnInteger("/PhysicsList/SwitchOpticalPhotons", this);
  switchOpticalPhotonsCmd->SetGuidance("Turn On (1) or off (0) Optical photon creation and propogation in the scintillators");
  switchOpticalPhotonsCmd->SetParameterName("OpticalPhotonSwitch",false);
  switchOpticalPhotonsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottPhysicsListMessenger::~MottPhysicsListMessenger() {

  delete physicsDir;
  delete switchOpticalPhotonsCmd;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MottPhysicsListMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {

  if( command == switchOpticalPhotonsCmd ) {
   
    myPhysicsList->SetOpticalPhotonSwitch(switchOpticalPhotonsCmd->GetNewIntValue(newValue));
    
    G4cout << "Set OpticalPhotonSwitch = " << newValue << G4endl; 
    
  }
  
}

