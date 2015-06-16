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
// $Id: MottEventActionMessenger.cc,v 1.1 2014/01/16 23:40:03 mjmchugh Exp mjmchugh $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "MottEventAction.hh"
#include "MottEventActionMessenger.hh"

#include "globals.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottEventActionMessenger::MottEventActionMessenger(MottEventAction* myMPGA)
{

  myEventAction = myMPGA;

  eventDir = new G4UIdirectory("/EventAction/");
  eventDir->SetGuidance("Select types of events to record");

  storeAllEventsCmd = new G4UIcmdWithAnInteger("/EventAction/StoreAllEvents", this);
  storeAllEventsCmd->SetGuidance("Store all events (1) or only those which hit the detectors (0)");
  storeAllEventsCmd->SetParameterName("StoreAllFlag",false);
  storeAllEventsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottEventActionMessenger::~MottEventActionMessenger() {

  delete eventDir;
  delete storeAllEventsCmd;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MottEventActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {

  if( command == storeAllEventsCmd ) {
   
    myEventAction->SetStoreAll(storeAllEventsCmd->GetNewIntValue(newValue));
    
    G4cout << "Set StoreAll = " << newValue << G4endl; 
    
  }
  
}

