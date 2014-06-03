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
// MottRunActionMessenger.cc
// Created 2014-04-30
// Martin McHugh
// mjmchugh@jlab.org
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "MottRunAction.hh"
#include "MottRunActionMessenger.hh"

#include "globals.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottRunActionMessenger::MottRunActionMessenger(MottRunAction* runAction)
:myRunAction(runAction)
{ 
  G4cout << "\tEntering MottRunActionMessenger::MottRunActionMessenger()" << G4endl;

  runActionDir = new G4UIdirectory("/Analysis/");
  runActionDir->SetGuidance("Name output files and other run things");
  
  rootFileNameCmd = new G4UIcmdWithAString("/Analysis/RootFileName",this);
  rootFileNameCmd->SetGuidance("Set Name of file with output ROOT tree");
  rootFileNameCmd->SetParameterName("RootFileName",false);
  rootFileNameCmd->SetDefaultValue("MottSim.root");
  rootFileNameCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  G4cout << "\tLeaving MottRunActionMessenger::MottRunActionMessenger()" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottRunActionMessenger::~MottRunActionMessenger()
{
  G4cout << "\tEntering MottRunActionMessenger::~MottRunActionMessenger()" << G4endl;
  
  if(rootFileNameCmd) delete rootFileNameCmd;
  if(runActionDir) delete runActionDir;
  
  G4cout << "\tLeaving MottRunActionMessenger::~MottRunActionMessenger()" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MottRunActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if(command == rootFileNameCmd) myRunAction->SetRootFileName(newValue);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
