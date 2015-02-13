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
// $Id: MottDetectorMessenger.cc,v 3.4 2014/01/16 23:40:03 mjmchugh Exp mjmchugh $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "MottDetectorMessenger.hh"

#include "MottDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottDetectorMessenger::MottDetectorMessenger(MottDetectorConstruction* myDet)
:myDetector(myDet)
{ 
  TargetDir = new G4UIdirectory("/Target/");
  TargetDir->SetGuidance("Set Target Foil Properties");
  
  TargetLengthCmd = new G4UIcmdWithADoubleAndUnit("/Target/SetTargetLength",this);
  TargetLengthCmd->SetGuidance("Select thickness of the target");
  TargetLengthCmd->SetParameterName("Thickness",false);
  TargetLengthCmd->SetUnitCategory("Length");
  TargetLengthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  TargetInCmd = new G4UIcmdWithoutParameter("/Target/TargetIn", this);
  TargetInCmd->SetGuidance("Places the target in the beam");
  TargetInCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  TargetOutCmd = new G4UIcmdWithoutParameter("/Target/TargetOut", this);
  TargetOutCmd->SetGuidance("Removes the target from the scattering chamber");
  TargetOutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  TargetMaterCmd = new G4UIcmdWithAString("/Target/SetTargetMaterial", this);
  TargetMaterCmd->SetGuidance("Enter the chemical symbol of the target foil");
  TargetMaterCmd->SetParameterName("Symbol",1);
  TargetMaterCmd->SetCandidates("Au Ag");
  TargetMaterCmd->SetDefaultValue("Au");
  TargetMaterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  SteppingDir = new G4UIdirectory("/Stepping/");
  SteppingDir->SetGuidance("Stepping Action settings for the world");      
        
  StepMaxCmd = new G4UIcmdWithADoubleAndUnit("/Stepping/stepMax",this); 
  StepMaxCmd->SetGuidance("Define a step max");
  StepMaxCmd->SetParameterName("stepMax",false);
  StepMaxCmd->SetUnitCategory("Length");
  StepMaxCmd->AvailableForStates(G4State_Idle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottDetectorMessenger::~MottDetectorMessenger()
{
  delete TargetLengthCmd;
  delete TargetInCmd;
  delete TargetOutCmd;
  delete StepMaxCmd;  
  delete SteppingDir;
  delete TargetDir;
  delete TargetMaterCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MottDetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == TargetLengthCmd )
   { myDetector->SetTargetLength(TargetLengthCmd->GetNewDoubleValue(newValue));}
      
  if( command == StepMaxCmd )
   { myDetector->SetMaxStep(StepMaxCmd->GetNewDoubleValue(newValue));}
   
  if( command == TargetInCmd ) 
   { myDetector->SetTargetIn(); }
   
  if( command == TargetOutCmd ) 
   { myDetector->SetTargetOut(); }

  if( command == TargetMaterCmd )
   { myDetector->SetTargetMater(newValue); }
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
