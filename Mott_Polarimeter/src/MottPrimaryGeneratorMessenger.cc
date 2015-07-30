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
// $Id: MottPrimaryGeneratorMessenger.cc,v 1.1 2014/01/16 23:40:03 mjmchugh Exp mjmchugh $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "MottPrimaryGeneratorAction.hh"
#include "MottPrimaryGeneratorMessenger.hh"

#include "globals.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottPrimaryGeneratorMessenger::MottPrimaryGeneratorMessenger(MottPrimaryGeneratorAction* myMPGA)
{

  myPrimaryGeneratorAction = myMPGA;

  beamDir = new G4UIdirectory("/Beam/");
  beamDir->SetGuidance("Set various beam properties");
  
  beamEnergyCmd = new G4UIcmdWithADoubleAndUnit("/Beam/SetBeamEnergy", this);
  beamEnergyCmd->SetGuidance("Select incident beam energy");
  beamEnergyCmd->SetParameterName("BeamEnergy", false);
  beamEnergyCmd->SetUnitCategory("Energy");
  beamEnergyCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  energySpreadCmd = new G4UIcmdWithADoubleAndUnit("/Beam/SetEnergySpread", this);
  energySpreadCmd->SetGuidance("Select incident beam energy spread");
  energySpreadCmd->SetParameterName("EnergySpread", false);
  energySpreadCmd->SetUnitCategory("Energy");
  energySpreadCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  beamDiameterCmd = new G4UIcmdWithADoubleAndUnit("/Beam/SetBeamDiameter", this);
  beamDiameterCmd->SetGuidance("Select incident beam diameter");
  beamDiameterCmd->SetParameterName("BeamDiameter", false);
  beamDiameterCmd->SetUnitCategory("Length");
  beamDiameterCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  primaryDir = new G4UIdirectory("/PrimaryGenerator/");
  primaryDir->SetGuidance("Primary Generator Commands");
  
  eventTypeCmd = new G4UIcmdWithAnInteger("/PrimaryGenerator/EventType",this);
  eventTypeCmd->SetGuidance("Choose type of primary event ");
  eventTypeCmd->SetGuidance(" 0 - Throw from upstream into the target ");
  eventTypeCmd->SetGuidance(" 1 - Throw single scattered e- at the detectors (default)");
  eventTypeCmd->SetGuidance(" 2 - Throw double scattered e- at the detectors");
  eventTypeCmd->SetGuidance(" 3 - Throw single scattered e- into user specified angular range");
  eventTypeCmd->SetParameterName("EventType", false);
  eventTypeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  thetaMinCmd = new G4UIcmdWithADoubleAndUnit("/PrimaryGenerator/SetThetaMin", this);
  thetaMinCmd->SetGuidance("Select Minimum Scattering Angle");
  thetaMinCmd->SetParameterName("ThetaMin", false);
  thetaMinCmd->SetUnitCategory("Angle");
  thetaMinCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  thetaMaxCmd = new G4UIcmdWithADoubleAndUnit("/PrimaryGenerator/SetThetaMax", this);
  thetaMaxCmd->SetGuidance("Select Maximum Scattering Angle");
  thetaMaxCmd->SetParameterName("ThetaMax", false);
  thetaMaxCmd->SetUnitCategory("Angle");
  thetaMaxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  phiMinCmd = new G4UIcmdWithADoubleAndUnit("/PrimaryGenerator/SetPhiMin", this);
  phiMinCmd->SetGuidance("Select Minimum Azimuthal Angle");
  phiMinCmd->SetParameterName("PhiMin", false);
  phiMinCmd->SetUnitCategory("Angle");
  phiMinCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  phiMaxCmd = new G4UIcmdWithADoubleAndUnit("/PrimaryGenerator/SetPhiMax", this);
  phiMaxCmd->SetGuidance("Select Maximum Azimuthal Angle");
  phiMaxCmd->SetParameterName("PhiMax", false);
  phiMaxCmd->SetUnitCategory("Angle");
  phiMaxCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottPrimaryGeneratorMessenger::~MottPrimaryGeneratorMessenger() {

  delete beamDir;
  delete beamEnergyCmd;
  delete energySpreadCmd;
  delete beamDiameterCmd;

  delete primaryDir;
  delete eventTypeCmd;
  delete thetaMinCmd;
  delete thetaMaxCmd;
  delete phiMinCmd;
  delete phiMaxCmd;
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MottPrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue) {

  if( command == beamEnergyCmd ) {
   
    myPrimaryGeneratorAction->SetBeamEnergy(beamEnergyCmd->GetNewDoubleValue(newValue));
    
    G4cout << "Set Beam Energy to: " 
           << myPrimaryGeneratorAction->GetBeamEnergy()/MeV
           << " MeV" << G4endl; 
    
  }

  if( command == energySpreadCmd ) {
  
    myPrimaryGeneratorAction->SetEnergySpread(energySpreadCmd->GetNewDoubleValue(newValue));
  
    G4cout << "Set Beam Energy Spread to: " 
           << myPrimaryGeneratorAction->GetEnergySpread()/MeV 
           << " MeV" << G4endl;
           
  }
  
  if( command == beamDiameterCmd ) {
  
    myPrimaryGeneratorAction->SetBeamDiameter(beamDiameterCmd->GetNewDoubleValue(newValue));
    
    G4cout << "Set Beam Diameter (FWHM) to: " 
           << myPrimaryGeneratorAction->GetBeamDiameter()/mm 
           << " mm" << G4endl;
  
  }

  if( command == eventTypeCmd ) {
    
    myPrimaryGeneratorAction->SetEventType(eventTypeCmd->GetNewIntValue(newValue));
  
    G4cout << "Set primary EventType to : "
           << myPrimaryGeneratorAction->GetEventType() << G4endl;

  }

  if( command == thetaMinCmd ) {

    myPrimaryGeneratorAction->SetThetaMin(thetaMinCmd->GetNewDoubleValue(newValue));

    G4cout << "Set Minimum Scattering Angle [0 deg, 180 deg] to: "
           << myPrimaryGeneratorAction->GetThetaMin()/deg
           << " deg" << G4endl;

  } 

  if( command == thetaMaxCmd ) {

    myPrimaryGeneratorAction->SetThetaMax(thetaMaxCmd->GetNewDoubleValue(newValue));

    G4cout << "Set Maximum Scattering Angle [0 deg, 180 deg] to: "
           << myPrimaryGeneratorAction->GetThetaMax()/deg
           << " deg" << G4endl;

  } 

  if( command == phiMinCmd ) {

    myPrimaryGeneratorAction->SetPhiMin(phiMinCmd->GetNewDoubleValue(newValue));

    G4cout << "Set Minimum Azimuthal Angle [-180 deg, 180 deg] to: "
           << myPrimaryGeneratorAction->GetPhiMin()/deg
           << " deg" << G4endl;

  } 

  if( command == phiMaxCmd ) {

    myPrimaryGeneratorAction->SetPhiMax(phiMaxCmd->GetNewDoubleValue(newValue));

    G4cout << "Set Maximum Azimuthal Angle [-180 deg, 180 deg] to: "
           << myPrimaryGeneratorAction->GetPhiMax()/deg
           << " deg" << G4endl;

  } 


}

