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
// $Id: MottDetectorMessenger.hh,v 3.3 2014/01/14 03:17:14 mjmchugh Exp $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef MottDetectorMessenger_h
#define MottDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class MottDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MottDetectorMessenger: public G4UImessenger
{
  public:
    MottDetectorMessenger(MottDetectorConstruction*);
   ~MottDetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    MottDetectorConstruction*   myDetector;
    
    // Mott Target
    G4UIdirectory*              TargetDir;    
    G4UIcmdWithoutParameter*	TargetInCmd;
    G4UIcmdWithoutParameter*	TargetOutCmd;
    G4UIcmdWithADoubleAndUnit*  TargetLengthCmd;
    
    // Mott Stepping stuff
    G4UIdirectory*              SteppingDir;
    G4UIcmdWithADoubleAndUnit*  StepMaxCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

