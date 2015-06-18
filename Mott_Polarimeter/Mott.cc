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
// $Id: Mott.cc,v 3.4 2013/05/30 18:18:23 mjmchugh Exp $
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "MottDetectorConstruction.hh"
#include "MottPhysicsList.hh"
#include "MottPrimaryGeneratorAction.hh"
#include "MottRunAction.hh"
#include "MottEventAction.hh"
#include "MottSteppingAction.hh"
#include "MottSteppingVerbose.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  // User Verbose output class
  //
  G4VSteppingVerbose* verbosity = new MottSteppingVerbose;
  G4VSteppingVerbose::SetInstance(verbosity);

  // Choose the random engine and set seed based on execution time.
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  G4long seed=time(0);
  CLHEP::HepRandom::setTheSeed(seed);
  
  // Run manager
  //
  G4RunManager * runManager = new G4RunManager;

  // User Initialization classes (mandatory)
  //
  MottDetectorConstruction* detector = new MottDetectorConstruction;
  runManager->SetUserInitialization(detector);
  //
  G4VUserPhysicsList* physics = new MottPhysicsList;
  runManager->SetUserInitialization(physics);
   
  // User Action classes
  //
  G4VUserPrimaryGeneratorAction* gen_action = new MottPrimaryGeneratorAction(detector);
  runManager->SetUserAction(gen_action);
  //
  G4UserRunAction* run_action = new MottRunAction;
  runManager->SetUserAction(run_action);
  //
  G4UserEventAction* event_action = new MottEventAction;
  runManager->SetUserAction(event_action);
  //
  G4UserSteppingAction* stepping_action = new MottSteppingAction;
  runManager->SetUserAction(stepping_action);

  // Initialize G4 kernel
  //
  runManager->Initialize();
      
#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif    
     
  // Get the pointer to the User Interface manager
  //
  G4UImanager * UImanager = G4UImanager::GetUIpointer();  

  if (argc!=1)   // batch mode  
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }
  else           // interactive mode : define UI session
    { 
#ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);
#ifdef G4VIS_USE
      UImanager->ApplyCommand("/control/execute vis.mac");     
#endif
      ui->SessionStart();
      delete ui;
#endif
     
#ifdef G4VIS_USE
      delete visManager;
#endif     
    }

  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

  delete runManager;
  delete verbosity;

  G4long time2 = time(0);
 
  G4cout << time2 << " - " << seed << " = " << time2-seed << G4endl; 

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

