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
// $Id: MottEventAction.hh,v 3.6 2013/12/02 21:53:35 mjmchugh Exp $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#ifndef MottEventAction_h
#define MottEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"
class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MottEventAction : public G4UserEventAction
{
  public:
    MottEventAction();
   ~MottEventAction();

  public:
    void BeginOfEventAction(const G4Event*);
    void EndOfEventAction(const G4Event*);
    
    inline void SetEdep(G4double edep, G4int i) { Edep[i] = edep; }; 
    inline void SetdEdep(G4double edep, G4int i) { dEdep[i] = edep; };
    inline void SetEtrackL(G4double dl, G4int i) { EtrackL[i] = dl; }; 
    inline void SetdEtrackL(G4double dl, G4int i) { dEtrackL[i] = dl; };
    inline void SetNumEPhotons(G4int nhits, G4int i) { NumEPhotons[i] = nhits; };
    inline void SetNumdEPhotons(G4int nhits, G4int i) { NumdEPhotons[i] = nhits; };
    
    /* Dump Plate
    inline void SetBeEnergyDeposited(G4double edep, G4int i) { BeEnergyDeposited[i] = edep; };
    inline void AddBeEnergyDeposited(G4double edep, G4int i) { BeEnergyDeposited[i] += edep; };
    G4double GetBeEnergyDeposited(G4int i) { return BeEnergyDeposited[i]; };
    inline void SetCuEnergyDeposited(G4double edep, G4int i) { CuEnergyDeposited[i] = edep; };
    inline void AddCuEnergyDeposited(G4double edep, G4int i) { CuEnergyDeposited[i] += edep; };
    G4double GetCuEnergyDeposited(G4int i) { return CuEnergyDeposited[i]; };
    */
    
    G4double GetEdep(G4int i) { return Edep[i]; };
    G4double GetdEdep(G4int i) { return dEdep[i]; };
    G4double GetEtrackL(G4int i) { return EtrackL[i]; };
    G4double GetdEtrackL(G4int i) { return dEtrackL[i]; };
    G4int GetNumEPhotons(G4int i) { return NumEPhotons[i]; };    
    G4int GetNumdEPhotons(G4int i) { return NumdEPhotons[i]; };  
    G4int GetNEPE(G4int i) { return NEPE[i]; };  
    G4int GetndEPE(G4int i) { return NdEPE[i]; };

  private:
  
    // Main Detectors
    G4double Edep[4];
    G4double dEdep[4];
    G4double EtrackL[4];
    G4double dEtrackL[4];
    G4int NumEPhotons[4]; 
    G4int NumdEPhotons[4];
    G4int NEPE[4];
    G4int NdEPE[4];     

    // Dump "Detectors"
    //G4double BeEnergyDeposited[6*20];
    //G4double CuEnergyDeposited[18*20];
     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
