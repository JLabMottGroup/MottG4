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
// $Id: MottDumpParameterisation.cc,v 1.1 2013/12/02 21:53:35 mjmchugh Exp $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "MottDumpParameterisation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottDumpParameterisation::MottDumpParameterisation(G4int NoLayers,
        					   G4double startZ,
        					   G4double spacingZ,
        					   G4double DepthChamber,
        					   G4int NoRings,
        					   G4double radialStep)
{

  fNoLayers   =  NoLayers;
  fNoRings    =  NoRings;
  fStartZ     =  startZ; 
  fHalfDepth  =  DepthChamber*0.5;
  fSpacing    =  spacingZ;
  fRadialStep =  radialStep;
  
  if( NoLayers > 0 ){
    if (spacingZ < DepthChamber) {
      G4Exception("MottDumpParameterisation::MottDumpParameterisation()",
                  "InvalidSetup", FatalException,
                  "Depth>Spacing");
    }
  }
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MottDumpParameterisation::~MottDumpParameterisation() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MottDumpParameterisation::ComputeTransformation
  (const G4int copyNo, G4VPhysicalVolume* physVol) const
{

  G4int Ring = copyNo%fNoRings;
  G4int Layer = (copyNo-Ring)/fNoRings;
  
  G4double Zposition = fStartZ + Layer*fSpacing;
  G4ThreeVector origin(0,0,Zposition);
  
  physVol->SetTranslation(origin);
  physVol->SetRotation(0);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MottDumpParameterisation::ComputeDimensions 
  (G4Tubs& solidVol, const G4int copyNo, const G4VPhysicalVolume*) const 
{

  G4int Ring = copyNo%fNoRings;

  G4double IR = Ring*fRadialStep;
  G4double OR = (Ring+1)*fRadialStep;
  
  solidVol.SetInnerRadius(IR);
  solidVol.SetOuterRadius(OR);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
