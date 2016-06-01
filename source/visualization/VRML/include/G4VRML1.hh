// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VRML1.hh,v 2.3 1998/11/09 19:30:45 allison Exp $
// GEANT4 tag $Name: geant4-00 $
//
// G4VRML1.hh
// Yasuhide Sawada and Satoshi Tanaka

#if defined (G4VIS_BUILD_VRML_DRIVER) || (G4VIS_USE_VRML)

#ifndef G4VRML1_HH
#define G4VRML1_HH

#include "G4VGraphicsSystem.hh"
#include "FRClient.h"

class G4VScene;

//#define FR_VRML_DEFAULT_PORT	40801
//#define FR_VRML_PORT_ENV	"FR_VRML_PORT"
//#define FR_VRML_HOST_NAME_ENV	"FR_VRML_HOST_NAME"

#include "G4VRMLNetConfig.hh"

class G4VRML1: public G4VGraphicsSystem {
public:
	G4VRML1(); 
	~G4VRML1();
	G4VScene* CreateScene(const G4String& name = "");
	G4VView*  CreateView(G4VScene&, const G4String& name = "");

public:
	const G4String& getHostName() { return fHostName; }
	G4int getPort() { return fPort; }

private:
	G4String	fHostName;
	G4int		fPort;

};

#endif //G4VRML1_HH
#endif //G4VIS_BUILD_VRML_DRIVER