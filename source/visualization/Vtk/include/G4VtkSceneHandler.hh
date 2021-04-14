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
//
// 
// John Allison  5th April 2001
// A template for a simplest possible graphics driver.
//?? Lines or sections marked like this require specialisation for your driver.

#ifndef G4VTKSCENEHANDLER_HH
#define G4VTKSCENEHANDLER_HH

#define G4VTKDEBUG  // Comment this out to suppress debug code.

#include "G4VSceneHandler.hh"

#include "G4VtkViewer.hh"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wextra-semi"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkLine.h"
#include "vtkNamedColors.h"
#include "vtkProperty.h"
#include "vtkTextProperty.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkTextActor.h"
#include "vtkBillboardTextActor3D.h"
#include "vtkRegularPolygonSource.h"
#include "vtkGlyph2D.h"
#include "vtkGlyph3D.h"
#include "vtkMatrix4x4.h"
#pragma GCC diagnostic pop

class G4VtkSceneHandler: public G4VSceneHandler {
  friend class G4VtkViewer;

public:
  G4VtkSceneHandler(G4VGraphicsSystem& system, const G4String& name);
  virtual ~G4VtkSceneHandler();

  ////////////////////////////////////////////////////////////////
  // Required implementation of pure virtual functions...

  void AddPrimitive(const G4Polyline&);
  void AddPrimitive(const G4Text&);
  void AddPrimitive(const G4Circle&);
  void AddPrimitive(const G4Square&);
  void AddPrimitive(const G4Polyhedron&);
  // Further optional AddPrimitive methods.  Explicitly invoke base
  // class methods if not otherwise defined to avoid warnings about
  // hiding of base class methods.
  void AddPrimitive(const G4Polymarker& polymarker)
  {G4VSceneHandler::AddPrimitive (polymarker);}
  void AddPrimitive(const G4Scale& scale)
  {G4VSceneHandler::AddPrimitive (scale);}

protected:
  static G4int         fSceneIdCount;  // Counter for Vtk scene handlers.

private:
#ifdef G4VTKDEBUG
  void PrintThings();
#endif

};

#endif
