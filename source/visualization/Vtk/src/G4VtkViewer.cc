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

#include "G4VtkViewer.hh"

#include "G4VSceneHandler.hh"
#include "G4VtkSceneHandler.hh"

G4VtkViewer::G4VtkViewer (G4VSceneHandler& sceneHandler, const G4String& name) :
G4VViewer(sceneHandler, sceneHandler.IncrementViewCount(), name) {

  G4cout << "G4VtkViewer::G4VtkViewer" << G4endl;
  renderWindow->SetSize(900, 900);
  renderWindow->SetWindowName("Vtk viewer");

  //renderer->SetBackground(1,1,1);
  renderWindow->AddRenderer(renderer);
  renderWindowInteractor->SetRenderWindow(renderWindow);
}

G4VtkViewer::~G4VtkViewer() {}

void G4VtkViewer::SetView() {
#ifdef G4VTKDEBUG
  G4cout << "G4VtkViewer::SetView() called." << G4endl;
#endif
}

void G4VtkViewer::ClearView() {
#ifdef G4VTKDEBUG
  G4cout << "G4VtkViewer::ClearView() called." << G4endl;
#endif
}

void G4VtkViewer::DrawView() {
#ifdef G4VTKDEBUG
  G4cout << "G4VtkViewer::DrawView() called." << G4endl;
#endif

  // First, a view should decide when to re-visit the G4 kernel.
  // Sometimes it might not be necessary, e.g., if the scene is stored
  // in a graphical database (e.g., OpenGL's display lists) and only
  // the viewing angle has changed.  But graphics systems without a
  // graphical database will always need to visit the G4 kernel.

  NeedKernelVisit();  // Default is - always visit G4 kernel.
  // Note: this routine sets the fNeedKernelVisit flag of *all* the
  // views of the scene.

  ProcessView();      // The basic logic is here.

  // Then a view may have more to do, e.g., display the graphical
  // database.  That code should come here...

  // ...before finally...
  FinishView();       // Flush streams and/or swap buffers.
}

void G4VtkViewer::ShowView() {
#ifdef G4VTKDEBUG
  G4cout << "G4VtkViewer::ShowView() called." << G4endl;
  // static_cast<G4VtkSceneHandler&>(fSceneHandler).PrintStores();
#endif

  renderWindow->Render();
  renderWindowInteractor->Start();
}
