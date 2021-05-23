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

#ifndef G4VTKVIEWER_HH
#define G4VTKVIEWER_HH

// #define G4VTKDEBUG  // Comment this out to suppress debug code.

#include "G4VViewer.hh"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wextra-semi"
#include "vtkObject.h"
#include "vtkAutoInit.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderWindow.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkInteractorStyleTerrain.h"
#include "vtkTextActor.h"
#pragma GCC diagnostic pop

// VTK_MODULE_INIT(vtkRenderingOpenGL2)
VTK_MODULE_INIT(vtkRenderingFreeType)

class vtkGeant4Callback : public vtkCommand
{
public:
  static vtkGeant4Callback* New() {return new vtkGeant4Callback;}

  vtkGeant4Callback() { fVP = nullptr; }
  void SetGeant4ViewParameters(G4ViewParameters *VP) { fVP = VP;}

  virtual void Execute(vtkObject *caller, unsigned long, void*)
  {
    vtkRenderer *ren = reinterpret_cast<vtkRenderer *>(caller);
    vtkCamera *cam = ren->GetActiveCamera();
    fVP->SetCurrentTargetPoint(G4Point3D(cam->GetFocalPoint()));
    fVP->SetViewpointDirection((G4Point3D(cam->GetPosition()) - G4Point3D(cam->GetFocalPoint())).unit());
  }

protected:
  G4ViewParameters *fVP;
};

class vtkInfoCallback : public vtkCommand
{
public:
  static vtkInfoCallback *New() { return new vtkInfoCallback; }

  vtkInfoCallback() {
    t1 = std::chrono::steady_clock::now();
    t2 = std::chrono::steady_clock::now();
  }
  void SetTextActor(vtkTextActor *txt) { this->TextActor = txt; }

  virtual void Execute(vtkObject *caller, unsigned long, void*)
  {
    vtkRenderer *ren = reinterpret_cast<vtkRenderer *>(caller);
    int      nActors = ren->GetActors()->GetNumberOfItems();
    vtkCamera   *cam = ren->GetActiveCamera();
    if(!cam) return;

    double      *pos     = cam->GetPosition();
    double viewAngle     = cam->GetViewAngle();
    double distance      = cam->GetDistance();
    double parallelScale = cam->GetParallelScale();

    if(!pos) return;

    // Get current time
    t2 = std::chrono::steady_clock::now();

    // Frame rate calculation
    std::chrono::duration<double> tdiff = t2-t1;
    t1 = t2;
    float fps = 1.0/tdiff.count();

    // String for display
    sprintf(this->TextBuff,"camera position : %.1f  %.1f %.1f \n"
                           "view angle      : %.1f\n"
                           "distance        : %.1f\n"
                           "parallel scale  : %.1f\n"
                           "number actors   : %i\n"
                           "fps             : %.1f",pos[0], pos[1], pos[2], viewAngle, distance, parallelScale, nActors, fps);
    if(this->TextActor) {
      this->TextActor->SetInput(this->TextBuff);
    }
  }
protected:
  vtkTextActor *TextActor;
  char TextBuff[256];
  std::chrono::time_point<std::chrono::steady_clock> t1;
  std::chrono::time_point<std::chrono::steady_clock> t2;
};

class G4VtkViewer: public G4VViewer {
public:
  G4VtkViewer(G4VSceneHandler &, const G4String &name);
  virtual ~G4VtkViewer();
  void SetView();
  void SetGeant4View() {};
  void ClearView();
  void DrawView();
  void ShowView();

  vtkNew<vtkTextActor>              infoTextActor;
  vtkNew<vtkInfoCallback>           infoCallback;
  vtkNew<vtkGeant4Callback>         geant4Callback;
  vtkNew<vtkLight>                  light;
  vtkNew<vtkCamera>                 camera;
  vtkNew<vtkRenderer>               renderer;
  vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
  vtkNew<vtkRenderWindow>           renderWindow;

private:

};

#endif
