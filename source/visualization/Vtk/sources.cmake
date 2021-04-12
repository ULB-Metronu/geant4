#------------------------------------------------------------------------------
# sources.cmake
# Module : G4visQt3D
# Package: Geant4.src.G4visualization.G4visQt3D
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
#
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
GEANT4_DEFINE_MODULE(NAME G4visVtk
    HEADERS
        G4Vtk.hh
        G4VtkSceneHandler.hh
        G4VtkViewer.hh
    SOURCES
        G4Vtk.cc
        G4VtkSceneHandler.cc
        G4VtkViewer.cc
    GRANULAR_DEPENDENCIES
        G4csg
        G4geometrymng
        G4globman
        G4graphics_reps
        G4hits
        G4intercoms
        G4interfaces
        G4modeling
        G4specsolids
        G4tracking
        G4vis_management
    GLOBAL_DEPENDENCIES
        G4digits_hits
        G4geometry
        G4global
        G4graphics_reps
        G4intercoms
        G4interfaces
        G4modeling
        G4tracking
        G4vis_management
    LINK_LIBRARIES
        /opt/local/lib/libvtkCommonCore-8.2.dylib
        /opt/local/lib/libvtkCommonDataModel-8.2.dylib
        /opt/local/lib/libvtkRenderingCore-8.2.dylib
        /opt/local/lib/libvtkRenderingOpenGL2-8.2.dylib
        )

# List any source specific properties here
geant4_module_include_directories(G4visVtk PUBLIC /opt/local/include/vtk-8.2/)

