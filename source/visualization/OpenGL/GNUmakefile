# $Id: GNUmakefile,v 2.3 1998/07/07 17:11:13 allison Exp $

name := G4OpenGL

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/G4VIS_BUILD.gmk
include $(G4INSTALL)/config/interactivity.gmk

CPPFLAGS += -I$(G4BASE)/visualization/management/include
CPPFLAGS += -I$(G4BASE)/visualization/modeling/include
CPPFLAGS += -I$(G4BASE)/global/management/include
CPPFLAGS += -I$(G4BASE)/global/HEPGeometry/include
CPPFLAGS += -I$(G4BASE)/geometry/solids/CSG/include
CPPFLAGS += -I$(G4BASE)/graphics_reps/include
CPPFLAGS += -I$(G4BASE)/geometry/management/include
CPPFLAGS += -I$(G4BASE)/geometry/volumes/include
CPPFLAGS += -I$(G4BASE)/intercoms/include
CPPFLAGS += -I$(G4BASE)/interfaces/common/include

ifdef G4VIS_BUILD_OPENGL_DRIVER
  CPPFLAGS += $(OGLFLAGS)
endif
ifdef G4VIS_BUILD_OPENGLXAW_DRIVER
  CPPFLAGS += $(XAWFLAGS)
endif
ifdef G4VIS_BUILD_OPENGLXM_DRIVER
  CPPFLAGS += $(XMFLAGS)
endif
ifdef G4VIS_BUILD_OPENGLX_DRIVER
  CPPFLAGS += $(X11FLAGS)
endif

include $(G4INSTALL)/config/common.gmk