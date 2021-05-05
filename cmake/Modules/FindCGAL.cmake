set(CGAL_PREFIX "" CACHE PATH "path ")

find_path(CGAL_INCLUDE_DIR CGAL
    PATHS ${CGAL_PREFIX}/include /usr/include /usr/local/include )

if(CGAL_INCLUDE_DIR)
    set(CGAL_FOUND TRUE)
endif()

if(CGAL_FOUND)
   set(CGAL_INCLUDE_DIRS ${CGAL_INCLUDE_DIR})

   if(NOT CGAL_FIND_QUIETLY)
      MESSAGE(STATUS "Found CGAL: ${CGAL_INCLUDE_DIR}")
   endif()
elseif(CGAL_FOUND)
   if(CGAL_FIND_REQUIRED)
      message(FATAL_ERROR "Could not find CGAL")
   endif()
endif()
