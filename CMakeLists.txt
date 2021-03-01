cmake_minimum_required(VERSION 3.10)
project(triplelens)

set(CMAKE_CXX_STANDARD 20)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")

SET(CMAKE_CXX_FLAGS_PROF
  "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} -pg -no-pie -fno-builtin"
   CACHE STRING "Flags used for cprof profiler build."
   FORCE)

set (CCC src/CCC)
set (IMAGES src/Images)
set (LENS src/LensCore)
set (LAGUERRE src/Laguerre)
set (LIGHTCURVE src/LightCurves)
set (UTILS src/Utils)


####################################
# main executable
####################################

set(DIRECTORIES ${LENS} ${UTILS} ${CCC} ${LAGUERRE} ${LIGHTCURVE})

set(SOURCE_FILES "../src/main.cc")


foreach(DIR ${DIRECTORIES})
  file(GLOB SOURCE_FILES_DIR ${DIR}/*.cc ${DIR}/*.h)
  list(APPEND SOURCE_FILES ${SOURCE_FILES_DIR})
endforeach()

list(APPEND SOURCE_FILES ${IMAGES}/imgpoint.cc
                         ${IMAGES}/imgpoint.h)

# Adding the executable
add_executable(triple_lens_test ${SOURCE_FILES})

target_include_directories(triple_lens_test PUBLIC ${DIRECTORIES} ${IMAGES})

###################################
# liblcirs
###################################

set(DIRECTORIES ${LENS} ${UTILS} ${CCC} ${LAGUERRE} ${LIGHTCURVE})

set(SOURCE_FILES "")

foreach(DIR ${DIRECTORIES})
  file(GLOB SOURCE_FILES_DIR ${DIR}/*.cc ${DIR}/*.h)
  list(APPEND SOURCE_FILES ${SOURCE_FILES_DIR})
endforeach()

# Images directory is added separately to avoid building object
# for imgpointcoeff.cc
add_library(lcirs SHARED ${SOURCE_FILES} ${IMAGES}/imgpoint.cc
                                         ${IMAGES}/imgpoint.h)
                                         

target_include_directories(lcirs PUBLIC ${DIRECTORIES} ${IMAGES})

###################################
# libimg
###################################

set(DIRECTORIES ${LENS} ${LAGUERRE})

set(SOURCE_FILES "")

foreach(DIR ${DIRECTORIES})
  file(GLOB SOURCE_FILES_DIR ${DIR}/*.cc ${DIR}/*.h)
  list(APPEND SOURCE_FILES ${SOURCE_FILES_DIR})
endforeach()

# Images directory is added separately to avoid building object
# for imgpointcoeff.cc
add_library(img SHARED ${SOURCE_FILES} ${IMAGES}/imgpoint.cc
                                       ${IMAGES}/imgpoint.h)

target_include_directories(img PUBLIC ${DIRECTORIES} ${IMAGES})

###################################
# libccc
###################################

set(DIRECTORIES ${LENS} ${LAGUERRE} ${CCC})

set(SOURCE_FILES "")

foreach(DIR ${DIRECTORIES})
  file(GLOB SOURCE_FILES_DIR ${DIR}/*.cc ${DIR}/*.h)
  list(APPEND SOURCE_FILES ${SOURCE_FILES_DIR})
endforeach()

add_library(ccc SHARED ${SOURCE_FILES} ${UTILS}/utils.h )

target_include_directories(ccc PUBLIC ${DIRECTORIES} ${UTILS})


#####################################
# lenslib
#####################################

set(HEADER_FILES src/LensCore/lens.h)

set(HEADER_DIRECTORY src/LensCore)
set(SOURCE_FILES src/LensCore/lens.cc)

add_library(lenslib STATIC 
            src/LensCore/lens.cc ${HEADER_FILES})

target_include_directories(lenslib PRIVATE ${HEADER_DIRECTORY})



