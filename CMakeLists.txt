cmake_minimum_required(VERSION 3.10)
project(triplelens)

set(CMAKE_CXX_STANDARD 17)

set (CCC src/CCC)
set (IMAGES src/Images)
set (LENS src/LensCore)
set (LAGUERRE src/Laguerre)
set (LIGHTCURVE src/LightCurves)
set (UTILS src/Utils)


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

#add_subdirectory(src)
#add_subdirectory(tests)
#add_subdirectory(tests/googletest)
