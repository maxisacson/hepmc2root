cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(hepmc2root)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(HEPMC_PATH "/usr/local" CACHE PATH "Path to HepMC2 installation")
set(ROOT_CONFIG "root-config" CACHE FILEPATH "Path to root-config")

execute_process(COMMAND ${ROOT_CONFIG} --prefix OUTPUT_VARIABLE ROOT_PATH)
list(APPEND CMAKE_PREFIX_PATH ${ROOT_PATH})
find_package(ROOT)

include(${ROOT_USE_FILE})
# ROOT_GENERATE_DICTIONARY(Dict LINKDEF LinkDef.h)
REFLEX_GENERATE_DICTIONARY(G__Classes classes.h SELECTION classes.xml)
add_library(G__ClassesDict SHARED G__Classes.cxx)
target_link_libraries(G__ClassesDict ${ROOT_LIBRARIES})

include_directories(${PROJECT_SOURCE_DIR}/include ${HEPMC_PATH}/include ${ROOT_INCLUDE_DIRS})

add_executable(hepmc2root hepmc2root.cxx G__Classes.cxx)
target_link_libraries(hepmc2root ${HEPMC_PATH}/lib/libHepMC.so ${ROOT_LIBRARIES})
