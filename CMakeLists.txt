cmake_minimum_required (VERSION 2.8)
set(PROYECTO "DWGReverb")
project (${PROYECTO})

#include_directories("C:/SupercolliderRepos/pluginfork/sc3-plugins/source/DWGUGens/dwglib")
include_directories(${SC_PATH}/include/plugin_interface)
include_directories(${SC_PATH}/include/common)
include_directories(${SC_PATH}/common)
include_directories(${SC_PATH}/external_libraries/libsndfile/)


set(CMAKE_SHARED_MODULE_PREFIX "")
if(APPLE OR WIN32)
set(CMAKE_SHARED_MODULE_SUFFIX ".scx")
endif()
	
if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse2 -mstackrealign -mfpmath=sse")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2 -mstackrealign -mfpmath=sse")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()
aux_source_directory(. src)
add_library(${PROYECTO} MODULE ${src})

add_library(${PROYECTO}_supernova MODULE ${src})
set_property(TARGET ${PROYECTO}_supernova APPEND PROPERTY COMPILE_DEFINITIONS  "SUPERNOVA")
include_directories(${SC_PATH}/external_libraries
					${SC_PATH}/external_libraries/boost-lockfree
					${SC_PATH}/external_libraries/nova-simd
					${SC_PATH}/external_libraries/nova-tt
					)
include_directories(${SC_PATH}/external_libraries/boost)
include_directories(${SC_PATH}/external_libraries/nova-simd)
add_definitions("-DNOVA_SIMD")

#INSTALL(TARGETS ${PROYECTO} ${PROYECTO}_supernova LIBRARY DESTINATION "C:/SupercolliderRepos/plugins/mingw48/source/")
INSTALL(TARGETS ${PROYECTO} ${PROYECTO}_supernova LIBRARY DESTINATION ${INSTALL_DIR})