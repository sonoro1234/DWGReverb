cmake_minimum_required (VERSION 2.8)
set(PROYECTO "DWGReverb")
project (${PROYECTO})


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

if(SUPERNOVA)
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
endif()



#install
if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    set(INSTALL_DESTINATION "lib/SuperCollider/plugins")
    if (QUARKS)
        set(INSTALL_DESTINATION_LADSPALIST
                "share/SuperCollider/SC3plugins/LadspaUGen")
        set(INSTALL_DESTINATION_DISTRO
                "share/SuperCollider/SC3plugins")
    else()
        set(INSTALL_DESTINATION_LADSPALIST
                "share/SuperCollider/Extensions/SC3plugins/LadspaUGen")
        set(INSTALL_DESTINATION_DISTRO
                "share/SuperCollider/Extensions/SC3plugins")
    endif(QUARKS)
else()
    set(INSTALL_DESTINATION "SC3plugins")
    set(INSTALL_DESTINATION_LADSPALIST "SC3plugins/LadspaUGen")
    set(INSTALL_DESTINATION_DISTRO "SC3plugins")
endif()


if (NOT APPLE)
    install(TARGETS ${PROYECTO} DESTINATION ${INSTALL_DESTINATION})
    if(SUPERNOVA)
        install(TARGETS ${PROYECTO}_supernova DESTINATION ${INSTALL_DESTINATION})
    endif()
else()
    install(TARGETS ${PROYECTO} DESTINATION ${INSTALL_DESTINATION}/${PROYECTO})
    if(SUPERNOVA)
        install(TARGETS ${PROYECTO}_supernova DESTINATION ${INSTALL_DESTINATION}/${PROYECTO})
    endif()
endif()

install(DIRECTORY "sclang/" DESTINATION ${INSTALL_DESTINATION_DISTRO}/${PROYECTO} PATTERN "*")
