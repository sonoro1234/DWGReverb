cmake_minimum_required (VERSION 3.13)
set(PROYECTO "DWGReverb")
project (${PROYECTO})


include_directories(${SC_PATH}/include/plugin_interface)
include_directories(${SC_PATH}/include/common)
include_directories(${SC_PATH}/common)
include_directories(${SC_PATH}/external_libraries/libsndfile/)

include_directories(external_libraries/nova-simd)
add_definitions(-DNOVA_SIMD)


set(CMAKE_SHARED_MODULE_PREFIX "")
if(APPLE OR WIN32)
	set(CMAKE_SHARED_MODULE_SUFFIX ".scx")
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
	set(CMAKE_COMPILER_IS_CLANG 1)
endif()
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_CLANG)
    #add_definitions(-fvisibility=hidden)

    include (CheckCCompilerFlag)
    include (CheckCXXCompilerFlag)

    CHECK_C_COMPILER_FLAG(-msse HAS_SSE)
    CHECK_CXX_COMPILER_FLAG(-msse HAS_CXX_SSE)

    if (HAS_SSE)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse")
    endif()
    if (HAS_CXX_SSE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse")
    endif()

    CHECK_C_COMPILER_FLAG(-msse2 HAS_SSE2)
    CHECK_CXX_COMPILER_FLAG(-msse2 HAS_CXX_SSE2)

    if (HAS_SSE2)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse2")
    endif()
    if (HAS_CXX_SSE2)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse2")
    endif()

    CHECK_C_COMPILER_FLAG(-mfpmath=sse HAS_FPMATH_SSE)
    CHECK_CXX_COMPILER_FLAG(-mfpmath=sse HAS_CXX_FPMATH_SSE)

    if (HAS_FPMATH_SSE)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mfpmath=sse")
    endif()
    if (HAS_CXX_FPMATH_SSE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfpmath=sse")
    endif()


	#if(CPP11)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
		if(CMAKE_COMPILER_IS_CLANG)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
		endif()
	#endif()
endif()
if(MINGW)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mstackrealign")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mstackrealign")
endif()

aux_source_directory(. src)

add_library(${PROYECTO} MODULE ${src})


if (SUPERNOVA)
    include_directories(${SC_PATH}/external_libraries/nova-tt)
    # actually just boost.atomic
    include_directories(${SC_PATH}/external_libraries/boost)
    include_directories(${SC_PATH}/external_libraries/boost_lockfree)
    include_directories(${SC_PATH}/external_libraries/boost-lockfree)

    add_library(${PROYECTO}_supernova MODULE ${src})
    set_property(TARGET ${PROYECTO}_supernova APPEND PROPERTY COMPILE_DEFINITIONS  "SUPERNOVA")
endif()



#install
if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    set(INSTALL_DESTINATION "lib/SuperCollider/plugins")
    set(INSTALL_DESTINATION_DISTRO "share/SuperCollider/Extensions/Myplugins")
else()
    set(INSTALL_DESTINATION "Myplugins")
    set(INSTALL_DESTINATION_DISTRO "Myplugins")
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

install(DIRECTORY "sc/" DESTINATION ${INSTALL_DESTINATION_DISTRO}/${PROYECTO} PATTERN "*")
install(FILES lua/${PROYECTO}.lua DESTINATION ${INSTALL_DESTINATION_DISTRO}/lua )
