cmake_minimum_required(VERSION 3.0)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

project (bigbwt)
# set( CMAKE_VERBOSE_MAKEFILE on )

include_directories(${PROJECT_SOURCE_DIR})

set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -Wall -Wextra")
set(CMAKE_C_FLAGS "-O3 -Wall -std=c99 -g")

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -ggdb -g -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -fstrict-aliasing -march=native -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -ggdb -Ofast -fstrict-aliasing -march=native")


add_executable(bwtparse bwtparse.c gsa/gsacak.c utils.c malloc_count.c)
target_link_libraries(bwtparse dl)

add_executable(bwtparse64 bwtparse.c gsa/gsacak.c utils.c malloc_count.c)
target_link_libraries(bwtparse64 dl)
target_compile_options(bwtparse64 PUBLIC -DM64)

add_executable(simplebwt simplebwt.c gsa/gsacak.c)

add_executable(simplebwt64 simplebwt.c gsa/gsacak.c)
target_compile_options(simplebwt64 PUBLIC -DM64)

add_executable(newscanNT.x newscan.cpp malloc_count.c utils.c)
target_link_libraries(newscanNT.x z dl)
target_compile_options(newscanNT.x PUBLIC -DNOTHREADS) 

add_executable(newscan.x newscan.cpp malloc_count.c utils.c xerrors.c)
target_link_libraries(newscan.x z dl pthread)
target_compile_options(newscan.x PUBLIC -DM64) 

add_executable(pscan.x pscan.cpp pscan.hpp malloc_count.c utils.c xerrors.c)
target_link_libraries(pscan.x dl pthread)

add_executable(pfbwt.x pfbwt.cpp pfthreads.hpp gsa/gsacak.c utils.c xerrors.c malloc_count.c)
target_link_libraries(pfbwt.x dl pthread)

add_executable(pfbwt64.x pfbwt.cpp pfthreads.hpp gsa/gsacak.c utils.c xerrors.c malloc_count.c)
target_link_libraries(pfbwt64.x dl pthread)
target_compile_options(pfbwt64.x PUBLIC -DM64)

add_executable(pfbwtNT.x pfbwt.cpp gsa/gsacak.c utils.c malloc_count.c)
target_link_libraries(pfbwtNT.x dl)
target_compile_options(pfbwtNT.x PUBLIC -DNOTHREADS)

add_executable(pfbwtNT64.x pfbwt.cpp gsa/gsacak.c utils.c malloc_count.c)
target_link_libraries(pfbwtNT64.x dl)
target_compile_options(pfbwtNT64.x PUBLIC -DM64 -DNOTHREADS)

add_executable(unparse unparse.c utils.c malloc_count.c)
target_link_libraries(unparse dl)

add_executable(remap remap.c)
target_link_libraries(remap m)

configure_file(${PROJECT_SOURCE_DIR}/bigbwt ${PROJECT_BINARY_DIR}/bigbwt)

install(TARGETS bwtparse bwtparse64 simplebwt simplebwt64 newscan.x pscan.x pfbwt.x pfbwt64.x newscanNT.x pfbwtNT.x pfbwtNT64.x unparse remap
        DESTINATION bin)
install(PROGRAMS bigbwt DESTINATION bin)
