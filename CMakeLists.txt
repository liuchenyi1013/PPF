cmake_minimum_required(VERSION 3.0)
project(pcd_write)
set(CMAKE_CXX_STANDARD 11)

SET(CMAKE_BUILD_TYPE "Debug")
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g2 -ggdb")
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")

find_package(PCL 1.8 REQUIRED)
include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

add_executable(pcd_write pcd_write.cpp t_hash_int.cpp   ppf3d_public.h hash_murmur.hpp hash_murmur64.hpp  icp_copy.hpp ppf3d_public.h stdafx.h  t_hash_int.hpp t_hash_int.cpp c_utils.h)

target_link_libraries (pcd_write ${PCL_LIBRARIES})

install(TARGETS pcd_write RUNTIME DESTINATION bin)
