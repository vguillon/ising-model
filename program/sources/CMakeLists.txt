cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

project(Ising)


add_executable(Ising main.cpp Ising/ising1d.cpp Ising/ising2d.cpp)
target_include_directories(Ising PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/Ising)


