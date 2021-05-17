cmake_minimum_required(VERSION 3.14)

project(inmarsatc LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(inmarsatc SHARED
  inmarsatc.cpp
  inmarsatc.h
)

target_compile_definitions(inmarsatc PRIVATE INMARSATC_LIBRARY)
