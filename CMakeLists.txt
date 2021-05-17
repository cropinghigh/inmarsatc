cmake_minimum_required(VERSION 3.14)

project(inmarsatc LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "/usr/")
endif()

add_library(inmarsatc SHARED
  inmarsatc.cpp
  inmarsatc.h
)

target_compile_definitions(inmarsatc PRIVATE INMARSATC_LIBRARY)

install(TARGETS inmarsatc DESTINATION lib)
install(FILES inmarsatc.h DESTINATION include)
