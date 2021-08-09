cmake_minimum_required(VERSION 3.14)

project(inmarsatc LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "/usr/")
endif()

add_library(inmarsatc_demodulator SHARED
  inmarsatc_demodulator.cpp
  inmarsatc_demodulator.h
)

add_library(inmarsatc_decoder SHARED
  inmarsatc_decoder.cpp
  inmarsatc_decoder.h
)

add_library(inmarsatc_parser SHARED
  inmarsatc_parser.cpp
  inmarsatc_parser.h
)

target_compile_definitions(inmarsatc_demodulator PRIVATE INMARSATC_LIBRARY)
target_compile_definitions(inmarsatc_decoder PRIVATE INMARSATC_LIBRARY)
target_compile_definitions(inmarsatc_parser PRIVATE INMARSATC_LIBRARY)

install(TARGETS inmarsatc_demodulator inmarsatc_decoder inmarsatc_parser DESTINATION lib)
install(FILES inmarsatc_demodulator.h inmarsatc_decoder.h inmarsatc_parser.h DESTINATION include)
