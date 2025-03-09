cmake_minimum_required(VERSION 3.24)
cmake_policy(SET CMP0135 NEW)

project(simplex
  HOMEPAGE_URL "https://github.com/ytsao/simplex-method.git"
  LANGUAGES CXX
)

add_executable(simplex main.cpp)
target_include_directories(simplex PRIVATE include)
target_compile_features(simplex PRIVATE cxx_std_17)
