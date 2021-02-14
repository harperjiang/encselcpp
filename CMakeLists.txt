cmake_minimum_required (VERSION 3.16)
project (encsel)



include_directories(sparsehash)

add_executable(encsel_featuretime encsel_featuretime.cc feature.cc)
add_executable(encsel_int encsel_int.cc feature.cc)
add_executable(encsel_str encsel_str.cc feature.cc)
add_compile_options("-march=skylake")