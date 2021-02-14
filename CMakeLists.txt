cmake_minimum_required (VERSION 3.16)
project (encsel)

add_executable(encsel_featuretime encsel_featuretime.cc feature.cc)
add_compile_options("-march=skylake")