cmake_minimum_required(VERSION 3.10)
project(ofdm_algo CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
  # Plain (a+bi)(c+di) complex multiply/divide without the C99 Annex G inf/NaN
  # recovery path (__muldc3), which is very slow in MinGW's libgcc. Safe here:
  # signals never contain inf/NaN.
  add_compile_options(-fcx-limited-range)
endif()

add_library(ofdm
  src/qpsk.cpp
  src/fft.cpp
  src/ofdm.cpp
  src/channel.cpp
  src/metrics.cpp
  src/pa.cpp
)
target_include_directories(ofdm PUBLIC include)

foreach(app ofdm_demo ber_sim papr_ccdf pa_sim)
  add_executable(${app} apps/${app}.cpp)
  target_link_libraries(${app} PRIVATE ofdm)
endforeach()

enable_testing()
add_executable(ofdm_tests tests/test_main.cpp)
target_link_libraries(ofdm_tests PRIVATE ofdm)
add_test(NAME ofdm_tests COMMAND ofdm_tests)
