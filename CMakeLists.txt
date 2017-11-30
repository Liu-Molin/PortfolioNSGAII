cmake_minimum_required(VERSION 3.8)
project(PortfolioNSGAII)

set(CMAKE_CXX_STANDARD 11)
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRES} ${GSLCBLAS_INCLUDE_DIRS} /usr/local/include)
link_libraries(${GSL_LIBRARIES} ${GSLBLAS_LIBRARIES} -fopenmp)

set(SOURCE_FILES main.cpp loader.h nsga2.h randG.h)
add_executable(PortfolioNSGAII  ${SOURCE_FILES})
target_link_libraries (PortfolioNSGAII ${GSL_LIBRARIES} ${GSLBLAS_LIBRARIES})