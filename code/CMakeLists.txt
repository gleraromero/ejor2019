cmake_minimum_required(VERSION 2.8.4)
project(ejor2019)

# Set compiler flags.
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++14 -m64")
if ("${CMAKE_BUILD_TYPE}" MATCHES release OR "${CMAKE_BUILD_TYPE}" MATCHES Release)
    message("-- Compiling in Release Mode")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -DNDEBUG")
elseif ("${CMAKE_BUILD_TYPE}" MATCHES debug OR "${CMAKE_BUILD_TYPE}" MATCHES Debug)
    message("-- Compiling in Debug Mode")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()

if ("${RUNNER}" MATCHES ON)
    message("-- Runner compilation")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DRUNNER")
endif()

# Add include directory.
include_directories("include/")

# Including GOC.
add_subdirectory(goc)
include_directories(goc/include)

# Create binaries.
add_executable(ejor2019 src/main.cpp src/vrp_instance.cpp src/preprocess_travel_times.cpp src/gcs.cpp src/preprocess_service_waiting.cpp src/preprocess_time_windows.cpp src/preprocess_capacity.cpp src/ttbf.cpp src/cttbf.cpp src/model.cpp src/weak_pi.cpp src/weak_sigma.cpp src/weak_pisigma.cpp src/tdfi.cpp src/tdcsi.cpp)
target_link_libraries(ejor2019 goc)
target_link_libraries(ejor2019)