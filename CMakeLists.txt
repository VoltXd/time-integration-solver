cmake_minimum_required(VERSION 3.9) 

project(Time_Integration_Solver LANGUAGES CUDA CXX)

# Source & Header location
file(GLOB_RECURSE SOURCE_FILES src/*.cpp)
file(GLOB_RECURSE HEADER_FILES include/*.hpp)
file(GLOB_RECURSE CUDA_SOURCE_FILES src/*.cu)
file(GLOB_RECURSE CUDA_HEADER_FILES include/*.cuh)

# Make the executable with the headers & sources
add_executable(${PROJECT_NAME} ${SOURCE_FILES} ${HEADER_FILES} ${CUDA_SOURCE_FILES} ${CUDA_HEADER_FILES})

# Set the directory where CMakeLists.txt is to be the working directory
set_property(TARGET ${PROJECT_NAME} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

# Target the include folder (used to make headers visible in visual studio)
target_include_directories(${PROJECT_NAME} PRIVATE include)