cmake_minimum_required(VERSION 3.12)
project(JelenProject)

set(CMAKE_CXX_STANDARD 17)

# --- SPOA SEKCIJA ---
set(spoa_build_exe OFF CACHE BOOL "" FORCE)
set(spoa_build_tests OFF CACHE BOOL "" FORCE)
add_subdirectory(vendor/spoa)

# --- BIOPARSER SEKCIJA ---
# Bioparser is a dependency of SPOA, located in build/_deps/bioparser-src
add_library(bioparser INTERFACE)
target_include_directories(bioparser INTERFACE ${CMAKE_SOURCE_DIR}/vendor/spoa/build/_deps/bioparser-src/include)

# --- TVOJ PROJEKT ---
# Add your own source files here
add_executable(jelen_analiza 
    main.cpp
    src/MyClass.cpp
)

# Include your headers
target_include_directories(jelen_analiza PRIVATE include)

# Link libraries (bioparser is included via spoa)
target_link_libraries(jelen_analiza spoa bioparser)

# Ensure SPOA is built first
add_dependencies(jelen_analiza spoa)