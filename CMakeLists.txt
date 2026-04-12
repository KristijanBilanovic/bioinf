cmake_minimum_required(VERSION 3.12)
project(JelenProject)

set(CMAKE_CXX_STANDARD 17)

# --- SPOA SEKCIJA ---
# Isključujemo nepotrebne stvari iz SPOA-e da se brže kompajlira
set(spoa_build_exe OFF CACHE BOOL "" FORCE)
set(spoa_build_tests OFF CACHE BOOL "" FORCE)

# add_subdirectory kaže CMake-u: "Uđi u ovaj folder i nađi njegov CMakeLists.txt"
# To će kompajlirati SPOA-u zajedno s tvojim projektom.
add_subdirectory(vendor/spoa)

# --- BIOPARSER SEKCIJA ---
# Bioparser je "header-only", što znači da nema .cpp datoteka, samo .hpp.
# Njega ne moramo kompajlirati, samo kažemo gdje su mu datoteke.
add_library(bioparser INTERFACE)
target_include_directories(bioparser INTERFACE vendor/bioparser/include)

# --- TVOJ PROJEKT ---
add_executable(jelen_analiza main.cpp)

# target_link_libraries je KLJUČNA naredba. 
# Ona kaže: "Uzmi moj kompajlirani main.cpp i spoji ga s kompajliranom SPOA bibliotekom."
target_link_libraries(jelen_analiza spoa bioparser)