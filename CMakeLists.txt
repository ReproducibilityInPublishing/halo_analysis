cmake_minimum_required (VERSION 2.8.11)
project (halo_analysis)

find_package(PkgConfig)

pkg_check_modules(libconfig++ REQUIRED libconfig++)
pkg_check_modules(ArgParse REQUIRED ArgParse)

include_directories(${ArgParse_INCLUDEDIR} ${libconfig++_INCLUDEDIR})

message(STATUS "ArgParse_FOUND: " ${ArgParse_FOUND})
message(STATUS "ArgParse_LIBRARIES: " ${ArgParse_LIBRARIES})
message(STATUS "ArgParse_LIBRARY_DIRS: " ${ArgParse_LIBRARY_DIRS})
message(STATUS "ArgParse_LDFLAGS: " ${ArgParse_LDFLAGS})
message(STATUS "ArgParse_LDFLAGS_OTHER: " ${ArgParse_LDFLAGS_OTHER})
message(STATUS "ArgParse_INCLUDE_DIRS: " ${ArgParse_INCLUDE_DIRS})
message(STATUS "ArgParse_CFLAGS: " ${ArgParse_CFLAGS})
message(STATUS "ArgParse_CFLAGS_OTHER: " ${ArgParse_CFLAGS_OTHER})
message(STATUS "ArgParse_VERSION: " ${ArgParse_VERSION})
message(STATUS "ArgParse_PREFIX: " ${ArgParse_PREFIX})
message(STATUS "ArgParse_INCLUDEDIR: " ${ArgParse_INCLUDEDIR})
message(STATUS "ArgParse_LIBDIR: " ${ArgParse_LIBDIR})
message(STATUS "libconfig++_INCLUDEDIR: " ${libconfig++_INCLUDEDIR})

add_executable(superhalo_finding superhalo_finding.cpp)

target_link_libraries(superhalo_finding ${libconfig++_LDFLAGS})
target_link_libraries(superhalo_finding ${ArgParse_LDFLAGS})

install (TARGETS superhalo_finding DESTINATION bin)
