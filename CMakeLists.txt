cmake_minimum_required (VERSION 2.8.11)
project (halo_analysis)

find_package(PkgConfig)

pkg_check_modules(libconfig++ REQUIRED libconfig++)
pkg_check_modules(ArgParse REQUIRED ArgParse)

message(STATUS "ArgParse_INCLUDEDIR: " ${ArgParse_INCLUDEDIR})
message(STATUS "ArgParse_LIBDIR: " ${ArgParse_LIBDIR})
message(STATUS "ArgParse_LDFLAGS: " ${ArgParse_LDFLAGS})
message(STATUS "libconfig++_INCLUDEDIR: " ${libconfig++_INCLUDEDIR})
message(STATUS "libconfig++_LIBDIR: " ${libconfig++_LIBDIR})
message(STATUS "libconfig++_LDFLAGS: " ${libconfig++_LDFLAGS})

include_directories(${ArgParse_INCLUDEDIR} ${libconfig++_INCLUDEDIR})

add_executable(superhalo_finding superhalo_finding.cpp)

target_link_libraries(superhalo_finding -L${libconfig++_LIBDIR} ${libconfig++_LDFLAGS})
target_link_libraries(superhalo_finding -L${ArgParse_LIBDIR} ${ArgParse_LDFLAGS})

install (TARGETS superhalo_finding DESTINATION bin)
