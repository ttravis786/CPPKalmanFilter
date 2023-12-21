cmake_minimum_required(VERSION 3.22)
project(cplusdemo)

set(CMAKE_CXX_STANDARD 17)

add_executable(cplusdemo main.cpp
        src/CombKalFilterBranch.cpp
        src/CombKalFilterBranch.h
        src/CombKalFilter.cpp
        src/CombKalFilter.h
        test/CombKalFilterTest.cpp
        main.cpp
        src/MathFunctions.cpp
        src/MathFunctions.h
        src/CellFinder.cpp
        src/CellFinder.h
        src/CSVFunctions.cpp
        src/CSVFunctions.h
)

set(EIGEN3_DIR /home/tt1020/CPP_LIBS/eigen-3.4.0)

include_directories(${EIGEN3_DIR})

