cmake_minimum_required(VERSION 3.27)
project(convex_hull_3d)

set(CMAKE_CXX_STANDARD 23)

add_executable(convex_hull_3d main.cpp
        convex_hull_3D.h
        Face.h
        TetrahedronCreator.h
        BipartiteGraph.h
        CustomUnorderedSet.h)
