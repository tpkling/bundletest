cmake_minimum_required(VERSION 2.8)
project(bundletet)

 set(SOURCE_FILES main.cpp KerrRayTDerivs.cpp KerrRayTFuncs.cpp KerrRayVDots.cpp Coords.h InitialValues.cpp
        InitialValues.h KerrRay.cpp KerrRay.h)

 add_executable(bundletest ${SOURCE_FILES})