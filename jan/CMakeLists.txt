cmake_minimum_required(VERSION 3.15)
project(Big_XBWT)

set(CMAKE_CXX_STANDARD 17)

add_executable(newscan
        newscan.cpp
        utils.c
        utils.h)

add_executable(bwtparse
        bwtparse.cpp
        utils.c
        )
