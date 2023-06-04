cmake_minimum_required(VERSION 3.20)
project(fgt C)

set(CMAKE_C_STANDARD 11)

add_executable(fgt main.c)
target_link_libraries(fgt m)
