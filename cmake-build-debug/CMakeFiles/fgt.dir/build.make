# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/siyuchen/Downloads/clion/clion-2021.2.3/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/siyuchen/Downloads/clion/clion-2021.2.3/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/siyuchen/CLionProjects/fgt

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/siyuchen/CLionProjects/fgt/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/fgt.dir/depend.make
# Include the progress variables for this target.
include CMakeFiles/fgt.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/fgt.dir/flags.make

CMakeFiles/fgt.dir/main.c.o: CMakeFiles/fgt.dir/flags.make
CMakeFiles/fgt.dir/main.c.o: ../main.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/siyuchen/CLionProjects/fgt/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building C object CMakeFiles/fgt.dir/main.c.o"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/fgt.dir/main.c.o -c /home/siyuchen/CLionProjects/fgt/main.c

CMakeFiles/fgt.dir/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/fgt.dir/main.c.i"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/siyuchen/CLionProjects/fgt/main.c > CMakeFiles/fgt.dir/main.c.i

CMakeFiles/fgt.dir/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/fgt.dir/main.c.s"
	/usr/bin/cc $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/siyuchen/CLionProjects/fgt/main.c -o CMakeFiles/fgt.dir/main.c.s

# Object files for target fgt
fgt_OBJECTS = \
"CMakeFiles/fgt.dir/main.c.o"

# External object files for target fgt
fgt_EXTERNAL_OBJECTS =

fgt: CMakeFiles/fgt.dir/main.c.o
fgt: CMakeFiles/fgt.dir/build.make
fgt: CMakeFiles/fgt.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/siyuchen/CLionProjects/fgt/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking C executable fgt"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/fgt.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/fgt.dir/build: fgt
.PHONY : CMakeFiles/fgt.dir/build

CMakeFiles/fgt.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/fgt.dir/cmake_clean.cmake
.PHONY : CMakeFiles/fgt.dir/clean

CMakeFiles/fgt.dir/depend:
	cd /home/siyuchen/CLionProjects/fgt/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/siyuchen/CLionProjects/fgt /home/siyuchen/CLionProjects/fgt /home/siyuchen/CLionProjects/fgt/cmake-build-debug /home/siyuchen/CLionProjects/fgt/cmake-build-debug /home/siyuchen/CLionProjects/fgt/cmake-build-debug/CMakeFiles/fgt.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/fgt.dir/depend
