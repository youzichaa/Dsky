# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/chenyu97/github.com/openfhe-development

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/chenyu97/github.com/openfhe-development/build

# Utility rule file for allmodules.

# Include the progress variables for this target.
include CMakeFiles/allmodules.dir/progress.make

CMakeFiles/allmodules: lib/libOPENFHEcore.so
CMakeFiles/allmodules: lib/libOPENFHEpke.so
CMakeFiles/allmodules: lib/libOPENFHEbinfhe.so


allmodules: CMakeFiles/allmodules
allmodules: CMakeFiles/allmodules.dir/build.make

.PHONY : allmodules

# Rule to build all files generated by this target.
CMakeFiles/allmodules.dir/build: allmodules

.PHONY : CMakeFiles/allmodules.dir/build

CMakeFiles/allmodules.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/allmodules.dir/cmake_clean.cmake
.PHONY : CMakeFiles/allmodules.dir/clean

CMakeFiles/allmodules.dir/depend:
	cd /home/chenyu97/github.com/openfhe-development/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chenyu97/github.com/openfhe-development /home/chenyu97/github.com/openfhe-development /home/chenyu97/github.com/openfhe-development/build /home/chenyu97/github.com/openfhe-development/build /home/chenyu97/github.com/openfhe-development/build/CMakeFiles/allmodules.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/allmodules.dir/depend

