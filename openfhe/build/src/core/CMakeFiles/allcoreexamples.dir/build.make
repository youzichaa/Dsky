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

# Utility rule file for allcoreexamples.

# Include the progress variables for this target.
include src/core/CMakeFiles/allcoreexamples.dir/progress.make

allcoreexamples: src/core/CMakeFiles/allcoreexamples.dir/build.make

.PHONY : allcoreexamples

# Rule to build all files generated by this target.
src/core/CMakeFiles/allcoreexamples.dir/build: allcoreexamples

.PHONY : src/core/CMakeFiles/allcoreexamples.dir/build

src/core/CMakeFiles/allcoreexamples.dir/clean:
	cd /home/chenyu97/github.com/openfhe-development/build/src/core && $(CMAKE_COMMAND) -P CMakeFiles/allcoreexamples.dir/cmake_clean.cmake
.PHONY : src/core/CMakeFiles/allcoreexamples.dir/clean

src/core/CMakeFiles/allcoreexamples.dir/depend:
	cd /home/chenyu97/github.com/openfhe-development/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chenyu97/github.com/openfhe-development /home/chenyu97/github.com/openfhe-development/src/core /home/chenyu97/github.com/openfhe-development/build /home/chenyu97/github.com/openfhe-development/build/src/core /home/chenyu97/github.com/openfhe-development/build/src/core/CMakeFiles/allcoreexamples.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/core/CMakeFiles/allcoreexamples.dir/depend

