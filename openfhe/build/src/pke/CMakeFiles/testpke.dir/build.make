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

# Utility rule file for testpke.

# Include the progress variables for this target.
include src/pke/CMakeFiles/testpke.dir/progress.make

src/pke/CMakeFiles/testpke: unittest/pke_tests
src/pke/CMakeFiles/testpke: src/pke/runpketests


src/pke/runpketests:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/chenyu97/github.com/openfhe-development/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating runpketests"
	/home/chenyu97/github.com/openfhe-development/build/unittest/pke_tests

testpke: src/pke/CMakeFiles/testpke
testpke: src/pke/runpketests
testpke: src/pke/CMakeFiles/testpke.dir/build.make

.PHONY : testpke

# Rule to build all files generated by this target.
src/pke/CMakeFiles/testpke.dir/build: testpke

.PHONY : src/pke/CMakeFiles/testpke.dir/build

src/pke/CMakeFiles/testpke.dir/clean:
	cd /home/chenyu97/github.com/openfhe-development/build/src/pke && $(CMAKE_COMMAND) -P CMakeFiles/testpke.dir/cmake_clean.cmake
.PHONY : src/pke/CMakeFiles/testpke.dir/clean

src/pke/CMakeFiles/testpke.dir/depend:
	cd /home/chenyu97/github.com/openfhe-development/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chenyu97/github.com/openfhe-development /home/chenyu97/github.com/openfhe-development/src/pke /home/chenyu97/github.com/openfhe-development/build /home/chenyu97/github.com/openfhe-development/build/src/pke /home/chenyu97/github.com/openfhe-development/build/src/pke/CMakeFiles/testpke.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/pke/CMakeFiles/testpke.dir/depend

