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

# Utility rule file for testcore.

# Include the progress variables for this target.
include src/core/CMakeFiles/testcore.dir/progress.make

src/core/CMakeFiles/testcore: unittest/core_tests
src/core/CMakeFiles/testcore: src/core/runcoretests


src/core/runcoretests:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/chenyu97/github.com/openfhe-development/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating runcoretests"
	/home/chenyu97/github.com/openfhe-development/build/unittest/core_tests

testcore: src/core/CMakeFiles/testcore
testcore: src/core/runcoretests
testcore: src/core/CMakeFiles/testcore.dir/build.make

.PHONY : testcore

# Rule to build all files generated by this target.
src/core/CMakeFiles/testcore.dir/build: testcore

.PHONY : src/core/CMakeFiles/testcore.dir/build

src/core/CMakeFiles/testcore.dir/clean:
	cd /home/chenyu97/github.com/openfhe-development/build/src/core && $(CMAKE_COMMAND) -P CMakeFiles/testcore.dir/cmake_clean.cmake
.PHONY : src/core/CMakeFiles/testcore.dir/clean

src/core/CMakeFiles/testcore.dir/depend:
	cd /home/chenyu97/github.com/openfhe-development/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chenyu97/github.com/openfhe-development /home/chenyu97/github.com/openfhe-development/src/core /home/chenyu97/github.com/openfhe-development/build /home/chenyu97/github.com/openfhe-development/build/src/core /home/chenyu97/github.com/openfhe-development/build/src/core/CMakeFiles/testcore.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/core/CMakeFiles/testcore.dir/depend
