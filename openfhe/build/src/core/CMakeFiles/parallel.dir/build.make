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

# Include any dependencies generated for this target.
include src/core/CMakeFiles/parallel.dir/depend.make

# Include the progress variables for this target.
include src/core/CMakeFiles/parallel.dir/progress.make

# Include the compile flags for this target's objects.
include src/core/CMakeFiles/parallel.dir/flags.make

src/core/CMakeFiles/parallel.dir/examples/parallel.cpp.o: src/core/CMakeFiles/parallel.dir/flags.make
src/core/CMakeFiles/parallel.dir/examples/parallel.cpp.o: ../src/core/examples/parallel.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chenyu97/github.com/openfhe-development/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/core/CMakeFiles/parallel.dir/examples/parallel.cpp.o"
	cd /home/chenyu97/github.com/openfhe-development/build/src/core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/parallel.dir/examples/parallel.cpp.o -c /home/chenyu97/github.com/openfhe-development/src/core/examples/parallel.cpp

src/core/CMakeFiles/parallel.dir/examples/parallel.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/parallel.dir/examples/parallel.cpp.i"
	cd /home/chenyu97/github.com/openfhe-development/build/src/core && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chenyu97/github.com/openfhe-development/src/core/examples/parallel.cpp > CMakeFiles/parallel.dir/examples/parallel.cpp.i

src/core/CMakeFiles/parallel.dir/examples/parallel.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/parallel.dir/examples/parallel.cpp.s"
	cd /home/chenyu97/github.com/openfhe-development/build/src/core && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chenyu97/github.com/openfhe-development/src/core/examples/parallel.cpp -o CMakeFiles/parallel.dir/examples/parallel.cpp.s

# Object files for target parallel
parallel_OBJECTS = \
"CMakeFiles/parallel.dir/examples/parallel.cpp.o"

# External object files for target parallel
parallel_EXTERNAL_OBJECTS =

bin/examples/core/parallel: src/core/CMakeFiles/parallel.dir/examples/parallel.cpp.o
bin/examples/core/parallel: src/core/CMakeFiles/parallel.dir/build.make
bin/examples/core/parallel: lib/libOPENFHEcore.so.1.2.1
bin/examples/core/parallel: src/core/CMakeFiles/parallel.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/chenyu97/github.com/openfhe-development/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../bin/examples/core/parallel"
	cd /home/chenyu97/github.com/openfhe-development/build/src/core && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/parallel.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/core/CMakeFiles/parallel.dir/build: bin/examples/core/parallel

.PHONY : src/core/CMakeFiles/parallel.dir/build

src/core/CMakeFiles/parallel.dir/clean:
	cd /home/chenyu97/github.com/openfhe-development/build/src/core && $(CMAKE_COMMAND) -P CMakeFiles/parallel.dir/cmake_clean.cmake
.PHONY : src/core/CMakeFiles/parallel.dir/clean

src/core/CMakeFiles/parallel.dir/depend:
	cd /home/chenyu97/github.com/openfhe-development/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chenyu97/github.com/openfhe-development /home/chenyu97/github.com/openfhe-development/src/core /home/chenyu97/github.com/openfhe-development/build /home/chenyu97/github.com/openfhe-development/build/src/core /home/chenyu97/github.com/openfhe-development/build/src/core/CMakeFiles/parallel.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/core/CMakeFiles/parallel.dir/depend

