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
include src/core/CMakeFiles/sampling.dir/depend.make

# Include the progress variables for this target.
include src/core/CMakeFiles/sampling.dir/progress.make

# Include the compile flags for this target's objects.
include src/core/CMakeFiles/sampling.dir/flags.make

src/core/CMakeFiles/sampling.dir/examples/sampling.cpp.o: src/core/CMakeFiles/sampling.dir/flags.make
src/core/CMakeFiles/sampling.dir/examples/sampling.cpp.o: ../src/core/examples/sampling.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chenyu97/github.com/openfhe-development/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/core/CMakeFiles/sampling.dir/examples/sampling.cpp.o"
	cd /home/chenyu97/github.com/openfhe-development/build/src/core && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/sampling.dir/examples/sampling.cpp.o -c /home/chenyu97/github.com/openfhe-development/src/core/examples/sampling.cpp

src/core/CMakeFiles/sampling.dir/examples/sampling.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/sampling.dir/examples/sampling.cpp.i"
	cd /home/chenyu97/github.com/openfhe-development/build/src/core && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chenyu97/github.com/openfhe-development/src/core/examples/sampling.cpp > CMakeFiles/sampling.dir/examples/sampling.cpp.i

src/core/CMakeFiles/sampling.dir/examples/sampling.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/sampling.dir/examples/sampling.cpp.s"
	cd /home/chenyu97/github.com/openfhe-development/build/src/core && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chenyu97/github.com/openfhe-development/src/core/examples/sampling.cpp -o CMakeFiles/sampling.dir/examples/sampling.cpp.s

# Object files for target sampling
sampling_OBJECTS = \
"CMakeFiles/sampling.dir/examples/sampling.cpp.o"

# External object files for target sampling
sampling_EXTERNAL_OBJECTS =

bin/examples/core/sampling: src/core/CMakeFiles/sampling.dir/examples/sampling.cpp.o
bin/examples/core/sampling: src/core/CMakeFiles/sampling.dir/build.make
bin/examples/core/sampling: lib/libOPENFHEcore.so.1.2.1
bin/examples/core/sampling: src/core/CMakeFiles/sampling.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/chenyu97/github.com/openfhe-development/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../bin/examples/core/sampling"
	cd /home/chenyu97/github.com/openfhe-development/build/src/core && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sampling.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/core/CMakeFiles/sampling.dir/build: bin/examples/core/sampling

.PHONY : src/core/CMakeFiles/sampling.dir/build

src/core/CMakeFiles/sampling.dir/clean:
	cd /home/chenyu97/github.com/openfhe-development/build/src/core && $(CMAKE_COMMAND) -P CMakeFiles/sampling.dir/cmake_clean.cmake
.PHONY : src/core/CMakeFiles/sampling.dir/clean

src/core/CMakeFiles/sampling.dir/depend:
	cd /home/chenyu97/github.com/openfhe-development/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chenyu97/github.com/openfhe-development /home/chenyu97/github.com/openfhe-development/src/core /home/chenyu97/github.com/openfhe-development/build /home/chenyu97/github.com/openfhe-development/build/src/core /home/chenyu97/github.com/openfhe-development/build/src/core/CMakeFiles/sampling.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/core/CMakeFiles/sampling.dir/depend

