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
include benchmark/CMakeFiles/binfhe-ginx.dir/depend.make

# Include the progress variables for this target.
include benchmark/CMakeFiles/binfhe-ginx.dir/progress.make

# Include the compile flags for this target's objects.
include benchmark/CMakeFiles/binfhe-ginx.dir/flags.make

benchmark/CMakeFiles/binfhe-ginx.dir/src/binfhe-ginx.cpp.o: benchmark/CMakeFiles/binfhe-ginx.dir/flags.make
benchmark/CMakeFiles/binfhe-ginx.dir/src/binfhe-ginx.cpp.o: ../benchmark/src/binfhe-ginx.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chenyu97/github.com/openfhe-development/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object benchmark/CMakeFiles/binfhe-ginx.dir/src/binfhe-ginx.cpp.o"
	cd /home/chenyu97/github.com/openfhe-development/build/benchmark && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/binfhe-ginx.dir/src/binfhe-ginx.cpp.o -c /home/chenyu97/github.com/openfhe-development/benchmark/src/binfhe-ginx.cpp

benchmark/CMakeFiles/binfhe-ginx.dir/src/binfhe-ginx.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/binfhe-ginx.dir/src/binfhe-ginx.cpp.i"
	cd /home/chenyu97/github.com/openfhe-development/build/benchmark && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chenyu97/github.com/openfhe-development/benchmark/src/binfhe-ginx.cpp > CMakeFiles/binfhe-ginx.dir/src/binfhe-ginx.cpp.i

benchmark/CMakeFiles/binfhe-ginx.dir/src/binfhe-ginx.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/binfhe-ginx.dir/src/binfhe-ginx.cpp.s"
	cd /home/chenyu97/github.com/openfhe-development/build/benchmark && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chenyu97/github.com/openfhe-development/benchmark/src/binfhe-ginx.cpp -o CMakeFiles/binfhe-ginx.dir/src/binfhe-ginx.cpp.s

# Object files for target binfhe-ginx
binfhe__ginx_OBJECTS = \
"CMakeFiles/binfhe-ginx.dir/src/binfhe-ginx.cpp.o"

# External object files for target binfhe-ginx
binfhe__ginx_EXTERNAL_OBJECTS =

bin/benchmark/binfhe-ginx: benchmark/CMakeFiles/binfhe-ginx.dir/src/binfhe-ginx.cpp.o
bin/benchmark/binfhe-ginx: benchmark/CMakeFiles/binfhe-ginx.dir/build.make
bin/benchmark/binfhe-ginx: lib/libOPENFHEpke.so.1.2.1
bin/benchmark/binfhe-ginx: lib/libOPENFHEbinfhe.so.1.2.1
bin/benchmark/binfhe-ginx: lib/libOPENFHEcore.so.1.2.1
bin/benchmark/binfhe-ginx: lib/libbenchmark.a
bin/benchmark/binfhe-ginx: /usr/lib/x86_64-linux-gnu/librt.so
bin/benchmark/binfhe-ginx: benchmark/CMakeFiles/binfhe-ginx.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/chenyu97/github.com/openfhe-development/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/benchmark/binfhe-ginx"
	cd /home/chenyu97/github.com/openfhe-development/build/benchmark && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/binfhe-ginx.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
benchmark/CMakeFiles/binfhe-ginx.dir/build: bin/benchmark/binfhe-ginx

.PHONY : benchmark/CMakeFiles/binfhe-ginx.dir/build

benchmark/CMakeFiles/binfhe-ginx.dir/clean:
	cd /home/chenyu97/github.com/openfhe-development/build/benchmark && $(CMAKE_COMMAND) -P CMakeFiles/binfhe-ginx.dir/cmake_clean.cmake
.PHONY : benchmark/CMakeFiles/binfhe-ginx.dir/clean

benchmark/CMakeFiles/binfhe-ginx.dir/depend:
	cd /home/chenyu97/github.com/openfhe-development/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chenyu97/github.com/openfhe-development /home/chenyu97/github.com/openfhe-development/benchmark /home/chenyu97/github.com/openfhe-development/build /home/chenyu97/github.com/openfhe-development/build/benchmark /home/chenyu97/github.com/openfhe-development/build/benchmark/CMakeFiles/binfhe-ginx.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : benchmark/CMakeFiles/binfhe-ginx.dir/depend
