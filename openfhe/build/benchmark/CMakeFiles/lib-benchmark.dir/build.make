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
include benchmark/CMakeFiles/lib-benchmark.dir/depend.make

# Include the progress variables for this target.
include benchmark/CMakeFiles/lib-benchmark.dir/progress.make

# Include the compile flags for this target's objects.
include benchmark/CMakeFiles/lib-benchmark.dir/flags.make

benchmark/CMakeFiles/lib-benchmark.dir/src/lib-benchmark.cpp.o: benchmark/CMakeFiles/lib-benchmark.dir/flags.make
benchmark/CMakeFiles/lib-benchmark.dir/src/lib-benchmark.cpp.o: ../benchmark/src/lib-benchmark.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chenyu97/github.com/openfhe-development/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object benchmark/CMakeFiles/lib-benchmark.dir/src/lib-benchmark.cpp.o"
	cd /home/chenyu97/github.com/openfhe-development/build/benchmark && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/lib-benchmark.dir/src/lib-benchmark.cpp.o -c /home/chenyu97/github.com/openfhe-development/benchmark/src/lib-benchmark.cpp

benchmark/CMakeFiles/lib-benchmark.dir/src/lib-benchmark.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/lib-benchmark.dir/src/lib-benchmark.cpp.i"
	cd /home/chenyu97/github.com/openfhe-development/build/benchmark && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chenyu97/github.com/openfhe-development/benchmark/src/lib-benchmark.cpp > CMakeFiles/lib-benchmark.dir/src/lib-benchmark.cpp.i

benchmark/CMakeFiles/lib-benchmark.dir/src/lib-benchmark.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/lib-benchmark.dir/src/lib-benchmark.cpp.s"
	cd /home/chenyu97/github.com/openfhe-development/build/benchmark && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chenyu97/github.com/openfhe-development/benchmark/src/lib-benchmark.cpp -o CMakeFiles/lib-benchmark.dir/src/lib-benchmark.cpp.s

# Object files for target lib-benchmark
lib__benchmark_OBJECTS = \
"CMakeFiles/lib-benchmark.dir/src/lib-benchmark.cpp.o"

# External object files for target lib-benchmark
lib__benchmark_EXTERNAL_OBJECTS =

bin/benchmark/lib-benchmark: benchmark/CMakeFiles/lib-benchmark.dir/src/lib-benchmark.cpp.o
bin/benchmark/lib-benchmark: benchmark/CMakeFiles/lib-benchmark.dir/build.make
bin/benchmark/lib-benchmark: lib/libOPENFHEpke.so.1.2.1
bin/benchmark/lib-benchmark: lib/libOPENFHEbinfhe.so.1.2.1
bin/benchmark/lib-benchmark: lib/libOPENFHEcore.so.1.2.1
bin/benchmark/lib-benchmark: lib/libbenchmark.a
bin/benchmark/lib-benchmark: /usr/lib/x86_64-linux-gnu/librt.so
bin/benchmark/lib-benchmark: benchmark/CMakeFiles/lib-benchmark.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/chenyu97/github.com/openfhe-development/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/benchmark/lib-benchmark"
	cd /home/chenyu97/github.com/openfhe-development/build/benchmark && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/lib-benchmark.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
benchmark/CMakeFiles/lib-benchmark.dir/build: bin/benchmark/lib-benchmark

.PHONY : benchmark/CMakeFiles/lib-benchmark.dir/build

benchmark/CMakeFiles/lib-benchmark.dir/clean:
	cd /home/chenyu97/github.com/openfhe-development/build/benchmark && $(CMAKE_COMMAND) -P CMakeFiles/lib-benchmark.dir/cmake_clean.cmake
.PHONY : benchmark/CMakeFiles/lib-benchmark.dir/clean

benchmark/CMakeFiles/lib-benchmark.dir/depend:
	cd /home/chenyu97/github.com/openfhe-development/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chenyu97/github.com/openfhe-development /home/chenyu97/github.com/openfhe-development/benchmark /home/chenyu97/github.com/openfhe-development/build /home/chenyu97/github.com/openfhe-development/build/benchmark /home/chenyu97/github.com/openfhe-development/build/benchmark/CMakeFiles/lib-benchmark.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : benchmark/CMakeFiles/lib-benchmark.dir/depend
