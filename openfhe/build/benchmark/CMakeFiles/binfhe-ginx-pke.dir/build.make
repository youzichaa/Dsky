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
include benchmark/CMakeFiles/binfhe-ginx-pke.dir/depend.make

# Include the progress variables for this target.
include benchmark/CMakeFiles/binfhe-ginx-pke.dir/progress.make

# Include the compile flags for this target's objects.
include benchmark/CMakeFiles/binfhe-ginx-pke.dir/flags.make

benchmark/CMakeFiles/binfhe-ginx-pke.dir/src/binfhe-ginx-pke.cpp.o: benchmark/CMakeFiles/binfhe-ginx-pke.dir/flags.make
benchmark/CMakeFiles/binfhe-ginx-pke.dir/src/binfhe-ginx-pke.cpp.o: ../benchmark/src/binfhe-ginx-pke.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chenyu97/github.com/openfhe-development/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object benchmark/CMakeFiles/binfhe-ginx-pke.dir/src/binfhe-ginx-pke.cpp.o"
	cd /home/chenyu97/github.com/openfhe-development/build/benchmark && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/binfhe-ginx-pke.dir/src/binfhe-ginx-pke.cpp.o -c /home/chenyu97/github.com/openfhe-development/benchmark/src/binfhe-ginx-pke.cpp

benchmark/CMakeFiles/binfhe-ginx-pke.dir/src/binfhe-ginx-pke.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/binfhe-ginx-pke.dir/src/binfhe-ginx-pke.cpp.i"
	cd /home/chenyu97/github.com/openfhe-development/build/benchmark && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chenyu97/github.com/openfhe-development/benchmark/src/binfhe-ginx-pke.cpp > CMakeFiles/binfhe-ginx-pke.dir/src/binfhe-ginx-pke.cpp.i

benchmark/CMakeFiles/binfhe-ginx-pke.dir/src/binfhe-ginx-pke.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/binfhe-ginx-pke.dir/src/binfhe-ginx-pke.cpp.s"
	cd /home/chenyu97/github.com/openfhe-development/build/benchmark && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chenyu97/github.com/openfhe-development/benchmark/src/binfhe-ginx-pke.cpp -o CMakeFiles/binfhe-ginx-pke.dir/src/binfhe-ginx-pke.cpp.s

# Object files for target binfhe-ginx-pke
binfhe__ginx__pke_OBJECTS = \
"CMakeFiles/binfhe-ginx-pke.dir/src/binfhe-ginx-pke.cpp.o"

# External object files for target binfhe-ginx-pke
binfhe__ginx__pke_EXTERNAL_OBJECTS =

bin/benchmark/binfhe-ginx-pke: benchmark/CMakeFiles/binfhe-ginx-pke.dir/src/binfhe-ginx-pke.cpp.o
bin/benchmark/binfhe-ginx-pke: benchmark/CMakeFiles/binfhe-ginx-pke.dir/build.make
bin/benchmark/binfhe-ginx-pke: lib/libOPENFHEpke.so.1.2.1
bin/benchmark/binfhe-ginx-pke: lib/libOPENFHEbinfhe.so.1.2.1
bin/benchmark/binfhe-ginx-pke: lib/libOPENFHEcore.so.1.2.1
bin/benchmark/binfhe-ginx-pke: lib/libbenchmark.a
bin/benchmark/binfhe-ginx-pke: /usr/lib/x86_64-linux-gnu/librt.so
bin/benchmark/binfhe-ginx-pke: benchmark/CMakeFiles/binfhe-ginx-pke.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/chenyu97/github.com/openfhe-development/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/benchmark/binfhe-ginx-pke"
	cd /home/chenyu97/github.com/openfhe-development/build/benchmark && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/binfhe-ginx-pke.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
benchmark/CMakeFiles/binfhe-ginx-pke.dir/build: bin/benchmark/binfhe-ginx-pke

.PHONY : benchmark/CMakeFiles/binfhe-ginx-pke.dir/build

benchmark/CMakeFiles/binfhe-ginx-pke.dir/clean:
	cd /home/chenyu97/github.com/openfhe-development/build/benchmark && $(CMAKE_COMMAND) -P CMakeFiles/binfhe-ginx-pke.dir/cmake_clean.cmake
.PHONY : benchmark/CMakeFiles/binfhe-ginx-pke.dir/clean

benchmark/CMakeFiles/binfhe-ginx-pke.dir/depend:
	cd /home/chenyu97/github.com/openfhe-development/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chenyu97/github.com/openfhe-development /home/chenyu97/github.com/openfhe-development/benchmark /home/chenyu97/github.com/openfhe-development/build /home/chenyu97/github.com/openfhe-development/build/benchmark /home/chenyu97/github.com/openfhe-development/build/benchmark/CMakeFiles/binfhe-ginx-pke.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : benchmark/CMakeFiles/binfhe-ginx-pke.dir/depend
