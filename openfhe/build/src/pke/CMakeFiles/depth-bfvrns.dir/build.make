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
include src/pke/CMakeFiles/depth-bfvrns.dir/depend.make

# Include the progress variables for this target.
include src/pke/CMakeFiles/depth-bfvrns.dir/progress.make

# Include the compile flags for this target's objects.
include src/pke/CMakeFiles/depth-bfvrns.dir/flags.make

src/pke/CMakeFiles/depth-bfvrns.dir/examples/depth-bfvrns.cpp.o: src/pke/CMakeFiles/depth-bfvrns.dir/flags.make
src/pke/CMakeFiles/depth-bfvrns.dir/examples/depth-bfvrns.cpp.o: ../src/pke/examples/depth-bfvrns.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chenyu97/github.com/openfhe-development/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/pke/CMakeFiles/depth-bfvrns.dir/examples/depth-bfvrns.cpp.o"
	cd /home/chenyu97/github.com/openfhe-development/build/src/pke && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/depth-bfvrns.dir/examples/depth-bfvrns.cpp.o -c /home/chenyu97/github.com/openfhe-development/src/pke/examples/depth-bfvrns.cpp

src/pke/CMakeFiles/depth-bfvrns.dir/examples/depth-bfvrns.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/depth-bfvrns.dir/examples/depth-bfvrns.cpp.i"
	cd /home/chenyu97/github.com/openfhe-development/build/src/pke && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chenyu97/github.com/openfhe-development/src/pke/examples/depth-bfvrns.cpp > CMakeFiles/depth-bfvrns.dir/examples/depth-bfvrns.cpp.i

src/pke/CMakeFiles/depth-bfvrns.dir/examples/depth-bfvrns.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/depth-bfvrns.dir/examples/depth-bfvrns.cpp.s"
	cd /home/chenyu97/github.com/openfhe-development/build/src/pke && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chenyu97/github.com/openfhe-development/src/pke/examples/depth-bfvrns.cpp -o CMakeFiles/depth-bfvrns.dir/examples/depth-bfvrns.cpp.s

# Object files for target depth-bfvrns
depth__bfvrns_OBJECTS = \
"CMakeFiles/depth-bfvrns.dir/examples/depth-bfvrns.cpp.o"

# External object files for target depth-bfvrns
depth__bfvrns_EXTERNAL_OBJECTS =

bin/examples/pke/depth-bfvrns: src/pke/CMakeFiles/depth-bfvrns.dir/examples/depth-bfvrns.cpp.o
bin/examples/pke/depth-bfvrns: src/pke/CMakeFiles/depth-bfvrns.dir/build.make
bin/examples/pke/depth-bfvrns: lib/libOPENFHEpke.so.1.2.1
bin/examples/pke/depth-bfvrns: lib/libOPENFHEbinfhe.so.1.2.1
bin/examples/pke/depth-bfvrns: lib/libOPENFHEcore.so.1.2.1
bin/examples/pke/depth-bfvrns: src/pke/CMakeFiles/depth-bfvrns.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/chenyu97/github.com/openfhe-development/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../bin/examples/pke/depth-bfvrns"
	cd /home/chenyu97/github.com/openfhe-development/build/src/pke && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/depth-bfvrns.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/pke/CMakeFiles/depth-bfvrns.dir/build: bin/examples/pke/depth-bfvrns

.PHONY : src/pke/CMakeFiles/depth-bfvrns.dir/build

src/pke/CMakeFiles/depth-bfvrns.dir/clean:
	cd /home/chenyu97/github.com/openfhe-development/build/src/pke && $(CMAKE_COMMAND) -P CMakeFiles/depth-bfvrns.dir/cmake_clean.cmake
.PHONY : src/pke/CMakeFiles/depth-bfvrns.dir/clean

src/pke/CMakeFiles/depth-bfvrns.dir/depend:
	cd /home/chenyu97/github.com/openfhe-development/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chenyu97/github.com/openfhe-development /home/chenyu97/github.com/openfhe-development/src/pke /home/chenyu97/github.com/openfhe-development/build /home/chenyu97/github.com/openfhe-development/build/src/pke /home/chenyu97/github.com/openfhe-development/build/src/pke/CMakeFiles/depth-bfvrns.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/pke/CMakeFiles/depth-bfvrns.dir/depend

