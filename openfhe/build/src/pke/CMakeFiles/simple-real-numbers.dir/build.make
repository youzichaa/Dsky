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
include src/pke/CMakeFiles/simple-real-numbers.dir/depend.make

# Include the progress variables for this target.
include src/pke/CMakeFiles/simple-real-numbers.dir/progress.make

# Include the compile flags for this target's objects.
include src/pke/CMakeFiles/simple-real-numbers.dir/flags.make

src/pke/CMakeFiles/simple-real-numbers.dir/examples/simple-real-numbers.cpp.o: src/pke/CMakeFiles/simple-real-numbers.dir/flags.make
src/pke/CMakeFiles/simple-real-numbers.dir/examples/simple-real-numbers.cpp.o: ../src/pke/examples/simple-real-numbers.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chenyu97/github.com/openfhe-development/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/pke/CMakeFiles/simple-real-numbers.dir/examples/simple-real-numbers.cpp.o"
	cd /home/chenyu97/github.com/openfhe-development/build/src/pke && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/simple-real-numbers.dir/examples/simple-real-numbers.cpp.o -c /home/chenyu97/github.com/openfhe-development/src/pke/examples/simple-real-numbers.cpp

src/pke/CMakeFiles/simple-real-numbers.dir/examples/simple-real-numbers.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/simple-real-numbers.dir/examples/simple-real-numbers.cpp.i"
	cd /home/chenyu97/github.com/openfhe-development/build/src/pke && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chenyu97/github.com/openfhe-development/src/pke/examples/simple-real-numbers.cpp > CMakeFiles/simple-real-numbers.dir/examples/simple-real-numbers.cpp.i

src/pke/CMakeFiles/simple-real-numbers.dir/examples/simple-real-numbers.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/simple-real-numbers.dir/examples/simple-real-numbers.cpp.s"
	cd /home/chenyu97/github.com/openfhe-development/build/src/pke && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chenyu97/github.com/openfhe-development/src/pke/examples/simple-real-numbers.cpp -o CMakeFiles/simple-real-numbers.dir/examples/simple-real-numbers.cpp.s

# Object files for target simple-real-numbers
simple__real__numbers_OBJECTS = \
"CMakeFiles/simple-real-numbers.dir/examples/simple-real-numbers.cpp.o"

# External object files for target simple-real-numbers
simple__real__numbers_EXTERNAL_OBJECTS =

bin/examples/pke/simple-real-numbers: src/pke/CMakeFiles/simple-real-numbers.dir/examples/simple-real-numbers.cpp.o
bin/examples/pke/simple-real-numbers: src/pke/CMakeFiles/simple-real-numbers.dir/build.make
bin/examples/pke/simple-real-numbers: lib/libOPENFHEpke.so.1.2.1
bin/examples/pke/simple-real-numbers: lib/libOPENFHEbinfhe.so.1.2.1
bin/examples/pke/simple-real-numbers: lib/libOPENFHEcore.so.1.2.1
bin/examples/pke/simple-real-numbers: src/pke/CMakeFiles/simple-real-numbers.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/chenyu97/github.com/openfhe-development/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../bin/examples/pke/simple-real-numbers"
	cd /home/chenyu97/github.com/openfhe-development/build/src/pke && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/simple-real-numbers.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/pke/CMakeFiles/simple-real-numbers.dir/build: bin/examples/pke/simple-real-numbers

.PHONY : src/pke/CMakeFiles/simple-real-numbers.dir/build

src/pke/CMakeFiles/simple-real-numbers.dir/clean:
	cd /home/chenyu97/github.com/openfhe-development/build/src/pke && $(CMAKE_COMMAND) -P CMakeFiles/simple-real-numbers.dir/cmake_clean.cmake
.PHONY : src/pke/CMakeFiles/simple-real-numbers.dir/clean

src/pke/CMakeFiles/simple-real-numbers.dir/depend:
	cd /home/chenyu97/github.com/openfhe-development/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chenyu97/github.com/openfhe-development /home/chenyu97/github.com/openfhe-development/src/pke /home/chenyu97/github.com/openfhe-development/build /home/chenyu97/github.com/openfhe-development/build/src/pke /home/chenyu97/github.com/openfhe-development/build/src/pke/CMakeFiles/simple-real-numbers.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/pke/CMakeFiles/simple-real-numbers.dir/depend

