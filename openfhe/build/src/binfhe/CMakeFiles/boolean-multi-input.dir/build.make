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
include src/binfhe/CMakeFiles/boolean-multi-input.dir/depend.make

# Include the progress variables for this target.
include src/binfhe/CMakeFiles/boolean-multi-input.dir/progress.make

# Include the compile flags for this target's objects.
include src/binfhe/CMakeFiles/boolean-multi-input.dir/flags.make

src/binfhe/CMakeFiles/boolean-multi-input.dir/examples/boolean-multi-input.cpp.o: src/binfhe/CMakeFiles/boolean-multi-input.dir/flags.make
src/binfhe/CMakeFiles/boolean-multi-input.dir/examples/boolean-multi-input.cpp.o: ../src/binfhe/examples/boolean-multi-input.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/chenyu97/github.com/openfhe-development/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object src/binfhe/CMakeFiles/boolean-multi-input.dir/examples/boolean-multi-input.cpp.o"
	cd /home/chenyu97/github.com/openfhe-development/build/src/binfhe && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/boolean-multi-input.dir/examples/boolean-multi-input.cpp.o -c /home/chenyu97/github.com/openfhe-development/src/binfhe/examples/boolean-multi-input.cpp

src/binfhe/CMakeFiles/boolean-multi-input.dir/examples/boolean-multi-input.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/boolean-multi-input.dir/examples/boolean-multi-input.cpp.i"
	cd /home/chenyu97/github.com/openfhe-development/build/src/binfhe && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/chenyu97/github.com/openfhe-development/src/binfhe/examples/boolean-multi-input.cpp > CMakeFiles/boolean-multi-input.dir/examples/boolean-multi-input.cpp.i

src/binfhe/CMakeFiles/boolean-multi-input.dir/examples/boolean-multi-input.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/boolean-multi-input.dir/examples/boolean-multi-input.cpp.s"
	cd /home/chenyu97/github.com/openfhe-development/build/src/binfhe && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/chenyu97/github.com/openfhe-development/src/binfhe/examples/boolean-multi-input.cpp -o CMakeFiles/boolean-multi-input.dir/examples/boolean-multi-input.cpp.s

# Object files for target boolean-multi-input
boolean__multi__input_OBJECTS = \
"CMakeFiles/boolean-multi-input.dir/examples/boolean-multi-input.cpp.o"

# External object files for target boolean-multi-input
boolean__multi__input_EXTERNAL_OBJECTS =

bin/examples/binfhe/boolean-multi-input: src/binfhe/CMakeFiles/boolean-multi-input.dir/examples/boolean-multi-input.cpp.o
bin/examples/binfhe/boolean-multi-input: src/binfhe/CMakeFiles/boolean-multi-input.dir/build.make
bin/examples/binfhe/boolean-multi-input: lib/libOPENFHEbinfhe.so.1.2.1
bin/examples/binfhe/boolean-multi-input: lib/libOPENFHEcore.so.1.2.1
bin/examples/binfhe/boolean-multi-input: src/binfhe/CMakeFiles/boolean-multi-input.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/chenyu97/github.com/openfhe-development/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../../bin/examples/binfhe/boolean-multi-input"
	cd /home/chenyu97/github.com/openfhe-development/build/src/binfhe && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/boolean-multi-input.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/binfhe/CMakeFiles/boolean-multi-input.dir/build: bin/examples/binfhe/boolean-multi-input

.PHONY : src/binfhe/CMakeFiles/boolean-multi-input.dir/build

src/binfhe/CMakeFiles/boolean-multi-input.dir/clean:
	cd /home/chenyu97/github.com/openfhe-development/build/src/binfhe && $(CMAKE_COMMAND) -P CMakeFiles/boolean-multi-input.dir/cmake_clean.cmake
.PHONY : src/binfhe/CMakeFiles/boolean-multi-input.dir/clean

src/binfhe/CMakeFiles/boolean-multi-input.dir/depend:
	cd /home/chenyu97/github.com/openfhe-development/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chenyu97/github.com/openfhe-development /home/chenyu97/github.com/openfhe-development/src/binfhe /home/chenyu97/github.com/openfhe-development/build /home/chenyu97/github.com/openfhe-development/build/src/binfhe /home/chenyu97/github.com/openfhe-development/build/src/binfhe/CMakeFiles/boolean-multi-input.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/binfhe/CMakeFiles/boolean-multi-input.dir/depend

