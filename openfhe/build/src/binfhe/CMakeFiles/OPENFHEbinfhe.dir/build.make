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
include src/binfhe/CMakeFiles/OPENFHEbinfhe.dir/depend.make

# Include the progress variables for this target.
include src/binfhe/CMakeFiles/OPENFHEbinfhe.dir/progress.make

# Include the compile flags for this target's objects.
include src/binfhe/CMakeFiles/OPENFHEbinfhe.dir/flags.make

# Object files for target OPENFHEbinfhe
OPENFHEbinfhe_OBJECTS =

# External object files for target OPENFHEbinfhe
OPENFHEbinfhe_EXTERNAL_OBJECTS = \
"/home/chenyu97/github.com/openfhe-development/build/src/binfhe/CMakeFiles/binfheobj.dir/lib/binfhe-base-scheme.cpp.o" \
"/home/chenyu97/github.com/openfhe-development/build/src/binfhe/CMakeFiles/binfheobj.dir/lib/binfhe-constants-impl.cpp.o" \
"/home/chenyu97/github.com/openfhe-development/build/src/binfhe/CMakeFiles/binfheobj.dir/lib/binfhecontext.cpp.o" \
"/home/chenyu97/github.com/openfhe-development/build/src/binfhe/CMakeFiles/binfheobj.dir/lib/lwe-pke.cpp.o" \
"/home/chenyu97/github.com/openfhe-development/build/src/binfhe/CMakeFiles/binfheobj.dir/lib/rgsw-acc-cggi.cpp.o" \
"/home/chenyu97/github.com/openfhe-development/build/src/binfhe/CMakeFiles/binfheobj.dir/lib/rgsw-acc-dm.cpp.o" \
"/home/chenyu97/github.com/openfhe-development/build/src/binfhe/CMakeFiles/binfheobj.dir/lib/rgsw-acc-lmkcdey.cpp.o" \
"/home/chenyu97/github.com/openfhe-development/build/src/binfhe/CMakeFiles/binfheobj.dir/lib/rgsw-acc.cpp.o" \
"/home/chenyu97/github.com/openfhe-development/build/src/binfhe/CMakeFiles/binfheobj.dir/lib/rgsw-cryptoparameters.cpp.o"

lib/libOPENFHEbinfhe.so.1.2.1: src/binfhe/CMakeFiles/binfheobj.dir/lib/binfhe-base-scheme.cpp.o
lib/libOPENFHEbinfhe.so.1.2.1: src/binfhe/CMakeFiles/binfheobj.dir/lib/binfhe-constants-impl.cpp.o
lib/libOPENFHEbinfhe.so.1.2.1: src/binfhe/CMakeFiles/binfheobj.dir/lib/binfhecontext.cpp.o
lib/libOPENFHEbinfhe.so.1.2.1: src/binfhe/CMakeFiles/binfheobj.dir/lib/lwe-pke.cpp.o
lib/libOPENFHEbinfhe.so.1.2.1: src/binfhe/CMakeFiles/binfheobj.dir/lib/rgsw-acc-cggi.cpp.o
lib/libOPENFHEbinfhe.so.1.2.1: src/binfhe/CMakeFiles/binfheobj.dir/lib/rgsw-acc-dm.cpp.o
lib/libOPENFHEbinfhe.so.1.2.1: src/binfhe/CMakeFiles/binfheobj.dir/lib/rgsw-acc-lmkcdey.cpp.o
lib/libOPENFHEbinfhe.so.1.2.1: src/binfhe/CMakeFiles/binfheobj.dir/lib/rgsw-acc.cpp.o
lib/libOPENFHEbinfhe.so.1.2.1: src/binfhe/CMakeFiles/binfheobj.dir/lib/rgsw-cryptoparameters.cpp.o
lib/libOPENFHEbinfhe.so.1.2.1: src/binfhe/CMakeFiles/OPENFHEbinfhe.dir/build.make
lib/libOPENFHEbinfhe.so.1.2.1: lib/libOPENFHEcore.so.1.2.1
lib/libOPENFHEbinfhe.so.1.2.1: src/binfhe/CMakeFiles/OPENFHEbinfhe.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/chenyu97/github.com/openfhe-development/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Linking CXX shared library ../../lib/libOPENFHEbinfhe.so"
	cd /home/chenyu97/github.com/openfhe-development/build/src/binfhe && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/OPENFHEbinfhe.dir/link.txt --verbose=$(VERBOSE)
	cd /home/chenyu97/github.com/openfhe-development/build/src/binfhe && $(CMAKE_COMMAND) -E cmake_symlink_library ../../lib/libOPENFHEbinfhe.so.1.2.1 ../../lib/libOPENFHEbinfhe.so.1 ../../lib/libOPENFHEbinfhe.so

lib/libOPENFHEbinfhe.so.1: lib/libOPENFHEbinfhe.so.1.2.1
	@$(CMAKE_COMMAND) -E touch_nocreate lib/libOPENFHEbinfhe.so.1

lib/libOPENFHEbinfhe.so: lib/libOPENFHEbinfhe.so.1.2.1
	@$(CMAKE_COMMAND) -E touch_nocreate lib/libOPENFHEbinfhe.so

# Rule to build all files generated by this target.
src/binfhe/CMakeFiles/OPENFHEbinfhe.dir/build: lib/libOPENFHEbinfhe.so

.PHONY : src/binfhe/CMakeFiles/OPENFHEbinfhe.dir/build

src/binfhe/CMakeFiles/OPENFHEbinfhe.dir/clean:
	cd /home/chenyu97/github.com/openfhe-development/build/src/binfhe && $(CMAKE_COMMAND) -P CMakeFiles/OPENFHEbinfhe.dir/cmake_clean.cmake
.PHONY : src/binfhe/CMakeFiles/OPENFHEbinfhe.dir/clean

src/binfhe/CMakeFiles/OPENFHEbinfhe.dir/depend:
	cd /home/chenyu97/github.com/openfhe-development/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chenyu97/github.com/openfhe-development /home/chenyu97/github.com/openfhe-development/src/binfhe /home/chenyu97/github.com/openfhe-development/build /home/chenyu97/github.com/openfhe-development/build/src/binfhe /home/chenyu97/github.com/openfhe-development/build/src/binfhe/CMakeFiles/OPENFHEbinfhe.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/binfhe/CMakeFiles/OPENFHEbinfhe.dir/depend

