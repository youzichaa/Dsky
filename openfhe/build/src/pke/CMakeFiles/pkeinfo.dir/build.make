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

# Utility rule file for pkeinfo.

# Include the progress variables for this target.
include src/pke/CMakeFiles/pkeinfo.dir/progress.make

src/pke/CMakeFiles/pkeinfo: src/pke/pkeinfocmd


src/pke/pkeinfocmd:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/chenyu97/github.com/openfhe-development/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating pkeinfocmd"
	cd /home/chenyu97/github.com/openfhe-development/build/src/pke && echo Builds OPENFHEpke and these apps: advanced-ckks-bootstrapping advanced-real-numbers-128 advanced-real-numbers ckks-noise-flooding depth-bfvrns-behz depth-bfvrns depth-bgvrns function-evaluation inner-product iterative-ckks-bootstrapping linearwsum-evaluation polynomial-evaluation pre-buffer pre-hra-secure rotation scheme-switching-serial scheme-switching simple-ckks-bootstrapping simple-integers-bgvrns simple-integers-serial-bgvrns simple-integers-serial simple-integers simple-real-numbers-serial simple-real-numbers skyline tckks-interactive-mp-bootstrapping-Chebyshev tckks-interactive-mp-bootstrapping threshold-fhe-5p threshold-fhe

pkeinfo: src/pke/CMakeFiles/pkeinfo
pkeinfo: src/pke/pkeinfocmd
pkeinfo: src/pke/CMakeFiles/pkeinfo.dir/build.make

.PHONY : pkeinfo

# Rule to build all files generated by this target.
src/pke/CMakeFiles/pkeinfo.dir/build: pkeinfo

.PHONY : src/pke/CMakeFiles/pkeinfo.dir/build

src/pke/CMakeFiles/pkeinfo.dir/clean:
	cd /home/chenyu97/github.com/openfhe-development/build/src/pke && $(CMAKE_COMMAND) -P CMakeFiles/pkeinfo.dir/cmake_clean.cmake
.PHONY : src/pke/CMakeFiles/pkeinfo.dir/clean

src/pke/CMakeFiles/pkeinfo.dir/depend:
	cd /home/chenyu97/github.com/openfhe-development/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chenyu97/github.com/openfhe-development /home/chenyu97/github.com/openfhe-development/src/pke /home/chenyu97/github.com/openfhe-development/build /home/chenyu97/github.com/openfhe-development/build/src/pke /home/chenyu97/github.com/openfhe-development/build/src/pke/CMakeFiles/pkeinfo.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : src/pke/CMakeFiles/pkeinfo.dir/depend
