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

# Utility rule file for benchmarkinfo.

# Include the progress variables for this target.
include benchmark/CMakeFiles/benchmarkinfo.dir/progress.make

benchmark/CMakeFiles/benchmarkinfo: benchmark/benchmarkinfocmd


benchmark/benchmarkinfocmd:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/chenyu97/github.com/openfhe-development/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating benchmarkinfocmd"
	cd /home/chenyu97/github.com/openfhe-development/build/benchmark && echo Builds benchmark lib and these apps: IntegerMath Lattice NbTheory VectorMath bfv-mult-method-benchmark binfhe-ap-pke binfhe-ap binfhe-ginx-pke binfhe-ginx binfhe-lmkcdey binfhe-paramsets compare-bfv-hps-leveled-vs-behz compare-bfvrns-vs-bgvrns fft-ckks-encoding lib-benchmark mult-vs-square poly-benchmark-16k poly-benchmark-1k poly-benchmark-4k poly-benchmark-64k poly-benchmark-8k serialize-ckks

benchmarkinfo: benchmark/CMakeFiles/benchmarkinfo
benchmarkinfo: benchmark/benchmarkinfocmd
benchmarkinfo: benchmark/CMakeFiles/benchmarkinfo.dir/build.make

.PHONY : benchmarkinfo

# Rule to build all files generated by this target.
benchmark/CMakeFiles/benchmarkinfo.dir/build: benchmarkinfo

.PHONY : benchmark/CMakeFiles/benchmarkinfo.dir/build

benchmark/CMakeFiles/benchmarkinfo.dir/clean:
	cd /home/chenyu97/github.com/openfhe-development/build/benchmark && $(CMAKE_COMMAND) -P CMakeFiles/benchmarkinfo.dir/cmake_clean.cmake
.PHONY : benchmark/CMakeFiles/benchmarkinfo.dir/clean

benchmark/CMakeFiles/benchmarkinfo.dir/depend:
	cd /home/chenyu97/github.com/openfhe-development/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/chenyu97/github.com/openfhe-development /home/chenyu97/github.com/openfhe-development/benchmark /home/chenyu97/github.com/openfhe-development/build /home/chenyu97/github.com/openfhe-development/build/benchmark /home/chenyu97/github.com/openfhe-development/build/benchmark/CMakeFiles/benchmarkinfo.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : benchmark/CMakeFiles/benchmarkinfo.dir/depend

