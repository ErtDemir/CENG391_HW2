# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.20

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2021.2.3\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2021.2.3\bin\cmake\win\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "C:\Users\Ertugrul Demir\Desktop\libceng391"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "C:\Users\Ertugrul Demir\Desktop\libceng391\cmake-build-debug"

# Include any dependencies generated for this target.
include app/CMakeFiles/image-test.dir/depend.make
# Include the progress variables for this target.
include app/CMakeFiles/image-test.dir/progress.make

# Include the compile flags for this target's objects.
include app/CMakeFiles/image-test.dir/flags.make

app/CMakeFiles/image-test.dir/image_test.cpp.obj: app/CMakeFiles/image-test.dir/flags.make
app/CMakeFiles/image-test.dir/image_test.cpp.obj: app/CMakeFiles/image-test.dir/includes_CXX.rsp
app/CMakeFiles/image-test.dir/image_test.cpp.obj: ../app/image_test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="C:\Users\Ertugrul Demir\Desktop\libceng391\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object app/CMakeFiles/image-test.dir/image_test.cpp.obj"
	cd /d C:\Users\ERTUGR~1\Desktop\LIBCEN~1\CMAKE-~1\app && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\image-test.dir\image_test.cpp.obj -c "C:\Users\Ertugrul Demir\Desktop\libceng391\app\image_test.cpp"

app/CMakeFiles/image-test.dir/image_test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/image-test.dir/image_test.cpp.i"
	cd /d C:\Users\ERTUGR~1\Desktop\LIBCEN~1\CMAKE-~1\app && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "C:\Users\Ertugrul Demir\Desktop\libceng391\app\image_test.cpp" > CMakeFiles\image-test.dir\image_test.cpp.i

app/CMakeFiles/image-test.dir/image_test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/image-test.dir/image_test.cpp.s"
	cd /d C:\Users\ERTUGR~1\Desktop\LIBCEN~1\CMAKE-~1\app && C:\MinGW\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "C:\Users\Ertugrul Demir\Desktop\libceng391\app\image_test.cpp" -o CMakeFiles\image-test.dir\image_test.cpp.s

# Object files for target image-test
image__test_OBJECTS = \
"CMakeFiles/image-test.dir/image_test.cpp.obj"

# External object files for target image-test
image__test_EXTERNAL_OBJECTS =

app/image-test.exe: app/CMakeFiles/image-test.dir/image_test.cpp.obj
app/image-test.exe: app/CMakeFiles/image-test.dir/build.make
app/image-test.exe: src/libceng391.a
app/image-test.exe: app/CMakeFiles/image-test.dir/linklibs.rsp
app/image-test.exe: app/CMakeFiles/image-test.dir/objects1.rsp
app/image-test.exe: app/CMakeFiles/image-test.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="C:\Users\Ertugrul Demir\Desktop\libceng391\cmake-build-debug\CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable image-test.exe"
	cd /d C:\Users\ERTUGR~1\Desktop\LIBCEN~1\CMAKE-~1\app && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\image-test.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
app/CMakeFiles/image-test.dir/build: app/image-test.exe
.PHONY : app/CMakeFiles/image-test.dir/build

app/CMakeFiles/image-test.dir/clean:
	cd /d C:\Users\ERTUGR~1\Desktop\LIBCEN~1\CMAKE-~1\app && $(CMAKE_COMMAND) -P CMakeFiles\image-test.dir\cmake_clean.cmake
.PHONY : app/CMakeFiles/image-test.dir/clean

app/CMakeFiles/image-test.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" "C:\Users\Ertugrul Demir\Desktop\libceng391" "C:\Users\Ertugrul Demir\Desktop\libceng391\app" "C:\Users\Ertugrul Demir\Desktop\libceng391\cmake-build-debug" "C:\Users\Ertugrul Demir\Desktop\libceng391\cmake-build-debug\app" "C:\Users\Ertugrul Demir\Desktop\libceng391\cmake-build-debug\app\CMakeFiles\image-test.dir\DependInfo.cmake" --color=$(COLOR)
.PHONY : app/CMakeFiles/image-test.dir/depend

