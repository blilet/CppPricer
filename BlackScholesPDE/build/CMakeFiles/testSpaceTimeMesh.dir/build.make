# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.30

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

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.30.5/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.30.5/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build

# Include any dependencies generated for this target.
include CMakeFiles/testSpaceTimeMesh.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/testSpaceTimeMesh.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/testSpaceTimeMesh.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/testSpaceTimeMesh.dir/flags.make

CMakeFiles/testSpaceTimeMesh.dir/tests/testSpaceTimeMesh.cpp.o: CMakeFiles/testSpaceTimeMesh.dir/flags.make
CMakeFiles/testSpaceTimeMesh.dir/tests/testSpaceTimeMesh.cpp.o: /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/tests/testSpaceTimeMesh.cpp
CMakeFiles/testSpaceTimeMesh.dir/tests/testSpaceTimeMesh.cpp.o: CMakeFiles/testSpaceTimeMesh.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/testSpaceTimeMesh.dir/tests/testSpaceTimeMesh.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testSpaceTimeMesh.dir/tests/testSpaceTimeMesh.cpp.o -MF CMakeFiles/testSpaceTimeMesh.dir/tests/testSpaceTimeMesh.cpp.o.d -o CMakeFiles/testSpaceTimeMesh.dir/tests/testSpaceTimeMesh.cpp.o -c /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/tests/testSpaceTimeMesh.cpp

CMakeFiles/testSpaceTimeMesh.dir/tests/testSpaceTimeMesh.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testSpaceTimeMesh.dir/tests/testSpaceTimeMesh.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/tests/testSpaceTimeMesh.cpp > CMakeFiles/testSpaceTimeMesh.dir/tests/testSpaceTimeMesh.cpp.i

CMakeFiles/testSpaceTimeMesh.dir/tests/testSpaceTimeMesh.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testSpaceTimeMesh.dir/tests/testSpaceTimeMesh.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/tests/testSpaceTimeMesh.cpp -o CMakeFiles/testSpaceTimeMesh.dir/tests/testSpaceTimeMesh.cpp.s

CMakeFiles/testSpaceTimeMesh.dir/src/Asset.cpp.o: CMakeFiles/testSpaceTimeMesh.dir/flags.make
CMakeFiles/testSpaceTimeMesh.dir/src/Asset.cpp.o: /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/Asset.cpp
CMakeFiles/testSpaceTimeMesh.dir/src/Asset.cpp.o: CMakeFiles/testSpaceTimeMesh.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/testSpaceTimeMesh.dir/src/Asset.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testSpaceTimeMesh.dir/src/Asset.cpp.o -MF CMakeFiles/testSpaceTimeMesh.dir/src/Asset.cpp.o.d -o CMakeFiles/testSpaceTimeMesh.dir/src/Asset.cpp.o -c /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/Asset.cpp

CMakeFiles/testSpaceTimeMesh.dir/src/Asset.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testSpaceTimeMesh.dir/src/Asset.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/Asset.cpp > CMakeFiles/testSpaceTimeMesh.dir/src/Asset.cpp.i

CMakeFiles/testSpaceTimeMesh.dir/src/Asset.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testSpaceTimeMesh.dir/src/Asset.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/Asset.cpp -o CMakeFiles/testSpaceTimeMesh.dir/src/Asset.cpp.s

CMakeFiles/testSpaceTimeMesh.dir/src/ItoProcess.cpp.o: CMakeFiles/testSpaceTimeMesh.dir/flags.make
CMakeFiles/testSpaceTimeMesh.dir/src/ItoProcess.cpp.o: /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/ItoProcess.cpp
CMakeFiles/testSpaceTimeMesh.dir/src/ItoProcess.cpp.o: CMakeFiles/testSpaceTimeMesh.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/testSpaceTimeMesh.dir/src/ItoProcess.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testSpaceTimeMesh.dir/src/ItoProcess.cpp.o -MF CMakeFiles/testSpaceTimeMesh.dir/src/ItoProcess.cpp.o.d -o CMakeFiles/testSpaceTimeMesh.dir/src/ItoProcess.cpp.o -c /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/ItoProcess.cpp

CMakeFiles/testSpaceTimeMesh.dir/src/ItoProcess.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testSpaceTimeMesh.dir/src/ItoProcess.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/ItoProcess.cpp > CMakeFiles/testSpaceTimeMesh.dir/src/ItoProcess.cpp.i

CMakeFiles/testSpaceTimeMesh.dir/src/ItoProcess.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testSpaceTimeMesh.dir/src/ItoProcess.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/ItoProcess.cpp -o CMakeFiles/testSpaceTimeMesh.dir/src/ItoProcess.cpp.s

CMakeFiles/testSpaceTimeMesh.dir/src/MeshUtils.cpp.o: CMakeFiles/testSpaceTimeMesh.dir/flags.make
CMakeFiles/testSpaceTimeMesh.dir/src/MeshUtils.cpp.o: /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/MeshUtils.cpp
CMakeFiles/testSpaceTimeMesh.dir/src/MeshUtils.cpp.o: CMakeFiles/testSpaceTimeMesh.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/testSpaceTimeMesh.dir/src/MeshUtils.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testSpaceTimeMesh.dir/src/MeshUtils.cpp.o -MF CMakeFiles/testSpaceTimeMesh.dir/src/MeshUtils.cpp.o.d -o CMakeFiles/testSpaceTimeMesh.dir/src/MeshUtils.cpp.o -c /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/MeshUtils.cpp

CMakeFiles/testSpaceTimeMesh.dir/src/MeshUtils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testSpaceTimeMesh.dir/src/MeshUtils.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/MeshUtils.cpp > CMakeFiles/testSpaceTimeMesh.dir/src/MeshUtils.cpp.i

CMakeFiles/testSpaceTimeMesh.dir/src/MeshUtils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testSpaceTimeMesh.dir/src/MeshUtils.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/MeshUtils.cpp -o CMakeFiles/testSpaceTimeMesh.dir/src/MeshUtils.cpp.s

CMakeFiles/testSpaceTimeMesh.dir/src/Pricers.cpp.o: CMakeFiles/testSpaceTimeMesh.dir/flags.make
CMakeFiles/testSpaceTimeMesh.dir/src/Pricers.cpp.o: /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/Pricers.cpp
CMakeFiles/testSpaceTimeMesh.dir/src/Pricers.cpp.o: CMakeFiles/testSpaceTimeMesh.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/testSpaceTimeMesh.dir/src/Pricers.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testSpaceTimeMesh.dir/src/Pricers.cpp.o -MF CMakeFiles/testSpaceTimeMesh.dir/src/Pricers.cpp.o.d -o CMakeFiles/testSpaceTimeMesh.dir/src/Pricers.cpp.o -c /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/Pricers.cpp

CMakeFiles/testSpaceTimeMesh.dir/src/Pricers.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testSpaceTimeMesh.dir/src/Pricers.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/Pricers.cpp > CMakeFiles/testSpaceTimeMesh.dir/src/Pricers.cpp.i

CMakeFiles/testSpaceTimeMesh.dir/src/Pricers.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testSpaceTimeMesh.dir/src/Pricers.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/Pricers.cpp -o CMakeFiles/testSpaceTimeMesh.dir/src/Pricers.cpp.s

CMakeFiles/testSpaceTimeMesh.dir/src/main.cpp.o: CMakeFiles/testSpaceTimeMesh.dir/flags.make
CMakeFiles/testSpaceTimeMesh.dir/src/main.cpp.o: /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/main.cpp
CMakeFiles/testSpaceTimeMesh.dir/src/main.cpp.o: CMakeFiles/testSpaceTimeMesh.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/testSpaceTimeMesh.dir/src/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testSpaceTimeMesh.dir/src/main.cpp.o -MF CMakeFiles/testSpaceTimeMesh.dir/src/main.cpp.o.d -o CMakeFiles/testSpaceTimeMesh.dir/src/main.cpp.o -c /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/main.cpp

CMakeFiles/testSpaceTimeMesh.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testSpaceTimeMesh.dir/src/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/main.cpp > CMakeFiles/testSpaceTimeMesh.dir/src/main.cpp.i

CMakeFiles/testSpaceTimeMesh.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testSpaceTimeMesh.dir/src/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/main.cpp -o CMakeFiles/testSpaceTimeMesh.dir/src/main.cpp.s

# Object files for target testSpaceTimeMesh
testSpaceTimeMesh_OBJECTS = \
"CMakeFiles/testSpaceTimeMesh.dir/tests/testSpaceTimeMesh.cpp.o" \
"CMakeFiles/testSpaceTimeMesh.dir/src/Asset.cpp.o" \
"CMakeFiles/testSpaceTimeMesh.dir/src/ItoProcess.cpp.o" \
"CMakeFiles/testSpaceTimeMesh.dir/src/MeshUtils.cpp.o" \
"CMakeFiles/testSpaceTimeMesh.dir/src/Pricers.cpp.o" \
"CMakeFiles/testSpaceTimeMesh.dir/src/main.cpp.o"

# External object files for target testSpaceTimeMesh
testSpaceTimeMesh_EXTERNAL_OBJECTS =

testSpaceTimeMesh: CMakeFiles/testSpaceTimeMesh.dir/tests/testSpaceTimeMesh.cpp.o
testSpaceTimeMesh: CMakeFiles/testSpaceTimeMesh.dir/src/Asset.cpp.o
testSpaceTimeMesh: CMakeFiles/testSpaceTimeMesh.dir/src/ItoProcess.cpp.o
testSpaceTimeMesh: CMakeFiles/testSpaceTimeMesh.dir/src/MeshUtils.cpp.o
testSpaceTimeMesh: CMakeFiles/testSpaceTimeMesh.dir/src/Pricers.cpp.o
testSpaceTimeMesh: CMakeFiles/testSpaceTimeMesh.dir/src/main.cpp.o
testSpaceTimeMesh: CMakeFiles/testSpaceTimeMesh.dir/build.make
testSpaceTimeMesh: CMakeFiles/testSpaceTimeMesh.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX executable testSpaceTimeMesh"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/testSpaceTimeMesh.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/testSpaceTimeMesh.dir/build: testSpaceTimeMesh
.PHONY : CMakeFiles/testSpaceTimeMesh.dir/build

CMakeFiles/testSpaceTimeMesh.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/testSpaceTimeMesh.dir/cmake_clean.cmake
.PHONY : CMakeFiles/testSpaceTimeMesh.dir/clean

CMakeFiles/testSpaceTimeMesh.dir/depend:
	cd /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build/CMakeFiles/testSpaceTimeMesh.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/testSpaceTimeMesh.dir/depend

