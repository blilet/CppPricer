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
include CMakeFiles/testItoProcess.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/testItoProcess.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/testItoProcess.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/testItoProcess.dir/flags.make

CMakeFiles/testItoProcess.dir/tests/testItoProcess.cpp.o: CMakeFiles/testItoProcess.dir/flags.make
CMakeFiles/testItoProcess.dir/tests/testItoProcess.cpp.o: /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/tests/testItoProcess.cpp
CMakeFiles/testItoProcess.dir/tests/testItoProcess.cpp.o: CMakeFiles/testItoProcess.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/testItoProcess.dir/tests/testItoProcess.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testItoProcess.dir/tests/testItoProcess.cpp.o -MF CMakeFiles/testItoProcess.dir/tests/testItoProcess.cpp.o.d -o CMakeFiles/testItoProcess.dir/tests/testItoProcess.cpp.o -c /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/tests/testItoProcess.cpp

CMakeFiles/testItoProcess.dir/tests/testItoProcess.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testItoProcess.dir/tests/testItoProcess.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/tests/testItoProcess.cpp > CMakeFiles/testItoProcess.dir/tests/testItoProcess.cpp.i

CMakeFiles/testItoProcess.dir/tests/testItoProcess.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testItoProcess.dir/tests/testItoProcess.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/tests/testItoProcess.cpp -o CMakeFiles/testItoProcess.dir/tests/testItoProcess.cpp.s

CMakeFiles/testItoProcess.dir/src/Asset.cpp.o: CMakeFiles/testItoProcess.dir/flags.make
CMakeFiles/testItoProcess.dir/src/Asset.cpp.o: /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/Asset.cpp
CMakeFiles/testItoProcess.dir/src/Asset.cpp.o: CMakeFiles/testItoProcess.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/testItoProcess.dir/src/Asset.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testItoProcess.dir/src/Asset.cpp.o -MF CMakeFiles/testItoProcess.dir/src/Asset.cpp.o.d -o CMakeFiles/testItoProcess.dir/src/Asset.cpp.o -c /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/Asset.cpp

CMakeFiles/testItoProcess.dir/src/Asset.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testItoProcess.dir/src/Asset.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/Asset.cpp > CMakeFiles/testItoProcess.dir/src/Asset.cpp.i

CMakeFiles/testItoProcess.dir/src/Asset.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testItoProcess.dir/src/Asset.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/Asset.cpp -o CMakeFiles/testItoProcess.dir/src/Asset.cpp.s

CMakeFiles/testItoProcess.dir/src/ItoProcess.cpp.o: CMakeFiles/testItoProcess.dir/flags.make
CMakeFiles/testItoProcess.dir/src/ItoProcess.cpp.o: /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/ItoProcess.cpp
CMakeFiles/testItoProcess.dir/src/ItoProcess.cpp.o: CMakeFiles/testItoProcess.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/testItoProcess.dir/src/ItoProcess.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testItoProcess.dir/src/ItoProcess.cpp.o -MF CMakeFiles/testItoProcess.dir/src/ItoProcess.cpp.o.d -o CMakeFiles/testItoProcess.dir/src/ItoProcess.cpp.o -c /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/ItoProcess.cpp

CMakeFiles/testItoProcess.dir/src/ItoProcess.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testItoProcess.dir/src/ItoProcess.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/ItoProcess.cpp > CMakeFiles/testItoProcess.dir/src/ItoProcess.cpp.i

CMakeFiles/testItoProcess.dir/src/ItoProcess.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testItoProcess.dir/src/ItoProcess.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/ItoProcess.cpp -o CMakeFiles/testItoProcess.dir/src/ItoProcess.cpp.s

CMakeFiles/testItoProcess.dir/src/MeshUtils.cpp.o: CMakeFiles/testItoProcess.dir/flags.make
CMakeFiles/testItoProcess.dir/src/MeshUtils.cpp.o: /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/MeshUtils.cpp
CMakeFiles/testItoProcess.dir/src/MeshUtils.cpp.o: CMakeFiles/testItoProcess.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/testItoProcess.dir/src/MeshUtils.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testItoProcess.dir/src/MeshUtils.cpp.o -MF CMakeFiles/testItoProcess.dir/src/MeshUtils.cpp.o.d -o CMakeFiles/testItoProcess.dir/src/MeshUtils.cpp.o -c /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/MeshUtils.cpp

CMakeFiles/testItoProcess.dir/src/MeshUtils.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testItoProcess.dir/src/MeshUtils.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/MeshUtils.cpp > CMakeFiles/testItoProcess.dir/src/MeshUtils.cpp.i

CMakeFiles/testItoProcess.dir/src/MeshUtils.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testItoProcess.dir/src/MeshUtils.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/MeshUtils.cpp -o CMakeFiles/testItoProcess.dir/src/MeshUtils.cpp.s

CMakeFiles/testItoProcess.dir/src/Pricers.cpp.o: CMakeFiles/testItoProcess.dir/flags.make
CMakeFiles/testItoProcess.dir/src/Pricers.cpp.o: /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/Pricers.cpp
CMakeFiles/testItoProcess.dir/src/Pricers.cpp.o: CMakeFiles/testItoProcess.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/testItoProcess.dir/src/Pricers.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testItoProcess.dir/src/Pricers.cpp.o -MF CMakeFiles/testItoProcess.dir/src/Pricers.cpp.o.d -o CMakeFiles/testItoProcess.dir/src/Pricers.cpp.o -c /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/Pricers.cpp

CMakeFiles/testItoProcess.dir/src/Pricers.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testItoProcess.dir/src/Pricers.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/Pricers.cpp > CMakeFiles/testItoProcess.dir/src/Pricers.cpp.i

CMakeFiles/testItoProcess.dir/src/Pricers.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testItoProcess.dir/src/Pricers.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/Pricers.cpp -o CMakeFiles/testItoProcess.dir/src/Pricers.cpp.s

CMakeFiles/testItoProcess.dir/src/main.cpp.o: CMakeFiles/testItoProcess.dir/flags.make
CMakeFiles/testItoProcess.dir/src/main.cpp.o: /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/main.cpp
CMakeFiles/testItoProcess.dir/src/main.cpp.o: CMakeFiles/testItoProcess.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/testItoProcess.dir/src/main.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testItoProcess.dir/src/main.cpp.o -MF CMakeFiles/testItoProcess.dir/src/main.cpp.o.d -o CMakeFiles/testItoProcess.dir/src/main.cpp.o -c /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/main.cpp

CMakeFiles/testItoProcess.dir/src/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/testItoProcess.dir/src/main.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/main.cpp > CMakeFiles/testItoProcess.dir/src/main.cpp.i

CMakeFiles/testItoProcess.dir/src/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/testItoProcess.dir/src/main.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/src/main.cpp -o CMakeFiles/testItoProcess.dir/src/main.cpp.s

# Object files for target testItoProcess
testItoProcess_OBJECTS = \
"CMakeFiles/testItoProcess.dir/tests/testItoProcess.cpp.o" \
"CMakeFiles/testItoProcess.dir/src/Asset.cpp.o" \
"CMakeFiles/testItoProcess.dir/src/ItoProcess.cpp.o" \
"CMakeFiles/testItoProcess.dir/src/MeshUtils.cpp.o" \
"CMakeFiles/testItoProcess.dir/src/Pricers.cpp.o" \
"CMakeFiles/testItoProcess.dir/src/main.cpp.o"

# External object files for target testItoProcess
testItoProcess_EXTERNAL_OBJECTS =

testItoProcess: CMakeFiles/testItoProcess.dir/tests/testItoProcess.cpp.o
testItoProcess: CMakeFiles/testItoProcess.dir/src/Asset.cpp.o
testItoProcess: CMakeFiles/testItoProcess.dir/src/ItoProcess.cpp.o
testItoProcess: CMakeFiles/testItoProcess.dir/src/MeshUtils.cpp.o
testItoProcess: CMakeFiles/testItoProcess.dir/src/Pricers.cpp.o
testItoProcess: CMakeFiles/testItoProcess.dir/src/main.cpp.o
testItoProcess: CMakeFiles/testItoProcess.dir/build.make
testItoProcess: CMakeFiles/testItoProcess.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX executable testItoProcess"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/testItoProcess.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/testItoProcess.dir/build: testItoProcess
.PHONY : CMakeFiles/testItoProcess.dir/build

CMakeFiles/testItoProcess.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/testItoProcess.dir/cmake_clean.cmake
.PHONY : CMakeFiles/testItoProcess.dir/clean

CMakeFiles/testItoProcess.dir/depend:
	cd /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build /Users/timou/desktop/BlackScholesPDE/BlackScholesPDE/build/CMakeFiles/testItoProcess.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/testItoProcess.dir/depend

