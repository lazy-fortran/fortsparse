# Cross build of fortsparse for WebAssembly with LFortran and emscripten, used
# by tools/wasm/test_wasm.sh:
#   cmake -S . -B build-wasm -G Ninja -DCMAKE_BUILD_TYPE=Release \
#       -DCMAKE_TOOLCHAIN_FILE=tools/wasm/emscripten-lfortran.cmake
# Needs LFortran built with WITH_TARGET_WASM (LFORTRAN or lfortran on PATH),
# the emscripten SDK (EMSDK_PATH) and node (NODE or node on PATH).
# FORTSPARSE_WASM_ENVIRONMENT selects how the test programs run: node
# (default, under ctest, with the host file system and environment), shell (a
# JavaScript shell such as JavaScriptCore's jsc) or web (a browser page); the
# latter two are single files with the fixtures of bench/fixtures embedded.
set(CMAKE_SYSTEM_NAME Emscripten)
set(CMAKE_SYSTEM_PROCESSOR wasm32)

if(NOT DEFINED ENV{EMSDK_PATH})
    message(FATAL_ERROR "set EMSDK_PATH to the emscripten SDK")
endif()
set(_emsdk "$ENV{EMSDK_PATH}")
set(CMAKE_C_COMPILER "${_emsdk}/upstream/emscripten/emcc")
set(CMAKE_CXX_COMPILER "${_emsdk}/upstream/emscripten/em++")
set(CMAKE_AR "${_emsdk}/upstream/emscripten/emar" CACHE FILEPATH "")
set(CMAKE_RANLIB "${_emsdk}/upstream/emscripten/emranlib" CACHE FILEPATH "")
if(DEFINED ENV{LFORTRAN})
    set(CMAKE_Fortran_COMPILER "$ENV{LFORTRAN}")
else()
    find_program(CMAKE_Fortran_COMPILER lfortran REQUIRED)
endif()
get_filename_component(_lf_bin "${CMAKE_Fortran_COMPILER}" DIRECTORY)
set(_runtime "${_lf_bin}/../lib/lfortran_runtime_wasm_emcc.o")

set(CMAKE_Fortran_FLAGS_INIT
    "--target=wasm32-unknown-emscripten --separate-compilation --realloc-lhs-arrays --legacy-array-sections")
set(CMAKE_C_FLAGS_INIT "-msimd128")

set(FORTSPARSE_WASM_ENVIRONMENT node CACHE STRING "Where the test programs run: node, shell or web")
set(_link "-sALLOW_MEMORY_GROWTH -sSTACK_SIZE=8mb -sEXIT_RUNTIME --pre-js ${CMAKE_CURRENT_LIST_DIR}/environment.js")
if(FORTSPARSE_WASM_ENVIRONMENT STREQUAL "node")
    string(APPEND _link " -sENVIRONMENT=node -sNODERAWFS")
else()
    # the test fixtures a test opens relative to the source directory
    get_filename_component(_fixtures "${CMAKE_CURRENT_LIST_DIR}/../../bench/fixtures" ABSOLUTE)
    string(APPEND _link " -sENVIRONMENT=${FORTSPARSE_WASM_ENVIRONMENT} -sSINGLE_FILE --embed-file ${_fixtures}@bench/fixtures")
endif()
# emcc links the Fortran objects with LFortran's runtime.
set(CMAKE_Fortran_LINK_EXECUTABLE
    "${CMAKE_C_COMPILER} -O2 <LINK_FLAGS> <OBJECTS> -o <TARGET> <LINK_LIBRARIES> ${_runtime} ${_link}")
# C and C++ programs that call fortsparse (the interface tests) need the runtime too.
set(CMAKE_C_STANDARD_LIBRARIES_INIT "${_runtime} ${_link}")
set(CMAKE_CXX_STANDARD_LIBRARIES_INIT "${_runtime} ${_link}")
set(CMAKE_EXECUTABLE_SUFFIX ".js")
if(DEFINED ENV{NODE})
    set(CMAKE_CROSSCOMPILING_EMULATOR "$ENV{NODE}")
else()
    find_program(_node node REQUIRED)
    set(CMAKE_CROSSCOMPILING_EMULATOR "${_node}")
endif()
