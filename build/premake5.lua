solution "Pneumo-ABM"
    configurations { "Debug", "Release" }
    platforms { "x64" }
    location ("./" .. _ACTION)
    configuration "not windows"
        buildoptions { 
            "-std=c++11",
            "-Wno-unknown-pragmas"
        }
    configuration "vs*"
        includedirs { "$(boost_prefix)" }
    configuration "Debug"
        flags { "Symbols" }
	optimize "Off"
    configuration "Release"
        optimize "Full"

project "Pneumo-ABM"
    kind "ConsoleApp"
    language "C++"
    targetname "pneumo-abm"
    files { 
       "../source/*.h",
       "../source/*.cpp"
    }
    flags { 
       "Unicode",
       "NoEditAndContinue",
       "NoManifest",
       "NoPCH"
    }
    debugargs { "inputs=../inputs outputs=../outputs" }
    configuration "Debug"
	targetdir "../bin/debug"
    configuration "Release"
        flags { "LinkTimeOptimization" }
	targetdir "../bin/release"
