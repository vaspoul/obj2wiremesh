-- all paths are relative to this script

workspace("obj2mesh")
	configurations( { "Debug", "Release" } )
	platforms({"Win64"})

targetdir("BUILD/binaries/%{cfg.buildcfg}")
objdir("BUILD/intermediate")
startproject("obj2mesh")
exceptionhandling ("Off")

filter "system:Windows"
	systemversion "latest"
	  
filter({"configurations:Debug*"})
	defines({ "DEBUG", "CONFIG_DEBUG" })
	optimize("Off")
	symbols("On")
	runtime("Release")
	staticruntime("off")

filter( {"configurations:Release*"} )
	defines({ "NDEBUG", "CONFIG_RELEASE" })
	optimize("Speed")
	symbols("On")
	runtime("Release")
	staticruntime("off")

filter { "platforms:Win64" }
    system "Windows"
    architecture "x86_64"
	defines({"CONFIG_WINDOWS", "CONFIG_WIN64", "CONFIG_DX12"})

-- deactivate previous filter
filter {}

flags({"FatalCompileWarnings", "MultiProcessorCompile"})
warnings("Extra")
disablewarnings ( 
{
	"4100", -- unreferenced formal parameter
	"4101", -- unreferenced local variable
	"4189", -- local variable is initialized but not referenced
	"4201", -- nonstandard extension used: nameless struct/union
	"4505", -- unreferenced local function has been removed 
	"4251", -- 'type' : class 'type1' needs to have dll-interface to be used by clients of class 'type2'
	"4127", -- conditional expression is constant
	"26812", -- warning about unscoped enums
	"4324", -- 'struct_name' : structure was padded due to __declspec(align())
} )

fatalwarnings  (
{
	"4296", -- 'operator' : expression is always false
})

project("obj2mesh")
	location("BUILD/projects")
	kind("WindowedApp")
	language("C++")
	cppdialect("C++17")
	characterset("Unicode")
	
	links({"d3d12.lib", "dxgi.lib", "d3dcompiler.lib", "xinput.lib"})

	debugdir(".")

	filter( {"platforms:Win64"} )

	filter {}

	files({ "**.h", "**.cpp", "**.resx" })
	
	links( {"comctl32.lib"} )

project("Premake")
	kind("None")
	location("BUILD/projects")
	files({ "*.lua" })

