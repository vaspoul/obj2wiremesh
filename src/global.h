#pragma once

#if defined(CONFIG_WINDOWS)

#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <windowsx.h>
#include <shellapi.h>
#include <oleidl.h>
#include <processthreadsapi.h>
#include <commdlg.h>
#endif

#include <stdint.h>
#include <cstdlib>
#include <vector>
#include <array>
#include <list>
#include <map>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <string_view>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <initializer_list>
#include <functional>
#include <stdlib.h>
#include <charconv>
#include <type_traits>
#include <time.h>
#include <assert.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include "hidusage.h"

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

#if !defined(_M_CEE)
// Some std stuff not supported when compiling with /clr
#include <thread>
#include <mutex>
#include <shared_mutex>
#include <future>
#endif

#include <queue>
#include <stack>
#include <atomic>
#include <cstdio>
#include <filesystem>
#include <random>
#include <typeinfo>
#include <typeindex>
#include <bitset>
#include <filesystem>
#include <chrono>
#include <regex>
#include <stdint.h>

using namespace std::chrono_literals;

#include <immintrin.h> // platforms?
//#include <intrin.h> // MSVC specific?

#include <d3d12.h>
#include <dxgi1_4.h>
#include <D3Dcompiler.h>
#include <Xinput.h>

#include <wrl/client.h>
using namespace Microsoft::WRL;

// Am I going to regret this?
typedef uint8_t ubyte;
typedef int8_t sbyte;

#include "float2.h"
#include "float3.h"
#include "float4.h"
#include "matrix44.h"
#include "ivector.h"
#include "vectorOps.h"
#include "matrixOps.h"

namespace std
{
	namespace fs = filesystem;

	template<class Key, class Value>
	using umap = unordered_map<Key, Value>;

	template<class Value>
	using uset = unordered_set<Value>;
}

#define SAFE_RELEASE(x)			{ if (x != 0) { x->Release(); x=0; } }
#define SAFE_DELETE(x)			{ if (x != 0) { delete x; x=0; } }
#define SAFE_FREE(x)			{ if (x != 0) { free(x); x=0; } }
#define ZeroMem(x)				memset(&x, 0, sizeof(x))


#define STATIC_ASSERT(e) typedef char __CompileTimeASSERT__##__LINE__[(e)?1:-1]

#define ALIGN(value, alignmentSize) (((value) + (decltype(value))(alignmentSize - 1)) & (~(decltype(value))(alignmentSize - 1)))

#ifndef ARRAYSIZE
#define ARRAYSIZE(x) sizeof((x))/sizeof(x[0])
#endif

#define V(x)		{ HRESULT hr = (x); assert(hr == S_OK); }
