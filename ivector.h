#pragma once

template<class T, unsigned int count>
struct TVector
{};

struct float2;
struct float3;
struct float4;

template<class T>
struct TVector<T, 2>
{
	static const TVector<T, 2> k00;
	static const TVector<T, 2> kFF;
	static const TVector<T, 2> kF0;
	static const TVector<T, 2> k0F;

	TVector()
		: x(0), y(0)
	{}

	TVector(const T& _x, const T& _y)
		: x(_x), y(_y)
	{}

	TVector(const float2& v);

	bool operator==(const TVector& other) const
	{
		return x == other.x && y == other.y;
	}

	bool operator!=(const TVector& other) const
	{
		return x != other.x || y != other.y;
	}

	void operator+=(const TVector& other)
	{
		x += other.x;
		y += other.y;
	}

	TVector operator+(const TVector& other) const
	{
		TVector result;
		result.x = x + other.x;
		result.y = y + other.y;
		return result;
	}

	void operator-=(const TVector& other)
	{
		x -= other.x;
		y -= other.y;
	}

	TVector operator-(const TVector& other) const
	{
		TVector result;
		result.x = x - other.x;
		result.y = y - other.y;
		return result;
	}

	void operator*=(const TVector& other)
	{
		x *= other.x;
		y *= other.y;
	}

	void operator*=(T value)
	{
		x *= value;
		y *= value;
	}

	void operator/=(T value)
	{
		x /= value;
		y /= value;
	}

	T& operator[](uint32_t index)
	{
		return values[index];
	}

	const T& operator[](uint32_t index) const
	{
		return values[index];
	}

	union
	{
		T values[2];
		struct { T x, y; };
		struct { T r, g; };
	};
};

template<typename Archive, class T>
Archive& operator<<(Archive& archive, TVector<T, 2>& v)
{
	archive << v.values;

	return archive;
}

template<class T> inline const TVector<T, 2> TVector<T, 2>::k00 = TVector<T, 2>(0,0);
template<class T> inline const TVector<T, 2> TVector<T, 2>::kFF = TVector<T, 2>(std::numeric_limits<T>::max(),std::numeric_limits<T>::max());
template<class T> inline const TVector<T, 2> TVector<T, 2>::kF0 = TVector<T, 2>(std::numeric_limits<T>::max(),0);
template<class T> inline const TVector<T, 2> TVector<T, 2>::k0F = TVector<T, 2>(0,std::numeric_limits<T>::max());

template<class T>
struct TVector<T, 3>
{
	static const TVector<T, 3> k000;
	static const TVector<T, 3> kFFF;
	static const TVector<T, 3> kF00;
	static const TVector<T, 3> k0F0;
	static const TVector<T, 3> k00F;
	static const TVector<T, 3> kFF0;
	static const TVector<T, 3> kF0F;
	static const TVector<T, 3> k0FF;

	TVector()
		: x(0), y(0), z(0)
	{}

	TVector(const TVector<T, 2>& xy, const T& z)
		: x(xy.x), y(xy.y), z(z)
	{}

	TVector(const T& _x, const T& _y, const T& _z)
		: x(_x), y(_y), z(_z)
	{}

	TVector(const float3& v);

	bool operator==(const TVector& other) const
	{
		return x == other.x && y == other.y && z == other.z;
	}

	bool operator!=(const TVector& other) const
	{
		return x != other.x || y != other.y || z != other.z;
	}

	void operator+=(const TVector& other)
	{
		x += other.x;
		y += other.y;
		z += other.z;
	}

	TVector operator+(const TVector& other) const
	{
		TVector result;
		result.x = x + other.x;
		result.y = y + other.y;
		result.z = z + other.z;
		return result;
	}

	void operator-=(const TVector& other)
	{
		x -= other.x;
		y -= other.y;
		z -= other.z;
	}

	TVector operator-(const TVector& other) const
	{
		TVector result;
		result.x = x - other.x;
		result.y = y - other.y;
		result.z = z - other.z;
		return result;
	}

	void operator*=(const TVector& other)
	{
		x *= other.x;
		y *= other.y;
		z *= other.z;
	}

	void operator*=(T value)
	{
		x *= value;
		y *= value;
		w *= value;
	}

	void operator/=(T value)
	{
		x /= value;
		y /= value;
		z /= value;
	}

	T& operator[](uint32_t index)
	{
		return values[index];
	}

	const T& operator[](uint32_t index) const
	{
		return values[index];
	}

	union
	{
		T values[3];
		struct { T x, y, z; };
		struct { T r, g, b; };
	};
};

template<typename Archive, class T>
Archive& operator<<(Archive& archive, TVector<T, 3>& v)
{
	archive << v.values;

	return archive;
}

#define X std::numeric_limits<T>::max()
template<class T> inline const TVector<T, 3> TVector<T, 3>::k000 = TVector<T, 3>(0,0,0);
template<class T> inline const TVector<T, 3> TVector<T, 3>::kFFF = TVector<T, 3>(X,X,X);
template<class T> inline const TVector<T, 3> TVector<T, 3>::kF00 = TVector<T, 3>(X,0,0);
template<class T> inline const TVector<T, 3> TVector<T, 3>::k0F0 = TVector<T, 3>(0,X,0);
template<class T> inline const TVector<T, 3> TVector<T, 3>::k00F = TVector<T, 3>(0,0,X);
template<class T> inline const TVector<T, 3> TVector<T, 3>::kFF0 = TVector<T, 3>(X,X,0);
template<class T> inline const TVector<T, 3> TVector<T, 3>::kF0F = TVector<T, 3>(X,0,X);
template<class T> inline const TVector<T, 3> TVector<T, 3>::k0FF = TVector<T, 3>(0,X,X);
#undef X

template<class T>
struct TVector<T, 4>
{
	static const TVector<T, 4> k0000;
	static const TVector<T, 4> kFFFF;
	
	static const TVector<T, 4> kF000;
	static const TVector<T, 4> k0F00;
	static const TVector<T, 4> k00F0;
	static const TVector<T, 4> k000F;

	static const TVector<T, 4> kFF00;
	static const TVector<T, 4> kF0F0;
	static const TVector<T, 4> kF00F;
	static const TVector<T, 4> k0FF0;
	static const TVector<T, 4> k0F0F;
	static const TVector<T, 4> k00FF;

	static const TVector<T, 4> kFFF0;
	static const TVector<T, 4> k0FFF;
	static const TVector<T, 4> kF0FF;
	static const TVector<T, 4> kFF0F;


	TVector()
		: x(0), y(0), z(0), w(0)
	{}

	TVector(const T& _x, const T& _y, const T& _z, const T& _w)
		: x(_x), y(_y), z(_z), w(_w)
	{}

	TVector(const TVector<T, 2>& xy, const TVector<T, 2>& zw)
		: x(xy.x), y(xy.y), z(zw.x), w(zw.y)
	{}

	TVector(const TVector<T, 3>& xyz, const T& w)
		: x(xyz.x), y(xyz.y), z(xyz.z), w(w)
	{}

	TVector(const float4& v);

	bool operator==(const TVector& other) const
	{
		return x == other.x && y == other.y && z == other.z && w == other.w;
	}

	bool operator!=(const TVector& other) const
	{
		return x != other.x || y != other.y || z != other.z || w != other.w;
	}

	void operator+=(const TVector& other)
	{
		x += other.x;
		y += other.y;
		z += other.z;
		w += other.w;
	}

	TVector operator+(const TVector& other) const
	{
		TVector result;
		result.x = x + other.x;
		result.y = y + other.y;
		result.z = z + other.z;
		result.w = w + other.w;
		return result;
	}

	void operator-=(const TVector& other)
	{
		x -= other.x;
		y -= other.y;
		z -= other.z;
		w -= other.w;
	}

	TVector operator-(const TVector& other) const
	{
		TVector result;
		result.x = x - other.x;
		result.y = y - other.y;
		result.z = z - other.z;
		result.w = w - other.w;
		return result;
	}

	void operator*=(const TVector& other)
	{
		x *= other.x;
		y *= other.y;
		z *= other.z;
		w *= other.w;
	}

	void operator*=(T value)
	{
		x *= value;
		y *= value;
		z *= value;
		w *= value;
	}

	void operator/=(T value)
	{
		x /= value;
		y /= value;
		z /= value;
		w /= value;
	}

	T& operator[](uint32_t index)
	{
		return values[index];
	}

	const T& operator[](uint32_t index) const
	{
		return values[index];
	}

	union
	{
		T values[4];
		struct { T x, y, z, w; };
		struct { T r, g, b, a; };
		struct  
		{
			TVector<T,3> xyz;
			T w;
		};

		struct  
		{
			TVector<T,2> xy;
			TVector<T,2> zw;
		};
	};
};

template<typename Archive, class T>
Archive& operator<<(Archive& archive, TVector<T, 4>& v)
{
	archive << v.values;

	return archive;
}

#define X std::numeric_limits<T>::max()
template<class T> inline const TVector<T, 4> TVector<T, 4>::k0000 = TVector<T, 4>(0,0,0,0);
template<class T> inline const TVector<T, 4> TVector<T, 4>::kFFFF = TVector<T, 4>(X,X,X,X);
template<class T> inline const TVector<T, 4> TVector<T, 4>::kF000 = TVector<T, 4>(X,0,0,0);
template<class T> inline const TVector<T, 4> TVector<T, 4>::k0F00 = TVector<T, 4>(0,X,0,0);
template<class T> inline const TVector<T, 4> TVector<T, 4>::k00F0 = TVector<T, 4>(0,0,X,0);
template<class T> inline const TVector<T, 4> TVector<T, 4>::k000F = TVector<T, 4>(0,0,0,X);
template<class T> inline const TVector<T, 4> TVector<T, 4>::kFF00 = TVector<T, 4>(X,X,0,0);
template<class T> inline const TVector<T, 4> TVector<T, 4>::kF0F0 = TVector<T, 4>(X,0,X,0);
template<class T> inline const TVector<T, 4> TVector<T, 4>::kF00F = TVector<T, 4>(X,0,0,X);
template<class T> inline const TVector<T, 4> TVector<T, 4>::k0FF0 = TVector<T, 4>(0,X,X,0);
template<class T> inline const TVector<T, 4> TVector<T, 4>::k0F0F = TVector<T, 4>(0,X,0,X);
template<class T> inline const TVector<T, 4> TVector<T, 4>::k00FF = TVector<T, 4>(0,0,X,X);
template<class T> inline const TVector<T, 4> TVector<T, 4>::kFFF0 = TVector<T, 4>(X,X,X,0);
template<class T> inline const TVector<T, 4> TVector<T, 4>::k0FFF = TVector<T, 4>(0,X,X,X);
template<class T> inline const TVector<T, 4> TVector<T, 4>::kF0FF = TVector<T, 4>(X,0,X,X);
template<class T> inline const TVector<T, 4> TVector<T, 4>::kFF0F = TVector<T, 4>(X,X,0,X);
#undef X


typedef uint8_t					ubyte1;
typedef TVector<uint8_t, 2>		ubyte2;
typedef TVector<uint8_t, 3>		ubyte3;
typedef TVector<uint8_t, 4>		ubyte4;

namespace Colors
{
	static const ubyte4 White	( 255, 255, 255, 255);
	static const ubyte4 Black	(   0,   0,   0, 255);
	static const ubyte4 Gray	( 128, 128,	128, 255);
	static const ubyte4 Red		( 255,   0,   0, 255);
	static const ubyte4 Green	(   0, 255,   0, 255);
	static const ubyte4 Blue	(   0,   0, 255, 255);
	static const ubyte4 Maroon	( 128,   0,   0, 255);
	static const ubyte4 Navy	(   0,   0, 128, 255);
	static const ubyte4 Teal	(   0, 128, 128, 255);
	static const ubyte4 Olive	( 128, 128,   0, 255);
	static const ubyte4 Purple	( 128,   0, 128, 255);
	static const ubyte4 Cyan	(   0, 255, 255, 255);
	static const ubyte4 Yellow	( 255, 255,   0, 255);
	static const ubyte4 Orange	( 255, 128,   0, 255);
}


//typedef int8_t					byte;
typedef int8_t					sbyte1;
typedef TVector<int8_t, 2>		sbyte2;
typedef TVector<int8_t, 3>		sbyte3;
typedef TVector<int8_t, 4>		sbyte4;

typedef int32_t					int1;
typedef TVector<int32_t, 2>		int2;
typedef TVector<int32_t, 3>		int3;
typedef TVector<int32_t, 4>		int4;

typedef uint32_t				uint1;
typedef TVector<uint32_t, 2>	uint2;
typedef TVector<uint32_t, 3>	uint3;
typedef TVector<uint32_t, 4>	uint4;

typedef TVector<int16_t, 2>		short2;
typedef TVector<int16_t, 3>		short3;
typedef TVector<int16_t, 4>		short4;

typedef TVector<uint16_t, 2>	ushort2;
typedef TVector<uint16_t, 3>	ushort3;
typedef TVector<uint16_t, 4>	ushort4;

typedef bool					bool1;
typedef TVector<bool, 2>		bool2;
typedef TVector<bool, 3>		bool3;
typedef TVector<bool, 4>		bool4;


namespace std
{
	template<class T>
	struct hash<TVector<T, 2> >
	{
		std::size_t operator()(const TVector<T, 2>& v) const noexcept
		{
			std::size_t result = 0;

			hash_combine(result, v.x);
			hash_combine(result, v.y);

			return result;
		}
	};

	template<class T>
	struct hash<TVector<T, 3> >
	{
		std::size_t operator()(const TVector<T, 3>& v) const noexcept
		{
			std::size_t result = 0;

			hash_combine(result, v.x);
			hash_combine(result, v.y);
			hash_combine(result, v.z);

			return result;
		}
	};

	template<class T>
	struct hash<TVector<T, 4> >
	{
		std::size_t operator()(const TVector<T, 4>& v) const noexcept
		{
			std::size_t result = 0;

			hash_combine(result, v.x);
			hash_combine(result, v.y);
			hash_combine(result, v.z);
			hash_combine(result, v.w);

			return result;
		}
	};
}
