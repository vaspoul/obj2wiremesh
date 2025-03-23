#pragma once

struct float4
{
	static const float4 k0000;
	static const float4 k1111;
	static const float4 k1000;
	static const float4 k0100;
	static const float4 k0010;
	static const float4 k0001;
	static const float4 k1100;
	static const float4 k1010;
	static const float4 k1001;
	static const float4 k0110;
	static const float4 k0101;
	static const float4 k0011;
	static const float4 k1110;
	static const float4 k0111;
	static const float4 k1011;
	static const float4 k1101;
	static const float4 kFloatMax;
	static const float4 kNegativeFloatMax;

	union
	{
		float values[4];

		struct
		{
			float x, y, z, w;
		};

		struct
		{
			float x, y, width, height; // for viewports
		};

		struct 
		{
			float3	xyz;
			float	w;
		};

		struct 
		{
			float2	xy;
			float2	zw;
		};

		struct 
		{
			float	x;
			float3	yzw;
		};

		struct
		{
			float r, g, b, a;
		};

		struct 
		{
			float3	rgb;
			float	a;
		};

		struct 
		{
			float2	rg;
			float2	ba;
		};

		struct 
		{
			float	r;
			float3	gba;
		};
	};

	float4()
		: x(0), y(0), z(0), w(0)
	{}

	float4(const float4& other) = default;

	explicit float4(float x1, float y1, float z1, float w1)
		: x(x1), y(y1), z(z1), w(w1)
	{}

	explicit float4(float x1, float3 yzw)
		: x(x1), y(yzw.x), z(yzw.y), w(yzw.z)
	{}

	explicit float4(float2 xy, float2 zw)
		: x(xy.x), y(xy.y), z(zw.x), w(zw.y)
	{}

	explicit float4(float2 xy, float z, float w)
		: x(xy.x), y(xy.y), z(z), w(w)
	{}

	explicit float4(const float* d)
		: x(d[0]), y(d[1]), z(d[2]), w(d[3])
	{}

	explicit float4(const float3& v3, float w_=1.0f)
		: x(v3.x), y(v3.y), z(v3.z), w(w_)
	{}
	    
	void Zero()
	{
		x = y = z = w = 0;
	}

	void operator=(float f)
	{
		x = f;
		y = f;
		z = f;
		w = f;
	}

	//void operator=(const float4& v)
	//{
	//	x = v.x;
	//	y = v.y;
	//	z = v.z;
	//	w = v.w;
	//}

    float& operator[](int index)
    {
		return *((&x) + index);
    }

    const float& operator[](int index) const
    {
		return *((&x) + index);
    }
 
	float length() const
	{
		return sqrtf(x*x + y*y + z*z + w*w);
	}

	float lengthSqr() const
	{
		return (x*x + y*y + z*z + w*w);
	}

	float3 normalize()
	{
		float mag = length();

		assert(fabs(mag) > 1e-5);
    
		(*this) *= (1.0f/mag);

		return (*this);
	}

	float4 unit() const
	{
		float mag = length();

		assert(fabs(mag) > 1e-5);
    
		return (*this) * (1.0f/mag);
	}

	float4 abs() const
	{
		return float4( fabs(x), fabs(y), fabs(z), fabs(w) );
	}

	float4 operator-() const
	{
		return float4(-x, -y, -z, -w);
	}

	const float4& operator+() const
	{
		return *this;
	}

	float4 operator-(const float4& v) const
	{
		return float4(	x - v.x,
						y - v.y,
						z - v.z,
						w - v.w);
	}

	float4 operator-(float v) const
	{
		return float4(	x - v,
						y - v,
						z - v,
						w - v);
	}

	float4 operator+(const float4& v) const
	{
		return float4(	x + v.x,
						y + v.y,
						z + v.z, 
						w + v.w);
	}

	float4 operator+(float v) const
	{
		return float4(	x + v,
						y + v,
						z + v, 
						w + v);
	}

	float4 operator*(const float4& v) const
	{
		return float4(	x * v.x,
						y * v.y,
						z * v.z,
						w * v.w);
	}

	float4 operator%(const float4& v) const
	{
		return float4(	 y * v.z - z * v.y,
	    				-x * v.z + z * v.x,
		    			 x * v.y - y * v.x,
						 1.0f);
	}

	float4 operator*(float v) const
	{
		return float4(x*v, y*v, z*v, w*v);
	}

	float4 operator/(const float4& v) const
	{
		return float4(x/v.x, y/v.y, z/v.z, w/v.w);
	}

	float4 operator/(float v) const
	{
		float one_over_v(1.0f/v);
		return float4(x * one_over_v, y * one_over_v, z * one_over_v, w * one_over_v);
	}

	bool operator==(const float4& v) const
	{
		return	fabs(v.x-x) < 1e-5 &&
				fabs(v.y-y) < 1e-5 &&
				fabs(v.z-z) < 1e-5 &&
				fabs(v.w-w) < 1e-5;
	}

	bool operator!=(const float4& v) const
	{
		return	fabs(v.x-x) > 1e-5 ||
				fabs(v.y-y) > 1e-5 ||
				fabs(v.z-z) > 1e-5 ||
				fabs(v.w-w) > 1e-5;
	}

	void operator-=(const float4& v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
		w -= v.w;
	}

	void operator-=(float v)
	{
		x -= v;
		y -= v;
		z -= v;
		w -= v;
	}

	void operator+=(const float4& v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
		w += v.w;
	}

	void operator+=(float v)
	{
		x += v;
		y += v;
		z += v;
		w += v;
	}

	void operator*=(const float4& v)
	{
		x *= v.x;
		y *= v.y;
		z *= v.z;
		w *= v.w;
	}

	void operator*=(float v)
	{
		x *= v;
		y *= v;
		z *= v;
		w *= v;
	}

	void operator/=(const float4& v)
	{
		*this = *this / v;
	}

	void operator/=(float v)
	{
		*this = *this / v;
	}

	float4 operator<(const float4& v) const
	{
		return	float4(	x < v.x ? 1.0f : 0,
						y < v.y ? 1.0f : 0,
						z < v.z ? 1.0f : 0,
						w < v.w ? 1.0f : 0);
	}

	float4 operator<=(const float4& v) const
	{
		return	float4(	x <= v.x ? 1.0f : 0,
						y <= v.y ? 1.0f : 0,
						z <= v.z ? 1.0f : 0,
						w <= v.w ? 1.0f : 0 );
	}

	float4 operator>(const float4& v) const
	{
		return	float4(	x > v.x ? 1.0f : 0,
						y > v.y ? 1.0f : 0,
						z > v.z ? 1.0f : 0,
						w > v.w ? 1.0f : 0 );
	}

	float4 operator>=(const float4& v) const
	{
		return	float4(	x >= v.x ? 1.0f : 0,
						y >= v.y ? 1.0f : 0,
						z >= v.z ? 1.0f : 0,
						w >= v.w ? 1.0f : 0 );
	}

	bool AllTrue() const
	{
		return x>0 && y>0 && z>0 && w>0;
	}

	bool AllFalse() const
	{
		return x==0 && y==0 && z==0 && z==0;
	}

	bool AnyTrue() const
	{
		return !AllFalse();
	}

	bool AnyFalse() const
	{
		return !AllTrue();
	}

	operator const float3& () const	{ return xyz;	}
	operator float3()				{ return xyz;	}
};

template<typename Archive>
Archive& operator<<(Archive& archive, float4& v)
{
	archive << v.values;

	return archive;
}


#if _MSC_VER >= 1920 // VS 2019 onwards
inline const float4 float4::k0000 = float4(0,0,0,0);
inline const float4 float4::k1111 = float4(1,1,1,1);
inline const float4 float4::k1000 = float4(1,0,0,0);
inline const float4 float4::k0100 = float4(0,1,0,0);
inline const float4 float4::k0010 = float4(0,0,1,0);
inline const float4 float4::k0001 = float4(0,0,0,1);
inline const float4 float4::k1100 = float4(1,1,0,0);
inline const float4 float4::k1010 = float4(1,0,1,0);
inline const float4 float4::k1001 = float4(1,0,0,1);
inline const float4 float4::k0110 = float4(0,1,1,0);
inline const float4 float4::k0101 = float4(0,1,0,1);
inline const float4 float4::k0011 = float4(0,0,1,1);
inline const float4 float4::k1110 = float4(1,1,1,0);
inline const float4 float4::k0111 = float4(0,1,1,1);
inline const float4 float4::k1011 = float4(1,0,1,1);
inline const float4 float4::k1101 = float4(1,1,0,1);
inline const float4 float4::kFloatMax = float4(FLT_MAX,FLT_MAX,FLT_MAX,FLT_MAX);
inline const float4 float4::kNegativeFloatMax = float4(-FLT_MAX,-FLT_MAX,-FLT_MAX,-FLT_MAX);
#endif

inline float4 operator*(const float& left, const float4& right)
{
	return right * left;
}

inline float4 operator/(float nominator, const float4& denominator)
{
	return float4(nominator/denominator.x, nominator/denominator.y, nominator/denominator.z, nominator/denominator.w);
}
