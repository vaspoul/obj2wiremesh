#pragma once

struct float3
{
	static const float3 k000;
	static const float3 k111;
	static const float3 k100;
	static const float3 k010;
	static const float3 k001;
	static const float3 k110;
	static const float3 k101;
	static const float3 k011;
	static const float3 kFloatMax;
	static const float3 kNegativeFloatMax;

	union
	{
		float values[3];

		struct
		{
			float x, y, z;
		};
	};

	float3()
		: x(0), y(0), z(0)
	{}

	float3(const float3& other) = default;

	explicit float3(float x1, float y1, float z1)
		: x(x1), y(y1), z(z1)
	{}

	explicit float3(const float* d)
		: x(d[0]), y(d[1]), z(d[2])
	{}

	void Zero()
	{
		x = y = z = 0;
	}

	void operator=(float f)
	{
		x = f;
		y = f;
		z = f;
	}

	//void operator=(const float3& a)
	//{
	//	x = a.x;
	//	y = a.y;
	//	z = a.z;
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
		return sqrtf(x*x + y*y + z*z);
	}

	float lengthSqr() const
	{
		return (x*x + y*y + z*z);
	}

	float3 normalize()
	{
		float mag = length();

		assert(fabs(mag) > 1e-5);
    
		(*this) *= (1.0f/mag);

		return (*this);
	}

	float3 unit() const
	{
		float mag = length();

		assert(fabs(mag) > FLT_MIN);
    
		return (*this) * (1.0f/mag);
	}

	float3 abs() const
	{
		return float3( fabs(x), fabs(y), fabs(z) );
	}

	float dot(const float3& v) const
	{
		return x * v.x + y * v.y + z * v.z;
	}

	float3 operator-() const
	{
		return float3(-x, -y, -z);
	}

	const float3& operator+() const
	{
		return *this;
	}

	float3 operator-(const float3& a) const
	{
		return float3( x - a.x,
						y - a.y,
						z - a.z);
	}

	float3 operator-(float v) const
	{
		return float3(  x - v,
						y - v,
						z - v);
	}

	float3 operator+(const float3& a) const
	{
		return float3(	x + a.x,
						y + a.y,
						z + a.z);
	}

	float3 operator+(float v) const
	{
		return float3(	x + v,
						y + v,
						z + v);
	}

	float3 operator*(const float3& a) const
	{
		return float3(	x * a.x,
						y * a.y,
						z * a.z);
	}

	float3 operator%(const float3& a) const
	{
		return float3(	 y * a.z - z * a.y,
	    				-x * a.z + z * a.x,
		    			 x * a.y - y * a.x);
	}

	float3 operator*(float a) const
	{
		return float3(x*a, y*a, z*a);
	}

	float3 operator/(const float3& a) const
	{
		return float3(x/a.x, y/a.y, z/a.z);
	}

	float3 operator/(float a) const
	{
		float one_over_a(1.0f/a);
		return float3(x * one_over_a, y * one_over_a, z * one_over_a);
	}

	bool operator==(const float3& a) const
	{
		return	fabs(a.x-x) < 1e-5 &&
				fabs(a.y-y) < 1e-5 &&
				fabs(a.z-z) < 1e-5;
	}

	bool operator!=(const float3& a) const
	{
		return	fabs(a.x-x) > 1e-5 ||
				fabs(a.y-y) > 1e-5 ||
				fabs(a.z-z) > 1e-5;
	}

	void operator-=(const float3& a)
	{
		x -= a.x;
		y -= a.y;
		z -= a.z;
	}

	void operator-=(float v)
	{
		x -= v;
		y -= v;
		z -= v;
	}

	void operator+=(const float3& a)
	{
		x += a.x;
		y += a.y;
		z += a.z;
	}

	void operator+=(float v)
	{
		x += v;
		y += v;
		z += v;
	}

	void operator*=(const float3& a)
	{
		x *= a.x;
		y *= a.y;
		z *= a.z;
	}

	void operator*=(float a)
	{
		x *= a;
		y *= a;
		z *= a;
	}

	void operator/=(const float3& a)
	{
		*this = *this / a;
	}

	void operator/=(float a)
	{
		*this = *this / a;
	}

	float3 operator<(const float3& a) const
	{
		return	float3(	x < a.x ? 1.0f : 0,
						y < a.y ? 1.0f : 0,
						z < a.z ? 1.0f : 0);
	}

	float3 operator<=(const float3& a) const
	{
		return	float3(	x <= a.x ? 1.0f : 0,
						y <= a.y ? 1.0f : 0,
						z <= a.z ? 1.0f : 0 );
	}

	float3 operator>(const float3& a) const
	{
		return	float3(	x > a.x ? 1.0f : 0,
						y > a.y ? 1.0f : 0,
						z > a.z ? 1.0f : 0 );
	}

	float3 operator>=(const float3& a) const
	{
		return	float3(	x >= a.x ? 1.0f : 0,
						y >= a.y ? 1.0f : 0,
						z >= a.z ? 1.0f : 0 );
	}

	bool AllTrue() const
	{
		return x>0 && y>0 && z>0;
	}

	bool AllFalse() const
	{
		return x==0 && y==0 && z==0;
	}

	bool AnyTrue() const
	{
		return !AllFalse();
	}

	bool AnyFalse() const
	{
		return !AllTrue();
	}
};

template<typename Archive>
Archive& operator<<(Archive& archive, float3& v)
{
	archive << v.values;

	return archive;
}

#if _MSC_VER >= 1920 // VS 2019 onwards
inline const float3 float3::k000 = float3(0,0,0);
inline const float3 float3::k111 = float3(1,1,1);
inline const float3 float3::k100 = float3(1,0,0);
inline const float3 float3::k010 = float3(0,1,0);
inline const float3 float3::k001 = float3(0,0,1);
inline const float3 float3::k110 = float3(1,1,0);
inline const float3 float3::k101 = float3(1,0,1);
inline const float3 float3::k011 = float3(0,1,1);
inline const float3 float3::kFloatMax = float3(FLT_MAX,FLT_MAX,FLT_MAX);
inline const float3 float3::kNegativeFloatMax = float3(-FLT_MAX,-FLT_MAX,-FLT_MAX);
#endif

inline float3 operator*(const float& left, const float3& right)
{
	return right * left;
}

inline float3 operator/(float nominator, const float3& denominator)
{
	return float3(nominator/denominator.x, nominator/denominator.y, nominator/denominator.z);
}
