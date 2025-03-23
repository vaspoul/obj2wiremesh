#pragma once

struct float2
{
	static const float2 k00;
	static const float2 k11;
	static const float2 k10;
	static const float2 k01;
	static const float2 kFloatMax;
	static const float2 kNegativeFloatMax;

	union 
	{ 
		float values[2];

		struct
		{
			float x, y; 
		};

		struct
		{
			float r, g; 
		};

		struct
		{
			float u, v; 
		};
	};

	float2()
		: x(0)
		, y(0)
	{}

	float2(const float2& other) = default;

	explicit float2(float x1, float y1)
		: x(x1), y(y1)
	{}

	explicit float2(const float* d)
		: x(d[0]), y(d[1])
	{}
	    
	void Zero()
	{
		x = y = 0;
	}

	void operator=(float f)
	{
		x = f;
		y = f;
	}

	//void operator=(const float2& a)
	//{
	//	x = a.x;
	//	y = a.y;
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
		return sqrtf(x*x + y*y);
	}

	float lengthSqr() const
	{
		return (x*x + y*y);
	}

	float2 normalize()
	{
		float mag = length();

		assert(fabs(mag) > 1e-5);
    
		(*this) *= (1.0f/mag);

		return (*this);
	}

	float2 unit() const
	{
		float mag = length();

		assert(fabs(mag) > 1e-5);
    
		return (*this) * (1.0f/mag);
	}

	float2 abs() const
	{
		return float2( fabs(x), fabs(y) );
	}

	float2 tangent() const
	{
		return float2(-y, x);
	}

	float dot(const float2& a) const
	{
		return x * a.x + y * a.y;
	}

	float2 operator-() const
	{
		return float2(-x, -y);
	}

	const float2& operator+() const
	{
		return *this;
	}

	float2 operator-(const float2& a) const
	{
		return float2(	x - a.x,
						y - a.y);
	}

	float2 operator-(float a) const
	{
		return float2(	x - a,
						y - a);
	}

	float2 operator+(const float2& a) const
	{
		return float2(	x + a.x,
						y + a.y);
	}

	float2 operator+(float a) const
	{
		return float2(	x + a,
						y + a);
	}

	float2 operator*(const float2& a) const
	{
		return float2(	x * a.x,
						y * a.y);
	}

	float2 operator*(float a) const
	{
		return float2(x*a, y*a);
	}

	float2 operator/(const float2& a) const
	{
		return float2(x/a.x, y/a.y);
	}

	float2 operator/(float a) const
	{
		float one_over_a(1.0f/a);
		return float2(x * one_over_a, y * one_over_a);
	}

	bool operator==(const float2& a) const
	{
		return	fabs(a.x-x) < 1e-5 &&
				fabs(a.y-y) < 1e-5;
	}

	bool operator!=(const float2& a) const
	{
		return	fabs(a.x-x) > 1e-5 ||
				fabs(a.y-y) > 1e-5;
	}

	void operator-=(const float2& a)
	{
		x -= a.x;
		y -= a.y;
	}

	void operator-=(float a)
	{
		x -= a;
		y -= a;
	}

	void operator+=(const float2& a)
	{
		x += a.x;
		y += a.y;
	}

	void operator+=(float a)
	{
		x += a;
		y += a;
	}

	void operator*=(const float2& a)
	{
		x *= a.x;
		y *= a.y;
	}

	void operator*=(float a)
	{
		x *= a;
		y *= a;
	}

	void operator/=(const float2& a)
	{
		*this = *this / a;
	}

	void operator/=(float a)
	{
		*this = *this / a;
	}

	float2 operator<(const float2& a) const
	{
		return	float2(	x < a.x ? 1.0f : 0,
						y < a.y ? 1.0f : 0 );
	}

	float2 operator<=(const float2& a) const
	{
		return	float2(	x <= a.x ? 1.0f : 0,
						y <= a.y ? 1.0f : 0 );
	}

	float2 operator>(const float2& a) const
	{
		return	float2(	x > a.x ? 1.0f : 0,
						y > a.y ? 1.0f : 0 );
	}

	float2 operator>=(const float2& a) const
	{
		return	float2(	x >= a.x ? 1.0f : 0,
						y >= a.y ? 1.0f : 0 );
	}
	
	bool AllTrue() const
	{
		return x>0 && y>0;
	}

	bool AllFalse() const
	{
		return x==0 && y==0;
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
Archive& operator<<(Archive& archive, float2& v)
{
	archive << v.values;

	return archive;
}

#if _MSC_VER >= 1920 // VS 2019 onwards
inline const float2 float2::k00 = float2(0,0);
inline const float2 float2::k11 = float2(1,1);
inline const float2 float2::k10 = float2(1,0);
inline const float2 float2::k01 = float2(0,1);
inline const float2 float2::kFloatMax = float2(FLT_MAX,FLT_MAX);
inline const float2 float2::kNegativeFloatMax = float2(-FLT_MAX,-FLT_MAX);
#endif

inline float2 operator*(const float& left, const float2& right)
{
	return right * left;
}

inline float2 operator/(float nominator, const float2& denominator)
{
	return float2(nominator/denominator.x, nominator/denominator.y);
}
