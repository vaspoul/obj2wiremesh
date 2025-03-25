#pragma once

#include "float2.h"
#include "float3.h"
#include "float4.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>

#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

// HLSL-like functions

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
// float
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

// Returns the [01] term of cubic hermite interpolation
// See http://en.wikipedia.org/wiki/Cubic_Hermite_spline
static inline float hermite01(const float& v)
{
	return v*v*(3.0f - 2.0f * v);
}

//Test if all components of x are nonzero.
static inline bool all(const float& v)
{
	return (v != 0);
}

//Test if any component of x is nonzero.
static inline bool any(const float& v)
{
	return (v != 0);
}

//Clamps x to the range [min, max].
static inline float clamp(const float& v, const float& vMin, const float& vMax)
{
	return (v<=vMin)?vMin:(v>=vMax)?vMax:v;
}

//Converts x from radians to degrees.
static inline float degrees(const float& v)
{
	return v * 180.0f / float(M_PI);
}

//Returns the distance between two points, a and b.
static inline float distance(const float& a, const float& b)
{
	return (b-a);
}

//Returns the square of the distance between two points, a and b.
static inline float distanceSqr(const float& a, const float& b)
{
	return (b - a) * (b - a);
}

//Returns the fractional part f of x, such that f is a value greater than or equal to 0, and less than 1.
static inline float frac(const float& v)
{
	return v - floor(v);
}

//Returns a + s(b - a). This linearly interpolates between a and b, such that the return value is a when s is 0, and b when s is 1.
static inline float lerp(const float& a, const float& b, const float& s)
{
	return a + (b-a) * s;
}

////Returns the base-2 logarithm of x. If x is negative, the function returns indefinite. If x is 0, the function returns +INF.
//static inline float log2(const float& v)
//{
//	static const float clog2 = logf(2.0f);
//	return logf(v) / clog2;
//}

//Converts x from degrees to radians.
static inline float radians(const float& v)
{
	return v / 180.0f * float(M_PI);
}

//Returns 1 / sqrt(x)
static inline float rsqrt(const float& v)
{
	return 1.0f/sqrtf(v);
}

//Clamps x to the range [0, 1]
static inline float saturate(const float& v)
{
	return clamp(v, 0, 1);
}

//Computes the sign of x. Returns -1 if x is less than 0, 0 if x equals 0, and 1 if x is greater than zero.
static inline float sign(const float& v)
{
	return v>0?1.0f:(v<0)?-1.0f:0;
}

//Returns 0 if x < min. Returns 1 if x > max. Returns a smooth Hermite interpolation between 0 and 1, if x is in the range [min, max].
static inline float smoothstep(const float& vMin, const float& vMax, const float& v)
{
	float interpolator(saturate(hermite01(v)));
	return lerp(vMin, vMax, interpolator);
}

//Returns (x >= a) ? 1 : 0
static inline float step(const float& a, const float& x)
{
	return x>=a?1.0f:0;
}

//Returns the sine and cosine of x. sin(x) is stored in the output parameter s. cos(x) is stored in the output parameter c.
static inline void sincos(const float& v, float& s, float& c)
{
	s = sinf(v);
	c = cosf(v);
}


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
// float2
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

//Absolute value (per component).
static inline float2 abs(const float2& v)
{
	return float2(fabs(v.x), fabs(v.y));
}

//Returns the arccosine of each component of x. Each component should be in the range [-1, 1].
static inline float2 acos(const float2& v)
{
	return float2(acosf(v.x), acosf(v.y));
}

//Test if all components of x are nonzero.
static inline bool all(const float2& v)
{
	return (v.x != 0) && (v.y != 0);
}

//Test if any component of x is nonzero.
static inline bool any(const float2& v)
{
	return (v.x != 0) || (v.y != 0);
}

//Returns the arcsine of each component of x. Each component should be in the range [-pi/2, pi/2].
static inline float2 asin(const float2& v)
{
	return float2(asinf(v.x), asinf(v.y));
}

//Returns the arctangent of x. The return values are in the range [-pi/2, pi/2].
static inline float2 atan(const float2& v)
{
	return float2(atanf(v.x), atanf(v.y));
}

//Returns the smallest integer which is greater than or equal to x.
static inline float2 ceil(const float2& v)
{
	return float2(ceilf(v.x), ceilf(v.y));
}

//Clamps x to the range [min, max].
static inline float2 clamp(const float2& v, const float2& vMin, const float2& vMax)
{
	return float2(clamp(v.x, vMin.x, vMax.x), clamp(v.y, vMin.y, vMax.y));
}

//Returns the cosine of x.
static inline float2 cos(const float2& v)
{
	return float2(cosf(v.x), cosf(v.y));
}

//Returns the hyperbolic cosine of x.
static inline float2 cosh(const float2& v)
{
	return float2(coshf(v.x), coshf(v.y));
}

//Returns the cross product of two 3D vectors a and b.
static inline float cross(const float2& a, const float2& b)
{
	return a.x * b.y - a.y * b.x;
}

//Converts x from radians to degrees.
static inline float2 degrees(const float2& v)
{
	return float2(v.x * 180.0f / float(M_PI), v.y * 180.0f / float(M_PI));
}

//Returns the distance between two points, a and b.
static inline float distance(const float2& a, const float2& b)
{
	return (b-a).length();
}

//Returns the square of the distance between two points, a and b.
static inline float distanceSqr(const float2& a, const float2& b)
{
	return (b - a).lengthSqr();
}

static inline float2 sub(const float2& left, const float2& right)
{
	return float2( left.x - right.x, left.y - right.y );
}

static inline float dot(const float2& left, const float2& right)
{
	return  left.x * right.x +
			left.y * right.y;
}

static inline float2 avg(const float2& left, const float2& right)
{
	return float2( (left.x + right.x) * 0.5f, (left.y + right.y) * 0.5f );
}

static inline float2 mad(const float2& a, const float2& b, const float2& c)
{
	return float2(  a.x * b.x + c.x, 
					a.y * b.y + c.y );
}

static inline float2 mad(const float2& a, const float b, const float2& c)
{
	return float2(  a.x * b + c.x, 
					a.y * b + c.y );
}

//Returns the base-e exponent.
static inline float2 exp(const float2& v)
{
	return float2(expf(v.x), expf(v.y));
}

//Base 2 Exp (per component).
static inline float2 exp2(const float2& v)
{
	return float2(powf(2.0f, v.x), powf(2.0f, v.y));
}

//Returns the greatest integer which is less than or equal to x.
static inline float2 floor(const float2& v)
{
	return float2(floorf(v.x), floorf(v.y));
}

//Returns the fractional part f of x, such that f is a value greater than or equal to 0, and less than 1.
static inline float2 frac(const float2& v)
{
	return v - floor(v);
}

//Returns the floating point remainder f of a / b such that a = i * b + f, where i is an integer, f has the same sign as x, and the absolute value of f is less than the absolute value of b.
static inline float2 fmod(const float2& a, const float2& b)
{
	return frac(a/b);
}

static inline float2 round(const float2& v)
{
	return float2(round(v.x), round(v.y));
}

//Returns the length of the vector v.
static inline float length(const float2& v)
{
	return v.length();
}

//Returns the length^2 of the vector v.
static inline float lengthSqr(const float2& v)
{
	return v.lengthSqr();
}

//Returns a + s(b - a). This linearly interpolates between a and b, such that the return value is a when s is 0, and b when s is 1.
static inline float2 lerp(const float2& a, const float2& b, const float2& s)
{
	return a + (b-a) * s;
}

//Returns a + s(b - a). This linearly interpolates between a and b, such that the return value is a when s is 0, and b when s is 1.
static inline float2 lerp(const float2& a, const float2& b, const float& s)
{
	return a + (b-a) * s;
}

//Returns the base-e logarithm of x. If x is negative, the function returns indefinite. If x is 0, the function returns +INF.
static inline float2 log(const float2& v)
{
	return float2(logf(v.x), logf(v.y));
}

//Returns the base-10 logarithm of x. If x is negative, the function returns indefinite. If x is 0, the function returns +INF.
static inline float2 log10(const float2& v)
{
	return float2(log10f(v.x), log10f(v.y));
}

//Returns the base-2 logarithm of x. If x is negative, the function returns indefinite. If x is 0, the function returns +INF.
static inline float2 log2(const float2& v)
{
	return float2(log2(v.x), log2(v.y));
}

//Selects the greater of a and b.
static inline float2 max(const float2& a, const float2& b)
{
	return float2(std::max(a.x, b.x), std::max(a.y, b.y));
}

//Selects the lesser of a and b.
static inline float2 min(const float2& a, const float2& b)
{
	return float2(std::min(a.x, b.x), std::min(a.y, b.y));
}

//Selects the greatest component
static inline float max(const float2& v)
{
	return std::max(v.x, v.y);
}

//Selects the least component
static inline float min(const float2& v)
{
	return std::min(v.x, v.y);
}

//Performs matrix multiplication between a and b. 
//If a is a vector, it is treated as a row vector. 
//If b is a vector, it is treated as a column vector. 
//The inner dimension acolumns and brows must be equal. The result has the dimension arows x bcolumns.
static inline float2 mul(const float2& a, const float2& b)
{
	return a * b;
}

//Returns the normalized vector v / length(v). If the length of v is 0, the result is indefinite.
static inline float2 normalize(const float2& v)
{
	float mag = v.length();

	assert(fabs(mag) > 1e-5);
    
    return v * (1.0f/mag);
}

// Returns the reciprocal of v
static inline float2 rcp(const float2& v)
{
    return float2(1.0f / v.x, 1.0f / v.y);
}

//Returns x^y.
static inline float2 pow(const float2& a, const float2& b)
{
	return float2(powf(a.x, b.x), powf(a.y, b.y));
}

static inline float2 pow(const float2& a, float b)
{
	return float2(powf(a.x, b), powf(a.y, b));
}

//Converts x from degrees to radians.
static inline float2 radians(const float2& v)
{
	return float2(v.x / 180.0f * float(M_PI), v.y / 180.0f * float(M_PI));
}

//Returns the reflection vector v, given the entering ray direction i, and the surface normal n, such that v = i - 2n * (i•n).
static inline float2 reflect(const float2& i, const float2& n)
{
	return i - 2.0f * dot(i, n) * n;
}

//Returns 1 / sqrt(x)
static inline float2 rsqrt(const float2& v)
{
	return float2(1.0f/sqrtf(v.x), 1.0f/sqrtf(v.y));
}

//Clamps x to the range [0, 1]
static inline float2 saturate(const float2& v)
{
	return float2(saturate(v.x), saturate(v.y));
}

//Computes the sign of x. Returns -1 if x is less than 0, 0 if x equals 0, and 1 if x is greater than zero.
static inline float2 sign(const float2& v)
{
	return float2(sign(v.x), sign(v.y));
}

//Returns the sine of x
static inline float2 sin(const float2& v)
{
	return float2(sinf(v.x), sinf(v.y));
}

//Returns the sine and cosine of x. sin(x) is stored in the output parameter s. cos(x) is stored in the output parameter c.
static inline void sincos(const float2& v, float2& s, float2& c)
{
	s = sin(v);
	c = cos(v);
}

//Returns the hyperbolic sine of x
static inline float2 sinh(const float2& v)
{
	return float2(sinhf(v.x), sinhf(v.y));
}


//Returns 0 if x < min. Returns 1 if x > max. Returns a smooth Hermite interpolation between 0 and 1, if x is in the range [min, max].
static inline float2 smoothstep(const float2& vMin, const float2& vMax, const float2& v)
{
	float2 interpolator(saturate(hermite01(v.x)), saturate(hermite01(v.y)));
	return lerp(vMin, vMax, interpolator);
}

//Square root (per component)
static inline float2 sqrt(const float2& v)
{
	return float2(sqrtf(v.x), sqrtf(v.y));
}

//Returns (x >= a) ? 1 : 0
static inline float2 step(const float2& a, const float2& x)
{
	return float2(x.x>=a.x?1.0f:0, x.y>=a.y?1.0f:0);
}

//Returns the tangent of x
static inline float2 tan(const float2& v)
{
	return float2(tanf(v.x), tanf(v.y));
}

//Returns the hyperbolic tangent of x
static inline float2 tanh(const float2& v)
{
	return float2(tanhf(v.x), tanhf(v.y));
}


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
// float3
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

//Absolute value (per component).
static inline float3 abs(const float3& v)
{
	return float3(fabs(v.x), fabs(v.y), fabs(v.z));
}

//Returns the arccosine of each component of x. Each component should be in the range [-1, 1].
static inline float3 acos(const float3& v)
{
	return float3(acosf(v.x), acosf(v.y), acosf(v.z));
}

//Test if all components of x are nonzero.
static inline bool all(const float3& v)
{
	return (v.x != 0) && (v.y != 0) && (v.z != 0);
}

//Test if any component of x is nonzero.
static inline bool any(const float3& v)
{
	return (v.x != 0) || (v.y != 0) || (v.z != 0);
}

//Returns the arcsine of each component of x. Each component should be in the range [-pi/2, pi/2].
static inline float3 asin(const float3& v)
{
	return float3(asinf(v.x), asinf(v.y), asinf(v.z));
}

//Returns the arctangent of x. The return values are in the range [-pi/2, pi/2].
static inline float3 atan(const float3& v)
{
	return float3(atanf(v.x), atanf(v.y), atanf(v.z));
}

//Returns the smallest integer which is greater than or equal to x.
static inline float3 ceil(const float3& v)
{
	return float3(ceilf(v.x), ceilf(v.y), ceilf(v.z));
}

//Clamps x to the range [min, max].
static inline float3 clamp(const float3& v, const float3& vMin, const float3& vMax)
{
	return float3(clamp(v.x, vMin.x, vMax.x), clamp(v.y, vMin.y, vMax.y), clamp(v.z, vMin.z, vMax.z));
}

//Returns the cosine of x.
static inline float3 cos(const float3& v)
{
	return float3(cosf(v.x), cosf(v.y), cosf(v.z));
}

//Returns the hyperbolic cosine of x.
static inline float3 cosh(const float3& v)
{
	return float3(coshf(v.x), coshf(v.y), coshf(v.z));
}

//Returns the cross product of two 3D vectors a and b.
static inline float3 cross(const float3& a, const float3& b)
{
	return float3(	 a.y * b.z - a.z * b.y,
	    			-a.x * b.z + a.z * b.x,
		    		 a.x * b.y - a.y * b.x);
}

//Converts x from radians to degrees.
static inline float3 degrees(const float3& v)
{
	return float3(v.x * 180.0f / float(M_PI), v.y * 180.0f / float(M_PI), v.z * 180.0f / float(M_PI));
}

//Returns the distance between two points, a and b.
static inline float distance(const float3& a, const float3& b)
{
	return (b-a).length();
}

//Returns the square of the distance between two points, a and b.
static inline float distanceSqr(const float3& a, const float3& b)
{
	return (b - a).lengthSqr();
}

static inline float dot(const float3& left, const float3& right)
{
	return  left.x * right.x +
			left.y * right.y +
			left.z * right.z;
}

static inline float3 sub(const float3& left, const float3& right)
{
	return float3( left.x - right.x, left.y - right.y, left.z - right.z );
}

static inline float3 avg(const float3& left, const float3& right)
{
	return float3( (left.x + right.x) * 0.5f, (left.y + right.y) * 0.5f, (left.z + right.z) * 0.5f );
}

static inline float3 mad(const float3& a, const float3& b, const float3& c)
{
	return float3(  a.x * b.x + c.x, 
					a.y * b.y + c.y, 
					a.z * b.z + c.z );
}

static inline float3 mad(const float3& a, const float b, const float3& c)
{
	return float3(  a.x * b + c.x, 
					a.y * b + c.y, 
					a.z * b + c.z );
}

//Returns the base-e exponent.
static inline float3 exp(const float3& v)
{
	return float3(expf(v.x), expf(v.y), expf(v.z));
}

//Base 2 Exp (per component).
static inline float3 exp2(const float3& v)
{
	return float3(powf(2.0f, v.x), powf(2.0f, v.y), powf(2.0f, v.z));
}

//Returns the greatest integer which is less than or equal to x.
static inline float3 floor(const float3& v)
{
	return float3(floorf(v.x), floorf(v.y), floorf(v.z));
}

//Returns the fractional part f of x, such that f is a value greater than or equal to 0, and less than 1.
static inline float3 frac(const float3& v)
{
	return v - floor(v);
}

//Returns the floating point remainder f of a / b such that a = i * b + f, where i is an integer, f has the same sign as x, and the absolute value of f is less than the absolute value of b.
static inline float3 fmod(const float3& a, const float3& b)
{
	return frac(a/b);
}

static inline float3 round(const float3& v)
{
	return float3(round(v.x), round(v.y), round(v.z));
}

//Returns the length of the vector v.
static inline float length(const float3& v)
{
	return v.length();
}

//Returns the length^2 of the vector v.
static inline float lengthSqr(const float3& v)
{
	return v.lengthSqr();
}

//Returns a + s(b - a). This linearly interpolates between a and b, such that the return value is a when s is 0, and b when s is 1.
static inline float3 lerp(const float3& a, const float3& b, const float3& s)
{
	return a + (b-a) * s;
}

//Returns a + s(b - a). This linearly interpolates between a and b, such that the return value is a when s is 0, and b when s is 1.
static inline float3 lerp(const float3& a, const float3& b, const float& s)
{
	return a + (b-a) * s;
}

//Returns the base-e logarithm of x. If x is negative, the function returns indefinite. If x is 0, the function returns +INF.
static inline float3 log(const float3& v)
{
	return float3(logf(v.x), logf(v.y), logf(v.z));
}

//Returns the base-10 logarithm of x. If x is negative, the function returns indefinite. If x is 0, the function returns +INF.
static inline float3 log10(const float3& v)
{
	return float3(log10f(v.x), log10f(v.y), log10f(v.z));
}

//Returns the base-2 logarithm of x. If x is negative, the function returns indefinite. If x is 0, the function returns +INF.
static inline float3 log2(const float3& v)
{
	return float3(log2(v.x), log2(v.y), log2(v.z));
}

//Selects the greater of a and b.
static inline float3 max(const float3& a, const float3& b)
{
	return float3(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));
}

//Selects the lesser of a and b.
static inline float3 min(const float3& a, const float3& b)
{
	return float3(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z));
}

//Selects the greatest component
static inline float max(const float3& v)
{
	return std::max(std::max(v.x, v.y), v.z);
}

//Selects the least component
static inline float min(const float3& v)
{
	return std::min(std::min(v.x, v.y), v.z);
}

//Performs matrix multiplication between a and b. 
//If a is a vector, it is treated as a row vector. 
//If b is a vector, it is treated as a column vector. 
//The inner dimension acolumns and brows must be equal. The result has the dimension arows x bcolumns.
static inline float3 mul(const float3& a, const float3& b)
{
	return a * b;
}

//Returns the normalized vector v / length(v). If the length of v is 0, the result is indefinite.
static inline float3 normalize(const float3& v)
{
	float mag = v.length();

	assert(fabs(mag) > 1e-5);
    
    return v * (1.0f/mag);
}

// Returns the reciprocal of v
static inline float3 rcp(const float3& v)
{
    return float3(1.0f / v.x, 1.0f / v.y, 1.0f / v.z);
}

//Returns x^y.
static inline float3 pow(const float3& a, const float3& b)
{
	return float3(powf(a.x, b.x), powf(a.y, b.y), powf(a.z, b.z));
}

static inline float3 pow(const float3& a, float b)
{
	return float3(powf(a.x, b), powf(a.y, b), powf(a.z, b));
}

//Converts x from degrees to radians.
static inline float3 radians(const float3& v)
{
	return float3(v.x / 180.0f * float(M_PI), v.y / 180.0f * float(M_PI), v.z / 180.0f * float(M_PI));
}

//Returns the reflection vector v, given the entering ray direction i, and the surface normal n, such that v = i - 2n * (i•n).
static inline float3 reflect(const float3& i, const float3& n)
{
	return i - 2.0f * dot(i, n) * n;
}

//Returns 1 / sqrt(x)
static inline float3 rsqrt(const float3& v)
{
	return float3(1.0f/sqrtf(v.x), 1.0f/sqrtf(v.y), 1.0f/sqrtf(v.z));
}

//Clamps x to the range [0, 1]
static inline float3 saturate(const float3& v)
{
	return float3(saturate(v.x), saturate(v.y), saturate(v.z));
}

//Computes the sign of x. Returns -1 if x is less than 0, 0 if x equals 0, and 1 if x is greater than zero.
static inline float3 sign(const float3& v)
{
	return float3(sign(v.x), sign(v.y), sign(v.z));
}

//Returns the sine of x
static inline float3 sin(const float3& v)
{
	return float3(sinf(v.x), sinf(v.y), sinf(v.z));
}

//Returns the sine and cosine of x. sin(x) is stored in the output parameter s. cos(x) is stored in the output parameter c.
static inline void sincos(const float3& v, float3& s, float3& c)
{
	s = sin(v);
	c = cos(v);
}

//Returns the hyperbolic sine of x
static inline float3 sinh(const float3& v)
{
	return float3(sinhf(v.x), sinhf(v.y), sinhf(v.z));
}


//Returns 0 if x < min. Returns 1 if x > max. Returns a smooth Hermite interpolation between 0 and 1, if x is in the range [min, max].
static inline float3 smoothstep(const float3& vMin, const float3& vMax, const float3& v)
{
	float3 interpolator(saturate(hermite01(v.x)), saturate(hermite01(v.y)), saturate(hermite01(v.z)));
	return lerp(vMin, vMax, interpolator);
}

//Square root (per component)
static inline float3 sqrt(const float3& v)
{
	return float3(sqrtf(v.x), sqrtf(v.y), sqrtf(v.z));
}

//Returns (x >= a) ? 1 : 0
static inline float3 step(const float3& a, const float3& x)
{
	return float3(x.x>=a.x?1.0f:0, x.y>=a.y?1.0f:0, x.z>=a.z?1.0f:0);
}

//Returns the tangent of x
static inline float3 tan(const float3& v)
{
	return float3(tanf(v.x), tanf(v.y), tanf(v.z));
}

//Returns the hyperbolic tangent of x
static inline float3 tanh(const float3& v)
{
	return float3(tanhf(v.x), tanhf(v.y), tanhf(v.z));
}


//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
// float4
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////

//Absolute value (per component).
static inline float4 abs(const float4& v)
{
	return float4(fabs(v.x), fabs(v.y), fabs(v.z), fabs(v.w));
}

//Returns the arccosine of each component of x. Each component should be in the range [-1, 1].
static inline float4 acos(const float4& v)
{
	return float4(acosf(v.x), acosf(v.y), acosf(v.z), acosf(v.w));
}

//Test if all components of x are nonzero.
static inline bool all(const float4& v)
{
	return (v.x != 0) && (v.y != 0) && (v.z != 0) && (v.w != 0);
}

//Test if any component of x is nonzero.
static inline bool any(const float4& v)
{
	return (v.x != 0) || (v.y != 0) || (v.z != 0) || (v.w != 0);
}

//Returns the arcsine of each component of x. Each component should be in the range [-pi/2, pi/2].
static inline float4 asin(const float4& v)
{
	return float4(asinf(v.x), asinf(v.y), asinf(v.z), asinf(v.w));
}

//Returns the arctangent of x. The return values are in the range [-pi/2, pi/2].
static inline float4 atan(const float4& v)
{
	return float4(atanf(v.x), atanf(v.y), atanf(v.z), atanf(v.w));
}

//Returns the smallest integer which is greater than or equal to x.
static inline float4 ceil(const float4& v)
{
	return float4(ceilf(v.x), ceilf(v.y), ceilf(v.z), ceilf(v.w));
}

//Clamps x to the range [min, max].
static inline float4 clamp(const float4& v, const float4& vMin, const float4& vMax)
{
	return float4(clamp(v.x, vMin.x, vMax.x), clamp(v.y, vMin.y, vMax.y), clamp(v.z, vMin.z, vMax.z), clamp(v.w, vMin.w, vMax.w));
}

//Returns the cosine of x.
static inline float4 cos(const float4& v)
{
	return float4(cosf(v.x), cosf(v.y), cosf(v.z), cosf(v.w));
}

//Returns the hyperbolic cosine of x.
static inline float4 cosh(const float4& v)
{
	return float4(coshf(v.x), coshf(v.y), coshf(v.z), coshf(v.w));
}

//Returns the cross product of two 3D vectors a and b.
static inline float4 cross(const float4& a, const float4& b)
{
	return float4(	 a.y * b.z - a.z * b.y,
	    			-a.x * b.z + a.z * b.x,
		    		 a.x * b.y - a.y * b.x,
					0);
}

//Converts x from radians to degrees.
static inline float4 degrees(const float4& v)
{
	return float4(v.x * 180.0f / float(M_PI), v.y * 180.0f / float(M_PI), v.z * 180.0f / float(M_PI), v.w * 180.0f / float(M_PI));
}

//Returns the distance between two points, a and b.
static inline float distance(const float4& a, const float4& b)
{
	return (b-a).xyz.length();
}

//Returns the square of the distance between two points, a and b.
static inline float distanceSqr(const float4& a, const float4& b)
{
	return (b - a).xyz.lengthSqr();
}

static inline float dot(const float4& left, const float4& right)
{
	return  left.x * right.x +
			left.y * right.y +
			left.z * right.z +
			left.w * right.w;
}

static inline float4 sub(const float4& left, const float4& right)
{
	return float4( left.x - right.x, left.y - right.y, left.z - right.z, left.w - right.w );
}

static inline float4 avg(const float4& left, const float4& right)
{
	return float4( (left.x + right.x) * 0.5f, (left.y + right.y) * 0.5f, (left.z + right.z) * 0.5f, (left.w + right.w) * 0.5f );
}

static inline float4 mad(const float4& a, const float4& b, const float4& c)
{
	return float4(  a.x * b.x + c.x, 
					a.y * b.y + c.y, 
					a.z * b.z + c.z, 
					a.w * b.w + c.w );
}

static inline float4 mad(const float4& a, float b, const float4& c)
{
	return float4(  a.x * b + c.x, 
					a.y * b + c.y, 
					a.z * b + c.z, 
					a.w * b + c.w );
}

//Returns the base-e exponent.
static inline float4 exp(const float4& v)
{
	return float4(expf(v.x), expf(v.y), expf(v.z), expf(v.w));
}

//Base 2 Exp (per component).
static inline float4 exp2(const float4& v)
{
	return float4(powf(2.0f, v.x), powf(2.0f, v.y), powf(2.0f, v.z), powf(2.0f, v.w));
}

//Returns the greatest integer which is less than or equal to x.
static inline float4 floor(const float4& v)
{
	return float4(floorf(v.x), floorf(v.y), floorf(v.z), floorf(v.w));
}

//Returns the fractional part f of x, such that f is a value greater than or equal to 0, and less than 1.
static inline float4 frac(const float4& v)
{
	return v - floor(v);
}

//Returns the floating point remainder f of a / b such that a = i * b + f, where i is an integer, f has the same sign as x, and the absolute value of f is less than the absolute value of b.
static inline float4 fmod(const float4& a, const float4& b)
{
	return frac(a/b);
}

static inline float4 round(const float4& v)
{
	return float4(round(v.x), round(v.y), round(v.z), round(v.w));
}

//Deadcode VGP 20/02/2008 //Returns the length of the vector v.
//Deadcode VGP 20/02/2008 static inline float length(const float4& v)
//Deadcode VGP 20/02/2008 {
//Deadcode VGP 20/02/2008 	return v.length();
//Deadcode VGP 20/02/2008 }
 
//Returns a + s(b - a). This linearly interpolates between a and b, such that the return value is a when s is 0, and b when s is 1.
static inline float4 lerp(const float4& a, const float4& b, const float4& s)
{
	return a + (b-a) * s;
}

//Returns a + s(b - a). This linearly interpolates between a and b, such that the return value is a when s is 0, and b when s is 1.
static inline float4 lerp(const float4& a, const float4& b, const float& s)
{
	return a + (b-a) * s;
}

//Returns the base-e logarithm of x. If x is negative, the function returns indefinite. If x is 0, the function returns +INF.
static inline float4 log(const float4& v)
{
	return float4(logf(v.x), logf(v.y), logf(v.z), logf(v.w));
}

//Returns the base-10 logarithm of x. If x is negative, the function returns indefinite. If x is 0, the function returns +INF.
static inline float4 log10(const float4& v)
{
	return float4(log10f(v.x), log10f(v.y), log10f(v.z), log10f(v.w));
}

//Returns the base-2 logarithm of x. If x is negative, the function returns indefinite. If x is 0, the function returns +INF.
static inline float4 log2(const float4& v)
{
	return float4(log2(v.x), log2(v.y), log2(v.z), log2(v.w));
}

//Selects the greater of a and b.
static inline float4 max(const float4& a, const float4& b)
{
	return float4(std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z), std::max(a.w, b.w));
}

//Selects the lesser of a and b.
static inline float4 min(const float4& a, const float4& b)
{
	return float4(std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z), std::min(a.w, b.w));
}


//Selects the greatest component
static inline float max(const float4& v)
{
	return std::max(std::max(std::max(v.x, v.y), v.z), v.w);
}

//Selects the least component
static inline float min(const float4& v)
{
	return std::min(std::min(std::min(v.x, v.y), v.z), v.w);
}

//Performs matrix multiplication between a and b. 
//If a is a vector, it is treated as a row vector. 
//If b is a vector, it is treated as a column vector. 
//The inner dimension acolumns and brows must be equal. The result has the dimension arows x bcolumns.
static inline float4 mul(const float4& a, const float4& b)
{
	return a * b;
}

//Returns the normalized vector v / length(v). If the length of v is 0, the result is indefinite.
static inline float4 normalize(const float4& v)
{
	float mag = v.length();

	assert(fabs(mag) > 1e-5);
    
    return v * (1.0f/mag);
}
 
// Returns the reciprocal of v
static inline float4 rcp(const float4& v)
{
    return float4(1.0f / v.x, 1.0f / v.y, 1.0f / v.z, 1.0f / v.w);
}

//Returns x^y.
static inline float4 pow(const float4& a, const float4& b)
{
	return float4(powf(a.x, b.x), powf(a.y, b.y), powf(a.z, b.z), powf(a.w, b.w));
}

static inline float4 pow(const float4& a, float b)
{
	return float4(powf(a.x, b), powf(a.y, b), powf(a.z, b), powf(a.w, b));
}

//Converts x from degrees to radians.
static inline float4 radians(const float4& v)
{
	return float4(v.x / 180.0f * float(M_PI), v.y / 180.0f * float(M_PI), v.z / 180.0f * float(M_PI), v.w / 180.0f * float(M_PI));
}

//Returns the reflection vector v, given the entering ray direction i, and the surface normal n, such that v = i - 2n * (i•n).
static inline float4 reflect(const float4& i, const float4& n)
{
	return i - 2.0f * dot(i, n) * n;
}

//Returns 1 / sqrt(x)
static inline float4 rsqrt(const float4& v)
{
	return float4(1.0f/sqrtf(v.x), 1.0f/sqrtf(v.y), 1.0f/sqrtf(v.z), 1.0f/sqrtf(v.w));
}

//Clamps x to the range [0, 1]
static inline float4 saturate(const float4& v)
{
	return float4(saturate(v.x), saturate(v.y), saturate(v.z), saturate(v.w));
}

//Computes the sign of x. Returns -1 if x is less than 0, 0 if x equals 0, and 1 if x is greater than zero.
static inline float4 sign(const float4& v)
{
	return float4(sign(v.x), sign(v.y), sign(v.z), sign(v.w));
}

//Returns the sine of x
static inline float4 sin(const float4& v)
{
	return float4(sinf(v.x), sinf(v.y), sinf(v.z), sinf(v.w));
}

//Returns the sine and cosine of x. sin(x) is stored in the output parameter s. cos(x) is stored in the output parameter c.
static inline void sincos(const float4& v, float4& s, float4& c)
{
	s = sin(v);
	c = cos(v);
}

//Returns the hyperbolic sine of x
static inline float4 sinh(const float4& v)
{
	return float4(sinhf(v.x), sinhf(v.y), sinhf(v.z), sinhf(v.w));
}


//Returns 0 if x < min. Returns 1 if x > max. Returns a smooth Hermite interpolation between 0 and 1, if x is in the range [min, max].
static inline float4 smoothstep(const float4& vMin, const float4& vMax, const float4& v)
{
	float4 interpolator(saturate(hermite01(v.x)), saturate(hermite01(v.y)), saturate(hermite01(v.z)), saturate(hermite01(v.w)));
	return lerp(vMin, vMax, interpolator);
}

//Square root (per component)
static inline float4 sqrt(const float4& v)
{
	return float4(sqrtf(v.x), sqrtf(v.y), sqrtf(v.z), sqrtf(v.w));
}

//Returns (x >= a) ? 1 : 0
static inline float4 step(const float4& a, const float4& x)
{
	return float4(x.x>=a.x?1.0f:0, x.y>=a.y?1.0f:0, x.z>=a.z?1.0f:0, x.w>=a.w?1.0f:0);
}

//Returns the tangent of x
static inline float4 tan(const float4& v)
{
	return float4(tanf(v.x), tanf(v.y), tanf(v.z), tanf(v.w));
}

//Returns the hyperbolic tangent of x
static inline float4 tanh(const float4& v)
{
	return float4(tanhf(v.x), tanhf(v.y), tanhf(v.z), tanhf(v.w));
}

template<class ByteVector>
ByteVector ConvertDirectionVector(const float2& floatVector)
{
	//#error Unsupported conversion
}

template<class ByteVector>
ByteVector ConvertDirectionVector(const float3& floatVector)
{
	//#error Unsupported conversion
}

template<class ByteVector>
ByteVector ConvertDirectionVector(const float4& floatVector)
{
	//#error Unsupported conversion
}

// https://math.stackexchange.com/questions/1036959/midpoint-of-the-shortest-distance-between-2-rays-in-3d
inline float4 MidPointBetweenRays(const float4& ray1Start, const float4& ray1Dir, const float4& ray2Start, const float4& ray2Dir, float* ray1_t = nullptr, float* ray2_t = nullptr)
{
	const float3 r1o = ray1Start.xyz;
	const float3 r1d = ray1Dir.xyz;
	const float3 r2o = ray2Start.xyz;
	const float3 r2d = ray2Dir.xyz;

	float t1 = ( dot((r2o - r1o), r1d) * dot(r2d, r2d) + ( dot( (r1o - r2o), r2d) ) * dot(r1d, r2d) ) / ( dot(r1d, r1d) * dot(r2d,r2d) - dot(r1d, r2d) * dot(r1d, r2d) );
	float t2 = ( dot((r1o - r2o), r2d) * dot(r1d, r1d) + ( dot( (r2o - r1o), r1d) ) * dot(r1d, r2d) ) / ( dot(r1d, r1d) * dot(r2d,r2d) - dot(r1d, r2d) * dot(r1d, r2d) );


	if (ray1_t)		(*ray1_t) = t1;
	if (ray2_t)		(*ray2_t) = t2;

	float3 m = (r1o + r1d * t1 + r2o + r2d * t2) * 0.5f;

	return float4(m,1);
}

inline float4 RayPlaneIntersection(const float4& rayStart, const float4& rayDir, const float4& planeEquation, float* tPtr = nullptr)
{
	float distanceToPlane = dot(rayStart, planeEquation);

	float t = distanceToPlane / -dot(rayDir.xyz, planeEquation.xyz);

	if (tPtr)
		(*tPtr) = t;

	float4 intersectionPos = rayStart + rayDir * t;

	return intersectionPos;
}

template<class Vector>
inline Vector ClosestPointOnRay(const Vector& rayStart, const Vector& rayDir, const Vector& p, float* tPtr = nullptr)
{
	Vector deltaPos = p-rayStart;

	float t = dot(deltaPos, rayDir) / dot(rayDir, rayDir);

	if (tPtr)
		(*tPtr) = t;

	Vector closestPoint = rayStart + rayDir * t;

	return closestPoint;
}

inline float4 MakePlaneEquation(const float3& planeDirection, const float3& pointOnPlane)
{
	float4 planeEq(normalize(planeDirection), 0);

	planeEq.w = -dot(planeEq.xyz, pointOnPlane);

	return planeEq;
}

// Returns number of points on the +ve side of the plane
inline int ClipEdgeAgainstPlane(const float4& planeEquation, float4& p0, float4& p1)
{
	float dist0 = dot(planeEquation, p0);
	float dist1 = dot(planeEquation, p1);

	// both inside
	if (dist0>=0 && dist1>=0)
		return 2;

	// both outside
	if (dist0<0 && dist1<0)
		return 0;

	// partially inside
	if (dist1>=0)
	{
		p0 += (p1-p0) * -dist0 / (dist1-dist0);
	}
	else
	{
		p1 += (p0-p1) * -dist1 / (dist0-dist1);
	}

	return 1;
}

// Returns number of points on the -ve side of the plane
inline int ClipEdgeAgainstPlaneNeg(const float4& planeEquation, float4& p0, float4& p1)
{
	float dist0 = dot(planeEquation, p0);
	float dist1 = dot(planeEquation, p1);

	// both outside
	if (dist0>0 && dist1>0)
		return 0;

	// both inside
	if (dist0<=0 && dist1<=0)
		return 2;

	// partially inside
	if (dist1>=0)
	{
		p1 += (p0-p1) * dist1 / (dist1-dist0);
	}
	else
	{
		p0 += (p1-p0) * dist0 / (dist0-dist1);
	}

	return 1;
}

inline float3 TriangleNormal(const float3& p0, const float3& p1, const float3& p2)
{
	float3 p01 = sub(p1, p0).unit();
	float3 p02 = sub(p2, p0).unit();
	float3 N = cross(p01, p02).unit();

	return N;
}

inline float DistanceToEdge(const float3& p0, const float3& p1, const float3& p)
{
	float edgeLengthSq = distanceSqr(p0, p1);

	float t = dot(p - p0, p1 - p0) / edgeLengthSq;

	t = saturate(t);

	float3 closestPoint = lerp(p0, p1, t);

	return distance(p, closestPoint);
}
