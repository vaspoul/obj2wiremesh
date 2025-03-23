#pragma once

#include "float3.h"
#include "float4.h"
#include "vectorOps.h"

struct matrix44
{
	static const matrix44 kIdentity;

	union
	{
		struct
		{
			float4 column[4];
		};

		struct
		{
			float	_00, _10, _20, _30;		// Column 0
			float	_01, _11, _21, _31;		// Column 1
			float	_02, _12, _22, _32;		// Column 2
			float	_03, _13, _23, _33;		// Column 3
		};

		struct
		{
			float	_a, _e, _i, _m;	// Column 0
			float	_b, _f, _j, _n;	// Column 1
			float	_c, _g, _k, _o;	// Column 2
			float	_d, _h, _l, _p;	// Column 3
		};

		struct
		{
			// Row-major access
			float	r0, r4,  r8, r12;	// Column 0
			float	r1, r5,  r9, r13;	// Column 1
			float	r2, r6, r10, r14;	// Column 2
			float	r3, r7, r11, r15;	// Column 3
		};
	
		struct
		{
			// Column-major access
			float	c0,   c1,  c2,  c3;	// Column 0
			float	c4,   c5,  c6,  c7;	// Column 1
			float	c8,   c9, c10, c11;	// Column 2
			float	c12, c13, c14, c15;	// Column 3
		};
	};

	matrix44()
		: _00(0), _10(0), _20(0), _30(0)
		, _01(0), _11(0), _21(0), _31(0)
		, _02(0), _12(0), _22(0), _32(0)
		, _03(0), _13(0), _23(0), _33(0)
	{}

	matrix44(	float _00_, float _01_, float _02_, float _03_,		// Row 0
				float _10_, float _11_, float _12_, float _13_,		// Row 1
				float _20_, float _21_, float _22_, float _23_,		// Row 2
				float _30_, float _31_, float _32_, float _33_)		// Row 3
		: _00(_00_), _10(_10_), _20(_20_), _30(_30_)
		, _01(_01_), _11(_11_), _21(_21_), _31(_31_)
		, _02(_02_), _12(_12_), _22(_22_), _32(_32_)
		, _03(_03_), _13(_13_), _23(_23_), _33(_33_)
	{}

	matrix44(const matrix44& M) = default;
	//{
	//	_00 = M._00;	_10 = M._10;	_20 = M._20;	_30 = M._30;
	//	_01 = M._01;	_11 = M._11;	_21 = M._21;	_31 = M._31;
	//	_02 = M._02;	_12 = M._12;	_22 = M._22;	_32 = M._32;
	//	_03 = M._03;	_13 = M._13;	_23 = M._23;	_33 = M._33;
	//}

	//void operator=(const matrix44& M)
	//{
	//	_00 = M._00;	_10 = M._10;	_20 = M._20;	_30 = M._30;
	//	_01 = M._01;	_11 = M._11;	_21 = M._21;	_31 = M._31;
	//	_02 = M._02;	_12 = M._12;	_22 = M._22;	_32 = M._32;
	//	_03 = M._03;	_13 = M._13;	_23 = M._23;	_33 = M._33;
	//}

	static matrix44 FromRows(const float4& row1, const float4& row2, const float4& row3, const float4& row4)
	{
		matrix44 temp;
		
		temp._00 = row1.x;
		temp._01 = row1.y;
		temp._02 = row1.z;
		temp._03 = row1.w;

		temp._10 = row2.x;
		temp._11 = row2.y;
		temp._12 = row2.z;
		temp._13 = row2.w;

		temp._20 = row3.x;
		temp._21 = row3.y;
		temp._22 = row3.z;
		temp._23 = row3.w;

		temp._30 = row4.x;
		temp._31 = row4.y;
		temp._32 = row4.z;
		temp._33 = row4.w;

		return temp;
	}

	static matrix44 FromColumns(const float4& c1, const float4& c2, const float4& c3, const float4& c4)
	{
		matrix44 temp;

		temp._00 = c1.x;
		temp._10 = c1.y;
		temp._20 = c1.z;
		temp._30 = c1.w;

		temp._01 = c2.x;
		temp._11 = c2.y;
		temp._21 = c2.z;
		temp._31 = c2.w;

		temp._02 = c3.x;
		temp._12 = c3.y;
		temp._22 = c3.z;
		temp._32 = c3.w;

		temp._03 = c4.x;
		temp._13 = c4.y;
		temp._23 = c4.z;
		temp._33 = c4.w;

		return temp;
	}

	void Zero()
	{
		_00 = _10 = _20 = _30 = 0;
		_01 = _11 = _21 = _31 = 0;
		_02 = _12 = _22 = _32 = 0;
		_03 = _13 = _23 = _33 = 0;
	}

	static matrix44 MakeDiagonal(float x, float y, float z, float w = 1)
	{
		return matrix44(	x,0,0,0,
							0,y,0,0,
							0,0,z,0,
							0,0,0,w);
	}

	static matrix44 MakeDiagonal(const float4& diag)
	{
		return matrix44(	diag.x,0,0,0,
							0,diag.y,0,0,
							0,0,diag.z,0,
							0,0,0,diag.w);
	}

	static matrix44 MakeDiagonal(const float3& diag)
	{
		return matrix44(	diag.x,0,0,0,
							0,diag.y,0,0,
							0,0,diag.z,0,
							0,0,0,1.0f);
	}

	matrix44& SetRotationX(float angleDegrees)
	{
		float angleRadians = radians(angleDegrees);

		float sf = sinf(angleRadians);
		float cf = cosf(angleRadians);

		_00 = 1;	_01 = 0;	_02 = 0;
		_10 = 0;	_11 = cf;	_12 = sf;
		_20 = 0;	_21 = -sf;	_22 = cf;

		return *this;
	}

	matrix44& SetRotationY(float angleDegrees)
	{
		float angleRadians = radians(angleDegrees);

		float sf = sinf(angleRadians);
		float cf = cosf(angleRadians);

		_00 = cf;	_01 = 0;	_02 = -sf;
		_10 = 0;	_11 = 1;	_12 = 0;
		_20 = sf;	_21 = 0;	_22 = cf;

		return *this;
	}

	matrix44& SetRotationZ(float angleDegrees)
	{
		float angleRadians = radians(angleDegrees);

		float sf = sinf(angleRadians);
		float cf = cosf(angleRadians);

		_00 = cf;	_01 = sf;	_02 = 0;
		_10 = -sf;	_11 = cf;	_12 = 0;
		_20 = 0;	_21 = 0;	_22 = 1;

		return *this;
	}

	matrix44& SetRotationXYZ(float angleXdegrees, float angleYdegrees, float angleZdegrees)
	{
		float radians_ax = radians(angleXdegrees);
		float radians_ay = radians(angleYdegrees);
		float radians_az = radians(angleZdegrees);

		float sx = sinf(radians_ax);
		float cx = cosf(radians_ax);

		float sy = sinf(radians_ay);
		float cy = cosf(radians_ay);

		float sz = sinf(radians_az);
		float cz = cosf(radians_az);

		// YXZ order (yaw, pitch, roll)

		_00	= cy * cz - sy * sx * sz;
		_01 = sy * sx * cz + sz * cy;
		_02 = -sy * cx;

		_10 = -sz * cx;
		_11 = cx * cz;
		_12 = sx;

		_20 = sy * cz + sx * sz * cy;
		_21 = sy * sz - sx * cy * cz;
		_22 = cy * cx;

		return *this;
	}

	static matrix44 MakeRotationXYZ(float3 angleXYZdegrees)
	{
		return MakeRotationXYZ(angleXYZdegrees.x, angleXYZdegrees.y, angleXYZdegrees.z);
	}

	static matrix44 MakeRotationXYZ(float angleXdegrees, float angleYdegrees, float angleZdegrees)
	{
		matrix44 mtx = kIdentity;
		mtx.SetRotationXYZ(angleXdegrees, angleYdegrees, angleZdegrees);
		return mtx;
	}

	static matrix44 MakeRotationX(float angleDegrees)
	{
		matrix44 mtx = kIdentity;
		mtx.SetRotationX(angleDegrees);
		return mtx;
	}

	static matrix44 MakeRotationY(float angleDegrees)
	{
		matrix44 mtx = kIdentity;
		mtx.SetRotationY(angleDegrees);
		return mtx;
	}

	static matrix44 MakeRotationZ(float angleDegrees)
	{
		matrix44 mtx = kIdentity;
		mtx.SetRotationZ(angleDegrees);
		return mtx;
	}

	static matrix44 MakeScale(float3 scales)
	{
		return MakeDiagonal(scales);
	}

	static matrix44 MakeScale(float sx, float sy, float sz)
	{
		return MakeDiagonal(sx, sy, sz);
	}

	static matrix44 MakeScaleBias(float3 scale, float3 offset)
	{
		return MakeDiagonal(scale).SetTranslation(offset);
	}

	static matrix44 MakeScaleBias(float sx, float sy, float sz, float ox, float oy, float oz)
	{
		return MakeDiagonal(sx, sy, sz).SetTranslation(ox, oy, oz);
	}

	matrix44& SetTranslation(float3 t)
	{
		return SetTranslation(t.x, t.y, t.z);
	}

	matrix44& SetTranslation(float tx, float ty, float tz)
	{
		_30 = tx;
		_31 = ty;
		_32 = tz;

		return *this;
	}

	static matrix44 MakeTranslation(float tx, float ty, float tz)
	{
		matrix44 mtx = kIdentity;
		mtx.SetTranslation(tx, ty, tz);
		return mtx;
	}

	static matrix44 MakeTranslation(float3 translation)
	{
		matrix44 mtx = kIdentity;
		mtx.SetTranslation(translation);
		return mtx;
	}

	matrix44& SetCameraView(const float3& cameraPos, const float3& cameraDirection, const float3& cameraUp)
	{
		float3 x = normalize(cameraUp % cameraDirection);
		float3 y = normalize(cameraDirection % x);

		_00 = x.x;					_01 = y.x;					_02 = cameraDirection.x;					_03 = 0;
		_10 = x.y;					_11 = y.y;					_12 = cameraDirection.y;					_13 = 0;
		_20 = x.z;					_21 = y.z;					_22 = cameraDirection.z;					_23 = 0;
		_30 = -dot(cameraPos, x);	_31 = -dot(cameraPos, y);	_32 = -dot(cameraPos, cameraDirection);		_33 = 1;

		return *this;
	}

	matrix44& LookAt(const float3& cameraPos, const float3& cameraTarget, const float3& cameraUp)
	{
		return SetCameraView(cameraPos, normalize(cameraTarget-cameraPos), cameraUp);
	}

	matrix44& LookAlong(const float3& cameraPos, const float3& cameraDirection, const float3& cameraUp)
	{
		return SetCameraView(cameraPos, normalize(cameraDirection), cameraUp);
	}

	static matrix44 MakeLookAt(const float3& cameraPos, const float3& cameraTarget, const float3& cameraUp)
	{
		matrix44 mtx;
		mtx.LookAt(cameraPos, cameraTarget, cameraUp);
		return mtx;
	}

	static matrix44 MakeLookAlong(const float3& cameraPos, const float3& cameraDirection, const float3& cameraUp)
	{
		matrix44 mtx;
		mtx.LookAlong(cameraPos, cameraDirection, cameraUp);
		return mtx;
	}

	matrix44& SetPerspective(float fovYdegrees, float aspect, float N, float F)
	{
		float yScale = 1.0f / ::tan(radians(fovYdegrees) * 0.5f);
		float xScale = yScale / aspect;

		column[0] = float4(xScale,    0,			0,				0			);
		column[1] = float4(  0,		yScale,		    0,				0			);
		column[2] = float4(  0,		  0,		F / (F-N),		-F*N / (F-N)	);
		column[3] = float4(  0,		  0,			1,				0			);

		return *this;
	}

	matrix44& SetPerspectiveFromTangents(float tanX, float tanY, float N, float F)
	{
		float xScale = 1.0f / tanX;
		float yScale = 1.0f / tanY;

		column[0] = float4(xScale,    0,			0,				0			);
		column[1] = float4(  0,		yScale,		    0,				0			);
		column[2] = float4(  0,		  0,		F / (F-N),		-F*N / (F-N)	);
		column[3] = float4(  0,		  0,			1,				0			);

		return *this;
	}

	matrix44& SetPerspective(float L, float R, float T, float B, float N, float F)
	{
		//	2*zn/(r-l)    0             0				0
		//	0             2*zn/(t-b)    0               0
		//	-(l+r)/(r-l)  -(t+b)/(t-b)  zf/(zf-zn)     1
		//	0             0             -zn*zf/(zf-zn)  0

		column[0] = float4(2*N/(R-L),			0,	   -(L+R)/(R-L),			   0	);
		column[1] = float4(        0,	2*N/(T-B),	   -(T+B)/(T-B),			   0	);
		column[2] = float4(        0,			0,		  F / (F-N),	-F*N / (F-N)	);
		column[3] = float4(        0,			0,				  1,			   0	);

		return *this;
	}

	static matrix44 MakePerspective(float fovYdegrees, float aspect, float N, float F)
	{
		matrix44 mtx;
		mtx.SetPerspective(fovYdegrees, aspect, N, F);
		return mtx;
	}

	static matrix44 MakePerspective(float L, float R, float T, float B, float N, float F)
	{
		matrix44 mtx;
		mtx.SetPerspective(L, R, T, B, N, F);
		return mtx;
	}

	matrix44& SetOrtho(float W, float H, float N, float F)
	{
		column[0] = float4(2.0f/W, 0, 0, 0);
		column[1] = float4(0, 2.0f/H, 0, 0);
		column[2] = float4(0, 0, 1.0f / (F-N), -N / (F-N));
		column[3] = float4(0, 0, 0, 1);
	
		return *this;
	}

	matrix44& SetOrtho(float L, float R, float T, float B, float N, float F)
	{
		column[0] = float4(2.0f/(R-L), 0, 0, (L+R)/(L-R));
		column[1] = float4(0, 2.0f/(T-B), 0, (T+B)/(B-T));
		column[2] = float4(0, 0, 1.0f / (F-N), -N / (F-N));
		column[3] = float4(0, 0, 0, 1);

		return *this;
	}

	static matrix44 MakeOrtho(float W, float H, float N, float F)
	{
		matrix44 mtx;
		mtx.SetOrtho(W, H, N, F);
		return mtx;
	}

	static matrix44 MakeOrtho(float L, float R, float T, float B, float N, float F)
	{
		matrix44 mtx;
		mtx.SetOrtho(L, R, T, B, N, F);
		return mtx;
	}

	matrix44& SetTransform(float3 position, float3 rotationAnglesDegrees, float3 scale)
	{
		matrix44 scaleMtx = MakeDiagonal(scale);
		matrix44 rotationMtx = MakeRotationXYZ(rotationAnglesDegrees);
		matrix44 translationMtx = MakeTranslation(position);

		(*this) = scaleMtx * rotationMtx * translationMtx;

		return *this;
	}

	static matrix44 MakeTransform(float3 position, float3 rotationAnglesDegrees, float3 scale)
	{
		matrix44 mtx = kIdentity;
		mtx.SetTransform(position, rotationAnglesDegrees, scale);
		return mtx;
	}

	matrix44& SetTransformFromDir(float3 position, float3 direction, float3 up)
	{
		float3 xAxis = normalize( cross(up, direction) );

		this->_00 = xAxis.x;
		this->_01 = xAxis.y;
		this->_02 = xAxis.z;

		this->_10 = up.x;
		this->_11 = up.y;
		this->_12 = up.z;

		this->_20 = direction.x;
		this->_21 = direction.y;
		this->_22 = direction.z;

		this->SetTranslation(position);

		return *this;
	}

	static matrix44 MakeTransformFromDir(float3 position, float3 direction, float3 up)
	{
		matrix44 mtx = kIdentity;
		mtx.SetTransformFromDir(position, direction, up);
		return mtx;
	}

	float4 Row(int index) const
	{
		switch(index)
		{
			case 0:	return float4(_00, _01, _02, _03); break;
			case 1:	return float4(_10, _11, _12, _13); break;
			case 2:	return float4(_20, _21, _22, _23); break;
			case 3:	return float4(_30, _31, _32, _33); break;
		}
		
		return float4();
	}

	const float4& Column(int index) const
	{
		return column[index];
	}

	float4& Column(int index)
	{
		return column[index];
	}

	matrix44 operator-() const
	{
		return matrix44(	-_00, -_01, -_02, -_03,
							-_10, -_11, -_12, -_13,
							-_20, -_21, -_22, -_23,
							-_30, -_31, -_32, -_33);

	}

	matrix44 operator+(const matrix44& M) const
	{
	    return matrix44(	_00 + M._00, _01 + M._01, _02 + M._02, _03 + M._03,
							_10 + M._10, _11 + M._11, _12 + M._12, _13 + M._13,
							_20 + M._20, _21 + M._21, _22 + M._22, _23 + M._23,
							_30 + M._30, _31 + M._31, _32 + M._32, _33 + M._33);
	}

	matrix44 operator-(const matrix44& M) const
	{
	    return matrix44(	_00 - M._00, _01 - M._01, _02 - M._02, _03 - M._03,
							_10 - M._10, _11 - M._11, _12 - M._12, _13 - M._13,
							_20 - M._20, _21 - M._21, _22 - M._22, _23 - M._23,
							_30 - M._30, _31 - M._31, _32 - M._32, _33 - M._33);
	}

	matrix44 operator*(const matrix44& M) const
	{
		matrix44 temp;
		temp._00 = _00 * M._00 + _01 * M._10 + _02 * M._20 + _03 * M._30;
		temp._01 = _00 * M._01 + _01 * M._11 + _02 * M._21 + _03 * M._31;
		temp._02 = _00 * M._02 + _01 * M._12 + _02 * M._22 + _03 * M._32;
		temp._03 = _00 * M._03 + _01 * M._13 + _02 * M._23 + _03 * M._33;

		temp._10 = _10 * M._00 + _11 * M._10 + _12 * M._20 + _13 * M._30;
		temp._11 = _10 * M._01 + _11 * M._11 + _12 * M._21 + _13 * M._31;
		temp._12 = _10 * M._02 + _11 * M._12 + _12 * M._22 + _13 * M._32;
		temp._13 = _10 * M._03 + _11 * M._13 + _12 * M._23 + _13 * M._33;

		temp._20 = _20 * M._00 + _21 * M._10 + _22 * M._20 + _23 * M._30;
		temp._21 = _20 * M._01 + _21 * M._11 + _22 * M._21 + _23 * M._31;
		temp._22 = _20 * M._02 + _21 * M._12 + _22 * M._22 + _23 * M._32;
		temp._23 = _20 * M._03 + _21 * M._13 + _22 * M._23 + _23 * M._33;

		temp._30 = _30 * M._00 + _31 * M._10 + _32 * M._20 + _33 * M._30;
		temp._31 = _30 * M._01 + _31 * M._11 + _32 * M._21 + _33 * M._31;
		temp._32 = _30 * M._02 + _31 * M._12 + _32 * M._22 + _33 * M._32;
		temp._33 = _30 * M._03 + _31 * M._13 + _32 * M._23 + _33 * M._33;

		return temp;
	}

	matrix44 operator*(float v) const
	{
		return matrix44(	v * _00, v * _01, v * _02, v * _03,
							v * _10, v * _11, v * _12, v * _13,
							v * _20, v * _21, v * _22, v * _23,
							v * _30, v * _31, v * _32, v * _33);
	}

	matrix44 operator/(float v) const
	{
		assert(fabs(v) > 1e-5);

		v = 1.0f / v;

		return matrix44(	v * _00, v * _01, v * _02, v * _03,
							v * _10, v * _11, v * _12, v * _13,
							v * _20, v * _21, v * _22, v * _23,
							v * _30, v * _31, v * _32, v * _33);
	}

	void operator-=(const matrix44& M)
	{
		*this = *this - M;
	}

	void operator+=(const matrix44& M)
	{
		*this = *this + M;
	}

	void operator*=(const matrix44& M)
	{
		*this = *this * M;
	}

	void operator*=(float v)
	{
		*this = *this * v;
	}

	void operator/=(float v)
	{
		*this = *this / v;
	}

	float4 GetXAxis() const
	{
		return float4(_00, _01, _02, 0);
	}

	float4 GetYAxis() const
	{
		return float4(_10, _11, _12, 0);
	}

	float4 GetZAxis() const
	{
		return float4(_20, _21, _22, 0);
	}

	float4 GetViewRight() const
	{
		return float4(_00, _10, _20, 0);
	}

	float4 GetViewUp() const
	{
		return float4(_01, _11, _21, 0);
	}

	float4 GetViewDir() const
	{
		return float4(_02, _12, _22, 0);
	}

	float3 GetTranslation() const
	{
		return float3(_30, _31, _32);
	}

	matrix44 GetViewAlignedWorldMtx(const float3& worldPos) const
	{
		float3 camPos = inverse43().GetTranslation();

		float3 zAxis = normalize(worldPos - camPos);
		float3 xAxis = normalize(cross(zAxis, GetViewUp().xyz));
		float3 yAxis = normalize(cross(zAxis, xAxis));

		return matrix44::FromRows( float4(xAxis,0), float4(yAxis,0), float4(zAxis,0), float4::k0001 );
	}

	void Decompose(float3* translation, float3* rotationDegrees, float3* scale)
	{
		if (translation)
		{
			translation->x = _30;
			translation->y = _31;
			translation->z = _32;
		}

		if (scale)
		{
			scale->x = column[0].xyz.length();
			scale->y = column[1].xyz.length();
			scale->z = column[2].xyz.length();
		}

		if (rotationDegrees)
		{
			float sy = sqrt(_12 * _12 + _22 * _22);

			rotationDegrees->x = atan2( _12, _22);
			rotationDegrees->y = atan2(-_02,  sy);
			rotationDegrees->z = atan2( _01, _00);

			*rotationDegrees = degrees(*rotationDegrees);
		}
	}

	void DecomposeProjection(float& FOVRadians, float& aspectRatio, float& N, float& F) const
	{
		//	2*zn/w    0         0				0
		//	0         2*zn/h    0               0
		//	0		  0			zf/(zf-zn)      1
		//	0         0         -zn*zf/(zf-zn)  0

		N = -_32 / _22;
		F = -_32 / (_22 - 1.0f);

		aspectRatio = _11 / _00;

		FOVRadians = atan(1.0f / _11) * 2.0f;
	}

	void DecomposeProjectionAssymetric(float& L, float& R, float& T, float& B, float& N, float& F) const
	{
		//	2*zn/(r-l)    0             0				0
		//	0             2*zn/(t-b)    0               0
		//	-(l+r)/(r-l)  -(t+b)/(t-b)  z f/(zf-zn)     1
		//	0             0             -zn*zf/(zf-zn)  0

		N = -_32 / _22;
		F = N * _22 / (_22 - 1.0f);

		float r_minus_l = 2.0f * N / _00;
		float r_plus_l = -_20 * r_minus_l;
		R = 0.5f * (r_plus_l + r_minus_l);
		L = r_plus_l - R;

		float t_minus_b = 2.0f * N / _11;
		float t_plus_b = -_21 * t_minus_b;
		T = 0.5f * (t_plus_b + t_minus_b);
		B = t_plus_b - T;
	}

	void ExtractProjectionParams(float4& projectionParams0, float4& projectionParams1) const
	{
		bool isOrtho = _23 == 0;
		float N = -_32 / _22;
		float F = isOrtho ? (1.0f - _32 / _22 ) : ( -_32 / (_22 - 1.0f) );

		projectionParams0 = float4( 1.0f / _00, 1.0f / _11, _22, _32 );
		projectionParams1 = float4( _30, _31, N, F );
	}

	matrix44 transpose() const
	{
		return matrix44(	_00, _10, _20, _30,
							_01, _11, _21, _31,
							_02, _12, _22, _32,
							_03, _13, _23, _33);
	}

	matrix44 abs() const
	{
		return matrix44(	fabs(_00), fabs(_10), fabs(_20), fabs(_30),
							fabs(_01), fabs(_11), fabs(_21), fabs(_31),
							fabs(_02), fabs(_12), fabs(_22), fabs(_32),
							fabs(_03), fabs(_13), fabs(_23), fabs(_33) );
	}

	// [a  b  c  d]
	// [e  f  g  h]
	// [i  j  k  l]
	// [m  n  o  p]
	inline matrix44 inverse() const
	{
		// 3x3 determinant
		// [a  b  c]
		// [d  e  f]
		// [g  h  i]
		#define determinant(a,b,c,d,e,f,g,h,i)	+a * (e*i - f*h)	\
												-b * (d*i - f*g)	\
												+c * (d*h - e*g)

		float deta = determinant(	_f, _g, _h, 
									_j, _k, _l, 
									_n, _o, _p);

		float detb = determinant(	_e, _g, _h, 
									_i, _k, _l, 
									_m, _o, _p);

		float detc = determinant(	_e, _f, _h, 
									_i, _j, _l, 
									_m, _n, _p);

		float detd = determinant(	_e, _f, _g, 
									_i, _j, _k, 
									_m, _n, _o);

		// ------

		float dete = determinant(	_b, _c, _d, 
									_j, _k, _l, 
									_n, _o, _p);

		float detf = determinant(	_a, _c, _d, 
									_i, _k, _l, 
									_m, _o, _p);

		float detg = determinant(	_a, _b, _d, 
									_i, _j, _l, 
									_m, _n, _p);

		float deth = determinant(	_a, _b, _c, 
									_i, _j, _k, 
									_m, _n, _o);

		// ------

		float deti = determinant(	_b, _c, _d, 
									_f, _g, _h, 
									_n, _o, _p);

		float detj = determinant(	_a, _c, _d, 
									_e, _g, _h, 
									_m, _o, _p);

		float detk = determinant(	_a, _b, _d, 
									_e, _f, _h, 
									_m, _n, _p);

		float detl = determinant(	_a, _b, _c,  
									_e, _f, _g,  
									_m, _n, _o);

		// ------

		float detm = determinant(	_b, _c, _d, 
									_f, _g, _h, 
									_j, _k, _l);

		float detn = determinant(	_a, _c, _d, 
									_e, _g, _h, 
									_i, _k, _l);

		float deto = determinant(	_a, _b, _d, 
									_e, _f, _h, 
									_i, _j, _l);

		float detp = determinant(	_a, _b, _c, 
									_e, _f, _g, 
									_i, _j, _k);

		// ------

		float det = +_a * deta
					-_b * detb
					+_c * detc
					-_d * detd;

		float ood = 1.0f / det;

		matrix44 out;
		
		// transpose as well
		out.c0 = +ood * deta;
		out.c1 = -ood * detb;
		out.c2 = +ood * detc;
		out.c3 = -ood * detd;
		out.c4 = -ood * dete;
		out.c5 = +ood * detf;
		out.c6 = -ood * detg;
		out.c7 = +ood * deth;
		out.c8 = +ood * deti;
		out.c9 = -ood * detj;
		out.c10 = +ood * detk;
		out.c11 = -ood * detl;
		out.c12 = -ood * detm;
		out.c13 = +ood * detn;
		out.c14 = -ood * deto;
		out.c15 = +ood * detp;

		return out;
	}

	matrix44 inverse43() const
	{
		return matrix44(	_00, _10, _20, 0,
							_01, _11, _21, 0,
							_02, _12, _22, 0,
							dot(-float3(_30, _31, _32), float3(_00, _01, _02)), dot(-float3(_30, _31, _32), float3(_10, _11, _12)), dot(-float3(_30, _31, _32), float3(_20, _21, _22)), 1);
	}

	bool operator==(const matrix44& M) const
	{
		return	fabs(M.c0 - c0) < 1e-5 &&
				fabs(M.c1 - c1) < 1e-5 &&
				fabs(M.c2 - c2) < 1e-5 &&
				fabs(M.c3 - c3) < 1e-5 &&
				fabs(M.c4 - c4) < 1e-5 &&
				fabs(M.c5 - c5) < 1e-5 &&
				fabs(M.c6 - c6) < 1e-5 &&
				fabs(M.c7 - c7) < 1e-5 &&
				fabs(M.c8 - c8) < 1e-5 &&
				fabs(M.c9 - c9) < 1e-5 &&
				fabs(M.c10 - c10) < 1e-5 &&
				fabs(M.c11 - c11) < 1e-5 &&
				fabs(M.c12 - c12) < 1e-5 &&
				fabs(M.c13 - c13) < 1e-5 &&
				fabs(M.c14 - c14) < 1e-5 &&
				fabs(M.c15 - c15) < 1e-5;
	}

	bool operator!=(const matrix44& M) const
	{
		return !((*this) == M);
	}

	template<typename Archive>
	void Serialize(Archive& archive)
	{
		if (archive.IsYAML())
		{
			float4 rows[4];

			if (archive.IsWriting())
			{
				rows[0] = Row(0);
				rows[1] = Row(1);
				rows[2] = Row(2);
				rows[3] = Row(3);
			}

			archive << MEMBER("Row0", rows[0]);
			archive << MEMBER("Row1", rows[1]);
			archive << MEMBER("Row2", rows[2]);
			archive << MEMBER("Row3", rows[3]);

			if (archive.IsReading())
			{
				(*this) = FromRows(rows[0], rows[1], rows[2], rows[3]);
			}
		}
		else
		{
			archive << MEMBER("Column0", column[0]);
			archive << MEMBER("Column1", column[1]);
			archive << MEMBER("Column2", column[2]);
			archive << MEMBER("Column3", column[3]);
		}
	}
};

typedef matrix44 float4x4;

#if _MSC_VER >= 1920 // VS 2019 onwards
inline const matrix44 matrix44::kIdentity = matrix44(	1,0,0,0, 
														0,1,0,0,
														0,0,1,0,
														0,0,0,1);
#endif

inline matrix44 transpose(const matrix44& m)
{
	return matrix44(	m._00, m._10, m._20, m._30,
						m._01, m._11, m._21, m._31,
						m._02, m._12, m._22, m._32,
						m._03, m._13, m._23, m._33);
}

inline matrix44 abs(const matrix44& m)
{
	matrix44 temp;
	temp.column[0] = abs(m.column[0]);
	temp.column[1] = abs(m.column[1]);
	temp.column[2] = abs(m.column[2]);
	temp.column[3] = abs(m.column[3]);

	return temp;
}

/// Matrix/float4 multiply.
inline float4 operator*(const float4& v, const matrix44& m)
{
    return float4(	m._00 * v.x + m._10 * v.y + m._20 * v.z + m._30 * v.w,
					m._01 * v.x + m._11 * v.y + m._21 * v.z + m._31 * v.w,
					m._02 * v.x + m._12 * v.y + m._22 * v.z + m._32 * v.w,
					m._03 * v.x + m._13 * v.y + m._23 * v.z + m._33 * v.w);
}

//inline float4 operator*(const float3& v, const matrix44& m)
//{
//    return float4(	m._00 * v.x + m._10 * v.y + m._20 * v.z + m._30,
//					m._01 * v.x + m._11 * v.y + m._21 * v.z + m._31,
//					m._02 * v.x + m._12 * v.y + m._22 * v.z + m._32,
//					m._03 * v.x + m._13 * v.y + m._23 * v.z + m._33);
//}

inline float4 operator*(const matrix44& m, const float4& v)
{
    return float4(	m._00 * v.x + m._01 * v.y + m._02 * v.z + m._03 * v.w,
					m._10 * v.x + m._11 * v.y + m._12 * v.z + m._13 * v.w,
					m._20 * v.x + m._21 * v.y + m._22 * v.z + m._23 * v.w,
					m._30 * v.x + m._31 * v.y + m._32 * v.z + m._33 * v.w);
}

inline matrix44 inverse(const matrix44& m)
{
	return m.inverse();
}
