
inline float4 mul(const float4& v, const float4x4& m)	{ return v * m; }
//inline float4 mul(const float3& v, const float4x4& m)	{ return v * m; }
inline float4 mul(const float4x4& m, const float4& v)	{ return m * v; }
inline float4x4 mul(const float4x4& a, const float4x4& b)	{ return a * b; }

inline matrix44 inverse43(const matrix44& m)
{
	float tx = dot(-float3(m._30, m._31, m._32), float3(m._00, m._01, m._02));
	float ty = dot(-float3(m._30, m._31, m._32), float3(m._10, m._11, m._12));
	float tz = dot(-float3(m._30, m._31, m._32), float3(m._20, m._21, m._22));

	return matrix44(	m._00, m._10, m._20, 0,
						m._01, m._11, m._21, 0, 
						m._02, m._12, m._22, 0,
						tx,    ty,    tz,    1);
}

inline float4 BackProject(const float2& uv, float depth, const matrix44& viewMtx, const matrix44& projMtx)
{
	float2 ndcPos = uv * float2(2, -2) + float2(-1, +1);

	float4 viewSpacePos;

	viewSpacePos.x = (ndcPos.x * depth - projMtx._30) / projMtx._00;
	viewSpacePos.y = (ndcPos.y * depth - projMtx._31) / projMtx._11;
	viewSpacePos.z = depth;
	viewSpacePos.w = 1;

	float4 worldSpacePos = mul(viewSpacePos, inverse43(viewMtx));

	return worldSpacePos;
}
