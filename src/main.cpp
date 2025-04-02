#include "global.h"
#include "AppInterface.h"

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

tinyobj::attrib_t					g_obj_attrib;
std::vector<tinyobj::shape_t>		g_obj_shapes;
std::vector<tinyobj::material_t>	g_obj_materials;
float3								g_modelCenter;
float3								g_modelSize;
float								g_avgTriangleArea = 0;
float								g_printbedSize = 210;

int									g_pointsPerTriangle = 20;
bool								g_stratifiedPoints = 1;
bool								g_irrationalPoints = 1;
std::vector<float3>					g_points;
std::vector<float3>					g_controlPoints;
std::vector<uint32_t>				g_triangleIndices;
float								g_processingTime = 0;

float								g_cellRadius = 5;
int									g_tubeSegmentCount = 20;
int									g_tubeAngleCount = 16;
float								g_tubeRadius = 0.5;

int									g_maxCellBoundaryIndex = -1;
bool								g_showJustOne = false;
bool								g_showControlPoints = false;
bool								g_showPoints = false;
bool								g_showCellsAxes = false;
bool								g_showCellsPoints = false;
bool								g_showCellsTubes = false;

float TriangleArea(const float3& v0, const float3& v1, const float3& v2)
{
	float3 crossProd = cross(v1-v0, v2-v0);
	return 0.5f * std::sqrt(crossProd.x * crossProd.x + crossProd.y * crossProd.y + crossProd.z * crossProd.z);
}

// Solve linear system using Gaussian elimination
void gaussianElimination(std::vector<std::vector<float>>& A, std::vector<float>& b)
{
	int n = (int)A.size();

	for (int i = 0; i < n; ++i)
	{
		// Partial pivoting
		int maxRow = i;
		for (int k = i + 1; k < n; ++k)
		{
			if (fabs(A[k][i]) > fabs(A[maxRow][i]))
			{
				maxRow = k;
			}
		}
		std::swap(A[i], A[maxRow]);
		std::swap(b[i], b[maxRow]);

		// Elimination
		for (int k = i + 1; k < n; ++k)
		{
			float factor = A[k][i] / A[i][i];
			for (int j = i; j < n; ++j)
			{
				A[k][j] -= factor * A[i][j];
			}
			b[k] -= factor * b[i];
		}
	}

	// Back substitution
	for (int i = n - 1; i >= 0; --i)
	{
		b[i] /= A[i][i];
		for (int j = 0; j < i; ++j)
		{
			b[j] -= A[j][i] * b[i];
		}
	}
}

std::vector<float> polynomialFit(const std::vector<float2>& v, int degree = 3)
{
	int n = (int)v.size();
	std::vector<std::vector<float>> A(degree + 1, std::vector<float>(degree + 1, 0.0f));
	std::vector<float> b(degree + 1, 0.0f);

	for (int i = 0; i <= degree; ++i)
	{
		for (int j = 0; j <= degree; ++j)
		{
			for (int k = 0; k < n; ++k)
			{
				A[i][j] += powf(v[k].x, (float)(i + j));
			}
		}

		for (int k = 0; k < n; ++k)
		{
			b[i] += v[k].y * powf(v[k].x, (float)i);
		}
	}

	gaussianElimination(A, b);

	return b;
}

float EvaluatePolynomial(const std::vector<float>& coeff, float x)
{
	float result = 0.0f;
	float power = 1.0f;

	for (float c : coeff)
	{
		result += c * power;
		power *= x;
	}

	return result;
}

float EvaluatePolynomialSlope(const std::vector<float>& coeff, float x)
{
	float result = 0.0f;
	float power = 1.0f;

	for (int i = 1; i < coeff.size(); ++i)
	{
		result += i * coeff[i] * power;
		power *= x;
	}

	return result;
}

float3 RandomPointInTriangle(const float3& v0, const float3& v1, const float3& v2, uint32_t i, uint32_t iMax)
{
	static std::random_device rd;
	static std::mt19937 gen(rd());
	static std::uniform_real_distribution<float> dist(0.0f, 1.0f);

	float u = dist(gen);
	float v = dist(gen);

	if (g_irrationalPoints)
	{
		//const float goldenRatio = (1 + sqrtf(5.0f)) / 2.0f;

		//// Fibonacci Spiral
		//float phi = (float)M_PI * 2.0f * i / goldenRatio;
		//float r = sqrtf((float)i / (float)iMax);
		//	
		//u = cosf(phi) * r * 0.5f + 0.5f;
		//v = sinf(phi) * r * 0.5f + 0.5f;

		auto radicalInverse = [](uint32_t bits) -> float
		{
			bits = (bits << 16) | (bits >> 16);
			bits = ((bits & 0x55555555) << 1) | ((bits & 0xAAAAAAAA) >> 1);
			bits = ((bits & 0x33333333) << 2) | ((bits & 0xCCCCCCCC) >> 2);
			bits = ((bits & 0x0F0F0F0F) << 4) | ((bits & 0xF0F0F0F0) >> 4);
			bits = ((bits & 0x00FF00FF) << 8) | ((bits & 0xFF00FF00) >> 8);

			return (float)(bits) * (float)2.3283064365386963e-10; // / 0x100000000
		};

		u = (float)i / (float)iMax;
		v = radicalInverse(i);
	}
	else if (g_stratifiedPoints)
	{
		uint32_t iMaxX = std::max(1u, (uint32_t)sqrtf((float)iMax));
	
		u = float(i % iMaxX) / float(iMaxX);
		v = float( (i / iMaxX) % iMaxX ) / float(iMaxX);

		u = saturate(u);
		v = saturate(v);
	}

	if (u + v > 1.0f)
	{
		u = 1.0f - u;
		v = 1.0f - v;
	}

	float w = 1.0f - u - v;
	
	return v0 * u + v1 * v + v2 * w;
}

float3 CircumCenter(const float3& a, const float3& b, const float3& c)
{
	float3 ab = b - a;
	float3 ac = c - a;

	float3 ab_cross_ac = cross(ab, ac);
	float denominator = 2.0f * (ab_cross_ac.x * ab_cross_ac.x + ab_cross_ac.y * ab_cross_ac.y + ab_cross_ac.z * ab_cross_ac.z);

	float3 ab_mid = avg(a, b);
	float3 ac_mid = avg(a, c);

	float t1 = ((ac.y * ab_cross_ac.z - ac.z * ab_cross_ac.y) * ab_mid.x + (ac.z * ab_cross_ac.x - ac.x * ab_cross_ac.z) * ab_mid.y + (ac.x * ab_cross_ac.y - ac.y * ab_cross_ac.x) * ab_mid.z) / denominator;
	float t2 = ((ab.y * ab_cross_ac.z - ab.z * ab_cross_ac.y) * ac_mid.x + (ab.z * ab_cross_ac.x - ab.x * ab_cross_ac.z) * ac_mid.y + (ab.x * ab_cross_ac.y - ab.y * ab_cross_ac.x) * ac_mid.z) / denominator;

	return a + ab * t1 + ac * t2;
}

bool LoadModel(const std::string& filename)
{
	std::string err;

	if (!tinyobj::LoadObj(&g_obj_attrib, &g_obj_shapes, &g_obj_materials, &err, filename.c_str()))
	{
		std::cerr << "Failed to load OBJ: " << err << std::endl;
		return false;
	}

	uint32_t	triangleCount = 0;
	float3		minPos = float3::kFloatMax;
	float3		maxPos = float3::kNegativeFloatMax;

	for (const tinyobj::shape_t& shape : g_obj_shapes)
	{
		for (size_t i = 0; i < shape.mesh.indices.size(); i += 3)
		{
			tinyobj::index_t idx0 = shape.mesh.indices[i];
			tinyobj::index_t idx1 = shape.mesh.indices[i + 1];
			tinyobj::index_t idx2 = shape.mesh.indices[i + 2];

			float3 v0 = float3(g_obj_attrib.vertices[3 * idx0.vertex_index], g_obj_attrib.vertices[3 * idx0.vertex_index + 1], g_obj_attrib.vertices[3 * idx0.vertex_index + 2]);
			float3 v1 = float3(g_obj_attrib.vertices[3 * idx1.vertex_index], g_obj_attrib.vertices[3 * idx1.vertex_index + 1], g_obj_attrib.vertices[3 * idx1.vertex_index + 2]);
			float3 v2 = float3(g_obj_attrib.vertices[3 * idx2.vertex_index], g_obj_attrib.vertices[3 * idx2.vertex_index + 1], g_obj_attrib.vertices[3 * idx2.vertex_index + 2]);

			minPos = min(minPos, v0);
			minPos = min(minPos, v1);
			minPos = min(minPos, v2);

			maxPos = max(maxPos, v0);
			maxPos = max(maxPos, v1);
			maxPos = max(maxPos, v2);

			float area = TriangleArea(v0, v1, v2);
			
			g_avgTriangleArea += area;

			++triangleCount;
		}
	}

	if (!triangleCount)
		return false;

	g_avgTriangleArea /= (float)triangleCount;

	g_modelCenter = avg(minPos, maxPos);
	g_modelSize = maxPos - minPos;

	// Scale to fit printbed
	{
		float uniformScale = g_printbedSize / std::max(g_modelSize.x, g_modelSize.z);

		for (tinyobj::real_t& v : g_obj_attrib.vertices)
		{
			v *= uniformScale;
		}
		
		g_avgTriangleArea	= 0;
		triangleCount		= 0;
		minPos				= float3::kFloatMax;
		maxPos				= float3::kNegativeFloatMax;

		for (const tinyobj::shape_t& shape : g_obj_shapes)
		{
			for (size_t i = 0; i < shape.mesh.indices.size(); i += 3)
			{
				tinyobj::index_t idx0 = shape.mesh.indices[i];
				tinyobj::index_t idx1 = shape.mesh.indices[i + 1];
				tinyobj::index_t idx2 = shape.mesh.indices[i + 2];

				float3 v0 = float3(g_obj_attrib.vertices[3 * idx0.vertex_index], g_obj_attrib.vertices[3 * idx0.vertex_index + 1], g_obj_attrib.vertices[3 * idx0.vertex_index + 2]);
				float3 v1 = float3(g_obj_attrib.vertices[3 * idx1.vertex_index], g_obj_attrib.vertices[3 * idx1.vertex_index + 1], g_obj_attrib.vertices[3 * idx1.vertex_index + 2]);
				float3 v2 = float3(g_obj_attrib.vertices[3 * idx2.vertex_index], g_obj_attrib.vertices[3 * idx2.vertex_index + 1], g_obj_attrib.vertices[3 * idx2.vertex_index + 2]);

				minPos = min(minPos, v0);
				minPos = min(minPos, v1);
				minPos = min(minPos, v2);

				maxPos = max(maxPos, v0);
				maxPos = max(maxPos, v1);
				maxPos = max(maxPos, v2);

				float area = TriangleArea(v0, v1, v2);
			
				g_avgTriangleArea += area;

				++triangleCount;
			}
		}

		g_avgTriangleArea /= (float)triangleCount;

		g_modelCenter = avg(minPos, maxPos);
		g_modelSize = maxPos - minPos;
	}

	return true;
}
struct CellBoundary
{
	int i0;
	int i1;
	float3 p0;
	float3 p1;
	float3 x,y,z;
	float3 size;
	float3 center;
	std::vector<float3> points;
	std::vector<float3>	boundaryPoints;
	std::vector<float3> tubePoints;
	std::vector<uint32_t> tubeIndices;
};

std::vector<CellBoundary> g_cellBoundaries;
std::unordered_map<uint64_t, uint32_t> g_cellBoundaryFromHash;	// hash to g_cellBoundaries index

uint64_t EdgeHash(uint32_t i0, uint32_t i1)
{
	return ((uint64_t)std::min(i0,i1) << 32) | ((uint64_t)std::max(i0,i1));
}

void GeneratePointCloudOld()
{
	auto startTime_ns = std::chrono::high_resolution_clock::now();
	
	g_cellBoundaries.clear();
	g_cellBoundaryFromHash.clear();
	g_points.clear();
	g_controlPoints.clear();
	g_triangleIndices.clear();

	// Generate point cloud
	{
		for (const tinyobj::shape_t& shape : g_obj_shapes)
		{
			for (size_t i = 0; i < shape.mesh.indices.size(); i += 3)
			{
				tinyobj::index_t idx0 = shape.mesh.indices[i];
				tinyobj::index_t idx1 = shape.mesh.indices[i + 1];
				tinyobj::index_t idx2 = shape.mesh.indices[i + 2];

				float3 v0 = float3(g_obj_attrib.vertices[3 * idx0.vertex_index], g_obj_attrib.vertices[3 * idx0.vertex_index + 1], g_obj_attrib.vertices[3 * idx0.vertex_index + 2]);
				float3 v1 = float3(g_obj_attrib.vertices[3 * idx1.vertex_index], g_obj_attrib.vertices[3 * idx1.vertex_index + 1], g_obj_attrib.vertices[3 * idx1.vertex_index + 2]);
				float3 v2 = float3(g_obj_attrib.vertices[3 * idx2.vertex_index], g_obj_attrib.vertices[3 * idx2.vertex_index + 1], g_obj_attrib.vertices[3 * idx2.vertex_index + 2]);

				float area = TriangleArea(v0, v1, v2);
			
				int numSamples = static_cast<int>(g_pointsPerTriangle * area / g_avgTriangleArea);

				for (int j = 0; j < numSamples; ++j)
				{
					g_points.push_back(RandomPointInTriangle(v0, v1, v2, j, numSamples));
				}
			}
		}
	}

	// Pick out control points and leave the rest as loose particles
	{
		std::vector<float3> reducedPoints;

		for (const float3& p : g_points)
		{
			bool keep = true;

			for (const float3& q : g_controlPoints) 
			{
				if (distance(p, q) < g_cellRadius) 
				{
					keep = false;
					break;
				}
			}

			if (keep)
				g_controlPoints.push_back(p);
			else
				reducedPoints.push_back(p);
		}

		g_points = reducedPoints;
	}

	// Assign particles to a boundary between its two closest control points
	{
		//for (int k=0; k!=g_simulationIterationCount; ++k)
		for (int i=0; i!=g_points.size(); ++i)
		{
			float3 p = g_points[i];

			float3 p0, p1;
			uint64_t edgeHash;

			// Find two closest points
			{
				uint32_t bestIndex0 = UINT32_MAX;
				uint32_t bestIndex1 = UINT32_MAX;
				float bestDistance0 = FLT_MAX;
				float bestDistance1 = FLT_MAX;

				for (int j=0; j!=g_controlPoints.size(); ++j)
				{
					float d = distance(g_controlPoints[j], p);

					if (d < g_cellRadius * 0.001f)
						continue;

					if (d <= bestDistance0)
					{
						bestIndex1 = bestIndex0;
						bestDistance1 = bestDistance0;

						bestIndex0 = j;
						bestDistance0 = d;
					}
					else if (d <= bestDistance1)
					{
						bestIndex1 = j;
						bestDistance1 = d;
					}
				}

				assert(bestIndex0 != UINT32_MAX);
				assert(bestIndex1 != UINT32_MAX);
				assert(bestIndex0 != bestIndex1);

				p0 = g_controlPoints[bestIndex0];
				p1 = g_controlPoints[bestIndex1];
				edgeHash = EdgeHash(bestIndex0, bestIndex1);
			}

			if (fabs(dot(normalize(p0-p1), normalize(p0-p))) >= 0.999f)
			{
				p = avg(p0,p1);
			}
			else
			{
				float3 edgeVec = p1 - p0;
				float3 edgeDir = normalize(p1 - p0);
				float3 planeNormal = TriangleNormal(p0, p1, p);
				float3 tangentDir = normalize( cross( planeNormal, edgeDir ) );

				float v = dot( p - p0, tangentDir );

				p = avg(p0,p1) + tangentDir * v;
			}

			if (g_cellBoundaryFromHash.find(edgeHash) == g_cellBoundaryFromHash.end())
			{
				g_cellBoundaryFromHash[edgeHash] = (uint32_t)g_cellBoundaries.size();
				g_cellBoundaries.push_back(CellBoundary());
			}

			CellBoundary& cellBoundary = g_cellBoundaries[ g_cellBoundaryFromHash[edgeHash] ];
			cellBoundary.p0 = p0;
			cellBoundary.p1 = p1;
			cellBoundary.points.push_back(p);
		}
	}

	// Approximate cell boundaries using polynomials and build tube geometry
	{
		for (int b=0; b != (int)g_cellBoundaries.size(); ++b)
		{
			CellBoundary& cellBoundary = g_cellBoundaries[b];

			const float3& p0 = cellBoundary.p0;
			const float3& p1 = cellBoundary.p1;
			float3 midP = avg(p0,p1);

			float3 avgN = float3::k000;
			int count = 0;

			for (const float3& p : cellBoundary.points)
			{
				if (DistanceToEdge(p0, p1, p) > 0.001f)
				{
					avgN += TriangleNormal(p0, p1, p);
					++count;
				}
			}

			if (count == 0)
				continue;

			avgN = normalize(avgN);

			float3 edgeZ = normalize(p1 - p0);
			float3 edgeY = avgN;
			float3 edgeX = normalize(cross(edgeY, edgeZ));

			cellBoundary.x = edgeX;
			cellBoundary.y = edgeY;
			cellBoundary.z = edgeZ;

			float pointsAvgX = 0;
			float pointsAvgY = 0;
			float pointsMinX = FLT_MAX;
			float pointsMinY = FLT_MAX;
			float pointsMaxX = -FLT_MAX;
			float pointsMaxY = -FLT_MAX;

			std::vector<float2> points2D;

			for (const float3& p : cellBoundary.points)
			{
				float pX = dot(p - p0, edgeX);
				float pY = dot(p - p0, edgeY);
				float pZ = dot(p - p0, edgeZ);

				pointsAvgX += pX;
				pointsAvgY += pY;
				pointsMinX = std::min(pointsMinX, pX);
				pointsMinY = std::min(pointsMinY, pY);
				pointsMaxX = std::max(pointsMaxX, pX);
				pointsMaxY = std::max(pointsMaxY, pY);

				points2D.push_back(float2(pX, pY));
			}

			cellBoundary.size.x = pointsMaxX - pointsMinX;
			cellBoundary.size.y = pointsMaxY - pointsMinY;
			cellBoundary.size.z = distance(p0, p1);

			cellBoundary.center.x = (pointsMaxX + pointsMinX) * 0.5f;
			cellBoundary.center.y = (pointsMaxY + pointsMinY) * 0.5f;
			cellBoundary.center.z = 0.5f * distance(p0, p1);

			if (cellBoundary.size.x < cellBoundary.size.y)
			{
				std::swap(pointsMinX, pointsMinY);
				std::swap(pointsMaxX, pointsMaxY);
				std::swap(cellBoundary.size.x, cellBoundary.size.y);
				std::swap(cellBoundary.center.x, cellBoundary.center.y);
				std::swap(cellBoundary.x, cellBoundary.y);
				std::swap(edgeX, edgeY);

				for (float2& p : points2D)
				{
					std::swap(p.x, p.y);
				}
			}

			std::sort(points2D.begin(), points2D.end(), [](const float2& a, const float2& b)->int
			{
				return a.x < b.x ? -1 : 1;
			});

			float aspectRatio = std::max(cellBoundary.size.x, cellBoundary.size.y) / std::min(cellBoundary.size.x, cellBoundary.size.y);

			int polynomialDegree = aspectRatio < 1.3f ? 5 : 1;

			std::vector<float> polynomialFactors = polynomialFit(points2D, std::max(1, std::min((int)points2D.size()-2, polynomialDegree)));

			bool valid = true;
			for (float c : polynomialFactors)
			{
				if (!std::isfinite(c))
				{
					valid = false;
					break;
				}
			}

			if (!valid)
				continue;

			float stepX = (pointsMaxX - pointsMinX) / (float)g_tubeSegmentCount;

			for (int s=0; s != g_tubeSegmentCount; ++s)
			{
				float x = pointsMinX + s * stepX;
				float y = EvaluatePolynomial(polynomialFactors, x);
				float dy_dx = EvaluatePolynomialSlope(polynomialFactors, x);

				float3 polyZ = normalize( edgeY * dy_dx + edgeX );
				float3 polyX = edgeZ;
				float3 polyY = normalize(cross(polyZ, polyX));

				float3 tubeP = midP + edgeX * x + edgeY * y;
				cellBoundary.tubePoints.push_back(tubeP);

				float angleStep = (float)M_PI * 2.0f / (float)g_tubeAngleCount;

				for (int a=0; a != g_tubeAngleCount; ++a)
				{
					float dx = cosf(angleStep * a);
					float dy = sinf(angleStep * a);

					float3 p = tubeP + polyX * dx  * g_tubeRadius + polyY * dy * g_tubeRadius;

					//cellBoundary.tubePoints.push_back(p);

					if (s>0)
					{
						if (a==0)
						{
							cellBoundary.tubeIndices.push_back( (s - 1) * g_tubeAngleCount + g_tubeAngleCount - 1 );
							cellBoundary.tubeIndices.push_back( (s - 0) * g_tubeAngleCount + g_tubeAngleCount - 1 );
							cellBoundary.tubeIndices.push_back( (s - 1) * g_tubeAngleCount + a - 0 );
							cellBoundary.tubeIndices.push_back( (s - 1) * g_tubeAngleCount + a - 0 );
							cellBoundary.tubeIndices.push_back( (s - 0) * g_tubeAngleCount + g_tubeAngleCount - 1 );
							cellBoundary.tubeIndices.push_back( (s - 0) * g_tubeAngleCount + a - 0 );
						}
						else
						{
							cellBoundary.tubeIndices.push_back( (s - 1) * g_tubeAngleCount + a - 1 );
							cellBoundary.tubeIndices.push_back( (s - 0) * g_tubeAngleCount + a - 1 );
							cellBoundary.tubeIndices.push_back( (s - 1) * g_tubeAngleCount + a - 0 );
							cellBoundary.tubeIndices.push_back( (s - 1) * g_tubeAngleCount + a - 0 );
							cellBoundary.tubeIndices.push_back( (s - 0) * g_tubeAngleCount + a - 1 );
							cellBoundary.tubeIndices.push_back( (s - 0) * g_tubeAngleCount + a - 0 );
						}
					}
				}
			}
		
			if (b == 2025)
				printf("");
		}
	}

	auto endTime_ns = std::chrono::high_resolution_clock::now();
	std::chrono::nanoseconds duration_ns = endTime_ns - startTime_ns;
	g_processingTime = duration_ns.count() / 1e9f;
}

void GeneratePointCloud()
{
	auto startTime_ns = std::chrono::high_resolution_clock::now();
	
	g_cellBoundaries.clear();
	g_cellBoundaryFromHash.clear();
	g_points.clear();
	g_controlPoints.clear();
	g_triangleIndices.clear();

	// Generate point cloud
	{
		for (const tinyobj::shape_t& shape : g_obj_shapes)
		{
			for (size_t i = 0; i < shape.mesh.indices.size(); i += 3)
			{
				tinyobj::index_t idx0 = shape.mesh.indices[i];
				tinyobj::index_t idx1 = shape.mesh.indices[i + 1];
				tinyobj::index_t idx2 = shape.mesh.indices[i + 2];

				float3 v0 = float3(g_obj_attrib.vertices[3 * idx0.vertex_index], g_obj_attrib.vertices[3 * idx0.vertex_index + 1], g_obj_attrib.vertices[3 * idx0.vertex_index + 2]);
				float3 v1 = float3(g_obj_attrib.vertices[3 * idx1.vertex_index], g_obj_attrib.vertices[3 * idx1.vertex_index + 1], g_obj_attrib.vertices[3 * idx1.vertex_index + 2]);
				float3 v2 = float3(g_obj_attrib.vertices[3 * idx2.vertex_index], g_obj_attrib.vertices[3 * idx2.vertex_index + 1], g_obj_attrib.vertices[3 * idx2.vertex_index + 2]);

				float area = TriangleArea(v0, v1, v2);
			
				int numSamples = static_cast<int>(g_pointsPerTriangle * area / g_avgTriangleArea);

				for (int j = 0; j < numSamples; ++j)
				{
					g_points.push_back(RandomPointInTriangle(v0, v1, v2, j, numSamples));
				}
			}
		}
	}

	// Pick out control points
	{
		for (const float3& p : g_points)
		{
			bool keep = true;

			for (const float3& q : g_controlPoints) 
			{
				if (distance(p, q) < g_cellRadius) 
				{
					keep = false;
					break;
				}
			}

			if (keep)
				g_controlPoints.push_back(p);
		}
	}

	{
		for (int i=0; i!=g_controlPoints.size(); ++i)
		{
			float3 p0 = g_controlPoints[i];

			float3 avgN		= float3::k000;
			float  avgCount	= 0;

			std::set<uint32_t> cellBoundaries;
			std::vector<float4> clippingPlanes;

			for (int j=0; j!=g_controlPoints.size(); ++j)
			{
				if (i==j)
					continue;

				float3 p1 = g_controlPoints[j];

				if (distance(p0, p1) > g_cellRadius * 2.0f)
					continue;

				uint64_t edgeHash01 = EdgeHash(i,j);

				if (g_cellBoundaryFromHash.find(edgeHash01) == g_cellBoundaryFromHash.end())
				{
					g_cellBoundaryFromHash[edgeHash01] = (uint32_t)g_cellBoundaries.size();
					g_cellBoundaries.push_back(CellBoundary());
				}

				uint32_t boundaryIndex01 = g_cellBoundaryFromHash[edgeHash01];
				cellBoundaries.insert(boundaryIndex01);

				float4 clipPlane = MakePlaneEquation(normalize(p1-p0), avg(p0,p1));
				clippingPlanes.push_back(clipPlane);

				for (int k=0; k!=g_controlPoints.size(); ++k)
				{
					if (i==k || j==k)
						continue;

					float3 p2 = g_controlPoints[k];

					if (distance(p0, p2) > g_cellRadius * 2.0f)
						continue;

					if (fabs(dot( normalize(p1-p0), normalize(p2-p0) )) > 0.999f)
						continue;

					uint64_t edgeHash02 = EdgeHash(i,k);

					if (g_cellBoundaryFromHash.find(edgeHash02) == g_cellBoundaryFromHash.end())
					{
						g_cellBoundaryFromHash[edgeHash02] = (uint32_t)g_cellBoundaries.size();
						g_cellBoundaries.push_back(CellBoundary());
					}

					uint32_t boundaryIndex02 = g_cellBoundaryFromHash[edgeHash02];
					cellBoundaries.insert(boundaryIndex02);

					float3 c = CircumCenter(p0, p1, p2);

					g_cellBoundaries[boundaryIndex01].boundaryPoints.push_back(c);
					g_cellBoundaries[boundaryIndex02].boundaryPoints.push_back(c);

					float3 N = TriangleNormal(p0, p1, p2);

					avgN += N;
					avgCount += 1;
				}
			}

			for (uint32_t index : cellBoundaries)
			{
				CellBoundary& cellBoundary = g_cellBoundaries[index];

				// Remove any boundary points outside any of the clipping planes
				for (auto itr = cellBoundary.boundaryPoints.begin(); itr != cellBoundary.boundaryPoints.end(); )
				{
					float4 p = float4( (*itr), 1 );
				
					bool removed = false;

					for (const float4& clippingPlane : clippingPlanes)
					{
						if (dot(clippingPlane, p) > 0)
						{
							itr = cellBoundary.boundaryPoints.erase(itr);
							removed = true;
							break;
						}
					}

					if (!removed)
					{
						++itr;
					}
				}
			}

			if (avgCount>0)
			{
				avgN /= avgCount;
			}
		}
	}

	auto endTime_ns = std::chrono::high_resolution_clock::now();
	std::chrono::nanoseconds duration_ns = endTime_ns - startTime_ns;
	g_processingTime = duration_ns.count() / 1e9f;
}

void SavePointCloud(const std::string& filename, const std::vector<float3>& points)
{
	std::ofstream outFile(filename);
	if (!outFile)
	{
		std::cerr << "Error: Unable to write file " << filename << std::endl;
		return;
	}
	for (const auto& p : points)
	{
		outFile << p.x << " " << p.y << " " << p.z << "\n";
	}
	std::cout << "Saved " << points.size() << " points to " << filename << std::endl;
}

void SaveToOBJ(const std::string& filename, const std::vector<float3>& points, const std::vector<uint32_t>& triangleIndices)
{
	std::ofstream file(filename);

	if (!file.is_open())
	{
		std::cerr << "Failed to open file: " << filename << std::endl;
		return;
	}

	for (const float3& v : points)
	{
		file << "v " << v.x << " " << v.y << " " << v.z << "\n";
	}

	for (size_t i = 0; i < triangleIndices.size(); i += 3)
	{
		file << "f " << triangleIndices[i+0]+1 << " " << triangleIndices[i+1]+1 << " " << triangleIndices[i+2]+1 << "\n";
	}

	file.close();
}

void AppInit(const std::vector<std::string>& args)
{
	if (args.size() < 2)
	{
		Quit();
		return;
	}

	const std::string inputFile = args[1];

	SetWindowTitle(inputFile);


	if (!LoadModel(inputFile))
	{
		Quit();
	}

	g_cameraController.SetTarget(float4(g_modelCenter, 1));
	g_cameraController.SetViewDirection(float4::k0010);
	g_cameraController.SetViewDistance(length(g_modelSize) * 2.0f);

	GeneratePointCloud();
}

void AppResize(uint32_t windowWidth, uint32_t windowHeight)
{}

void AppRender(uint32_t windowWidth, uint32_t windowHeight)
{
	static auto startTime_ns = std::chrono::high_resolution_clock::now();
	auto currentTime_ns = std::chrono::high_resolution_clock::now();
	
	std::chrono::nanoseconds duration_ns = currentTime_ns - startTime_ns;
	float duration_ms = duration_ns.count() / 1e6f;

	static uint32_t triIndexCount = 0;

	if (duration_ms >= 1000)
	{
		startTime_ns = currentTime_ns;
	
		if (triIndexCount < g_triangleIndices.size())
			triIndexCount += 3;
	}

	float modelRadius = length(g_modelSize) * 0.5f;
	float zNear = std::max(0.01f, std::min(0.1f, g_cameraController.GetViewDistance() * 0.5f));
	float zFar = std::max(1000.0f, std::max(zNear + modelRadius, distance(g_cameraController.GetViewPosition(), g_modelCenter) + modelRadius * 1.5f ));
	float viewDistance = std::max(g_cameraController.GetViewDistance(), zNear);
	assert(std::isfinite(zNear));
	assert(std::isfinite(zFar));
	assert(zNear<zFar);
	assert(std::isfinite(viewDistance));
	g_cameraController.SetViewDistance(viewDistance);
	g_cameraController.SetProjectionMtx( matrix44::MakePerspective( 60, (float)windowWidth / (float)windowHeight, zNear, zFar ) );

	// printbed
	{
		float3 printbedVerts[] = 
		{
			g_modelCenter + float3( -g_printbedSize * 0.5f, -g_modelSize.y * 0.5f, -g_printbedSize * 0.5f ),
			g_modelCenter + float3( +g_printbedSize * 0.5f, -g_modelSize.y * 0.5f, -g_printbedSize * 0.5f ),
			g_modelCenter + float3( +g_printbedSize * 0.5f, -g_modelSize.y * 0.5f, +g_printbedSize * 0.5f ),
			g_modelCenter + float3( -g_printbedSize * 0.5f, -g_modelSize.y * 0.5f, +g_printbedSize * 0.5f ),
			g_modelCenter + float3( -g_printbedSize * 0.5f, -g_modelSize.y * 0.5f, -g_printbedSize * 0.5f ),
		};

		uint32_t printbedIndices[] =
		{
			0, 3, 1,
			1, 3, 2,
		};

		DrawMesh(printbedVerts, 4, printbedIndices, 6, matrix44::kIdentity, float4(0.5f, 0.5f, 0.5f, 1.0f));
		DrawLineStrip(printbedVerts, 5, matrix44::kIdentity, float4(1,0,0,1));
	}

	{
		int count = 0;

		for (CellBoundary& cellBoundary : g_cellBoundaries)
		{
			if (g_showJustOne && g_maxCellBoundaryIndex>=0 && count != g_maxCellBoundaryIndex)
			{
				++count;
				continue;
			}
			else if (g_maxCellBoundaryIndex>=0 && count > g_maxCellBoundaryIndex)
			{
				break;
			}

			if (g_showCellsPoints)
			{
				//DrawPoints(cellBoundary.points.data(), (uint32_t)cellBoundary.points.size(), matrix44::kIdentity, float4::k1111 * 0.8f);
				DrawPoints(cellBoundary.boundaryPoints.data(), (uint32_t)cellBoundary.boundaryPoints.size(), matrix44::kIdentity, float4::k1111 * 0.8f);
			}

			if (g_showCellsTubes)
				DrawLineStrip(cellBoundary.tubePoints.data(), (int)cellBoundary.tubePoints.size(), matrix44::kIdentity, float4(235, 137, 52, 255)/255.0f);

			//DrawPoints(cellBoundary.tubePoints.data(), (uint32_t)cellBoundary.tubePoints.size(), matrix44::kIdentity, float4(1,1,1,1));

			//DrawMesh(cellBoundary.tubePoints.data(), (uint32_t)cellBoundary.tubePoints.size(), cellBoundary.tubeIndices.data(), (uint32_t)cellBoundary.tubeIndices.size(), matrix44::kIdentity, float4(235, 137, 52, 255)/255.0f);

			//DrawPoints(&cellBoundary.p0, 1, matrix44::kIdentity, float4(1,0,0,1));
			//DrawPoints(&cellBoundary.p1, 1, matrix44::kIdentity, float4(1,0,0,1));
			
			if (g_showCellsAxes)
			{
				float3 m = avg(cellBoundary.p0, cellBoundary.p1);

				float3 axisPoints[] = 
				{
					cellBoundary.p0 + cellBoundary.x * (cellBoundary.center.x - cellBoundary.size.x * 0.5f) + cellBoundary.y * cellBoundary.center.y + cellBoundary.z * cellBoundary.center.z,
					cellBoundary.p0 + cellBoundary.x * (cellBoundary.center.x + cellBoundary.size.x * 0.5f) + cellBoundary.y * cellBoundary.center.y + cellBoundary.z * cellBoundary.center.z, 

					cellBoundary.p0 + cellBoundary.y * (cellBoundary.center.y - cellBoundary.size.y * 0.5f) + cellBoundary.x * cellBoundary.center.x + cellBoundary.z * cellBoundary.center.z,
					cellBoundary.p0 + cellBoundary.y * (cellBoundary.center.y + cellBoundary.size.y * 0.5f) + cellBoundary.x * cellBoundary.center.x + cellBoundary.z * cellBoundary.center.z, 

					cellBoundary.p0 + cellBoundary.z * (cellBoundary.center.z - cellBoundary.size.z * 0.5f) + cellBoundary.x * cellBoundary.center.x + cellBoundary.y * cellBoundary.center.y,
					cellBoundary.p0 + cellBoundary.z * (cellBoundary.center.z + cellBoundary.size.z * 0.5f) + cellBoundary.x * cellBoundary.center.x + cellBoundary.y * cellBoundary.center.y, 
				};

				DrawLines(axisPoints + 0, 2, matrix44::kIdentity, float4(1.0f, 0.2f, 0.2f, 1));
				DrawLines(axisPoints + 2, 2, matrix44::kIdentity, float4(0.2f, 1.0f, 0.2f, 1));
				DrawLines(axisPoints + 4, 2, matrix44::kIdentity, float4(0.1f, 0.2f, 1.0f, 1));
			}

			++count;
		}
	}

	if (g_showControlPoints)
		DrawPoints(g_controlPoints.data(), (uint32_t)g_controlPoints.size(), matrix44::kIdentity, float4(1,0,0,1));

	if (g_showPoints)
	{
		//DrawMesh(g_points.data(), (uint32_t)g_points.size(), g_triangleIndices.data(), g_maxTriangleIndex>=0 ? g_maxTriangleIndex*3 : (uint32_t)g_triangleIndices.size(), matrix44::kIdentity, float4::k1111, float4::k1001);
		DrawPoints(g_points.data(), (uint32_t)g_points.size(), matrix44::kIdentity, float4(1,1,1,1));
	}

	ImGui::SetNextWindowSize(ImVec2(550, 450));
	ImGui::Begin("Controls");
	{
		if (ImGui::Button("Load Model"))
		{
			OPENFILENAME ofn;
			char szFile[260] = { 0 };

			ZeroMemory(&ofn, sizeof(ofn));
			ofn.lStructSize = sizeof(ofn);
			//ofn.hwndOwner = hwnd;
			ofn.lpstrFile = szFile;
			ofn.nMaxFile = sizeof(szFile);
			ofn.lpstrFilter = "*.OBJ\0*.obj\0";
			ofn.nFilterIndex = 0;
			ofn.lpstrFileTitle = NULL;
			ofn.nMaxFileTitle = 0;
			ofn.lpstrInitialDir = NULL;
			ofn.Flags = OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;

			if (GetOpenFileName(&ofn) == TRUE) 
			{
				LoadModel(ofn.lpstrFile);
				GeneratePointCloud();
			}
		}

		ImGui::SliderInt("Cell Boundaries", &g_maxCellBoundaryIndex, -1, (uint32_t)g_cellBoundaries.size());
		ImGui::SameLine();
		ImGui::Checkbox("Just one", &g_showJustOne);
		
		ImGui::SliderInt("Points per triangle", &g_pointsPerTriangle, 1, 100);
		ImGui::Checkbox("Stratified Points", &g_stratifiedPoints);
		ImGui::Checkbox("Irrational Points", &g_irrationalPoints);
		ImGui::SliderFloat("Cell Radius", &g_cellRadius, 0, 20);
		ImGui::SliderInt("Tube Segments", &g_tubeSegmentCount, 1, 30);
		ImGui::SliderInt("Tube Angle Count", &g_tubeAngleCount, 1, 32);
		ImGui::SliderFloat("Tube Radius", &g_tubeRadius, 0, 5);

		ImGui::Checkbox("Show Control Points", &g_showControlPoints);
		ImGui::Checkbox("Show Points", &g_showPoints);
		ImGui::Checkbox("Show Cell Axes", &g_showCellsAxes);
		ImGui::Checkbox("Show Cell Points", &g_showCellsPoints);
		ImGui::Checkbox("Show Cell Tubes", &g_showCellsTubes);

		if (ImGui::Button("Generate Points"))
		{
			GeneratePointCloud();
		}

		if (ImGui::Button("Teleport") || ImGui::IsKeyPressed(ImGuiKey_F))
		{
			int index = std::min(g_maxCellBoundaryIndex, (int)g_cellBoundaries.size());

			CellBoundary& cellBoundary = g_cellBoundaries[index];

			g_cameraController.SetTarget( float4(cellBoundary.p0 + cellBoundary.z * cellBoundary.center.z + cellBoundary.x * cellBoundary.center.x + cellBoundary.y * cellBoundary.center.y, 1) );
		}

		ImGui::Text("Processing took: %3.2fs", g_processingTime);
		ImGui::Text("Cell Boundaries: %d", g_cellBoundaries.size());
		ImGui::Text("Hover Vertex: %u", GetHoverVertexIndex());
	}
	ImGui::End();

	if (ImGui::IsKeyPressed(ImGuiKey_RightArrow))
		g_maxCellBoundaryIndex = std::min(g_maxCellBoundaryIndex+1, (int)g_cellBoundaries.size());

	if (ImGui::IsKeyPressed(ImGuiKey_UpArrow))
		g_maxCellBoundaryIndex = std::min(g_maxCellBoundaryIndex+10, (int)g_cellBoundaries.size());

	if (ImGui::IsKeyPressed(ImGuiKey_LeftArrow))
		g_maxCellBoundaryIndex = std::max(g_maxCellBoundaryIndex-1, -1);

	if (ImGui::IsKeyPressed(ImGuiKey_DownArrow))
		g_maxCellBoundaryIndex = std::max(g_maxCellBoundaryIndex-10, -1);

	if (ImGui::IsKeyPressed(ImGuiKey_F5))
		GeneratePointCloud();
}

void AppExit()
{
	SavePointCloud("points.xyz", g_points);
	SaveToOBJ("mesh.obj", g_points, g_triangleIndices); 
}