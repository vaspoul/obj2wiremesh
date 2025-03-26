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
bool								g_stratifiedPoints = false;
std::vector<float3>					g_points;
std::vector<float3>					g_controlPoints;
std::vector<uint32_t>				g_triangleIndices;

float								g_minDistance;
float								g_minPointDistanceFactor = 1.5f;
int									g_tubeSegments = 10;

int									g_maxCellBoundaryIndex = -1;
bool								g_showControlPoints = false;

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

float3 RandomPointInTriangle(const float3& v0, const float3& v1, const float3& v2, uint32_t i, uint32_t iMax)
{
	static std::random_device rd;
	static std::mt19937 gen(rd());
	static std::uniform_real_distribution<float> dist(0.0f, 1.0f);

	float u = dist(gen);
	float v = dist(gen);

	if (g_stratifiedPoints)
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
	float3 p0;
	float3 p1;
	std::vector<float3> points;
	std::vector<float3> tubePoints;
};

std::unordered_map<uint64_t, CellBoundary> g_cellBoundaries;

uint64_t EdgeHash(uint32_t i0, uint32_t i1)
{
	return ((uint64_t)std::min(i0,i1) << 32) | ((uint64_t)std::max(i0,i1));
}

void GeneratePointCloud(std::vector<float3>& points, float& minDistance)
{
	points.clear();
	g_cellBoundaries.clear();

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
				points.push_back(RandomPointInTriangle(v0, v1, v2, j, numSamples));
			}
		}
	}

	g_minDistance = sqrtf(g_avgTriangleArea * 2.0f) * g_minPointDistanceFactor;

	std::vector<float3> controlPoints;
	std::vector<float3> reducedPoints;

	for (const float3& p : points)
	{
		bool keep = true;

		for (const float3& q : controlPoints) 
		{
			if (distance(p, q) < g_minDistance) 
			{
				keep = false;
				break;
			}
		}

		if (keep)
			controlPoints.push_back(p);
		else
			reducedPoints.push_back(p);
	}

	points = reducedPoints;
	g_controlPoints = controlPoints;

	//for (int k=0; k!=g_simulationIterationCount; ++k)
	for (int i=0; i!=points.size(); ++i)
	{
		float3& p = points[i];

		float3 p0, p1;
		uint64_t edgeHash;

		// Find two closest points
		{
			uint32_t bestIndex0 = UINT32_MAX;
			uint32_t bestIndex1 = UINT32_MAX;
			float bestDistance0 = FLT_MAX;
			float bestDistance1 = FLT_MAX;

			for (int j=0; j!=controlPoints.size(); ++j)
			{
				float d = distance(controlPoints[j], p);

				if (d < g_minDistance * 0.001f)
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

			p0 = controlPoints[bestIndex0];
			p1 = controlPoints[bestIndex1];
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

		//assert( fabs(distance(p0, p) - distance(p1,p)) <= 1e-6 );

		CellBoundary& cellBoundary = g_cellBoundaries[edgeHash];
		cellBoundary.p0 = p0;
		cellBoundary.p1 = p1;
		cellBoundary.points.push_back(p);
	}
	

	points.clear();

	for (auto itr = g_cellBoundaries.begin(); itr != g_cellBoundaries.end(); ++itr)
	{
		CellBoundary& cellBoundary = (*itr).second;

		if (cellBoundary.points.size() < 5)
			continue;

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

		std::vector<float> polynomialFactors = polynomialFit(points2D, 6);

		float stepX = (pointsMaxX - pointsMinX) / (float)g_tubeSegments;

		for (int s=0; s != g_tubeSegments; ++s)
		{
			float x = pointsMinX + s * stepX;
			float y = EvaluatePolynomial(polynomialFactors, x);

			float3 tubeP = midP + edgeX * x + edgeY * y;

			cellBoundary.tubePoints.push_back(tubeP);
		}
	}
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

	GeneratePointCloud(g_points, g_minDistance);
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

	//DrawMesh(g_points.data(), (uint32_t)g_points.size(), g_triangleIndices.data(), g_maxTriangleIndex>=0 ? g_maxTriangleIndex*3 : (uint32_t)g_triangleIndices.size(), matrix44::kIdentity, float4::k1111, float4::k1001);
	DrawPoints(g_points.data(), (uint32_t)g_points.size(), matrix44::kIdentity, float4(1,1,1,1));
	
	{
		int count = 0;

		for (auto itr = g_cellBoundaries.begin(); itr != g_cellBoundaries.end(); ++itr)
		{
			CellBoundary& cellBoundary = (*itr).second;

			if (g_maxCellBoundaryIndex>=0 && count != g_maxCellBoundaryIndex)
			{
				++count;
				continue;
			}

			std::vector<float3> points = cellBoundary.tubePoints;
			points.insert(points.begin(), cellBoundary.p0);
			points.push_back(cellBoundary.p1);

			//DrawPoints(cellBoundary.points.data(), (uint32_t)cellBoundary.points.size(), matrix44::kIdentity, float4(0.5f, 0.5f, 0.5f, 1));
			//DrawLineStrip(points.data(), (int)points.size(), matrix44::kIdentity, count%2 ? float4(0.8f, 0.8f, 1.0f, 1.0f) : float4(0.8f, 1.0f, 0.8f, 1.0f) );
			//DrawLineStrip(points.data(), (int)points.size(), matrix44::kIdentity, float4::k1111);
			DrawPoints(cellBoundary.tubePoints.data(), (uint32_t)cellBoundary.tubePoints.size(), matrix44::kIdentity, float4(1,1,1,1));
			//DrawPoints(&cellBoundary.p0, 1, matrix44::kIdentity, float4(1,0,0,1));
			//DrawPoints(&cellBoundary.p1, 1, matrix44::kIdentity, float4(1,0,0,1));

			++count;
		}
	}

	if (g_showControlPoints)
		DrawPoints(g_controlPoints.data(), (uint32_t)g_controlPoints.size(), matrix44::kIdentity, float4(1,0,0,1));

	ImGui::SetNextWindowSize(ImVec2(550, 400));
	ImGui::Begin("Controls");
	{
		ImGui::SliderInt("Cell Boundaries", &g_maxCellBoundaryIndex, -1, std::min(2000u, (uint32_t)g_cellBoundaries.size()));

		ImGui::SameLine();
		if (ImGui::Button("-"))
		{
			g_maxCellBoundaryIndex = std::max(g_maxCellBoundaryIndex-1, -1);
		}

		ImGui::SameLine();
		if (ImGui::Button("+"))
		{
			g_maxCellBoundaryIndex = std::min(g_maxCellBoundaryIndex+1, (int)g_cellBoundaries.size());
		}

		
		ImGui::SliderInt("Points per triangle", &g_pointsPerTriangle, 1, 100);
		ImGui::Checkbox("Stratified Points", &g_stratifiedPoints);
		ImGui::SliderFloat("Min Point Distance Factor", &g_minPointDistanceFactor, 0, 5);
		ImGui::SliderInt("Tube Segments", &g_tubeSegments, 1, 30);

		ImGui::Checkbox("Show Control Points", &g_showControlPoints);

		if (ImGui::Button("Generate Points"))
		{
			GeneratePointCloud(g_points, g_minDistance);
		}

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
		GeneratePointCloud(g_points, g_minDistance);
}

void AppExit()
{
	SavePointCloud("points.xyz", g_points);
	SaveToOBJ("mesh.obj", g_points, g_triangleIndices); 
}