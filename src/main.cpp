#include "global.h"
#include "AppInterface.h"

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

std::vector<float3>				g_points;
std::vector<uint32_t>			g_triangleIndices;
float3							g_modelCenter;
float3							g_modelSize;
float							g_minDistance;
int								g_maxTriangleIndex = -1;

float TriangleArea(const float3& v0, const float3& v1, const float3& v2)
{
	float3 crossProd = cross(v1-v0, v2-v0);
	return 0.5f * std::sqrt(crossProd.x * crossProd.x + crossProd.y * crossProd.y + crossProd.z * crossProd.z);
}

float3 RandomPointInTriangle(const float3& v0, const float3& v1, const float3& v2, uint32_t i, uint32_t iMax)
{
	static std::random_device rd;
	static std::mt19937 gen(rd());
	static std::uniform_real_distribution<float> dist(0.0f, 1.0f);

	float u = dist(gen);
	float v = dist(gen);

	uint32_t iMaxX = std::max(1u, (uint32_t)sqrtf((float)iMax));
	
	u = float(i % iMaxX) / float(iMaxX);
	v = float(i / iMaxX) / float(iMaxX);

	u = saturate(u);
	v = saturate(v);

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

bool GeneratePointCloud(const std::string& filename, int pointsPerTriangle, std::vector<float3>& pointCloud, float& minDistance)
{
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string err;

	if (!tinyobj::LoadObj(&attrib,& shapes,& materials,& err, filename.c_str()))
	{
		std::cerr << "Failed to load OBJ: " << err << std::endl;
		return false;
	}

	float		avgTriangleArea = 0;
	uint32_t	triangleCount = 0;
	float3		minPos = float3::kFloatMax;
	float3		maxPos = float3::kNegativeFloatMax;

	for (const tinyobj::shape_t& shape : shapes)
	{
		for (size_t i = 0; i < shape.mesh.indices.size(); i += 3)
		{
			tinyobj::index_t idx0 = shape.mesh.indices[i];
			tinyobj::index_t idx1 = shape.mesh.indices[i + 1];
			tinyobj::index_t idx2 = shape.mesh.indices[i + 2];

			float3 v0 = float3(attrib.vertices[3 * idx0.vertex_index], attrib.vertices[3 * idx0.vertex_index + 1], attrib.vertices[3 * idx0.vertex_index + 2]);
			float3 v1 = float3(attrib.vertices[3 * idx1.vertex_index], attrib.vertices[3 * idx1.vertex_index + 1], attrib.vertices[3 * idx1.vertex_index + 2]);
			float3 v2 = float3(attrib.vertices[3 * idx2.vertex_index], attrib.vertices[3 * idx2.vertex_index + 1], attrib.vertices[3 * idx2.vertex_index + 2]);

			minPos = min(minPos, v0);
			minPos = min(minPos, v1);
			minPos = min(minPos, v2);

			maxPos = max(maxPos, v0);
			maxPos = max(maxPos, v1);
			maxPos = max(maxPos, v2);

			float area = TriangleArea(v0, v1, v2);
			
			avgTriangleArea += area;

			++triangleCount;
		}
	}

	if (!triangleCount)
		return false;

	g_modelCenter = avg(minPos, maxPos);
	g_modelSize = maxPos - minPos;

	avgTriangleArea /= (float)triangleCount;

	for (const tinyobj::shape_t& shape : shapes)
	{
		for (size_t i = 0; i < shape.mesh.indices.size(); i += 3)
		{
			tinyobj::index_t idx0 = shape.mesh.indices[i];
			tinyobj::index_t idx1 = shape.mesh.indices[i + 1];
			tinyobj::index_t idx2 = shape.mesh.indices[i + 2];

			float3 v0 = float3(attrib.vertices[3 * idx0.vertex_index], attrib.vertices[3 * idx0.vertex_index + 1], attrib.vertices[3 * idx0.vertex_index + 2]);
			float3 v1 = float3(attrib.vertices[3 * idx1.vertex_index], attrib.vertices[3 * idx1.vertex_index + 1], attrib.vertices[3 * idx1.vertex_index + 2]);
			float3 v2 = float3(attrib.vertices[3 * idx2.vertex_index], attrib.vertices[3 * idx2.vertex_index + 1], attrib.vertices[3 * idx2.vertex_index + 2]);

			float area = TriangleArea(v0, v1, v2);
			
			int numSamples = static_cast<int>(pointsPerTriangle * area / avgTriangleArea);

			for (int j = 0; j < numSamples; ++j)
			{
				pointCloud.push_back(RandomPointInTriangle(v0, v1, v2, j, numSamples));
			}
		}
	}

	g_minDistance = sqrtf(avgTriangleArea * 2.0f) * 1.0f;

	std::vector<float3> reducedPointCloud;

	for (const float3& p : pointCloud)
	{
		bool keep = true;

		for (const float3& q : reducedPointCloud) 
		{
			if (distance(p, q) < g_minDistance) 
			{
				keep = false;
				break;
			}
		}

		if (keep) 
			reducedPointCloud.push_back(p);
	}

	pointCloud = reducedPointCloud;

	return true;
}

std::vector<uint32_t>			g_edgeIndices;
std::vector<uint32_t>			g_edgeOppositeIndex;
std::unordered_set<uint64_t>	g_edgeHashes;
std::unordered_set<uint64_t>	g_fullEdges;

uint64_t EdgeHash(uint32_t i0, uint32_t i1)
{
	return ((uint64_t)std::min(i0,i1) << 32) | ((uint64_t)std::max(i0,i1));
}

void RollingBallTriangulationStart(const std::vector<float3>& points, std::vector<uint32_t>& triangleIndices, float ballRadius)
{
	const uint32_t numPoints = (uint32_t)points.size();

	// Find initial triangle
	{
		uint32_t index0 = 0;
		uint32_t index1 = 1;
		uint32_t index2 = 2;
		float bestDistance = 1e6;

		for (uint32_t i = 0; i != numPoints; ++i) 
		{
			if (i == index0)
				continue;

			float currentDist = distance(points[index0], points[i]);
			
			if (currentDist < bestDistance)
			{
				index1 = i;
				bestDistance = currentDist;
			}
		}

		float3 midPoint = avg(points[index0], points[index1]);
		float3 edgeDir01 = normalize(points[index1] - points[index0]);
		bestDistance = 1e6;

		for (uint32_t i = 0; i != numPoints; ++i) 
		{
			if (i == index0 || i == index1)
				continue;

			float currentDist = distance(midPoint, points[i]);
			float3 currentDir = normalize(points[i] - points[index0]);
			
			if ( currentDist < bestDistance && fabs(dot(edgeDir01, currentDir)) < 0.999f )
			{
				index2 = i;
				bestDistance = currentDist;
			}
		}

		std::swap(index0, index1);
		triangleIndices.push_back(index0);
		triangleIndices.push_back(index1);
		triangleIndices.push_back(index2);

		assert(TriangleArea(points[index0], points[index1], points[index2])>0);
	}

	g_edgeIndices.push_back(triangleIndices[0]);
	g_edgeIndices.push_back(triangleIndices[1]);
	g_edgeOppositeIndex.push_back(triangleIndices[2]);
	g_edgeHashes.insert(EdgeHash(triangleIndices[0], triangleIndices[1]));

	g_edgeIndices.push_back(triangleIndices[1]);
	g_edgeIndices.push_back(triangleIndices[2]);
	g_edgeOppositeIndex.push_back(triangleIndices[0]);
	g_edgeHashes.insert(EdgeHash(triangleIndices[1], triangleIndices[2]));

	g_edgeIndices.push_back(triangleIndices[2]);
	g_edgeIndices.push_back(triangleIndices[0]);
	g_edgeOppositeIndex.push_back(triangleIndices[1]);
	g_edgeHashes.insert(EdgeHash(triangleIndices[2], triangleIndices[0]));
}

void RollingBallTriangulation(const std::vector<float3>& points, std::vector<uint32_t>& triangleIndices, float ballRadius)
{
	triangleIndices.clear();
	g_edgeIndices.clear();
	g_edgeOppositeIndex.clear();
	g_edgeHashes.clear();
	g_fullEdges.clear();

	RollingBallTriangulationStart(g_points, triangleIndices, ballRadius);

	auto Range = [](float x, float xMin, float xMax) -> float
	{
		return saturate( (x - xMin) / (xMax - xMin) );
	};

	for (uint32_t i=0; i!=g_edgeIndices.size(); i += 2)
	{
		const uint32_t edgeIndex0 = g_edgeIndices[i+0];
		const uint32_t edgeIndex1 = g_edgeIndices[i+1];
		const uint32_t oppositeIndex = g_edgeOppositeIndex[i/2];

		if (g_fullEdges.count(EdgeHash(edgeIndex0, edgeIndex1)))
			continue;

		if (triangleIndices.size()/3 == 13)
			printf("");

		uint32_t bestIndex = UINT32_MAX;
		float bestScore = FLT_MIN;

		const float3 edgeDir = normalize(points[edgeIndex1] - points[edgeIndex0]);
		const float3 midPoint = avg(points[edgeIndex0], points[edgeIndex1]);
		const float3 triNormal0 = TriangleNormal( points[edgeIndex0], points[edgeIndex1], points[oppositeIndex] );
		//const float3 triTangent = normalize( points[oppositeIndex] - midPoint );
		const float3 triTangent = normalize( cross( triNormal0, edgeDir ) );

		for (uint32_t j=0; j!=points.size(); ++j)
		{
			if ( j == edgeIndex0 || j == edgeIndex1 || j == oppositeIndex )
				continue;

			if (g_fullEdges.count(EdgeHash(edgeIndex0, j)) || g_fullEdges.count(EdgeHash(edgeIndex1, j)))
				continue;

			float3 candidateDir = normalize(points[j] - points[edgeIndex0]);

			// Exclude colinear points
			if (fabs(dot(edgeDir, candidateDir)) >= 0.999f)
				continue;

			//float candidateDist = DistanceToEdge(points[edgeIndex0], points[edgeIndex1], points[j]);
			float candidateDist = distance(midPoint, points[j]);
			//float candidateDist = std::min( distance(points[edgeIndex0], points[j]), distance(points[edgeIndex1], points[j]) );
			float3 candidateNormal = TriangleNormal( points[edgeIndex1], points[edgeIndex0], points[j] );
			//float3 candidateTangent = normalize( points[j] - midPoint );
			float3 candidateTangent = normalize( cross( candidateNormal, -edgeDir ) );

			float candidatePlaneDist = dot( (points[j] - points[edgeIndex0]), triTangent );

			float candidateNormalDot = dot(candidateNormal, triNormal0);
			float candidateTangentDot = dot(candidateTangent, triTangent);

			float candidateAngle1 = fabs( dot( normalize(points[j] - points[edgeIndex0]), normalize(points[edgeIndex1] - points[edgeIndex0]) ) );
			float candidateAngle2 = fabs( dot( normalize(points[j] - points[edgeIndex1]), normalize(points[edgeIndex1] - points[edgeIndex0]) ) );
			float candidateMinAngle = std::min(candidateAngle1, candidateAngle2);

			//float hingeAngle = degrees( atan2f(candidateNormalDot, -candidateTangentDot) );
			float dy = dot(points[j] - midPoint, triNormal0);
			float dx = dot(points[j] - midPoint, triTangent);

			dy /= fabs(dy);
			dx /= fabs(dy);

			float hingeAngle = degrees( atan2f(dy, dx) );

			if (hingeAngle<0) hingeAngle += 360;

			float score = 0;

			//score += 6.0f * Range(candidateDist / ballRadius, 1, 0);
			//score += 4.0f * Range( fabs(candidateNormalDot), 0.2f, 1.0f );
			//score += 5.0f * Range( candidateTangentDot, 0.2f, -1.0f );
			score += 5.0f * Range( hingeAngle, 360, 0 );
			//score += 4.0f * Range( candidateMinAngle, 1.0f, 0.5f );

			score *= candidateDist <= ballRadius ? 1 : 0;
			score *= candidatePlaneDist <= 0 ? 1 : 0;
			//score *= fabs(candidateNormalDot) > 0.2f ? 1 : 0;

			if (score >= bestScore)
			{
				bestScore = score;
				bestIndex = j;
			}
		}

		uint32_t futureEdgeIndex = UINT32_MAX;

		//if (bestIndex != UINT32_MAX)
		//{
		//	for (uint32_t k=edgeIndex+2; k<g_edgeIndices.size(); k += 2)
		//	{
		//		if ( bestIndex == g_edgeIndices[k+0] || bestIndex == g_edgeIndices[k+1] || bestIndex == g_edgeOppositeIndex[k/2] )
		//			continue;

		//		float currentDist = DistanceToEdge(points[g_edgeIndices[k+0]], points[g_edgeIndices[k+1]], points[bestIndex]);

		//		if (currentDist < bestDistance)
		//		{
		//			bestIndex = UINT32_MAX;
		//			futureEdgeIndex = k;
		//			break;
		//		}
		//	}
		//}

		if (bestIndex != UINT32_MAX)
		{
			if (triangleIndices.size() == 17*3)
				printf("");

			float3 triNormal1 = TriangleNormal( points[edgeIndex1], points[edgeIndex0], points[bestIndex] );

			bool sameWinding = dot(triNormal0, triNormal1)>0;

			uint32_t newTriIndex0 = sameWinding ? edgeIndex1 : edgeIndex0;
			uint32_t newTriIndex1 = sameWinding ? edgeIndex0 : edgeIndex1;

			triangleIndices.push_back(newTriIndex0);
			triangleIndices.push_back(newTriIndex1);
			triangleIndices.push_back(bestIndex);

			assert(TriangleArea(points[newTriIndex0], points[newTriIndex1], points[bestIndex])>0);
			assert(g_edgeHashes.count(EdgeHash(newTriIndex0, newTriIndex1)) == 1);

			g_fullEdges.insert(EdgeHash(newTriIndex0, newTriIndex1));

			if (g_edgeHashes.count(EdgeHash(newTriIndex1, bestIndex)) == 0)
			{
				//g_edgeIndices.insert(g_edgeIndices.begin() + i + 2, bestIndex);
				//g_edgeIndices.insert(g_edgeIndices.begin() + i + 2, newTriIndex1);
				//g_edgeOppositeIndex.insert(g_edgeOppositeIndex.begin() + i/2 + 1, newTriIndex0);
				g_edgeIndices.push_back(newTriIndex1);
				g_edgeIndices.push_back(bestIndex);
				g_edgeOppositeIndex.push_back(newTriIndex0);
				g_edgeHashes.insert(EdgeHash(newTriIndex1, bestIndex));
			}
			else
			{
				g_fullEdges.insert(EdgeHash(newTriIndex1, bestIndex));
			}

			if (g_edgeHashes.count(EdgeHash(bestIndex, newTriIndex0)) == 0)
			{
				//g_edgeIndices.insert(g_edgeIndices.begin() + i + 2, newTriIndex0);
				//g_edgeIndices.insert(g_edgeIndices.begin() + i + 2, bestIndex);
				//g_edgeOppositeIndex.insert(g_edgeOppositeIndex.begin() + i/2 + 1, newTriIndex1);
				g_edgeIndices.push_back(bestIndex);
				g_edgeIndices.push_back(newTriIndex0);
				g_edgeOppositeIndex.push_back(newTriIndex1);
				g_edgeHashes.insert(EdgeHash(bestIndex, newTriIndex0));
			}
			else
			{
				g_fullEdges.insert(EdgeHash(bestIndex, newTriIndex0));
			}
		}

		if (triangleIndices.size() >= 3*10000)
			break;
	}
}

void SavePointCloud(const std::string& filename, const std::vector<float3>& pointCloud)
{
	std::ofstream outFile(filename);
	if (!outFile)
	{
		std::cerr << "Error: Unable to write file " << filename << std::endl;
		return;
	}
	for (const auto& p : pointCloud)
	{
		outFile << p.x << " " << p.y << " " << p.z << "\n";
	}
	std::cout << "Saved " << pointCloud.size() << " points to " << filename << std::endl;
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

	const int pointsPerTriangle = 20;
	
	if (!GeneratePointCloud(inputFile, pointsPerTriangle, g_points, g_minDistance))
	{
		Quit();
	}

	g_cameraController.SetTarget(float4(g_modelCenter, 1));
	g_cameraController.SetViewDirection(float4::k0010);
	g_cameraController.SetViewDistance(length(g_modelSize) * 2.0f);

	RollingBallTriangulation(g_points, g_triangleIndices, g_minDistance * 1.5f);
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

	DrawMesh(g_points.data(), (uint32_t)g_points.size(), g_triangleIndices.data(), g_maxTriangleIndex>=0 ? g_maxTriangleIndex*3 : (uint32_t)g_triangleIndices.size(), matrix44::kIdentity, float4::k1111, float4::k1001);
	DrawPoints(g_points.data(), (uint32_t)g_points.size(), matrix44::kIdentity, float4(1,1,1,1));

	ImGui::SetNextWindowSize(ImVec2(450, 200));
	ImGui::Begin("Controls");
	{
		ImGui::SliderInt("Triangles [1000]", &g_maxTriangleIndex, -1, std::min(1000u, (uint32_t)g_triangleIndices.size()/3));

		ImGui::SliderInt("Triangles", &g_maxTriangleIndex, -1, (uint32_t)g_triangleIndices.size()/3);
		
		ImGui::SameLine();
		if (ImGui::Button("-"))
		{
			g_maxTriangleIndex = std::max(g_maxTriangleIndex-1, -1);
		}

		ImGui::SameLine();
		if (ImGui::Button("+"))
		{
			g_maxTriangleIndex = std::min(g_maxTriangleIndex+1, (int)g_triangleIndices.size()/3);
		}

		if (ImGui::Button("Retriangulate"))
		{
			RollingBallTriangulation(g_points, g_triangleIndices, g_minDistance * 1.5f);
		}

		ImGui::Text("Hover Vertex: %u", GetHoverVertexIndex());
	}
	ImGui::End();

	if (ImGui::IsKeyPressed(ImGuiKey_RightArrow))
		g_maxTriangleIndex = std::min(g_maxTriangleIndex+1, (int)g_triangleIndices.size()/3);

	if (ImGui::IsKeyPressed(ImGuiKey_UpArrow))
		g_maxTriangleIndex = std::min(g_maxTriangleIndex+10, (int)g_triangleIndices.size()/3);

	if (ImGui::IsKeyPressed(ImGuiKey_LeftArrow))
		g_maxTriangleIndex = std::max(g_maxTriangleIndex-1, -1);

	if (ImGui::IsKeyPressed(ImGuiKey_DownArrow))
		g_maxTriangleIndex = std::max(g_maxTriangleIndex-10, -1);
}

void AppExit()
{
	SavePointCloud("points.xyz", g_points);
	SaveToOBJ("mesh.obj", g_points, g_triangleIndices); 
}