#include "global.h"
#include "AppInterface.h"

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

std::vector<float3>				g_points;
std::vector<uint32_t>			g_triangleIndices;
float3							g_modelCenter;
float3							g_modelSize;

float TriangleArea(const float3& v0, const float3& v1, const float3& v2)
{
	float3 crossProd = cross(v1-v0, v2-v0);
	return 0.5f * std::sqrt(crossProd.x * crossProd.x + crossProd.y * crossProd.y + crossProd.z * crossProd.z);
}

float3 RandomPointInTriangle(const float3& v0, const float3& v1, const float3& v2)
{
	static std::random_device rd;
	static std::mt19937 gen(rd());
	static std::uniform_real_distribution<float> dist(0.0f, 1.0f);

	float u = dist(gen);
	float v = dist(gen);
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
				pointCloud.push_back(RandomPointInTriangle(v0, v1, v2));
			}
		}
	}

	minDistance = sqrtf(avgTriangleArea) * 1.0f;

	std::vector<float3> reducedPointCloud;

	for (const float3& p : pointCloud)
	{
		bool keep = true;

		for (const float3& q : reducedPointCloud) 
		{
			if (distance(p, q) < minDistance) 
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

void RollingBallTriangulation(const std::vector<float3>& points, std::vector<uint32_t>& triangleIndices, float ballRadius)
{
	const uint32_t numPoints = (uint32_t)points.size();

	// Find initial triangle
	{
		uint32_t index0 = 0;
		uint32_t index1 = 1;
		uint32_t index2 = 2;
		float bestDistance1 = 1e6;
		float bestDistance2 = 1e6;

		for (uint32_t i = 0; i != numPoints; ++i) 
		{
			if (i == index0)
				continue;

			float currentDist = distance(points[index0], points[i]);
			
			if (currentDist < bestDistance1)
			{
				index2 = index1;
				bestDistance2 = bestDistance1;

				index1 = i;
				bestDistance1 = currentDist;
			}
			else if (currentDist < bestDistance2)
			{
				index2 = i;
				bestDistance2 = currentDist;
			}
		}

		triangleIndices.push_back(index0);
		triangleIndices.push_back(index1);
		triangleIndices.push_back(index2);
	}

	std::vector<uint32_t> edgeIndices;
	std::vector<uint32_t> edgeOppositeIndex;
	std::unordered_set<uint64_t> edgeHashes;

	auto EdgeHash = [](uint32_t i0, uint32_t i1) -> uint64_t
	{
		return ((uint64_t)std::min(i0,i1) << 32) | ((uint64_t)std::max(i0,i1));
	};

	edgeIndices.push_back(triangleIndices[0]);
	edgeIndices.push_back(triangleIndices[1]);
	edgeOppositeIndex.push_back(triangleIndices[2]);
	edgeHashes.insert(EdgeHash(triangleIndices[0], triangleIndices[1]));

	edgeIndices.push_back(triangleIndices[1]);
	edgeIndices.push_back(triangleIndices[2]);
	edgeOppositeIndex.push_back(triangleIndices[0]);
	edgeHashes.insert(EdgeHash(triangleIndices[1], triangleIndices[2]));

	edgeIndices.push_back(triangleIndices[2]);
	edgeIndices.push_back(triangleIndices[0]);
	edgeOppositeIndex.push_back(triangleIndices[1]);
	edgeHashes.insert(EdgeHash(triangleIndices[2], triangleIndices[0]));

	for (uint32_t i=0; i!=edgeIndices.size(); i += 2)
	{
		float bestDistance = 1e6;
		uint32_t bestIndex = UINT32_MAX;

		float3 midPoint = avg( points[edgeIndices[i+0]], points[edgeIndices[i+1]] );

		for (uint32_t j=0; j!=points.size(); ++j)
		{
			if ( j == edgeIndices[i+0] || j == edgeIndices[i+1] || j == edgeOppositeIndex[i/2] )
				continue;

			float currentDist = distance(midPoint, points[j]);

			if (currentDist < bestDistance)
			{
				bestDistance = currentDist;
				bestIndex = j;
			}
		}

		if (bestIndex != UINT32_MAX)
		{
			for (uint32_t k=i+2; k<edgeIndices.size(); k += 2)
			{
				if ( bestIndex == edgeIndices[k+0] || bestIndex == edgeIndices[k+1] || bestIndex == edgeOppositeIndex[k/2] )
					continue;

				float3 midPointK = avg( points[edgeIndices[k+0]], points[edgeIndices[k+1]] );
			
				float currentDist = distance(midPointK, points[bestIndex]);

				if (currentDist < bestDistance)
				{
					bestIndex = UINT32_MAX;
					break;
				}
			}
		}

		if (bestIndex != UINT32_MAX)
		{
			float3 triNormal0 = TriangleNormal( points[edgeIndices[i+0]], points[edgeIndices[i+1]], points[edgeOppositeIndex[i/2]] );
			float3 triNormal1 = TriangleNormal( points[edgeIndices[i+1]], points[edgeIndices[i+0]], points[bestIndex] );

			if (dot(triNormal0, triNormal1)>0)
			{
				triangleIndices.push_back(edgeIndices[i+1]);
				triangleIndices.push_back(edgeIndices[i+0]);
				triangleIndices.push_back(bestIndex);
			}
			else
			{
				triangleIndices.push_back(edgeIndices[i+0]);
				triangleIndices.push_back(edgeIndices[i+1]);
				triangleIndices.push_back(bestIndex);
			}

			if (edgeHashes.count(EdgeHash(edgeIndices[i+0], bestIndex)) == 0)
			{
				edgeIndices.push_back(edgeIndices[i+0]);
				edgeIndices.push_back(bestIndex);
				edgeOppositeIndex.push_back(edgeIndices[i+1]);
				edgeHashes.insert(EdgeHash(edgeIndices[i+0], bestIndex));
			}

			if (edgeHashes.count(EdgeHash(edgeIndices[i+1], bestIndex)) == 0)
			{
				edgeIndices.push_back(edgeIndices[i+1]);
				edgeIndices.push_back(bestIndex);
				edgeOppositeIndex.push_back(edgeIndices[i+0]);
				edgeHashes.insert(EdgeHash(edgeIndices[i+1], bestIndex));
			}
		}
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
	const std::string inputFile = args[1];

	SetWindowTitle(inputFile);

	const int pointsPerTriangle = 10;
	float minDistance;

	if (!GeneratePointCloud(inputFile, pointsPerTriangle, g_points, minDistance))
	{
		Quit();
	}

	g_cameraController.SetTarget(float4(g_modelCenter, 1));
	g_cameraController.SetViewDirection(float4::k0010);
	g_cameraController.SetViewDistance(length(g_modelSize) * 2.0f);

	RollingBallTriangulation(g_points, g_triangleIndices, minDistance * 1.5f);
}

void AppResize(uint32_t windowWidth, uint32_t windowHeight)
{}

void AppRender(uint32_t windowWidth, uint32_t windowHeight)
{
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

	DrawMesh(g_points, g_triangleIndices, matrix44::kIdentity, float4::k1111, float4::k1001);
}

void AppExit()
{
	SavePointCloud("points.xyz", g_points);
	SaveToOBJ("mesh.obj", g_points, g_triangleIndices); 
}