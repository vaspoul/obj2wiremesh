#pragma once
#include "ArcBallCameraController.h"
#include "imgui/imgui.h"

extern ArcBallCameraController g_cameraController;

// Functions the application needs to provide
void AppInit(const std::vector<std::string>& args);
void AppResize(uint32_t windowWidth, uint32_t windowHeight);
void AppRender(uint32_t windowWidth, uint32_t windowHeight);
void AppExit();

void DrawMesh(const float3* positions, uint32_t pointCount, uint32_t* indices, uint32_t indexCount, const matrix44& worldMtx, float4 fillColor, float4 wireColor = float4::k0000);
void DrawLines(const float3* positions, uint32_t pointCount, const matrix44& worldMtx, float4 color);
void DrawLineStrip(const float3* positions, uint32_t pointCount, const matrix44& worldMtx, float4 color);
void DrawPoints(const float3* positions, uint32_t pointCount, const matrix44& worldMtx, float4 color);

// Functions the application can call
void Quit();
void SetWindowTitle(const std::string& title);
uint32_t GetHoverVertexIndex();