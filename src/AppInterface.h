#pragma once
#include "ArcBallCameraController.h"

extern ArcBallCameraController g_cameraController;

// Functions the application needs to provide
void AppInit(const std::vector<std::string>& args);
void AppResize(uint32_t windowWidth, uint32_t windowHeight);
void AppRender(uint32_t windowWidth, uint32_t windowHeight);
void AppExit();

void DrawMesh(const std::vector<float3>& positions, const std::vector<uint32_t> indices, const matrix44& worldMtx, float4 fillColor, float4 wireColor = float4::k0000);
void DrawLines(const std::vector<float3>& positions, const matrix44& worldMtx, float4 color);

// Functions the application can call
void Quit();
void SetWindowTitle(const std::string& title);