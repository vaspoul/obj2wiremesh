#include "global.h"
#include "ArcBallCameraController.h"

ArcBallCameraController::ArcBallCameraController()
{
	Reset();
	m_previousTime_ns = std::chrono::high_resolution_clock::now();
}

void ArcBallCameraController::Reset()
{
	m_lookTarget					= float4(0,0,0,1);
	m_targetDistance				= 1.0f;
	m_viewPosition					= float4(0,0,-1,1);
	m_viewDirection					= float4(0,0,1,0);
	m_viewUp						= float4(0,1,0,0);
	m_moveSpeed						= 15.0f;
	m_moveVelocity					= float4(0,0,0,0);
	m_mouseSensitivityRotation		= 0.5f;
	m_mouseSensitivityWheelZoom		= 0.1f;
	m_mouseSensitivityDragZoom		= 0.002f;
	m_viewMatrix					= float4x4::kIdentity;
		
	UpdateViewMatrix();
}

bool ArcBallCameraController::ProcessInput(const InputCollector& inputState)
{
	auto currentTime_ns = std::chrono::high_resolution_clock::now();
	
	std::chrono::nanoseconds duration_ns = currentTime_ns - m_previousTime_ns;
	float frameDuration_s = duration_ns.count() / 1e9f;

	m_previousTime_ns = currentTime_ns;

	m_moveVelocity.Zero();

	const bool isPerspective = m_projectionMtx._23 > 0;
	const bool isOrtho = !isPerspective;

	if (inputState.IsMouseButtonDown(Mouse::Right))
	{
			 if (inputState.IsKeyDown(Keyboard::W))		m_moveVelocity.z =  m_moveSpeed;
		else if (inputState.IsKeyDown(Keyboard::S))		m_moveVelocity.z = -m_moveSpeed;

			 if (inputState.IsKeyDown(Keyboard::D))		m_moveVelocity.x =  m_moveSpeed;
		else if (inputState.IsKeyDown(Keyboard::A))		m_moveVelocity.x = -m_moveSpeed;

			 if (inputState.IsKeyDown(Keyboard::E))		m_moveVelocity.y =  m_moveSpeed;
		else if (inputState.IsKeyDown(Keyboard::Q))		m_moveVelocity.y = -m_moveSpeed;

			 if (inputState.IsKeyDown(Keyboard::Space))	m_moveVelocity.y =  m_moveSpeed;
		else if (inputState.IsKeyDown(Keyboard::C))		m_moveVelocity.y = -m_moveSpeed;
	}

	if (inputState.IsKeyDown(Keyboard::ShiftLeft))		m_moveVelocity *= 3.0f;

	//m_moveVelocity *= powf(m_targetDistance, 3);

	if (m_targetDistance < 150)
		m_moveVelocity *= powf(1.1f, m_targetDistance * 0.213f) - 1;
	else
		m_moveVelocity *= (m_targetDistance - 150) * (2.0f / 5.0f) + 20.0f;

	// Camera rotation
	if (isPerspective)
	{
		if ( (inputState.IsKeyDown(Keyboard::AltLeft) && inputState.IsMouseButtonDown(Mouse::Left)) || (!inputState.IsKeyDown(Keyboard::AltLeft) && inputState.IsMouseButtonDown(Mouse::Right)))
		{
			float4x4 rotationMtx = float4x4::MakeRotationXYZ(-inputState.GetMouseRelative().y * m_mouseSensitivityRotation, -inputState.GetMouseRelative().x * m_mouseSensitivityRotation, 0);

			float4 newDirectionViewSpace = float4(rotationMtx.Column(2).xyz, 0);

			m_viewDirection = float4(mul(m_viewMatrix, newDirectionViewSpace).xyz, 0);

			if (abs(dot(m_viewDirection, float4::k0100)) >= 0.999f)
			{
				m_viewUp = float4::k0010;
			}
			else
			{
				m_viewUp = float4::k0100;
			}

			float4 rightVector = normalize(cross(m_viewUp, m_viewDirection));

			m_viewUp = normalize(cross(m_viewDirection, rightVector));

			// FPS style movement
			if (inputState.IsMouseButtonDown(Mouse::Right))
			{
				m_lookTarget = m_viewPosition + m_viewDirection * m_targetDistance;
			}

			m_dirty = true;
		}
	}

	// Camera pan
	if ( /*(inputState.IsKeyDown(Keyboard::AltLeft) || isOrtho) &&*/ inputState.IsMouseDragging(Mouse::Middle))
	{
		if (!m_draggingMiddle)
		{
			m_draggingMiddle = true;
			m_lookTarget_dragStart = m_lookTarget;
		}

		float2 mousePosStartNDC = inputState.GetMouseDragStartPos(Mouse::Middle) * float2(2.0f, -2.0f) + float2(-1.0f, 1.0f);
		float2 mousePosNDC = inputState.GetMousePos() * float2(2.0f, -2.0f) + float2(-1.0f, 1.0f);

		float N = -m_projectionMtx._32 / m_projectionMtx._22;

		float viewSpaceZ = mul(m_lookTarget, m_viewMatrix).z;

		float3 rightVector = normalize(cross(m_viewUp, m_viewDirection));

		m_lookTarget = m_lookTarget_dragStart;
		m_lookTarget.xyz -= rightVector * (mousePosNDC.x - mousePosStartNDC.x) * viewSpaceZ / m_projectionMtx._00;
		m_lookTarget.xyz -= m_viewUp * (mousePosNDC.y - mousePosStartNDC.y) * viewSpaceZ / m_projectionMtx._11;

		m_dirty = true;
	}
	else if (inputState.IsMouseButtonUp(Mouse::Middle) && m_draggingMiddle)
	{
		m_draggingMiddle = false;
	}

	// Zoom in/out
	if (inputState.IsKeyDown(Keyboard::AltLeft) && inputState.IsMouseButtonDown(Mouse::Right))
	{
		if (!m_draggingRight)
		{
			m_draggingRight = true;
			m_targetDistance_dragStart = m_targetDistance;
		}

		float2 mousePosStartNDC = inputState.GetMouseDragStartPos(Mouse::Right) * float2(2.0f, -2.0f) + float2(-1.0f, 1.0f);
		float2 mousePosNDC = inputState.GetMousePos() * float2(2.0f, -2.0f) + float2(-1.0f, 1.0f);

		m_targetDistance = m_targetDistance_dragStart * powf(std::max(0.0f, 1 + (mousePosNDC.x - mousePosStartNDC.x) * 0.5f), 2.0f );

		float N = -m_projectionMtx._32 / m_projectionMtx._22;

		m_targetDistance = std::min( 2000000.0f, std::max(N * 1.5f, m_targetDistance) );

		assert(std::isfinite(m_targetDistance));

		m_orthoZoom *= 1.0f + inputState.GetMouseRelative().x * m_mouseSensitivityDragZoom;

		m_dirty = true;
	}
	else if (inputState.IsMouseButtonUp(Mouse::Middle) && m_draggingRight)
	{
		m_draggingRight = false;
	}

	// Zoom in/out
	if (inputState.GetMouseWheelDelta() != 0)
	{
		m_targetDistance += -inputState.GetMouseWheelDelta() * m_mouseSensitivityWheelZoom;

		m_targetDistance = std::max(0.0f, m_targetDistance);

		m_orthoZoom *= 1.0f + inputState.GetMouseWheelDelta() * m_mouseSensitivityWheelZoom;

		m_dirty = true;
	}

	if (m_moveVelocity.lengthSqr()>0)
	{
		m_lookTarget.xyz += mul(m_viewMatrix, m_moveVelocity * frameDuration_s).xyz;
		m_dirty = true;
	}

	//if (inputState.WasKeyReleased(Keyboard::Space))
	//{
	//	float4 absDirection = abs(m_viewDirection);

	//		 if (absDirection.x >= absDirection.y && absDirection.x >= absDirection.z)		{ m_viewDirection = float4(sign(m_viewDirection.x), 0, 0, 0); m_viewUp = float4::k0100; }
	//	else if (absDirection.y >= absDirection.x && absDirection.y >= absDirection.z)		{ m_viewDirection = float4(0, sign(m_viewDirection.y), 0, 0); m_viewUp = float4::k0010; }
	//	else if (absDirection.z >= absDirection.x && absDirection.z >= absDirection.y)		{ m_viewDirection = float4(0, 0, sign(m_viewDirection.z), 0); m_viewUp = float4::k0100; }

	//	m_lookTarget = float4::k0000;

	//	m_dirty = true;
	//}

	if (m_dirty)
	{
		UpdateViewMatrix();
	}

	return false;
}

void ArcBallCameraController::UpdateViewMatrix()
{
	if (m_dirty)
	{
		m_viewPosition = m_lookTarget - m_viewDirection * m_targetDistance;
		m_viewMatrix.LookAlong(m_viewPosition, m_viewDirection, m_viewUp);
		m_dirty = false;
	}
}

void ArcBallCameraController::SetViewDistance(float distance)
{
	m_targetDistance = distance;
	m_dirty = true;
}

