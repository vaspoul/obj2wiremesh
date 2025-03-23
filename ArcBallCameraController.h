#pragma once

#include "InputCollector.h"

class ArcBallCameraController
{
public:
	ArcBallCameraController();

	bool				ProcessInput(const InputCollector& inputState);

	void				SetTarget(const float4& targetPos);
	void				SetViewPosition(const float4& viewPos);
	void				SetViewDistance(float distance);
	void				SetViewDirection(const float4& viewDirection, const float4& upVector = float4::k0100);
	void				SetMoveSpeed(float speed);

	const float4x4&		GetViewMatrix() const;
	const float4&		GetViewPosition() const;
	const float4&		GetViewDirection() const;
	const float4&		GetTarget() const;
	float				GetViewDistance() const;
	float				GetMoveSpeed() const;
	void				Reset();

	void				SetOrthoZoom(float zoom);
	float				GetOrthoZoom() const;

	void				SetProjectionMtx(const matrix44& projMtx);
	const float4x4&		GetProjectionMtx() const;

private:
	void UpdateViewMatrix();

private:
	float4		m_lookTarget;
	float		m_targetDistance;
	float4		m_viewPosition;
	float4		m_viewDirection;
	float4		m_viewUp;
	float		m_moveSpeed = 15;
	float4		m_moveVelocity;
	float		m_mouseSensitivityRotation		= 0.5f;
	float		m_mouseSensitivityWheelZoom		= 0.1f;
	float		m_mouseSensitivityDragZoom		= 0.002f;
	matrix44	m_viewMatrix;
	matrix44	m_projectionMtx;
	float		m_orthoZoom = 1;
	bool		m_dirty = true;

	bool		m_draggingMiddle = false;
	float4		m_lookTarget_dragStart;

	bool		m_draggingRight = false;
	float		m_targetDistance_dragStart;

	std::chrono::high_resolution_clock::time_point m_previousTime_ns;
};

inline const float4x4& ArcBallCameraController::GetViewMatrix() const
{
	return m_viewMatrix;
}

inline const float4& ArcBallCameraController::GetViewPosition() const
{
	return m_viewPosition;
}

inline void	ArcBallCameraController::SetViewPosition(const float4& viewPos)
{
	m_lookTarget = viewPos + m_viewDirection * m_targetDistance;
	m_dirty = true;
}

inline void	ArcBallCameraController::SetViewDirection(const float4& viewDirection, const float4& upVector)
{
	m_viewDirection = float4(normalize(viewDirection.xyz), 0);

	float4 rightVector = normalize(cross(upVector, m_viewDirection));

	m_viewUp = normalize(cross(m_viewDirection, rightVector));


	m_dirty = true;
}

inline const float4& ArcBallCameraController::GetViewDirection() const
{
	return m_viewDirection;
}

inline const float4& ArcBallCameraController::GetTarget() const
{
	return m_lookTarget;
}

inline void	ArcBallCameraController::SetTarget(const float4& targetPos)
{
	m_lookTarget = targetPos;
	m_dirty = true;
}

inline float ArcBallCameraController::GetViewDistance() const
{
	return m_targetDistance;
}

inline float ArcBallCameraController::GetMoveSpeed() const
{
	return m_moveSpeed;
}

inline void ArcBallCameraController::SetMoveSpeed(float speed)
{
	m_moveSpeed = speed;
}

inline void ArcBallCameraController::SetProjectionMtx(const matrix44& projMtx)
{
	m_projectionMtx = projMtx;
}

inline const float4x4& ArcBallCameraController::GetProjectionMtx() const
{
	return m_projectionMtx;
}

inline void ArcBallCameraController::SetOrthoZoom(float zoom)
{
	m_orthoZoom = zoom;
}

inline float ArcBallCameraController::GetOrthoZoom() const
{
	return m_orthoZoom;
}

