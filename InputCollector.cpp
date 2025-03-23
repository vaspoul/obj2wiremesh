#include "global.h"
#include "InputCollector.h"

const uint64_t InputCollector::m_mouseDoubleClickThreshold_us = 300 * 1000;

#define ZeroMem(x)				memset(&x, 0, sizeof(x))

InputCollector::InputCollector()
	: m_hasInputFocus(false)
{
	ZeroMem(m_mouseStateCurrent);
	ZeroMem(m_keyboardStateCurrent);

	ZeroMem(m_mouseStatePrevious);
	ZeroMem(m_keyboardStatePrevious);
}

InputCollector::~InputCollector()
{}

void InputCollector::ResetInput()
{
	m_mouseStatePrevious	= m_mouseStateCurrent;
	m_keyboardStatePrevious = m_keyboardStateCurrent;

	m_mouseStateCurrent.positionDelta = int2(0,0);
	m_mouseStateCurrent.wheelDelta = 0;
}

