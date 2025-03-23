#pragma once

struct Mouse
{
	enum Button : uint8_t
	{
		Left,
		Right,
		Middle,
		Extra1,
		Extra2,

		ButtonCount
	};

	enum ButtonState : uint8_t
	{
		Released = 0,
		Pressed = 1,
	};
};

struct Keyboard
{
	enum Key : uint8_t
	{
		// Numeric values match * from WinUser.h
		Backspace		= 0x08,
		Tab				= 0x09,
		//CLEAR			= 0x0C,
		Return			= 0x0D,
		SHIFT			= 0x10,
		CONTROL			= 0x11,
		ALT				= 0x12,
		Pause			= 0x13,
		CapsLock		= 0x14,
		Escape			= 0x1B,
		Space			= 0x20,
		PageUp			= 0x21,
		PageDown		= 0x22,
		End				= 0x23,
		Home			= 0x24,
		ArrowLeft		= 0x25,
		ArrowUp			= 0x26,
		ArrowRight		= 0x27,
		ArrowDown		= 0x28,
		//SELECT		= 0x29,
		Print			= 0x2A,
		//EXECUTE		= 0x2B,
		//SNAPSHOT		= 0x2C,
		Insert			= 0x2D,
		Delete			= 0x2E,
		//HELP			= 0x2F,
		Number0			= 0x30,			// 0 - 9 are the same as ASCII '0' - '9' (0x30 - 0x39)
		Number1			= 0x31,
		Number2			= 0x32,
		Number3			= 0x33,
		Number4			= 0x34,
		Number5			= 0x35,
		Number6			= 0x36,
		Number7			= 0x37,
		Number8			= 0x38,
		Number9			= 0x39,
		A				= 0x41,	// A - Z are the same as ASCII 'A' - 'Z' (0x41 - 0x5A)
		B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V, W, X, Y, Z,
		WinKeyLeft		= 0x5B,
		WinKeyRight		= 0x5C,
		//APPS			= 0x5D,
		//SLEEP			= 0x5F,
		NumPad0			= 0x60,
		NumPad1			= 0x61,
		NumPad2			= 0x62,
		NumPad3			= 0x63,
		NumPad4			= 0x64,
		NumPad5			= 0x65,
		NumPad6			= 0x66,
		NumPad7			= 0x67,
		NumPad8			= 0x68,
		NumPad9			= 0x69,
		NumPadMultiply	= 0x6A,
		NumPadAdd		= 0x6B,
		//SEPARATOR	= 0x6C,
		NumPadSubtract	= 0x6D,
		NumPadDecimal	= 0x6E,
		NumPadDivide	= 0x6F,
		F1				= 0x70,
		F2				= 0x71,
		F3				= 0x72,
		F4				= 0x73,
		F5				= 0x74,
		F6				= 0x75,
		F7				= 0x76,
		F8				= 0x77,
		F9				= 0x78,
		F10				= 0x79,
		F11				= 0x7A,
		F12				= 0x7B,
		F13				= 0x7C,
		F14				= 0x7D,
		F15				= 0x7E,
		F16				= 0x7F,
		F17				= 0x80,
		F18				= 0x81,
		F19				= 0x82,
		F20				= 0x83,
		F21				= 0x84,
		F22				= 0x85,
		F23				= 0x86,
		F24				= 0x87,
		Numlock			= 0x90,
		Scroll			= 0x91,
		ShiftLeft		= 0xA0,
		ShiftRight		= 0xA1,
		ControlLeft		= 0xA2,
		ControlRight	= 0xA3,
		AltLeft			= 0xA4,
		AltRight		= 0xA5,
		NumPadEnter,

		COUNT = 0xFF
	};

	enum KeyState : uint8_t
	{
		Released = 0,
		Pressed = 1,
	};
};

class InputCollector
{
public:
	InputCollector();
	~InputCollector();

	void			ResetInput();

	bool			GetHasInputFocus() const;

	bool			IsKeyboardConnected() const;
	bool			IsKeyDown(Keyboard::Key key) const;
	bool			IsKeyUp(Keyboard::Key key) const;
	bool			IsShiftDown() const;
	bool			IsControlDown() const;
	bool			IsAltDown() const;
	bool			WasKeyPressed(Keyboard::Key key) const;
	bool			WasKeyReleased(Keyboard::Key key) const;

	bool			IsMouseConnected() const;
	const float2&	GetMousePos() const;
	const int2&		GetMousePosPixels() const;
	const int2&		GetMouseRelative() const;
	float			GetMouseWheelDelta() const;
	bool			IsAnyMouseButtonDown() const;
	bool			IsMouseButtonDown(Mouse::Button button) const;
	bool			IsMouseButtonUp(Mouse::Button button) const;
	bool			WasMouseButtonPressed(Mouse::Button button) const;
	bool			WasMouseButtonReleased(Mouse::Button button) const;
	const float2&	GetWindowSize() const;

	bool			IsMouseDragging(Mouse::Button button) const;
	int2			GetMouseDragDelta(Mouse::Button button) const;
	float2			GetMouseDragStartPos(Mouse::Button button) const;
	int2			GetMouseDragStartPosPixels(Mouse::Button button) const;

	void			OnInputFocus(bool focus);
	void			OnKeyboardConnected(bool connected);
	void			OnKeyboardKey(Keyboard::Key key, Keyboard::KeyState state);
	void			OnKeyboardModifiers(bool shiftDown, bool controlDown, bool altDown);
	void			OnCharInput(const uint16_t c);

	void			OnWindowSize(const int2& windowSize);
	void			OnMouseConnected(bool connected);
	void			OnMouseMove(const int2& pos);
	void			OnMouseMoveDelta(const int2& delta);
	void			OnMouseButton(Mouse::Button button, Mouse::ButtonState state);
	void			OnMouseWheel(float delta);

private:


	struct KeyboardState
	{
		bool	connected;
		bool	keyDown[Keyboard::COUNT];
	};

	struct MouseState
	{
		int2				positionPixels;
		float2				positionNormalized;
		int2				positionDelta;
		float				wheelDelta;
		Mouse::ButtonState	buttonDown[Mouse::ButtonCount];
		bool				connected;
	};

private:
	float2					m_windowSize;
	bool					m_hasInputFocus;

	MouseState				m_mouseStateCurrent;
	KeyboardState			m_keyboardStateCurrent;

	MouseState				m_mouseStatePrevious;
	KeyboardState			m_keyboardStatePrevious;

	int2					m_mouseDownPosPixels[Mouse::ButtonCount];
	float2					m_mouseDownPos[Mouse::ButtonCount];
	static const int		m_mouseDragDistanceThreshold = 5;	// in pixels
	static const uint64_t	m_mouseDoubleClickThreshold_us;
};

inline bool	InputCollector::GetHasInputFocus() const
{
	return m_hasInputFocus;
}

inline bool InputCollector::IsKeyboardConnected() const
{
	return m_keyboardStateCurrent.connected;
}

inline bool InputCollector::IsKeyDown(Keyboard::Key key) const
{
	return m_keyboardStateCurrent.keyDown[key];
}

inline bool InputCollector::IsKeyUp(Keyboard::Key key) const
{
	return !m_keyboardStateCurrent.keyDown[key];
}

inline bool	InputCollector::IsShiftDown() const
{
	return m_keyboardStateCurrent.keyDown[Keyboard::SHIFT];
}

inline bool	InputCollector::IsControlDown() const
{
	return m_keyboardStateCurrent.keyDown[Keyboard::CONTROL];
}

inline bool	InputCollector::IsAltDown() const
{
	return m_keyboardStateCurrent.keyDown[Keyboard::ALT];
}

inline bool InputCollector::WasKeyPressed(Keyboard::Key key) const
{
	return m_keyboardStateCurrent.keyDown[key] && !m_keyboardStatePrevious.keyDown[key];
}

inline bool InputCollector::WasKeyReleased(Keyboard::Key key) const
{
	return !m_keyboardStateCurrent.keyDown[key] && m_keyboardStatePrevious.keyDown[key];
}

inline bool InputCollector::IsMouseConnected() const
{
	return m_mouseStateCurrent.connected;
}

inline const float2& InputCollector::GetMousePos() const
{
	return m_mouseStateCurrent.positionNormalized;
}

inline const int2& InputCollector::GetMousePosPixels() const
{
	return m_mouseStateCurrent.positionPixels;
}

inline const int2& InputCollector::GetMouseRelative() const
{
	return m_mouseStateCurrent.positionDelta;
}

inline float InputCollector::GetMouseWheelDelta() const
{
	return m_mouseStateCurrent.wheelDelta;
}

inline bool InputCollector::IsAnyMouseButtonDown() const
{
	for (int i=0; i!=Mouse::ButtonCount; ++i)
		if (m_mouseStateCurrent.buttonDown[i])
			return true;

	return false;
}

inline bool InputCollector::IsMouseButtonDown(Mouse::Button button) const
{
	return m_mouseStateCurrent.buttonDown[button];
}

inline bool InputCollector::IsMouseButtonUp(Mouse::Button button) const
{
	return !m_mouseStateCurrent.buttonDown[button];
}

inline bool InputCollector::WasMouseButtonPressed(Mouse::Button button) const
{
	return m_mouseStateCurrent.buttonDown[button] && !m_mouseStatePrevious.buttonDown[button];
}

inline bool InputCollector::WasMouseButtonReleased(Mouse::Button button) const
{
	return !m_mouseStateCurrent.buttonDown[button] && m_mouseStatePrevious.buttonDown[button];
}

inline const float2& InputCollector::GetWindowSize() const
{
	return m_windowSize;
}


inline bool InputCollector::IsMouseDragging(Mouse::Button button) const
{
	// If it's down or just released
	if (m_mouseStateCurrent.buttonDown[button])// || m_mouseStatePrevious.buttonDown[button])
	{
		int2 delta = GetMouseDragDelta(button);

		if (std::max(std::abs(delta.x), std::abs(delta.y)) >= m_mouseDragDistanceThreshold)
		{
			return true;
		}
	}

	return false;
}

inline int2 InputCollector::GetMouseDragDelta(Mouse::Button button) const
{
	return GetMousePosPixels() - GetMouseDragStartPosPixels(button);
}

inline float2 InputCollector::GetMouseDragStartPos(Mouse::Button button) const
{
	return m_mouseDownPos[button];
}

inline int2 InputCollector::GetMouseDragStartPosPixels(Mouse::Button button) const
{
	return m_mouseDownPosPixels[button];
}



inline void	InputCollector::OnInputFocus(bool focus)
{
	if (!m_hasInputFocus && focus)
	{
		ZeroMem(m_keyboardStateCurrent.keyDown);
		ZeroMem(m_keyboardStatePrevious.keyDown);
	}

	m_hasInputFocus = focus;
}

inline void InputCollector::OnKeyboardConnected(bool connected)
{
	m_keyboardStateCurrent.connected = connected;
}

inline void InputCollector::OnKeyboardKey(Keyboard::Key key, Keyboard::KeyState state)
{
	m_keyboardStateCurrent.keyDown[key] = state == Keyboard::Pressed;
}

inline void InputCollector::OnKeyboardModifiers(bool shiftDown, bool controlDown, bool altDown)
{
	m_keyboardStateCurrent.keyDown[Keyboard::SHIFT] = shiftDown;
	m_keyboardStateCurrent.keyDown[Keyboard::CONTROL] = controlDown;
	m_keyboardStateCurrent.keyDown[Keyboard::ALT] = altDown;
}

inline void InputCollector::OnCharInput(const uint16_t c)
{
}

inline void InputCollector::OnWindowSize(const int2& windowSize)
{
	m_windowSize = float2((float)windowSize.x, (float)windowSize.y);
}

inline void InputCollector::OnMouseConnected(bool connected)
{
	m_mouseStateCurrent.connected = connected;
}

inline void InputCollector::OnMouseMove(const int2& pos)
{
	m_mouseStateCurrent.positionPixels = pos;
	m_mouseStateCurrent.positionNormalized = float2((float)pos.x, (float)pos.y) / m_windowSize;
}

inline void InputCollector::OnMouseMoveDelta(const int2& delta)
{
	m_mouseStateCurrent.positionDelta += delta;
}

inline void InputCollector::OnMouseButton(Mouse::Button button, Mouse::ButtonState state)
{
	m_mouseStateCurrent.buttonDown[button] = state;

	if (m_mouseStateCurrent.buttonDown[button] == Mouse::Pressed && m_mouseStatePrevious.buttonDown[button] == Mouse::Released)
	{
		m_mouseDownPos[button] = GetMousePos();
		m_mouseDownPosPixels[button] = GetMousePosPixels();
	}
}

inline void InputCollector::OnMouseWheel(float delta)
{
	m_mouseStateCurrent.wheelDelta += delta;
}
