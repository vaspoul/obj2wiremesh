#include "global.h"
#include "d3d.h"
#include "ArcBallCameraController.h"
#include "InputCollector.h"
#include "AppInterface.h"
#include "imgui/imgui.h"
#include "imgui/imgui_impl_win32.h"

HINSTANCE						g_hInst;
HWND							g_hWND;
uint32_t						g_windowWidth;
uint32_t						g_windowHeight;
bool							g_quit = false;
InputCollector					g_inputCollector;
ArcBallCameraController			g_cameraController;

extern IMGUI_IMPL_API LRESULT ImGui_ImplWin32_WndProcHandler(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);

void Resize()
{
	RECT clientRect;
	::GetClientRect(g_hWND, &clientRect);

	g_windowWidth = clientRect.right - clientRect.left;
	g_windowHeight = clientRect.bottom - clientRect.top;

	g_inputCollector.OnWindowSize(int2(g_windowWidth, g_windowHeight));

	ResizeD3D(g_hWND, g_windowWidth, g_windowHeight);

	AppResize(g_windowWidth, g_windowHeight);
}

LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam) 
{
	ImGui_ImplWin32_WndProcHandler(hWnd, message, wParam, lParam);

	switch (message) 
	{
		case WM_CREATE:
			g_hWND = hWnd;
			break;

		case WM_CLOSE:
			g_quit = true;
			DestroyWindow(hWnd);
			break;

		case WM_DESTROY:
			PostQuitMessage(0);
			break;

		case WM_SIZE:
		{
			Resize();

			break;
		}

		case WM_SIZING:
		{
			RECT clientRect;
			::GetClientRect(g_hWND, &clientRect);

			g_windowWidth = clientRect.right - clientRect.left;
			g_windowHeight = clientRect.bottom - clientRect.top;

			RenderD3D(g_windowWidth, g_windowHeight, g_inputCollector.GetMousePos());

			break;
		}

		case WM_INPUT:
		{
			bool foreground = GET_RAWINPUT_CODE_WPARAM(wParam) == RIM_INPUT;

			RAWINPUT rawInput;
			uint32_t rawInputSize = sizeof(rawInput);
			GetRawInputData((HRAWINPUT)lParam, RID_INPUT, &rawInput, &rawInputSize, sizeof(RAWINPUTHEADER));

			if (rawInput.header.dwType == RIM_TYPEMOUSE)
			{
				if (!(rawInput.data.mouse.usFlags & MOUSE_MOVE_ABSOLUTE))
				{
					g_inputCollector.OnMouseMoveDelta(int2( rawInput.data.mouse.lLastX, rawInput.data.mouse.lLastY ));
				}
			}
			else if (rawInput.header.dwType == RIM_TYPEKEYBOARD)
			{
				// See https://blog.molecular-matters.com/2011/09/05/properly-handling-keyboard-input/

				USHORT virtualKey	= rawInput.data.keyboard.VKey;
				USHORT scanCode		= rawInput.data.keyboard.MakeCode;
				USHORT flags		= rawInput.data.keyboard.Flags;
				const bool isE0		= ((flags & RI_KEY_E0) != 0);
				const bool isE1		= ((flags & RI_KEY_E1) != 0);
				bool pressed		= (flags & RI_KEY_BREAK) == 0;

				if (virtualKey == 255)
				{
					// discard "fake keys" which are part of an escaped sequence
					return false;
				}
				else if (virtualKey == VK_SHIFT)
				{
					// correct left-hand / right-hand SHIFT
					virtualKey = (USHORT)MapVirtualKey(scanCode, MAPVK_VSC_TO_VK_EX);
				}
				else if (virtualKey == VK_NUMLOCK)
				{
					// correct PAUSE/BREAK and NUM LOCK silliness, and set the extended bit
					scanCode = ((USHORT)MapVirtualKey(virtualKey, MAPVK_VK_TO_VSC) | 0x100);
				}

				if (isE1)
				{
					// for escaped sequences, turn the virtual key into the correct scan code using MapVirtualKey.
					// however, MapVirtualKey is unable to map VK_PAUSE (this is a known bug), hence we map that by hand.
					if (virtualKey == VK_PAUSE)
						scanCode = 0x45;
					else
						scanCode = (USHORT)MapVirtualKey(virtualKey, MAPVK_VK_TO_VSC);
				}

				switch (virtualKey)
				{
					case VK_CONTROL:	virtualKey = isE0 ? Keyboard::ControlRight : virtualKey = Keyboard::ControlLeft;	break;
					case VK_MENU:		virtualKey = isE0 ? Keyboard::AltRight : virtualKey = Keyboard::AltLeft;			break;
					case VK_RETURN:		if ( isE0)	virtualKey = Keyboard::NumPadEnter;										break;
					case VK_INSERT:		if (!isE0)	virtualKey = Keyboard::NumPad0;											break;
					case VK_DELETE:		if (!isE0)	virtualKey = Keyboard::NumPadDecimal;									break;
					case VK_HOME:		if (!isE0)	virtualKey = Keyboard::NumPad7;											break;
					case VK_END:		if (!isE0)	virtualKey = Keyboard::NumPad1;											break;
					case VK_PRIOR:		if (!isE0)	virtualKey = Keyboard::NumPad9;											break;
					case VK_NEXT:		if (!isE0)	virtualKey = Keyboard::NumPad3;											break;
					case VK_LEFT:		if (!isE0)	virtualKey = Keyboard::NumPad4;											break;
 					case VK_RIGHT:		if (!isE0)	virtualKey = Keyboard::NumPad6;											break;
 					case VK_UP:			if (!isE0)	virtualKey = Keyboard::NumPad8;											break;
					case VK_DOWN:		if (!isE0)	virtualKey = Keyboard::NumPad2;											break;
					case VK_CLEAR:		if (!isE0)	virtualKey = Keyboard::NumPad5;											break;
				}

				g_inputCollector.OnKeyboardKey((Keyboard::Key)virtualKey, pressed ? Keyboard::Pressed : Keyboard::Released);
			}
		}
		break;

		case WM_LBUTTONDOWN:
		{
			SetCapture(g_hWND);
			g_inputCollector.OnMouseButton(Mouse::Left, Mouse::Pressed);
			break;
		}

		case WM_RBUTTONDOWN:
		{
			SetCapture(g_hWND);
			g_inputCollector.OnMouseButton(Mouse::Right, Mouse::Pressed);
			break;
		}

		case WM_MBUTTONDOWN:
		{
			SetCapture(g_hWND);
			g_inputCollector.OnMouseButton(Mouse::Middle, Mouse::Pressed);
			break;
		}

		case WM_MOUSEMOVE:
		{
			int mouseX = GET_X_LPARAM(lParam); 
			int mouseY = GET_Y_LPARAM(lParam);
			g_inputCollector.OnMouseMove(int2(mouseX, mouseY));
			break;
		}

		case WM_LBUTTONUP:
		{
			ReleaseCapture();
			g_inputCollector.OnMouseButton(Mouse::Left, Mouse::Released);
			break;
		}

		case WM_RBUTTONUP:
		{
			ReleaseCapture();
			g_inputCollector.OnMouseButton(Mouse::Right, Mouse::Released);
			break;
		}

		case WM_MBUTTONUP:
		{
			ReleaseCapture();
			g_inputCollector.OnMouseButton(Mouse::Middle, Mouse::Released);
			break;
		}

		case WM_MOUSEWHEEL:
		{
			g_inputCollector.OnMouseWheel((float)GET_WHEEL_DELTA_WPARAM(wParam) / (float)WHEEL_DELTA );
		}
		break;

		case WM_LBUTTONDBLCLK:
		{
			g_inputCollector.OnMouseButton(Mouse::Left, Mouse::Pressed);
		}
		break;

		case WM_RBUTTONDBLCLK:
		{
			g_inputCollector.OnMouseButton(Mouse::Right, Mouse::Pressed);
		}
		break;

		case WM_MBUTTONDBLCLK:
		{
			g_inputCollector.OnMouseButton(Mouse::Middle, Mouse::Pressed);
		}
		break;

		case WM_CHAR:
		{
			g_inputCollector.OnCharInput((uint16_t)wParam);
		}
		break;

		case WM_SETFOCUS:
		{
			g_inputCollector.OnInputFocus(true);
		}
		break;

		case WM_KILLFOCUS:
		{
			g_inputCollector.OnInputFocus(false);
		}
		break;
	}

    return DefWindowProc(hWnd, message, wParam, lParam);
}


int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE, LPSTR, int nCmdShow)
{
    g_hInst = hInstance;
    
	WNDCLASS wc = {};
	wc.lpfnWndProc		= WndProc;
	wc.hInstance		= g_hInst;
	wc.lpszClassName	= L"obj2mesh";
	wc.hbrBackground	= (HBRUSH)GetStockObject(DKGRAY_BRUSH);
	wc.hCursor			= LoadCursor(NULL, IDC_ARROW);
	RegisterClass(&wc);
    
	RECT windowRect = {0, 0, 1280, 720};
	AdjustWindowRect(&windowRect, WS_OVERLAPPEDWINDOW, FALSE); // WS_OVERLAPPEDWINDOW includes title bar and borders

	g_hWND = CreateWindow(	L"obj2mesh", 
							L"obj2mesh", 
							WS_OVERLAPPEDWINDOW,
							CW_USEDEFAULT, CW_USEDEFAULT, 
							windowRect.right - windowRect.left,
							windowRect.bottom - windowRect.top,
							nullptr, 
							nullptr, 
							g_hInst, 
							nullptr);
    
	// Register for WM_INPUT messages
	{
		RAWINPUTDEVICE devices[2];

		// Keyboard
		devices[0].usUsagePage	= HID_USAGE_PAGE_GENERIC;
		devices[0].usUsage		= HID_USAGE_GENERIC_KEYBOARD;
		devices[0].dwFlags		= 0;//RIDEV_NOLEGACY;
		devices[0].hwndTarget	= NULL; // follow keyboard focus

		// Mouse
		devices[1].usUsagePage	= HID_USAGE_PAGE_GENERIC;
		devices[1].usUsage		= HID_USAGE_GENERIC_MOUSE;
		devices[1].dwFlags		= 0;//RIDEV_NOLEGACY;
		devices[1].hwndTarget	= NULL; // follow keyboard focus

		RegisterRawInputDevices(devices, ARRAYSIZE(devices), sizeof(RAWINPUTDEVICE));
	}

	// imgui
	{
		// Setup Dear ImGui context
		IMGUI_CHECKVERSION();
		ImGui::CreateContext();
		ImGuiIO& io = ImGui::GetIO(); (void)io;
		io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
		io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

		// Setup Dear ImGui style
		ImGui::StyleColorsDark();

		ImGui_ImplWin32_Init(g_hWND);

	}

	InitD3D(g_hWND);

	::ShowWindow(g_hWND, SW_MAXIMIZE);//nCmdShow);
	::UpdateWindow(g_hWND);

	std::vector<std::string> args;

	{
		const char* commandLine = ::GetCommandLineA();
		
		std::stringstream ss(commandLine);

		std::string token;

		while (ss >> token) 
		{
			args.push_back(token);
		}
	}

	AppInit(args);

	if (g_quit)
		return 0;
	

	while (!g_quit)
	{
		g_inputCollector.ResetInput();

		MSG msg;
		while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
		{
			if (msg.message == WM_QUIT)
				break;

			TranslateMessage(&msg);
			DispatchMessage(&msg);
		}

        ImGui_ImplWin32_NewFrame();

		if (!ImGui::GetIO().WantCaptureKeyboard)
		{
			g_cameraController.ProcessInput(g_inputCollector);
		}

		RenderD3D(g_windowWidth, g_windowHeight, g_inputCollector.GetMousePos());
	}

	AppExit();

	return 0;
}

void Quit()
{
	g_quit = true;
}

void SetWindowTitle(const std::string& title)
{
	::SetWindowTextA(g_hWND, title.c_str());
}
