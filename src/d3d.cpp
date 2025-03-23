#include "global.h"
#include "d3d.h"
#include "AppInterface.h"

#define D3DDEBUG	1

IDXGIFactory4*					g_dxgiFactory;
ID3D12Device*					g_device;
ID3D12DebugDevice*				g_debugDevice;
IDXGISwapChain3*				g_swapChain;
ID3D12Resource*					g_backBuffer[2];
HANDLE							g_swapChainWaitObject;
ID3D12CommandQueue*				g_commandQueue;
ID3D12DescriptorHeap*			g_rtvHeap;
ID3D12DescriptorHeap*			g_dsvHeap;
D3D12_CPU_DESCRIPTOR_HANDLE		g_rtvHandle[2];
ID3D12Resource*					g_depthStencilBuffer;
D3D12_CPU_DESCRIPTOR_HANDLE		g_dsvHandle;
ID3D12CommandAllocator*			g_commandAllocator[2];
ID3D12GraphicsCommandList*		g_commandList[2];
ID3D12GraphicsCommandList*		g_currentCommandList;
ID3D12Fence*					g_fence;
HANDLE							g_fenceEvent;
UINT64							g_fenceValue = 0;
ID3D12Resource*					g_vertexBuffer;
D3D12_VERTEX_BUFFER_VIEW		g_vertexBufferView;
ID3D12Resource*					g_indexBuffer;
D3D12_INDEX_BUFFER_VIEW			g_indexBufferView;
ID3DBlob*						g_vertexShaderBlob;
ID3DBlob*						g_pixelShaderBlob;
ID3D12RootSignature*			g_rootSignatureObject;
ID3D12PipelineState*			g_pso_Tri_Fill;
ID3D12PipelineState*			g_pso_Tri_Wire;
ID3D12PipelineState*			g_pso_Lines;
ID3D12RootSignature*			g_rootSignature;
ID3D12Resource*					g_cbCamera;
ID3D12Resource*					g_sbModelInstance;
bool							g_insideAppRender = false;
const uint32_t					g_maxVertexCount = 2*1024*1024;
const uint32_t					g_maxIndexCount = 2*1024*1024;
const uint32_t					g_maxModelInstanceCount = 2*1024;
uint32_t						g_currentVertexCount = 0;
uint32_t						g_currentIndexCount = 0;
uint32_t						g_currentModelInstanceCount = 0;

const char* g_vertexShaderCode = R"(
    cbuffer Camera : register(b0)
	{
        float4x4 g_viewProjMtx;
    };

	struct ModelInstanceCB
	{
		float4x4 modelMtx;
		float4	 color;
	};

    StructuredBuffer<ModelInstanceCB> g_ModelInstanceData : register(t0);

	cbuffer PushConstants : register(b1)
	{
        uint4 g_pushConstants;
    };

    struct VSInput 
	{
        float3 position : POSITION;
    };

    struct PSInput 
	{
        float4 position : SV_POSITION;
    };

    PSInput main(VSInput input) 
	{
        PSInput output;

        float4 pos = float4(input.position, 1.0);

        output.position = mul( mul(pos, g_ModelInstanceData[g_pushConstants.x].modelMtx), g_viewProjMtx);

        return output;
    }
)";

const char* g_pixelShaderCode = R"(
    cbuffer Camera : register(b0)
	{
        float4x4 g_viewProjMtx;
    };

	struct ModelInstanceCB
	{
		float4x4 modelMtx;
		float4	 color;
	};

    StructuredBuffer<ModelInstanceCB> g_ModelInstanceData : register(t0);

	cbuffer PushConstants : register(b1)
	{
        uint4 g_pushConstants;
    };

    struct PSInput 
	{
        float4 position : SV_POSITION;
    };

    float4 main(PSInput input) : SV_TARGET 
	{
        return g_ModelInstanceData[g_pushConstants.x].color;
    }
)";

struct CameraCB
{
	float4x4 viewProjMtx;
};

struct ModelInstanceCB
{
	float4x4 modelMtx;
	float4	 color;
};

void InitD3D(HWND hwnd) 
{
	UINT dxgiFactoryFlags = 0;

#if D3DDEBUG
	ComPtr<ID3D12Debug> debugController;

	if (SUCCEEDED(D3D12GetDebugInterface(IID_PPV_ARGS(&debugController))))
	{
		debugController->EnableDebugLayer();

		dxgiFactoryFlags |= DXGI_CREATE_FACTORY_DEBUG;
	}
#endif

	CreateDXGIFactory2(dxgiFactoryFlags, IID_PPV_ARGS(&g_dxgiFactory));

	D3D12CreateDevice(nullptr, D3D_FEATURE_LEVEL_11_0, IID_PPV_ARGS(&g_device));

#if D3DDEBUG
	{
		g_device->QueryInterface(IID_PPV_ARGS(&g_debugDevice));

		ID3D12InfoQueue* pInfoQueue = nullptr;

		if (SUCCEEDED(g_device->QueryInterface(IID_PPV_ARGS(&pInfoQueue))))
		{
			pInfoQueue->SetBreakOnSeverity(D3D12_MESSAGE_SEVERITY_CORRUPTION, true);
			pInfoQueue->SetBreakOnSeverity(D3D12_MESSAGE_SEVERITY_ERROR, true);
			pInfoQueue->SetBreakOnSeverity(D3D12_MESSAGE_SEVERITY_WARNING, true);

			D3D12_MESSAGE_ID ignoreIDs[] =
			{
				D3D12_MESSAGE_ID_UNKNOWN,
				// Workarounds for debug layer issues on hybrid-graphics systems
				//D3D12_MESSAGE_ID_EXECUTECOMMANDLISTS_WRONGSWAPCHAINBUFFERREFERENCE,
				//D3D12_MESSAGE_ID_RESOURCE_BARRIER_MISMATCHING_COMMAND_LIST_TYPE,
				//D3D12_MESSAGE_ID_DEVICE_REMOVAL_PROCESS_AT_FAULT,
				//D3D12_MESSAGE_ID_LIVE_DEVICE
			};

			//D3D12_MESSAGE_CATEGORY ignoreCategories[] =
			//{
			//};

			D3D12_MESSAGE_SEVERITY ignoreSeverities[] = 
			{
				D3D12_MESSAGE_SEVERITY_INFO,
			};

			D3D12_INFO_QUEUE_FILTER filter = {};
			filter.DenyList.NumIDs = ARRAYSIZE(ignoreIDs);
			filter.DenyList.pIDList = ignoreIDs;
			filter.DenyList.NumSeverities = ARRAYSIZE(ignoreSeverities);
			filter.DenyList.pSeverityList = ignoreSeverities;

			pInfoQueue->PushStorageFilter(&filter);
			pInfoQueue->Release();
		}

		ID3D12Debug* pDebug = nullptr;

		if (SUCCEEDED(D3D12GetDebugInterface(IID_PPV_ARGS(&pDebug))))
		{
			pDebug->EnableDebugLayer();
		}
	}
#endif

	D3D12_COMMAND_QUEUE_DESC queueDesc = {};
	queueDesc.Flags = D3D12_COMMAND_QUEUE_FLAG_NONE;
	queueDesc.Type = D3D12_COMMAND_LIST_TYPE_DIRECT;
	g_device->CreateCommandQueue(&queueDesc, IID_PPV_ARGS(&g_commandQueue));

	g_device->CreateCommandAllocator(D3D12_COMMAND_LIST_TYPE_DIRECT, IID_PPV_ARGS(&g_commandAllocator[0]));
	g_device->CreateCommandAllocator(D3D12_COMMAND_LIST_TYPE_DIRECT, IID_PPV_ARGS(&g_commandAllocator[1]));
	g_device->CreateCommandList(0, D3D12_COMMAND_LIST_TYPE_DIRECT, g_commandAllocator[0], nullptr, IID_PPV_ARGS(&g_commandList[0]));
	g_device->CreateCommandList(0, D3D12_COMMAND_LIST_TYPE_DIRECT, g_commandAllocator[1], nullptr, IID_PPV_ARGS(&g_commandList[1]));
	g_commandList[0]->Close();
	g_commandList[1]->Close();

	g_device->CreateFence(0, D3D12_FENCE_FLAG_NONE, IID_PPV_ARGS(&g_fence));
	g_fenceEvent = CreateEvent(nullptr, FALSE, FALSE, nullptr);

	// Descriptor heaps
	{
		D3D12_DESCRIPTOR_HEAP_DESC rtvHeapDesc = {};
		rtvHeapDesc.NumDescriptors = 2;
		rtvHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_RTV;
		g_device->CreateDescriptorHeap(&rtvHeapDesc, IID_PPV_ARGS(&g_rtvHeap));

		D3D12_DESCRIPTOR_HEAP_DESC dsvHeapDesc = {};
		dsvHeapDesc.NumDescriptors = 1;
		dsvHeapDesc.Type = D3D12_DESCRIPTOR_HEAP_TYPE_DSV;
		g_device->CreateDescriptorHeap(&dsvHeapDesc, IID_PPV_ARGS(&g_dsvHeap));
	}

	// Vertex & Index buffer
	{
		D3D12_HEAP_PROPERTIES resourceHeapProps;
		resourceHeapProps.Type					= D3D12_HEAP_TYPE_UPLOAD;
		resourceHeapProps.CPUPageProperty		= D3D12_CPU_PAGE_PROPERTY_UNKNOWN;
		resourceHeapProps.MemoryPoolPreference	= D3D12_MEMORY_POOL_UNKNOWN;
		resourceHeapProps.CreationNodeMask		= 0;
		resourceHeapProps.VisibleNodeMask		= 0;

		D3D12_RESOURCE_DESC resourceDesc;
		resourceDesc.Dimension			= D3D12_RESOURCE_DIMENSION_BUFFER;
		resourceDesc.Alignment			= 0;
		resourceDesc.Height				= 1;
		resourceDesc.DepthOrArraySize	= 1;
		resourceDesc.MipLevels			= 1;
		resourceDesc.Format				= DXGI_FORMAT_UNKNOWN;
		resourceDesc.SampleDesc			= {1,0};
		resourceDesc.Layout				= D3D12_TEXTURE_LAYOUT_ROW_MAJOR;
		resourceDesc.Flags				= D3D12_RESOURCE_FLAG_NONE;

		resourceDesc.Width = sizeof(float3) * g_maxVertexCount;

		g_device->CreateCommittedResource(	&resourceHeapProps,
											D3D12_HEAP_FLAG_NONE,
											&resourceDesc,
											D3D12_RESOURCE_STATE_GENERIC_READ,
											nullptr,
											IID_PPV_ARGS(&g_vertexBuffer) );

		resourceDesc.Width = sizeof(uint32_t) * g_maxIndexCount;

		g_device->CreateCommittedResource(	&resourceHeapProps,
											D3D12_HEAP_FLAG_NONE,
											&resourceDesc,
											D3D12_RESOURCE_STATE_GENERIC_READ,
											nullptr,
											IID_PPV_ARGS(&g_indexBuffer) );

		g_vertexBufferView.BufferLocation = g_vertexBuffer->GetGPUVirtualAddress();
		g_vertexBufferView.SizeInBytes = sizeof(float3) * g_maxVertexCount;
		g_vertexBufferView.StrideInBytes = sizeof(float3);

		g_indexBufferView.BufferLocation = g_indexBuffer->GetGPUVirtualAddress();
		g_indexBufferView.SizeInBytes = sizeof(uint32_t) * g_maxIndexCount;
		g_indexBufferView.Format = DXGI_FORMAT_R32_UINT;
	}

	// Constant Buffers
	{
		D3D12_HEAP_PROPERTIES resourceHeapProps;
		resourceHeapProps.Type					= D3D12_HEAP_TYPE_UPLOAD;
		resourceHeapProps.CPUPageProperty		= D3D12_CPU_PAGE_PROPERTY_UNKNOWN;
		resourceHeapProps.MemoryPoolPreference	= D3D12_MEMORY_POOL_UNKNOWN;
		resourceHeapProps.CreationNodeMask		= 0;
		resourceHeapProps.VisibleNodeMask		= 0;

		D3D12_RESOURCE_DESC resourceDesc;
		resourceDesc.Dimension			= D3D12_RESOURCE_DIMENSION_BUFFER;
		resourceDesc.Alignment			= 0;
		resourceDesc.Height				= 1;
		resourceDesc.DepthOrArraySize	= 1;
		resourceDesc.MipLevels			= 1;
		resourceDesc.Format				= DXGI_FORMAT_UNKNOWN;
		resourceDesc.SampleDesc			= {1,0};
		resourceDesc.Layout				= D3D12_TEXTURE_LAYOUT_ROW_MAJOR;
		resourceDesc.Flags				= D3D12_RESOURCE_FLAG_NONE;
		
		resourceDesc.Width = sizeof(CameraCB);

		g_device->CreateCommittedResource(	&resourceHeapProps,
											D3D12_HEAP_FLAG_NONE,
											&resourceDesc,
											D3D12_RESOURCE_STATE_GENERIC_READ,
											nullptr,
											IID_PPV_ARGS(&g_cbCamera) );
	}

	// Structured Buffers
	{
		D3D12_HEAP_PROPERTIES resourceHeapProps;
		resourceHeapProps.Type					= D3D12_HEAP_TYPE_UPLOAD;
		resourceHeapProps.CPUPageProperty		= D3D12_CPU_PAGE_PROPERTY_UNKNOWN;
		resourceHeapProps.MemoryPoolPreference	= D3D12_MEMORY_POOL_UNKNOWN;
		resourceHeapProps.CreationNodeMask		= 0;
		resourceHeapProps.VisibleNodeMask		= 0;

		D3D12_RESOURCE_DESC resourceDesc;
		resourceDesc.Dimension			= D3D12_RESOURCE_DIMENSION_BUFFER;
		resourceDesc.Alignment			= 0;
		resourceDesc.Height				= 1;
		resourceDesc.DepthOrArraySize	= 1;
		resourceDesc.MipLevels			= 1;
		resourceDesc.Format				= DXGI_FORMAT_UNKNOWN;
		resourceDesc.SampleDesc			= {1,0};
		resourceDesc.Layout				= D3D12_TEXTURE_LAYOUT_ROW_MAJOR;
		resourceDesc.Flags				= D3D12_RESOURCE_FLAG_NONE;

		resourceDesc.Width = sizeof(ModelInstanceCB) * g_maxModelInstanceCount;

		g_device->CreateCommittedResource(	&resourceHeapProps,
											D3D12_HEAP_FLAG_NONE,
											&resourceDesc,
											D3D12_RESOURCE_STATE_GENERIC_READ,
											nullptr,
											IID_PPV_ARGS(&g_sbModelInstance) );
	}

	// Root signature
	{
		D3D12_ROOT_PARAMETER rootParameters[3];

		rootParameters[0].ParameterType = D3D12_ROOT_PARAMETER_TYPE_CBV;
		rootParameters[0].ShaderVisibility = D3D12_SHADER_VISIBILITY_ALL;
		rootParameters[0].Descriptor.RegisterSpace = 0;
		rootParameters[0].Descriptor.ShaderRegister = 0;

		rootParameters[1].ParameterType = D3D12_ROOT_PARAMETER_TYPE_SRV;
		rootParameters[1].ShaderVisibility = D3D12_SHADER_VISIBILITY_ALL;
		rootParameters[1].Descriptor.RegisterSpace = 0;
		rootParameters[1].Descriptor.ShaderRegister = 0;

		rootParameters[2].ParameterType = D3D12_ROOT_PARAMETER_TYPE_32BIT_CONSTANTS;
		rootParameters[2].ShaderVisibility = D3D12_SHADER_VISIBILITY_ALL;
		rootParameters[2].Constants.ShaderRegister = 1;
		rootParameters[2].Constants.RegisterSpace = 0;
		rootParameters[2].Constants.Num32BitValues = 4;

		D3D12_ROOT_SIGNATURE_DESC rootSignatureDesc = {};
		rootSignatureDesc.NumParameters		= 3;
		rootSignatureDesc.pParameters		= rootParameters;
		rootSignatureDesc.NumStaticSamplers = 0;
		rootSignatureDesc.pStaticSamplers	= nullptr;
		rootSignatureDesc.Flags				= D3D12_ROOT_SIGNATURE_FLAG_ALLOW_INPUT_ASSEMBLER_INPUT_LAYOUT;

		ComPtr<ID3DBlob> signature;
		ComPtr<ID3DBlob> errorBlob;
		D3D12SerializeRootSignature(&rootSignatureDesc, D3D_ROOT_SIGNATURE_VERSION_1, signature.GetAddressOf(), errorBlob.GetAddressOf());

		g_device->CreateRootSignature(0, signature->GetBufferPointer(), signature->GetBufferSize(), IID_PPV_ARGS(&g_rootSignatureObject));
	}
   
	// Shaders & PSOs
	{
		ComPtr<ID3DBlob> errorBlob;

		V( D3DCompile(g_vertexShaderCode, strlen(g_vertexShaderCode), nullptr, nullptr, nullptr, "main", "vs_5_0", 0, 0, &g_vertexShaderBlob, errorBlob.GetAddressOf()) );
		const char* vsErrorStr = errorBlob ? (const char*)errorBlob->GetBufferPointer() : nullptr;

		V( D3DCompile(g_pixelShaderCode, strlen(g_pixelShaderCode), nullptr, nullptr, nullptr, "main", "ps_5_0", 0, 0, &g_pixelShaderBlob, errorBlob.GetAddressOf()) );
		const char* psErrorStr = errorBlob ? (const char*)errorBlob->GetBufferPointer() : nullptr;

		D3D12_GRAPHICS_PIPELINE_STATE_DESC defaultPSODesc = {};
	
		defaultPSODesc.pRootSignature			= g_rootSignatureObject;
		defaultPSODesc.VS						= { g_vertexShaderBlob->GetBufferPointer(), g_vertexShaderBlob->GetBufferSize() };
		defaultPSODesc.PS						= { g_pixelShaderBlob->GetBufferPointer(), g_pixelShaderBlob->GetBufferSize() };
		defaultPSODesc.PrimitiveTopologyType	= D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;
		defaultPSODesc.RTVFormats[0]			= DXGI_FORMAT_B8G8R8A8_UNORM;
		defaultPSODesc.DSVFormat				= DXGI_FORMAT_D32_FLOAT;
		defaultPSODesc.NumRenderTargets			= 1;
		defaultPSODesc.SampleDesc.Count			= 1;
		defaultPSODesc.SampleMask				= UINT_MAX;

		D3D12_INPUT_ELEMENT_DESC inputLayoutDesc[] =
		{
			{
				"position",									// LPCSTR SemanticName;
				0,											// UINT SemanticIndex;
				DXGI_FORMAT_R32G32B32_FLOAT,				// DXGI_FORMAT Format;
				0,											// UINT InputSlot;
				0,											// UINT AlignedByteOffset;
				D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, // D3D12_INPUT_CLASSIFICATION InputSlotClass;
				0,											// UINT InstanceDataStepRate;
			},
		};

		defaultPSODesc.InputLayout.pInputElementDescs	= inputLayoutDesc;
		defaultPSODesc.InputLayout.NumElements			= 1;


		defaultPSODesc.BlendState.AlphaToCoverageEnable					= FALSE;
		defaultPSODesc.BlendState.IndependentBlendEnable				= FALSE;
		defaultPSODesc.BlendState.RenderTarget[0].BlendEnable			= FALSE;
		defaultPSODesc.BlendState.RenderTarget[0].LogicOpEnable			= FALSE;
		defaultPSODesc.BlendState.RenderTarget[0].SrcBlend				= D3D12_BLEND_ONE;
		defaultPSODesc.BlendState.RenderTarget[0].DestBlend				= D3D12_BLEND_ZERO;
		defaultPSODesc.BlendState.RenderTarget[0].BlendOp				= D3D12_BLEND_OP_ADD;
		defaultPSODesc.BlendState.RenderTarget[0].SrcBlendAlpha			= D3D12_BLEND_ONE;
		defaultPSODesc.BlendState.RenderTarget[0].DestBlendAlpha		= D3D12_BLEND_ZERO;
		defaultPSODesc.BlendState.RenderTarget[0].BlendOpAlpha			= D3D12_BLEND_OP_ADD;
		defaultPSODesc.BlendState.RenderTarget[0].LogicOp				= D3D12_LOGIC_OP_NOOP;
		defaultPSODesc.BlendState.RenderTarget[0].RenderTargetWriteMask	= (1<<4)-1;

		defaultPSODesc.RasterizerState.FillMode					= D3D12_FILL_MODE_SOLID;
		defaultPSODesc.RasterizerState.CullMode					= D3D12_CULL_MODE_BACK;
		defaultPSODesc.RasterizerState.FrontCounterClockwise	= false;
		defaultPSODesc.RasterizerState.DepthBias				= 0;
		defaultPSODesc.RasterizerState.DepthBiasClamp			= 0;
		defaultPSODesc.RasterizerState.SlopeScaledDepthBias		= 0;
		defaultPSODesc.RasterizerState.DepthClipEnable			= TRUE;
		defaultPSODesc.RasterizerState.MultisampleEnable		= FALSE;
		defaultPSODesc.RasterizerState.AntialiasedLineEnable	= FALSE;
		defaultPSODesc.RasterizerState.ForcedSampleCount		= 0;
		defaultPSODesc.RasterizerState.ConservativeRaster		= D3D12_CONSERVATIVE_RASTERIZATION_MODE_OFF;

		defaultPSODesc.DepthStencilState.DepthEnable			= TRUE;
		defaultPSODesc.DepthStencilState.DepthWriteMask			= D3D12_DEPTH_WRITE_MASK_ALL;
		defaultPSODesc.DepthStencilState.DepthFunc				= D3D12_COMPARISON_FUNC_LESS_EQUAL;
		defaultPSODesc.DepthStencilState.StencilEnable			= FALSE;

		D3D12_GRAPHICS_PIPELINE_STATE_DESC psoDesc = defaultPSODesc;
		g_device->CreateGraphicsPipelineState(&psoDesc, IID_PPV_ARGS(&g_pso_Tri_Fill));

		psoDesc = defaultPSODesc;
		psoDesc.RasterizerState.FillMode = D3D12_FILL_MODE_WIREFRAME;
		psoDesc.RasterizerState.DepthBias = -5;
		g_device->CreateGraphicsPipelineState(&psoDesc, IID_PPV_ARGS(&g_pso_Tri_Wire));

		psoDesc = defaultPSODesc;
		psoDesc.PrimitiveTopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_LINE;
		g_device->CreateGraphicsPipelineState(&psoDesc, IID_PPV_ARGS(&g_pso_Lines));
	}
}

void ResizeD3D(HWND hwnd, uint32_t width, uint32_t height)
{
	// Color buffer
	{
		const uint32_t swapChainFlags = DXGI_SWAP_CHAIN_FLAG_FRAME_LATENCY_WAITABLE_OBJECT | DXGI_SWAP_CHAIN_FLAG_ALLOW_TEARING;

		if (g_swapChain == nullptr)
		{
			DXGI_SWAP_CHAIN_DESC1 swapChainDesc = {};

			swapChainDesc.Width				= width;
			swapChainDesc.Height			= height;
			swapChainDesc.Format			= DXGI_FORMAT_B8G8R8A8_UNORM;
			swapChainDesc.Stereo			= FALSE;
			swapChainDesc.SampleDesc.Count	= 1;
			swapChainDesc.BufferUsage		= DXGI_USAGE_RENDER_TARGET_OUTPUT;
			swapChainDesc.BufferCount		= 2;
			swapChainDesc.Scaling			= DXGI_SCALING_STRETCH;
			swapChainDesc.SwapEffect		= DXGI_SWAP_EFFECT_FLIP_DISCARD;
			swapChainDesc.AlphaMode			= DXGI_ALPHA_MODE_IGNORE;
			swapChainDesc.Flags				= swapChainFlags;

			Microsoft::WRL::ComPtr<IDXGISwapChain1> tempSwapChain;
			g_dxgiFactory->CreateSwapChainForHwnd(g_commandQueue, hwnd, &swapChainDesc, nullptr, nullptr, &tempSwapChain);
	
			Microsoft::WRL::ComPtr<IDXGISwapChain3> swapChain3;
			tempSwapChain.As(&swapChain3);

			g_swapChain = swapChain3.Get();
			g_swapChain->AddRef();

			uint32_t rtvIncrementSize = g_device->GetDescriptorHandleIncrementSize(D3D12_DESCRIPTOR_HEAP_TYPE_RTV);

			g_rtvHandle[0] = g_rtvHeap->GetCPUDescriptorHandleForHeapStart();
			g_rtvHandle[1] = g_rtvHandle[0];
			g_rtvHandle[1].ptr += rtvIncrementSize;

			g_swapChainWaitObject = g_swapChain->GetFrameLatencyWaitableObject();
		}
		else
		{
			//WaitForSingleObjectEx(g_swapChainWaitObject, 1000, true);

			SAFE_RELEASE(g_backBuffer[0]);
			SAFE_RELEASE(g_backBuffer[1]);

			V( g_swapChain->ResizeBuffers(	2, 
											width, 
											height, 
											DXGI_FORMAT_B8G8R8A8_UNORM,
											swapChainFlags ) );
		}

		g_swapChain->SetMaximumFrameLatency(1);

		g_swapChain->GetBuffer(0, IID_PPV_ARGS(&g_backBuffer[0]));
		g_swapChain->GetBuffer(1, IID_PPV_ARGS(&g_backBuffer[1]));

		D3D12_RENDER_TARGET_VIEW_DESC rtvDesc;
		rtvDesc.Format					= DXGI_FORMAT_B8G8R8A8_UNORM;
		rtvDesc.ViewDimension			= D3D12_RTV_DIMENSION_TEXTURE2D;
		rtvDesc.Texture2D.MipSlice		= 0;
		rtvDesc.Texture2D.PlaneSlice	= 0;

		g_device->CreateRenderTargetView(g_backBuffer[0], &rtvDesc, g_rtvHandle[0]);
		g_device->CreateRenderTargetView(g_backBuffer[1], &rtvDesc, g_rtvHandle[1]);
	}

	// Depth buffer
	{
		D3D12_RESOURCE_DESC depthStencilDesc = {};
		depthStencilDesc.Dimension			= D3D12_RESOURCE_DIMENSION_TEXTURE2D;
		depthStencilDesc.Width				= width;
		depthStencilDesc.Height				= height;
		depthStencilDesc.DepthOrArraySize	= 1;
		depthStencilDesc.MipLevels			= 1;
		depthStencilDesc.Format				= DXGI_FORMAT_D32_FLOAT;
		depthStencilDesc.SampleDesc.Count	= 1;
		depthStencilDesc.Layout				= D3D12_TEXTURE_LAYOUT_UNKNOWN;
		depthStencilDesc.Flags				= D3D12_RESOURCE_FLAG_ALLOW_DEPTH_STENCIL;

		D3D12_CLEAR_VALUE clearValue = {};
		clearValue.Format = DXGI_FORMAT_D32_FLOAT;
		clearValue.DepthStencil.Depth = 1.0f;
		clearValue.DepthStencil.Stencil = 0;

		D3D12_HEAP_PROPERTIES resourceHeapProps;

		resourceHeapProps.Type					= D3D12_HEAP_TYPE_DEFAULT;
		resourceHeapProps.CPUPageProperty		= D3D12_CPU_PAGE_PROPERTY_UNKNOWN;
		resourceHeapProps.MemoryPoolPreference	= D3D12_MEMORY_POOL_UNKNOWN;
		resourceHeapProps.CreationNodeMask		= 0;
		resourceHeapProps.VisibleNodeMask		= 0;

		g_device->CreateCommittedResource(	&resourceHeapProps,
											D3D12_HEAP_FLAG_NONE,
											&depthStencilDesc,
											D3D12_RESOURCE_STATE_DEPTH_WRITE,
											&clearValue,
											IID_PPV_ARGS(&g_depthStencilBuffer) );
		
		D3D12_DEPTH_STENCIL_VIEW_DESC dsvDesc;

		dsvDesc.Format				= DXGI_FORMAT_D32_FLOAT;
		dsvDesc.ViewDimension		= D3D12_DSV_DIMENSION_TEXTURE2D;
		dsvDesc.Flags				= D3D12_DSV_FLAG_NONE;
		dsvDesc.Texture2D.MipSlice	= 0;

		g_dsvHandle = g_dsvHeap->GetCPUDescriptorHandleForHeapStart();
		g_device->CreateDepthStencilView(g_depthStencilBuffer, &dsvDesc, g_dsvHandle);
	}
}

void RenderD3D(uint32_t windowWidth, uint32_t windowHeight) 
{
	WaitForSingleObjectEx(g_swapChainWaitObject, 1000, true);

	uint32_t currentBufferIndex = g_swapChain->GetCurrentBackBufferIndex();
	assert(currentBufferIndex<2);

	g_currentVertexCount = 0;
	g_currentIndexCount = 0;
	g_currentModelInstanceCount = 0;

	g_commandAllocator[currentBufferIndex]->Reset();
	g_commandList[currentBufferIndex]->Reset(g_commandAllocator[currentBufferIndex], nullptr);
	g_currentCommandList = g_commandList[currentBufferIndex];

	D3D12_RESOURCE_BARRIER barrier = {};
	barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
	barrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
	barrier.Transition.pResource = g_backBuffer[currentBufferIndex];
	barrier.Transition.Subresource = 0;
	barrier.Transition.StateBefore = D3D12_RESOURCE_STATE_PRESENT;
	barrier.Transition.StateAfter = D3D12_RESOURCE_STATE_RENDER_TARGET;

	g_currentCommandList->ResourceBarrier(1, &barrier);

	g_currentCommandList->OMSetRenderTargets(1, &g_rtvHandle[currentBufferIndex], FALSE, &g_dsvHandle);


	D3D12_VIEWPORT d3dViewport;
	d3dViewport.TopLeftX	= 0;
	d3dViewport.TopLeftY	= 0;
	d3dViewport.Width		= (float)windowWidth;
	d3dViewport.Height		= (float)windowHeight;
	d3dViewport.MinDepth	= 0;
	d3dViewport.MaxDepth	= 1;

	g_currentCommandList->RSSetViewports(1, &d3dViewport);

	D3D12_RECT d3dScissor;

	d3dScissor.left		= 0;
	d3dScissor.top		= 0;
	d3dScissor.right	= windowWidth;
	d3dScissor.bottom	= windowHeight;

	g_currentCommandList->RSSetScissorRects(1, &d3dScissor);

	float clearColor[] = { 0.2f, 0.4f, 0.6f, 1.0f };
	g_currentCommandList->ClearRenderTargetView(g_rtvHandle[currentBufferIndex], clearColor, 0, nullptr);
	g_currentCommandList->ClearDepthStencilView(g_dsvHandle, D3D12_CLEAR_FLAG_DEPTH, 1.0f, 0, 0, nullptr);

	void* pData;
	g_cbCamera->Map(0, nullptr, &pData);
	{
		CameraCB* cb = (CameraCB*)pData;

		float4x4 viewMtx = g_cameraController.GetViewMatrix();
		float4x4 projMtx = g_cameraController.GetProjectionMtx();

		cb->viewProjMtx = mul(viewMtx, projMtx);
	}
	g_cbCamera->Unmap(0, nullptr);
	
	g_currentCommandList->SetGraphicsRootSignature(g_rootSignatureObject);
    g_currentCommandList->IASetVertexBuffers(0, 1, &g_vertexBufferView);
	g_currentCommandList->IASetIndexBuffer(&g_indexBufferView);
	g_currentCommandList->SetGraphicsRootConstantBufferView(0, g_cbCamera->GetGPUVirtualAddress());
	g_currentCommandList->SetGraphicsRootShaderResourceView(1, g_sbModelInstance->GetGPUVirtualAddress() );

	{
		g_insideAppRender = true;

		AppRender(windowWidth, windowHeight);

		{
			std::vector<float3> axisPoints = 
			{
				float3(-1,0,0), float3(1,0,0),
				float3(0,-1,0), float3(0,1,0),
				float3(0,0,-1), float3(0,0,1),
			};

			DrawLines(axisPoints, matrix44::MakeTranslation(g_cameraController.GetTarget().xyz), float4::k1001);
		}

		g_insideAppRender = false;
	}

	barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
	barrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
	barrier.Transition.pResource = g_backBuffer[currentBufferIndex];
	barrier.Transition.Subresource = 0;
	barrier.Transition.StateBefore = D3D12_RESOURCE_STATE_RENDER_TARGET;
	barrier.Transition.StateAfter = D3D12_RESOURCE_STATE_PRESENT;

	g_currentCommandList->ResourceBarrier(1, &barrier);

	g_currentCommandList->Close();

	ID3D12CommandList* listsToExecute[] = { g_currentCommandList };
	g_commandQueue->ExecuteCommandLists(_countof(listsToExecute), listsToExecute);

	g_swapChain->Present(1, 0);

	const UINT64 currentFence = g_fenceValue;
	
	g_commandQueue->Signal(g_fence, currentFence);
	
	++g_fenceValue;

	if (g_fence->GetCompletedValue() < currentFence) 
	{
		g_fence->SetEventOnCompletion(currentFence, g_fenceEvent);
		WaitForSingleObject(g_fenceEvent, INFINITE);
	}
}

void DrawMesh(const std::vector<float3>& positions, const std::vector<uint32_t> indices, const matrix44& worldMtx, float4 fillColor, float4 wireColor)
{
	bool needsWire = wireColor != float4::k0000;

	assert(g_insideAppRender);
	assert( g_currentVertexCount + positions.size() <= g_maxVertexCount );
	assert( g_currentIndexCount + indices.size() <= g_maxIndexCount );
	assert( g_currentModelInstanceCount + needsWire ? 2 : 1 <= g_maxModelInstanceCount );

	void* pData;
	
	g_vertexBuffer->Map(0, nullptr, &pData);
	{
		memcpy( (float3*)pData + g_currentVertexCount, positions.data(), positions.size() * sizeof(float3));
	}
	g_vertexBuffer->Unmap(0, nullptr);

	g_indexBuffer->Map(0, nullptr, &pData);
	{
		memcpy( (uint32_t*)pData + g_currentIndexCount, indices.data(), indices.size() * sizeof(uint32_t));
	}
	g_indexBuffer->Unmap(0, nullptr);

	g_sbModelInstance->Map(0, nullptr, &pData);
	{
		ModelInstanceCB* cb = (ModelInstanceCB*)pData + g_currentModelInstanceCount;

		cb->modelMtx = worldMtx;
		cb->color = fillColor;
	}
	g_sbModelInstance->Unmap(0, nullptr);

	g_currentCommandList->SetPipelineState(g_pso_Tri_Fill);
	
	g_currentCommandList->IASetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
	
	uint4 instanceData(g_currentModelInstanceCount, 0, 0, 0);
	g_currentCommandList->SetGraphicsRoot32BitConstants(2, 4, &instanceData, 0);
    g_currentCommandList->DrawIndexedInstanced( (uint32_t)indices.size(), 1, g_currentIndexCount, g_currentVertexCount, 0);

	if (needsWire)
	{
		++g_currentModelInstanceCount;

		g_sbModelInstance->Map(0, nullptr, &pData);
		{
			ModelInstanceCB* cb = (ModelInstanceCB*)pData + g_currentModelInstanceCount;

			cb->modelMtx = worldMtx;
			cb->color = wireColor;
		}
		g_sbModelInstance->Unmap(0, nullptr);
	
		g_currentCommandList->SetPipelineState(g_pso_Tri_Wire);

		instanceData = uint4(g_currentModelInstanceCount, 0, 0, 0);
		g_currentCommandList->SetGraphicsRoot32BitConstants(2, 4, &instanceData, 0);

		g_currentCommandList->DrawIndexedInstanced( (uint32_t)indices.size(), 1, g_currentIndexCount, g_currentVertexCount, 0);
	}

	g_currentVertexCount += (uint32_t)positions.size();
	g_currentIndexCount += (uint32_t)indices.size();
	g_currentModelInstanceCount += 1;
}

void DrawLines(const std::vector<float3>& positions, const matrix44& worldMtx, float4 color)
{
	assert(g_insideAppRender);
	assert( g_currentVertexCount + positions.size() <= g_maxVertexCount );
	assert( g_currentModelInstanceCount + 1 <= g_maxModelInstanceCount );

	void* pData;
	
	g_vertexBuffer->Map(0, nullptr, &pData);
	{
		memcpy( (float3*)pData + g_currentVertexCount, positions.data(), positions.size() * sizeof(float3));
	}
	g_vertexBuffer->Unmap(0, nullptr);

	g_sbModelInstance->Map(0, nullptr, &pData);
	{
		ModelInstanceCB* cb = (ModelInstanceCB*)pData + g_currentModelInstanceCount;

		cb->modelMtx = worldMtx;
		cb->color = color;
	}
	g_sbModelInstance->Unmap(0, nullptr);

	g_currentCommandList->SetPipelineState(g_pso_Lines);

	g_currentCommandList->IASetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_LINELIST);

	uint4 instanceData(g_currentModelInstanceCount, 0, 0, 0);
	g_currentCommandList->SetGraphicsRoot32BitConstants(2, 4, &instanceData, 0);

	g_currentCommandList->DrawInstanced( (uint32_t)positions.size(), 1, g_currentVertexCount, 0);

	g_currentVertexCount += (uint32_t)positions.size();
	g_currentModelInstanceCount += 1;
}
