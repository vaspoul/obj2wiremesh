#include "global.h"
#include "d3d.h"
#include "AppInterface.h"
#include "imgui/imgui.h"

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
ID3D12DescriptorHeap*			g_srvHeap;
D3D12_CPU_DESCRIPTOR_HANDLE		g_rtvHandle[2];
ID3D12Resource*					g_depthStencilBuffer;
D3D12_CPU_DESCRIPTOR_HANDLE		g_dsvHandle;
ID3D12CommandAllocator*			g_commandAllocator[2];
ID3D12GraphicsCommandList*		g_commandList[2];
ID3D12GraphicsCommandList*		g_currentCommandList;
ID3D12Fence*					g_fence;
HANDLE							g_fenceEvent;
UINT64							g_fenceValue = 0;

ID3D12Resource*					g_vertexBuffer[2];
D3D12_VERTEX_BUFFER_VIEW		g_vertexBufferView[2];
ID3D12Resource*					g_indexBuffer[2];
D3D12_INDEX_BUFFER_VIEW			g_indexBufferView[2];

ID3D12Resource*					g_imguiVertexBuffer[2];
D3D12_VERTEX_BUFFER_VIEW		g_imguiVertexBufferView[2];
ID3D12Resource*					g_imguiIndexBuffer[2];
D3D12_INDEX_BUFFER_VIEW			g_imguiIndexBufferView[2];

ID3D12Resource*					g_imguiFontTex;
ID3D12RootSignature*			g_rootSignature;
ID3D12PipelineState*			g_pso_Tri_Fill;
ID3D12PipelineState*			g_pso_Tri_Wire;
ID3D12PipelineState*			g_pso_Lines;
ID3D12PipelineState*			g_pso_Points;
ID3D12PipelineState*			g_pso_imgui;
ID3D12Resource*					g_cbCamera;
ID3D12Resource*					g_sbModelInstance;
bool							g_insideAppRender = false;
const uint32_t					g_maxVertexCount = 2*1024*1024;
const uint32_t					g_maxIndexCount = 2*1024*1024;
const uint32_t					g_maxModelInstanceCount = 2*1024;
uint32_t						g_currentVertexCount = 0;
uint32_t						g_currentIndexCount = 0;
uint32_t						g_currentModelInstanceCount = 0;
uint32_t						g_currentBufferIndex = 0;

struct CameraCB
{
	float4x4 viewProjMtx;
	float4x4 orthoProjMtx;
};

struct ModelInstanceCB
{
	float4x4 modelMtx;
	float4	 color;
};

void RenderImGui();

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
				//D3D12_MESSAGE_ID_RESOURCE_BARRIER_MISMATCHING_g_currentCommandList_TYPE,
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

		D3D12_DESCRIPTOR_HEAP_DESC srvHeapDesc = {};
		srvHeapDesc.NumDescriptors = 1;
		srvHeapDesc.Type	= D3D12_DESCRIPTOR_HEAP_TYPE_CBV_SRV_UAV;
		srvHeapDesc.Flags	= D3D12_DESCRIPTOR_HEAP_FLAG_SHADER_VISIBLE;
		g_device->CreateDescriptorHeap(&srvHeapDesc, IID_PPV_ARGS(&g_srvHeap));
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

		for (int i=0; i!=2; ++i)
		{
			resourceDesc.Width = sizeof(float3) * g_maxVertexCount;

			g_device->CreateCommittedResource(	&resourceHeapProps,
												D3D12_HEAP_FLAG_NONE,
												&resourceDesc,
												D3D12_RESOURCE_STATE_GENERIC_READ,
												nullptr,
												IID_PPV_ARGS(&g_vertexBuffer[i]) );

			resourceDesc.Width = sizeof(uint32_t) * g_maxIndexCount;

			g_device->CreateCommittedResource(	&resourceHeapProps,
												D3D12_HEAP_FLAG_NONE,
												&resourceDesc,
												D3D12_RESOURCE_STATE_GENERIC_READ,
												nullptr,
												IID_PPV_ARGS(&g_indexBuffer[i]) );

			g_vertexBufferView[i].BufferLocation	= g_vertexBuffer[i]->GetGPUVirtualAddress();
			g_vertexBufferView[i].SizeInBytes		= sizeof(float3) * g_maxVertexCount;
			g_vertexBufferView[i].StrideInBytes		= sizeof(float3);

			g_indexBufferView[i].BufferLocation		= g_indexBuffer[i]->GetGPUVirtualAddress();
			g_indexBufferView[i].SizeInBytes		= sizeof(uint32_t) * g_maxIndexCount;
			g_indexBufferView[i].Format				= DXGI_FORMAT_R32_UINT;
		}
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
		D3D12_ROOT_PARAMETER rootParameters[4];

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

        D3D12_DESCRIPTOR_RANGE descRange = {};
        descRange.RangeType = D3D12_DESCRIPTOR_RANGE_TYPE_SRV;
        descRange.NumDescriptors = 1;
        descRange.BaseShaderRegister = 1;
        descRange.RegisterSpace = 0;
        descRange.OffsetInDescriptorsFromTableStart = 0;

        rootParameters[3].ParameterType = D3D12_ROOT_PARAMETER_TYPE_DESCRIPTOR_TABLE;
        rootParameters[3].DescriptorTable.NumDescriptorRanges = 1;
        rootParameters[3].DescriptorTable.pDescriptorRanges = &descRange;
        rootParameters[3].ShaderVisibility = D3D12_SHADER_VISIBILITY_PIXEL;

		D3D12_STATIC_SAMPLER_DESC samplers[1];

		samplers[0].Filter = D3D12_FILTER_MIN_MAG_MIP_LINEAR;
		samplers[0].AddressU = D3D12_TEXTURE_ADDRESS_MODE_CLAMP;
		samplers[0].AddressV = D3D12_TEXTURE_ADDRESS_MODE_CLAMP;
		samplers[0].AddressW = D3D12_TEXTURE_ADDRESS_MODE_CLAMP;
		samplers[0].MipLODBias = 0.f;
		samplers[0].MaxAnisotropy = 0;
		samplers[0].ComparisonFunc = D3D12_COMPARISON_FUNC_ALWAYS;
		samplers[0].BorderColor = D3D12_STATIC_BORDER_COLOR_TRANSPARENT_BLACK;
		samplers[0].MinLOD = 0.f;
		samplers[0].MaxLOD = 0.f;
		samplers[0].ShaderRegister = 0;
		samplers[0].RegisterSpace = 0;
		samplers[0].ShaderVisibility = D3D12_SHADER_VISIBILITY_PIXEL;

		D3D12_ROOT_SIGNATURE_DESC rootSignatureDesc = {};
		rootSignatureDesc.NumParameters		= ARRAYSIZE(rootParameters);
		rootSignatureDesc.pParameters		= rootParameters;
		rootSignatureDesc.NumStaticSamplers = ARRAYSIZE(samplers);
		rootSignatureDesc.pStaticSamplers	= samplers;
		rootSignatureDesc.Flags				= D3D12_ROOT_SIGNATURE_FLAG_ALLOW_INPUT_ASSEMBLER_INPUT_LAYOUT;

		ComPtr<ID3DBlob> signature;
		ComPtr<ID3DBlob> errorBlob;
		D3D12SerializeRootSignature(&rootSignatureDesc, D3D_ROOT_SIGNATURE_VERSION_1, signature.GetAddressOf(), errorBlob.GetAddressOf());

		g_device->CreateRootSignature(0, signature->GetBufferPointer(), signature->GetBufferSize(), IID_PPV_ARGS(&g_rootSignature));
	}
   
	// Shaders & PSOs
	{
		const char* vertexShaderCode = R"(
			cbuffer Camera : register(b0)
			{
				float4x4 g_viewProjMtx;
				float4x4 g_orthoProjMtx;
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

		const char* pixelShaderCode = R"(
			cbuffer Camera : register(b0)
			{
				float4x4 g_viewProjMtx;
				float4x4 g_orthoProjMtx;
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

			float4 main(PSInput input, bool isFrontFace : SV_IsFrontFace) : SV_TARGET 
			{
				return g_ModelInstanceData[g_pushConstants.x].color * (isFrontFace ? 1 : 0.5);
			}
		)";

		ComPtr<ID3DBlob> errorBlob;
		ID3DBlob* vertexShaderBlob;
		ID3DBlob* pixelShaderBlob;

		V( D3DCompile(vertexShaderCode, strlen(vertexShaderCode), nullptr, nullptr, nullptr, "main", "vs_5_0", 0, 0, &vertexShaderBlob, errorBlob.GetAddressOf()) );
		const char* vsErrorStr = errorBlob ? (const char*)errorBlob->GetBufferPointer() : nullptr;

		V( D3DCompile(pixelShaderCode, strlen(pixelShaderCode), nullptr, nullptr, nullptr, "main", "ps_5_0", 0, 0, &pixelShaderBlob, errorBlob.GetAddressOf()) );
		const char* psErrorStr = errorBlob ? (const char*)errorBlob->GetBufferPointer() : nullptr;

		D3D12_GRAPHICS_PIPELINE_STATE_DESC defaultPSODesc = {};
	
		defaultPSODesc.pRootSignature			= g_rootSignature;
		defaultPSODesc.VS						= { vertexShaderBlob->GetBufferPointer(), vertexShaderBlob->GetBufferSize() };
		defaultPSODesc.PS						= { pixelShaderBlob->GetBufferPointer(), pixelShaderBlob->GetBufferSize() };
		defaultPSODesc.PrimitiveTopologyType	= D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;
		defaultPSODesc.RTVFormats[0]			= DXGI_FORMAT_B8G8R8A8_UNORM;
		defaultPSODesc.DSVFormat				= DXGI_FORMAT_D24_UNORM_S8_UINT;
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
		defaultPSODesc.RasterizerState.CullMode					= D3D12_CULL_MODE_NONE;
		//defaultPSODesc.RasterizerState.CullMode					= D3D12_CULL_MODE_BACK;
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
		psoDesc.RasterizerState.DepthBias = -10000;
		g_device->CreateGraphicsPipelineState(&psoDesc, IID_PPV_ARGS(&g_pso_Tri_Wire));

		psoDesc = defaultPSODesc;
		psoDesc.PrimitiveTopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_LINE;
		psoDesc.RasterizerState.DepthBias = -10000;
		g_device->CreateGraphicsPipelineState(&psoDesc, IID_PPV_ARGS(&g_pso_Lines));

		psoDesc = defaultPSODesc;
		psoDesc.PrimitiveTopologyType = D3D12_PRIMITIVE_TOPOLOGY_TYPE_POINT;
		psoDesc.RasterizerState.DepthBias = -10000;
		g_device->CreateGraphicsPipelineState(&psoDesc, IID_PPV_ARGS(&g_pso_Points));
	}

	// imgui
	{
		// shader
		{
			const char* vertexShaderCode = R"(
				cbuffer Camera : register(b0)
				{
					float4x4 g_viewProjMtx;
					float4x4 g_orthoProjMtx;
				};

				struct VSInput 
				{
					float3 position : POSITION;
					float4 col		: COLOR0;
					float2 uv		: TEXCOORD0;
				};

				struct PSInput 
				{
					float4 position : SV_POSITION;
					float4 col		: COLOR0;
					float2 uv		: TEXCOORD0;
				};

				PSInput main(VSInput input) 
				{
					PSInput output;

					float4 pos = float4(input.position, 1.0);

					output.position = mul(pos, g_orthoProjMtx);
					output.col = input.col;
					output.uv  = input.uv;
					return output;
				}
			)";

			const char* pixelShaderCode = R"(
				SamplerState sampler0 : register(s0);
				Texture2D texture0 : register(t1);

				struct PSInput 
				{
					float4 position : SV_POSITION;
					float4 col		: COLOR0;
					float2 uv		: TEXCOORD0;
				};

				float4 main(PSInput input) : SV_TARGET 
				{
					return input.col * texture0.Sample(sampler0, input.uv);
				}
			)";

			ComPtr<ID3DBlob> errorBlob;
			ID3DBlob* vertexShaderBlob;
			ID3DBlob* pixelShaderBlob;


			V( D3DCompile(vertexShaderCode, strlen(vertexShaderCode), nullptr, nullptr, nullptr, "main", "vs_5_0", 0, 0, &vertexShaderBlob, errorBlob.GetAddressOf()) );
			const char* vsErrorStr = errorBlob ? (const char*)errorBlob->GetBufferPointer() : nullptr;

			V( D3DCompile(pixelShaderCode, strlen(pixelShaderCode), nullptr, nullptr, nullptr, "main", "ps_5_0", 0, 0, &pixelShaderBlob, errorBlob.GetAddressOf()) );
			const char* psErrorStr = errorBlob ? (const char*)errorBlob->GetBufferPointer() : nullptr;

			D3D12_GRAPHICS_PIPELINE_STATE_DESC psoDesc = {};
	
			psoDesc.pRootSignature			= g_rootSignature;
			psoDesc.VS						= { vertexShaderBlob->GetBufferPointer(), vertexShaderBlob->GetBufferSize() };
			psoDesc.PS						= { pixelShaderBlob->GetBufferPointer(), pixelShaderBlob->GetBufferSize() };
			psoDesc.PrimitiveTopologyType	= D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;
			psoDesc.RTVFormats[0]			= DXGI_FORMAT_B8G8R8A8_UNORM;
			psoDesc.DSVFormat				= DXGI_FORMAT_D24_UNORM_S8_UINT;
			psoDesc.NumRenderTargets		= 1;
			psoDesc.SampleDesc.Count		= 1;
			psoDesc.SampleMask				= UINT_MAX;

			D3D12_INPUT_ELEMENT_DESC inputLayoutDesc[] =
			{
				{ "POSITION", 0, DXGI_FORMAT_R32G32_FLOAT,   0, (UINT)offsetof(ImDrawVert, pos), D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
				{ "TEXCOORD", 0, DXGI_FORMAT_R32G32_FLOAT,   0, (UINT)offsetof(ImDrawVert, uv),  D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
				{ "COLOR",    0, DXGI_FORMAT_R8G8B8A8_UNORM, 0, (UINT)offsetof(ImDrawVert, col), D3D12_INPUT_CLASSIFICATION_PER_VERTEX_DATA, 0 },
			};

			psoDesc.InputLayout.pInputElementDescs	= inputLayoutDesc;
			psoDesc.InputLayout.NumElements			= 3;

			psoDesc.BlendState.AlphaToCoverageEnable					= false;
			psoDesc.BlendState.IndependentBlendEnable					= FALSE;
			psoDesc.BlendState.RenderTarget[0].BlendEnable				= true;
			psoDesc.BlendState.RenderTarget[0].SrcBlend					= D3D12_BLEND_SRC_ALPHA;
			psoDesc.BlendState.RenderTarget[0].DestBlend				= D3D12_BLEND_INV_SRC_ALPHA;
			psoDesc.BlendState.RenderTarget[0].BlendOp					= D3D12_BLEND_OP_ADD;
			psoDesc.BlendState.RenderTarget[0].SrcBlendAlpha			= D3D12_BLEND_ONE;
			psoDesc.BlendState.RenderTarget[0].DestBlendAlpha			= D3D12_BLEND_INV_SRC_ALPHA;
			psoDesc.BlendState.RenderTarget[0].BlendOpAlpha				= D3D12_BLEND_OP_ADD;
			psoDesc.BlendState.RenderTarget[0].RenderTargetWriteMask	= D3D12_COLOR_WRITE_ENABLE_ALL;


			psoDesc.RasterizerState.FillMode				= D3D12_FILL_MODE_SOLID;
			psoDesc.RasterizerState.CullMode				= D3D12_CULL_MODE_NONE;
			psoDesc.RasterizerState.FrontCounterClockwise	= false;
			psoDesc.RasterizerState.DepthBias				= 0;
			psoDesc.RasterizerState.DepthBiasClamp			= 0;
			psoDesc.RasterizerState.SlopeScaledDepthBias	= 0;
			psoDesc.RasterizerState.DepthClipEnable			= TRUE;
			psoDesc.RasterizerState.MultisampleEnable		= FALSE;
			psoDesc.RasterizerState.AntialiasedLineEnable	= FALSE;
			psoDesc.RasterizerState.ForcedSampleCount		= 0;
			psoDesc.RasterizerState.ConservativeRaster		= D3D12_CONSERVATIVE_RASTERIZATION_MODE_OFF;

			psoDesc.DepthStencilState.DepthEnable			= FALSE;
			psoDesc.DepthStencilState.DepthWriteMask		= D3D12_DEPTH_WRITE_MASK_ALL;
			psoDesc.DepthStencilState.DepthFunc				= D3D12_COMPARISON_FUNC_ALWAYS;
			psoDesc.DepthStencilState.StencilEnable			= FALSE;

			g_device->CreateGraphicsPipelineState(&psoDesc, IID_PPV_ARGS(&g_pso_imgui));
		}

		// Vertex & Index buffers
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

			for (int i=0; i!=2; ++i)
			{
				resourceDesc.Width = sizeof(ImDrawVert) * g_maxVertexCount;

				g_device->CreateCommittedResource(	&resourceHeapProps,
													D3D12_HEAP_FLAG_NONE,
													&resourceDesc,
													D3D12_RESOURCE_STATE_GENERIC_READ,
													nullptr,
													IID_PPV_ARGS(&g_imguiVertexBuffer[i]) );

				resourceDesc.Width = sizeof(uint16_t) * g_maxIndexCount;

				g_device->CreateCommittedResource(	&resourceHeapProps,
													D3D12_HEAP_FLAG_NONE,
													&resourceDesc,
													D3D12_RESOURCE_STATE_GENERIC_READ,
													nullptr,
													IID_PPV_ARGS(&g_imguiIndexBuffer[i]) );

				g_imguiVertexBufferView[i].BufferLocation	= g_imguiVertexBuffer[i]->GetGPUVirtualAddress();
				g_imguiVertexBufferView[i].SizeInBytes		= sizeof(ImDrawVert) * g_maxVertexCount;
				g_imguiVertexBufferView[i].StrideInBytes	= sizeof(ImDrawVert);

				g_imguiIndexBufferView[i].BufferLocation	= g_imguiIndexBuffer[i]->GetGPUVirtualAddress();
				g_imguiIndexBufferView[i].SizeInBytes		= sizeof(uint16_t) * g_maxIndexCount;
				g_imguiIndexBufferView[i].Format			= DXGI_FORMAT_R16_UINT;
			}
		}

		// Font texture
		{
			ImGuiIO& io = ImGui::GetIO();
			unsigned char* pixels;
			int width, height;
			io.Fonts->GetTexDataAsRGBA32(&pixels, &width, &height);


			D3D12_RESOURCE_DESC fontTexDesc = {};
			fontTexDesc.Dimension			= D3D12_RESOURCE_DIMENSION_TEXTURE2D;
			fontTexDesc.Width				= width;
			fontTexDesc.Height				= height;
			fontTexDesc.DepthOrArraySize	= 1;
			fontTexDesc.MipLevels			= 1;
			fontTexDesc.Format				= DXGI_FORMAT_R8G8B8A8_UNORM;
			fontTexDesc.SampleDesc.Count	= 1;
			fontTexDesc.Layout				= D3D12_TEXTURE_LAYOUT_UNKNOWN;
			fontTexDesc.Flags				= D3D12_RESOURCE_FLAG_NONE;

			D3D12_HEAP_PROPERTIES resourceHeapProps;
			resourceHeapProps.Type					= D3D12_HEAP_TYPE_DEFAULT;
			resourceHeapProps.CPUPageProperty		= D3D12_CPU_PAGE_PROPERTY_UNKNOWN;
			resourceHeapProps.MemoryPoolPreference	= D3D12_MEMORY_POOL_UNKNOWN;
			resourceHeapProps.CreationNodeMask		= 0;
			resourceHeapProps.VisibleNodeMask		= 0;

			g_device->CreateCommittedResource(	&resourceHeapProps,
												D3D12_HEAP_FLAG_NONE,
												&fontTexDesc,
												D3D12_RESOURCE_STATE_COPY_DEST,
												nullptr,
												IID_PPV_ARGS(&g_imguiFontTex) );
		

			UINT upload_pitch = (width * 4 + D3D12_TEXTURE_DATA_PITCH_ALIGNMENT - 1u) & ~(D3D12_TEXTURE_DATA_PITCH_ALIGNMENT - 1u);
			UINT upload_size = height * upload_pitch;
			D3D12_RESOURCE_DESC uploadBufferDesc;
			uploadBufferDesc.Dimension			= D3D12_RESOURCE_DIMENSION_BUFFER;
			uploadBufferDesc.Alignment			= 0;
			uploadBufferDesc.Width				= upload_size;
			uploadBufferDesc.Height				= 1;
			uploadBufferDesc.DepthOrArraySize	= 1;
			uploadBufferDesc.MipLevels			= 1;
			uploadBufferDesc.Format				= DXGI_FORMAT_UNKNOWN;
			uploadBufferDesc.SampleDesc.Count	= 1;
			uploadBufferDesc.SampleDesc.Quality = 0;
			uploadBufferDesc.Layout				= D3D12_TEXTURE_LAYOUT_ROW_MAJOR;
			uploadBufferDesc.Flags				= D3D12_RESOURCE_FLAG_NONE;

			resourceHeapProps.Type					= D3D12_HEAP_TYPE_UPLOAD;
			resourceHeapProps.CPUPageProperty		= D3D12_CPU_PAGE_PROPERTY_UNKNOWN;
			resourceHeapProps.MemoryPoolPreference	= D3D12_MEMORY_POOL_UNKNOWN;

			ID3D12Resource* uploadBuffer = nullptr;

			g_device->CreateCommittedResource(	&resourceHeapProps,
												D3D12_HEAP_FLAG_NONE,
												&uploadBufferDesc,
												D3D12_RESOURCE_STATE_GENERIC_READ,
												nullptr,
												IID_PPV_ARGS(&uploadBuffer) );

			void* pData = nullptr;
			uploadBuffer->Map(0, 0, &pData);
			{
				for (int y = 0; y < height; y++)
				{
					memcpy((void*) ((uintptr_t) pData + y * upload_pitch), pixels + y * width * 4, width * 4);
				}
			}
			uploadBuffer->Unmap(0, 0);

			D3D12_TEXTURE_COPY_LOCATION srcLocation = {};
			D3D12_TEXTURE_COPY_LOCATION dstLocation = {};
			{
				srcLocation.pResource							= uploadBuffer;
				srcLocation.Type								= D3D12_TEXTURE_COPY_TYPE_PLACED_FOOTPRINT;
				srcLocation.PlacedFootprint.Footprint.Format	= DXGI_FORMAT_R8G8B8A8_UNORM;
				srcLocation.PlacedFootprint.Footprint.Width		= width;
				srcLocation.PlacedFootprint.Footprint.Height	= height;
				srcLocation.PlacedFootprint.Footprint.Depth		= 1;
				srcLocation.PlacedFootprint.Footprint.RowPitch	= upload_pitch;

				dstLocation.pResource							= g_imguiFontTex;
				dstLocation.Type								= D3D12_TEXTURE_COPY_TYPE_SUBRESOURCE_INDEX;
				dstLocation.SubresourceIndex					= 0;
			}

			D3D12_RESOURCE_BARRIER barrier = {};
			barrier.Type					= D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
			barrier.Flags					= D3D12_RESOURCE_BARRIER_FLAG_NONE;
			barrier.Transition.pResource   = g_imguiFontTex;
			barrier.Transition.Subresource = D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES;
			barrier.Transition.StateBefore = D3D12_RESOURCE_STATE_COPY_DEST;
			barrier.Transition.StateAfter  = D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE;

			g_commandAllocator[0]->Reset();
			g_commandList[0]->Reset(g_commandAllocator[0], nullptr);
			g_commandList[0]->CopyTextureRegion(&dstLocation, 0, 0, 0, &srcLocation, nullptr);
			g_commandList[0]->ResourceBarrier(1, &barrier);
			g_commandList[0]->Close();

			ID3D12CommandList* listsToExecute[] = { g_commandList[0] };
			g_commandQueue->ExecuteCommandLists(_countof(listsToExecute), listsToExecute);

			g_commandQueue->Signal(g_fence, 1);

			g_fence->SetEventOnCompletion(1, g_fenceEvent);
			WaitForSingleObject(g_fenceEvent, INFINITE);

			uploadBuffer->Release();

			D3D12_SHADER_RESOURCE_VIEW_DESC srvDesc;
			ZeroMemory(&srvDesc, sizeof(srvDesc));
			srvDesc.Format						= DXGI_FORMAT_R8G8B8A8_UNORM;
			srvDesc.ViewDimension				= D3D12_SRV_DIMENSION_TEXTURE2D;
			srvDesc.Texture2D.MipLevels			= fontTexDesc.MipLevels;
			srvDesc.Texture2D.MostDetailedMip	= 0;
			srvDesc.Shader4ComponentMapping		= D3D12_DEFAULT_SHADER_4_COMPONENT_MAPPING;
			g_device->CreateShaderResourceView(g_imguiFontTex, &srvDesc, g_srvHeap->GetCPUDescriptorHandleForHeapStart());

			io.Fonts->SetTexID((ImTextureID)g_srvHeap->GetGPUDescriptorHandleForHeapStart().ptr);
		}
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
		depthStencilDesc.Format				= DXGI_FORMAT_D24_UNORM_S8_UINT;
		depthStencilDesc.SampleDesc.Count	= 1;
		depthStencilDesc.Layout				= D3D12_TEXTURE_LAYOUT_UNKNOWN;
		depthStencilDesc.Flags				= D3D12_RESOURCE_FLAG_ALLOW_DEPTH_STENCIL;

		D3D12_CLEAR_VALUE clearValue = {};
		clearValue.Format = DXGI_FORMAT_D24_UNORM_S8_UINT;
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

		dsvDesc.Format				= DXGI_FORMAT_D24_UNORM_S8_UINT;
		dsvDesc.ViewDimension		= D3D12_DSV_DIMENSION_TEXTURE2D;
		dsvDesc.Flags				= D3D12_DSV_FLAG_NONE;
		dsvDesc.Texture2D.MipSlice	= 0;

		g_dsvHandle = g_dsvHeap->GetCPUDescriptorHandleForHeapStart();
		g_device->CreateDepthStencilView(g_depthStencilBuffer, &dsvDesc, g_dsvHandle);
	}
}

void RenderD3D(uint32_t windowWidth, uint32_t windowHeight) 
{
	++g_fenceValue;

	if (g_fenceValue>=2 && g_fence->GetCompletedValue() < g_fenceValue-2) 
	{
		g_fence->SetEventOnCompletion(g_fenceValue-2, g_fenceEvent);
		WaitForSingleObject(g_fenceEvent, INFINITE);
	}

	//WaitForSingleObjectEx(g_swapChainWaitObject, 1000, true);

	g_currentBufferIndex = g_swapChain->GetCurrentBackBufferIndex();
	assert(g_currentBufferIndex<2);

	g_currentVertexCount = 0;
	g_currentIndexCount = 0;
	g_currentModelInstanceCount = 0;

	g_commandAllocator[g_currentBufferIndex]->Reset();
	g_commandList[g_currentBufferIndex]->Reset(g_commandAllocator[g_currentBufferIndex], nullptr);
	g_currentCommandList = g_commandList[g_currentBufferIndex];

	D3D12_RESOURCE_BARRIER barrier = {};
	barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
	barrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
	barrier.Transition.pResource = g_backBuffer[g_currentBufferIndex];
	barrier.Transition.Subresource = 0;
	barrier.Transition.StateBefore = D3D12_RESOURCE_STATE_PRESENT;
	barrier.Transition.StateAfter = D3D12_RESOURCE_STATE_RENDER_TARGET;

	g_currentCommandList->ResourceBarrier(1, &barrier);

	g_currentCommandList->SetDescriptorHeaps(1, &g_srvHeap);

	g_currentCommandList->OMSetRenderTargets(1, &g_rtvHandle[g_currentBufferIndex], FALSE, &g_dsvHandle);


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
	g_currentCommandList->ClearRenderTargetView(g_rtvHandle[g_currentBufferIndex], clearColor, 0, nullptr);
	g_currentCommandList->ClearDepthStencilView(g_dsvHandle, D3D12_CLEAR_FLAG_DEPTH, 1.0f, 0, 0, nullptr);

	void* pData;
	g_cbCamera->Map(0, nullptr, &pData);
	{
		CameraCB* cb = (CameraCB*)pData;

		float4x4 viewMtx = g_cameraController.GetViewMatrix();
		float4x4 projMtx = g_cameraController.GetProjectionMtx();

		cb->viewProjMtx = mul(viewMtx, projMtx);

		cb->orthoProjMtx = matrix44::MakeOrtho(0, (float)windowWidth, 0, (float)windowHeight, 0, 1);
	}
	g_cbCamera->Unmap(0, nullptr);
	
	g_currentCommandList->SetGraphicsRootSignature(g_rootSignature);
    g_currentCommandList->IASetVertexBuffers(0, 1, &g_vertexBufferView[g_currentBufferIndex]);
	g_currentCommandList->IASetIndexBuffer(&g_indexBufferView[g_currentBufferIndex]);
	g_currentCommandList->SetGraphicsRootConstantBufferView(0, g_cbCamera->GetGPUVirtualAddress());

	ImGui::NewFrame();

	// Render app content
	{
		g_insideAppRender = true;

		g_currentCommandList->SetGraphicsRootShaderResourceView(1, g_sbModelInstance->GetGPUVirtualAddress() );

		AppRender(windowWidth, windowHeight);

		{
			float3 axisPoints[] = 
			{
				float3(-1,0,0), float3(1,0,0),
				float3(0,-1,0), float3(0,1,0),
				float3(0,0,-1), float3(0,0,1),
			};

			DrawLines(axisPoints, 6, matrix44::MakeTranslation(g_cameraController.GetTarget().xyz), float4(0.5f, 0.5f, 0.5f, 1));
		}

		g_insideAppRender = false;
	}

	RenderImGui();

	barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
	barrier.Flags = D3D12_RESOURCE_BARRIER_FLAG_NONE;
	barrier.Transition.pResource = g_backBuffer[g_currentBufferIndex];
	barrier.Transition.Subresource = 0;
	barrier.Transition.StateBefore = D3D12_RESOURCE_STATE_RENDER_TARGET;
	barrier.Transition.StateAfter = D3D12_RESOURCE_STATE_PRESENT;

	g_currentCommandList->ResourceBarrier(1, &barrier);

	g_currentCommandList->Close();

	ID3D12CommandList* listsToExecute[] = { g_currentCommandList };
	g_commandQueue->ExecuteCommandLists(_countof(listsToExecute), listsToExecute);

	g_swapChain->Present(0, DXGI_PRESENT_ALLOW_TEARING);

	g_commandQueue->Signal(g_fence, g_fenceValue);
}

void RenderImGui()
{
	ImGui::Render();

	ImDrawData* draw_data = ImGui::GetDrawData();

	// Avoid rendering when minimized
	if (draw_data->DisplaySize.x <= 0.0f || draw_data->DisplaySize.y <= 0.0f)
		return;

	assert(draw_data->TotalVtxCount < g_maxVertexCount);
	assert(draw_data->TotalIdxCount < g_maxIndexCount);

	// Upload vertex/index data into a single contiguous GPU buffer
	{
		void* pVertexData;
		void* pIndexData;

		g_imguiVertexBuffer[g_currentBufferIndex]->Map(0, nullptr, &pVertexData);
		g_imguiIndexBuffer[g_currentBufferIndex]->Map(0, nullptr, &pIndexData);
		{
			ImDrawVert* vtx_dst = (ImDrawVert*)pVertexData;
			ImDrawIdx* idx_dst = (ImDrawIdx*)pIndexData;
		
			for (int n = 0; n < draw_data->CmdListsCount; n++)
			{
				const ImDrawList* draw_list = draw_data->CmdLists[n];

				memcpy(vtx_dst, draw_list->VtxBuffer.Data, draw_list->VtxBuffer.Size * sizeof(ImDrawVert));
				memcpy(idx_dst, draw_list->IdxBuffer.Data, draw_list->IdxBuffer.Size * sizeof(ImDrawIdx));

				vtx_dst += draw_list->VtxBuffer.Size;
				idx_dst += draw_list->IdxBuffer.Size;
			}
		}
		g_imguiVertexBuffer[g_currentBufferIndex]->Unmap(0, nullptr);
		g_imguiIndexBuffer[g_currentBufferIndex]->Unmap(0, nullptr);
	}

	// Bind state
	{
		g_currentCommandList->IASetVertexBuffers(0, 1, &g_imguiVertexBufferView[g_currentBufferIndex]);
		g_currentCommandList->IASetIndexBuffer(&g_imguiIndexBufferView[g_currentBufferIndex]);
		g_currentCommandList->IASetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
		g_currentCommandList->SetPipelineState(g_pso_imgui);
	}

	// Render command lists
	{
		int global_vtx_offset = 0;
		int global_idx_offset = 0;

		ImVec2 clip_off = draw_data->DisplayPos;

		for (int n = 0; n < draw_data->CmdListsCount; n++)
		{
			const ImDrawList* draw_list = draw_data->CmdLists[n];

			for (int cmd_i = 0; cmd_i < draw_list->CmdBuffer.Size; cmd_i++)
			{
				const ImDrawCmd* pcmd = &draw_list->CmdBuffer[cmd_i];
				if (pcmd->UserCallback != nullptr)
				{
					pcmd->UserCallback(draw_list, pcmd);
				}
				else
				{
					// Project scissor/clipping rectangles into framebuffer space
					ImVec2 clip_min(pcmd->ClipRect.x - clip_off.x, pcmd->ClipRect.y - clip_off.y);
					ImVec2 clip_max(pcmd->ClipRect.z - clip_off.x, pcmd->ClipRect.w - clip_off.y);

					if (clip_max.x <= clip_min.x || clip_max.y <= clip_min.y)
						continue;

					// Apply scissor/clipping rectangle
					const D3D12_RECT r = { (LONG)clip_min.x, (LONG)clip_min.y, (LONG)clip_max.x, (LONG)clip_max.y };
					
					g_currentCommandList->RSSetScissorRects(1, &r);

					// Bind texture, Draw
					D3D12_GPU_DESCRIPTOR_HANDLE texture_handle = {};
					texture_handle.ptr = (UINT64)pcmd->GetTexID();
					g_currentCommandList->SetGraphicsRootDescriptorTable(3, texture_handle);
					g_currentCommandList->DrawIndexedInstanced(pcmd->ElemCount, 1, pcmd->IdxOffset + global_idx_offset, pcmd->VtxOffset + global_vtx_offset, 0);
				}
			}

			global_idx_offset += draw_list->IdxBuffer.Size;
			global_vtx_offset += draw_list->VtxBuffer.Size;
		}
	}
}

void DrawMesh(const float3* positions, uint32_t pointCount, uint32_t* indices, uint32_t indexCount, const matrix44& worldMtx, float4 fillColor, float4 wireColor)
{
	bool needsWire = wireColor != float4::k0000;

	assert(g_insideAppRender);
	assert( g_currentVertexCount + pointCount <= g_maxVertexCount );
	assert( g_currentIndexCount + indexCount <= g_maxIndexCount );
	assert( g_currentModelInstanceCount + needsWire ? 2 : 1 <= g_maxModelInstanceCount );

	void* pData;
	
	g_vertexBuffer[g_currentBufferIndex]->Map(0, nullptr, &pData);
	{
		memcpy( (float3*)pData + g_currentVertexCount, positions, pointCount * sizeof(float3));
	}
	g_vertexBuffer[g_currentBufferIndex]->Unmap(0, nullptr);

	g_indexBuffer[g_currentBufferIndex]->Map(0, nullptr, &pData);
	{
		memcpy( (uint32_t*)pData + g_currentIndexCount, indices, indexCount * sizeof(uint32_t));
	}
	g_indexBuffer[g_currentBufferIndex]->Unmap(0, nullptr);

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
    g_currentCommandList->DrawIndexedInstanced( (uint32_t)indexCount, 1, g_currentIndexCount, g_currentVertexCount, 0);

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

		g_currentCommandList->DrawIndexedInstanced( (uint32_t)indexCount, 1, g_currentIndexCount, g_currentVertexCount, 0);
	}

	g_currentVertexCount += (uint32_t)pointCount;
	g_currentIndexCount += (uint32_t)indexCount;
	g_currentModelInstanceCount += 1;
}

void DrawLines(const float3* positions, uint32_t pointCount, const matrix44& worldMtx, float4 color)
{
	assert(g_insideAppRender);
	assert( g_currentVertexCount + pointCount <= g_maxVertexCount );
	assert( g_currentModelInstanceCount + 1 <= g_maxModelInstanceCount );

	void* pData;
	
	g_vertexBuffer[g_currentBufferIndex]->Map(0, nullptr, &pData);
	{
		memcpy( (float3*)pData + g_currentVertexCount, positions, pointCount * sizeof(float3));
	}
	g_vertexBuffer[g_currentBufferIndex]->Unmap(0, nullptr);

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

	g_currentCommandList->DrawInstanced( (uint32_t)pointCount, 1, g_currentVertexCount, 0);

	g_currentVertexCount += (uint32_t)pointCount;
	g_currentModelInstanceCount += 1;
}

void DrawPoints(const float3* positions, uint32_t pointCount, const matrix44& worldMtx, float4 color)
{
	assert(g_insideAppRender);
	assert( g_currentVertexCount + pointCount <= g_maxVertexCount );
	assert( g_currentModelInstanceCount + 1 <= g_maxModelInstanceCount );

	void* pData;
	
	g_vertexBuffer[g_currentBufferIndex]->Map(0, nullptr, &pData);
	{
		memcpy( (float3*)pData + g_currentVertexCount, positions, pointCount * sizeof(float3));
	}
	g_vertexBuffer[g_currentBufferIndex]->Unmap(0, nullptr);

	g_sbModelInstance->Map(0, nullptr, &pData);
	{
		ModelInstanceCB* cb = (ModelInstanceCB*)pData + g_currentModelInstanceCount;

		cb->modelMtx = worldMtx;
		cb->color = color;
	}
	g_sbModelInstance->Unmap(0, nullptr);

	g_currentCommandList->SetPipelineState(g_pso_Points);

	g_currentCommandList->IASetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_POINTLIST);

	uint4 instanceData(g_currentModelInstanceCount, 0, 0, 0);
	g_currentCommandList->SetGraphicsRoot32BitConstants(2, 4, &instanceData, 0);

	g_currentCommandList->DrawInstanced( (uint32_t)pointCount, 1, g_currentVertexCount, 0);

	g_currentVertexCount += (uint32_t)pointCount;
	g_currentModelInstanceCount += 1;
}