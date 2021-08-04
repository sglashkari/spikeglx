//****************************************************************************************************************************************************************************************
//	NetComClient.h
//****************************************************************************************************************************************************************************************
#pragma once
#include "Nlx_DataTypes.h"
#include <vector>
#include <string>

//These forward declares allow us to only distribute this
//header file, instead of the entire header tree.  The
//necessary header files are included in the cpp file
class NetworkNodeClient;

namespace NetCom {


	class NetComClient
	{
	public:
		__declspec(dllexport) NetComClient(void);
		__declspec(dllexport) virtual ~NetComClient(void);

		


		//client functions
		__declspec(dllexport) bool ConnectToServer(const wchar_t* const serverName, bool attemptRouterConnection = true);
		__declspec(dllexport) bool DisconnectFromServer();
		__declspec(dllexport) std::wstring GetClientVersionString() 
		{
			std::wstring clientVersionString(L"");
			wchar_t* pClientVersionString = NULL;

			GetClientVersionStringInternal(pClientVersionString);

			clientVersionString = pClientVersionString;

			FreeArray(pClientVersionString);

			return clientVersionString;
		}
		__declspec(dllexport) bool OpenStream(const wchar_t* const cheetahObjectName);
		__declspec(dllexport) bool CloseStream(const wchar_t* const cheetahObjectName);
		__declspec(dllexport) bool SendCommand(const wchar_t* const command, char*& reply, int& numBytesReturned);

		__declspec(dllexport) void SetCallbackFunctionSE( void (*myFunctionPtr)(void* myClassPtr, SERec* records, int numRecords, const wchar_t* const objectName), void* myClassPtr);
		__declspec(dllexport) void SetCallbackFunctionST( void (*myFunctionPtr)(void* myClassPtr, STRec* records, int numRecords, const wchar_t* const objectName), void* myClassPtr);
		__declspec(dllexport) void SetCallbackFunctionTT( void (*myFunctionPtr)(void* myClassPtr, TTRec* records, int numRecords, const wchar_t* const objectName), void* myClassPtr);
		__declspec(dllexport) void SetCallbackFunctionCSC( void (*myFunctionPtr)(void* myClassPtr, CRRec* records, int numRecords, const wchar_t objectName []), void* myClassPtr);
		__declspec(dllexport) void SetCallbackFunctionEV( void (*myFunctionPtr)(void* myClassPtr, EventRec* records, int numRecords, const wchar_t* const objectName), void* myClassPtr);
		__declspec(dllexport) void SetCallbackFunctionVT( void (*myFunctionPtr)(void* myClassPtr, VideoRec* records, int numRecords, const wchar_t* const objectName), void* myClassPtr);
		//setter functions
		__declspec(dllexport) bool SetApplicationName(const wchar_t* const myApplicationName);
		__declspec(dllexport) bool SetLogFileName(const wchar_t* const filename);

		//getter functions
		__declspec(dllexport) bool GetRingBufferSize(const wchar_t* const objectName, int& numRecords, int& recordSize);
		__declspec(dllexport) bool GetRingBufferData(const wchar_t* const objectName, char* dataBuffer, int numBytes, int& numBytesReturned);
		__declspec(dllexport) bool GetSamplingFrequency(const wchar_t* const objectName, int& samplingFrequency);
		__declspec(dllexport) bool LogMessageCheetah(const wchar_t* const message, int messageType = 0);
		__declspec(dllexport) bool GetADRange(const wchar_t* const objectName, int& maxVal, int& minVal );

		//getter status functions

		__declspec(dllexport) bool AreWeConnected();

		//C++ STL wrappers for C API functions
		//This saves having to worry about memory management and the DLL
		//boundary in an application linking against this library.
		//Since STL/CRT objects cannot pass the DLL boundary (the STL/CRT
		//versions may differ), these functions must be implemented in
		//the header so that the application that links to this library
		//compiles them using its version of the STL/CRT.

		//*******************************************************************
		bool GetCheetahObjectsAndTypes(std::vector<std::wstring>& cheetahObjects, std::vector<std::wstring>& cheetahTypes)
		{
			cheetahObjects.clear();
			cheetahTypes.clear();

			wchar_t** pCheetahObjects = NULL;
			wchar_t** pCheetahTypes = NULL;
			int numObjectsReturned = -1;

			if(GetCheetahObjectsAndTypesInternal(pCheetahObjects, pCheetahTypes, numObjectsReturned)) {
				for(int i = 0; i < numObjectsReturned; ++i) {
					cheetahObjects.push_back(std::wstring(pCheetahObjects[i]));
					FreeArray((void*)(pCheetahObjects[i]));
					cheetahTypes.push_back(std::wstring(pCheetahTypes[i]));
					FreeArray((void*)(pCheetahTypes[i]));
				}

				FreeArray(pCheetahObjects);
				FreeArray(pCheetahTypes);
				return true;
			} else {
				return false;
			}
		}

		//*******************************************************************
		bool SendCommand(const wchar_t command [], std::wstring& reply)
		{
			wchar_t* pReply = NULL;
			reply = L"";

			if(SendCommandInternal(command, pReply) == true) {
				reply = pReply;
				FreeArray(pReply);
				return true;
			} else {
				return false;
			}
		}

		//*****************************************************************************
		bool SendCommand (const wchar_t* const command, int& commandSucceeded, std::vector<std::wstring>& replyValues)
		{
			wchar_t** pReply = NULL;
			int numReplyStringsReturned = 0;
			commandSucceeded = -1;
			replyValues.clear();

			if(SendCommandInternal(command, commandSucceeded, pReply, numReplyStringsReturned)) {
				for(int i = 0; i < numReplyStringsReturned; ++i) {
					replyValues.push_back(std::wstring(pReply[i]));
					FreeArray((void*)(pReply[i]));
				}

				FreeArray(pReply);
				return true;
			} else {
				return false;
			}

		}

		//*******************************************************************
		std::wstring GetServerPCName(void)
		{
			std::wstring pcName(L"");
			wchar_t* pPCName = NULL;

			GetServerPCNameInternal(pPCName);

			pcName = pPCName;

			FreeArray(pPCName);

			return pcName;
		}

		//*******************************************************************
		std::wstring GetServerIPAddress(void)
		{
			std::wstring ipAddress(L"");
			wchar_t* pIPAddress = NULL;

			GetServerIPAddressInternal(pIPAddress);

			ipAddress = pIPAddress;

			FreeArray(pIPAddress);

			return ipAddress;
		}

		//*******************************************************************
		std::wstring GetServerApplicationName(void)
		{
			std::wstring appName(L"");
			wchar_t* pAppName = NULL;

			GetServerApplicationNameInternal(pAppName);

			appName = pAppName;

			FreeArray(pAppName);

			return appName;
		}

		//*******************************************************************
		bool GetFeature(const wchar_t* const objectName, const int featureNumber, std::wstring& feature)
		{
			feature = L"";
			wchar_t* pFeature = NULL;

			if(GetFeatureInternal(objectName, featureNumber, pFeature) == true) {
				feature = pFeature;
				FreeArray(pFeature);
				return true;
			}else {
				return false;
			}
		}

	private:

		bool SetComputerName();

		//C API general functions
		__declspec(dllexport) void GetClientVersionStringInternal(wchar_t*& clientVersionString);

		//C API cheetah AE functions
		__declspec(dllexport) bool GetCheetahObjectsAndTypesInternal(wchar_t**& cheetahObjects, wchar_t**& cheetahTypes, int& numObjectsReturned);
		__declspec(dllexport) bool GetFeatureInternal(const wchar_t* const objectName, const int featureNumber, wchar_t*& feature);

		//C Command functions
		__declspec(dllexport) bool SendCommandInternal(const wchar_t command [], wchar_t*& reply);
		__declspec(dllexport) bool SendCommandInternal(const wchar_t* const command, int& commandSucceeded, wchar_t**& pReplyStrings, int& numReplyStringsReturned);

		//C API server information functions
		__declspec(dllexport) void GetServerApplicationNameInternal(wchar_t*& appName);
		__declspec(dllexport) void GetServerPCNameInternal(wchar_t*& pcName);
		__declspec(dllexport) void GetServerIPAddressInternal(wchar_t*& ipAddress);

		//C API helper functions
		__declspec(dllexport) void FreeArray(void* arrayPtr);


		NetworkNodeClient* mNetworkConnection; //implementation object
	};

	//NetComClient factory methods so that we can allocate and deallocate
	//pointers to NetComClient objects on the DLL's heap. If a NetComClient*
	//is newed in some application, the memory for the NetComClient is in
	//the heap of that other application. However, anything newed by NetComClient
	//is in the heap of this DLL. So when ~NetComClient is called from the other
	//application, we get a heap corruption because we can't delete things in
	//this DLL's heap from the other application.
	//TODO: it might be a good idea to make the NetComClient class a factory with
	//private constructor.
	__declspec(dllexport) NetComClient* GetNewNetComClient(void);
	__declspec(dllexport) void DeleteNetComClient(NetComClient*);


};
