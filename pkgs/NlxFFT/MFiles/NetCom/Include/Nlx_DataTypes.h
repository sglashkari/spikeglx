// Created 1/8/01 Steve Franks
// This file 'works' with no modification under MS Visual C++ 6.0

// avoid duplicate definitions
#ifndef _NLX_DATATYPE_H
#define _NLX_DATATYPE_H

//This forces all records to be byte-packed for most efficient use of disk space & lack of confusion
#pragma pack(push, before_nlx_datatypes)
#pragma pack(1)

#ifndef Nlx_Odl
	#define Nlx_Odl
#endif


#include <string>

#define NlxRetVal int
const int NlxOK = 0;
const int NlxError = -1;
const int NlxNotYetInitialized = -1;


//hwss types need to be stored someone, not sure on this yet, talk with bri



enum NLX_SPIKE_SUB_CHANNEL_COUNT { NSSCC_UNKNOWN = 0, NSSCC_SE = 1, NSSCC_ST = 2, NSSCC_TT = 4 };


//this does not fit with the other data types that have been setup as it is a general type for all spikes
//this is probably what we want to move to in the future.  This is however the current header type so it must appear somewhere.
const wchar_t NlxFileHeaderDataTypeSpike [] = {L"Spike"};

const int DATAFILETYPE_INVALID = -1;
const int DATAFILETYPE_RESERVED = 0;
const int DATAFILETYPE_SESPIKE = 1;
const int DATAFILETYPE_STSPIKE = 2;
const int DATAFILETYPE_TIMESTAMP = 3;
const int DATAFILETYPE_TTSPIKE = 4;
const int DATAFILETYPE_CSC = 5;
const int DATAFILETYPE_EVENT = 6;
const int DATAFILETYPE_VIDEO = 7;
/**/

const int MAX_NUMELECTRODES = 4; //maximum number of electrodes in a spike AE

enum NLX_DATA_TYPE { NDT_INVALID, NDT_SINGLE_ELECTRODE, NDT_STEREOTRODE, NDT_TIMESTAMP, NDT_TETRODE, 
					 NDT_CSC, NDT_EVENT, NDT_VIDEO, NDT_RAW, NDT_MCLUST_TS, NDT_MEF };

//****************************************************************************************************************************
//****************************************************************************************************************************
static std::wstring GetDataTypeLabel(NLX_DATA_TYPE dataType)
{
	switch (dataType)
	{
	case NDT_INVALID:
		return(L"Invalid");
	case NDT_SINGLE_ELECTRODE:
		return(L"SingleElectrode");
	case NDT_STEREOTRODE:
		return(L"Stereotrode");
	case NDT_TIMESTAMP:
		return(L"Timestamp");
	case NDT_TETRODE:
		return(L"Tetrode");
	case NDT_CSC:
		return(L"CSC");
	case NDT_EVENT:
		return(L"Event");
	case NDT_VIDEO:
		return(L"Video");
	case NDT_RAW:
		return(L"Raw");
	case NDT_MCLUST_TS:
		return(L"MClustTS");
	case NDT_MEF:
		return(L"MEF");
	}
	return(L"Invalid");
}

//****************************************************************************************************************************
//****************************************************************************************************************************
static std::wstring GetFileExtension(NLX_DATA_TYPE dataType)
{
	switch (dataType)
	{
	case NDT_INVALID:
		return(L".error");
	case NDT_SINGLE_ELECTRODE:
		return(L".nse");
	case NDT_STEREOTRODE:
		return(L".nst");
	case NDT_TIMESTAMP:
		return(L".nts");
	case NDT_TETRODE:
		return(L".ntt");
	case NDT_CSC:
		return(L".ncs");
	case NDT_EVENT:
		return(L".nev");
	case NDT_VIDEO:
		return(L".nvt");
	case NDT_RAW:
		return(L".nrd");
	case NDT_MCLUST_TS:
		return(L".t");
	case NDT_MEF:
		return(L".mef");
	}
	return(L".error");
}



//////////////////
// AD(Raw) Datatype
//////////////////
//******************************************************************************************************************************************************************************************
//  AD(Raw) Record Format:
//  STX (or SOP)    2048
//  Packet ID       1
//  Size            0x0000002A   hex for (# A/D data wds + #extra wds = 32+10)
//  TimeStamp Hi    1 32 bit word
//  TimeStamp Low   1 32 bit word
//  CPU_Status_wd   1 32 bit word
//  Parallel_in     1 32 bit word
//  10 extras      10 32 bit words
//  A/D data       32 32 bit words
//  CRC             1 32 bit XOR of the entire packet including STX
//******************************************************************************************************************************************************************************************
const unsigned __int16 ADRecordStx = 2048;
const unsigned __int16 ADRecordID = 1;
const int ADRecordNumFieldsNoData = 18;

//since ADRecs are interleaved, this struct is used to simulate individual AD Records in an array
//similar to a TTRec*.  It should only be used where arrays of records are expected, and should never
//be part of an array itself.
struct ADRec {

	//constructor 
	ADRec(unsigned __int64* timestampArray, __int32* dataArray, int index = 0):
	qwTimeStamp(-1),
	snSample(-1)
	{
		//the timestamps and samples are contained in separate arrays because of AD interleaving
		mTimestampPtr = timestampArray;
		mDataPointer = dataArray;
		mIndex = index;
		
		//set the public members to be the array's value at the current index.
		if(mTimestampPtr != NULL) {
			qwTimeStamp = mTimestampPtr[mIndex];
		}
		if(mDataPointer != NULL) {
			snSample = mDataPointer[mIndex];
		}
	};

	//generic constructor, you can't actually use one of these structs for
	//anything useful other than assignment.
	ADRec():
	mTimestampPtr(NULL),
	mDataPointer(NULL),
	mIndex(0),
	qwTimeStamp(-1),
	snSample(-1)
	{

	};

	//emulate an array index off of this object.  This will not
	//work for ADRec* since the * is actually part of the core
	//language and cannot be overridden.  Using [] notation off of
	//an ADRec* will yeild unpredictable results
	ADRec operator[] (unsigned int index) {
		return ADRec(mTimestampPtr, mDataPointer, index);
	}

	unsigned __int64 qwTimeStamp; //the timestamp at the current index
	__int32 snSample; //the sample at the current index
private:
	unsigned __int64* mTimestampPtr; //internal pointer to the timestamp array
	__int32* mDataPointer; //internal pointer to the samples array
	int mIndex; //current offset for the timestamp and sample pointers
};


//these are used for NetCom, they were also used in cheetah 4 to denote filetype in a data files header
//currently used strings include (CSC, Spike, Event, Nothing for Video or TS)
#ifdef _UNICODE
	#define NetComSEDataType L"SEScAcqEnt"
	#define NetComSTDataType L"STScAcqEnt"
	#define NetComTTDataType L"TTScAcqEnt"
	#define NetComCSCDataType L"CscAcqEnt"
	#define NetComEventDataType L"EventAcqEnt"
	#define NetComVTDataType L"VTAcqEnt"
#else
	#define NetComSEDataType "SEScAcqEnt"
	#define NetComSTDataType "STScAcqEnt"
	#define NetComTTDataType "TTScAcqEnt"
	#define NetComCSCDataType "CscAcqEnt"
	#define NetComEventDataType "EventAcqEnt"
	#define NetComVTDataType "VTAcqEnt"
#endif




//////////////////
// Spike Datatypes
//////////////////


const int TT_NUMELECTRODES = 4;
const int ST_NUMELECTRODES = 2; 
const int SE_NUMELECTRODES = 1; 
const int MAX_PARAMS = 8;
const int SPIKE_NUMPOINTS = 32;

//IMPORTANT!  If you are getting a compile error on Nlx_Odl, you just need to add the line:
//IMPORTANT!																				#define Nlx_Odl
//IMPORTANT!		right before you '#include' this file.  This will make Nlx_Odl
//IMPORTANT!		equate out to nothing in your compiler, which will fix the error.
//IMPORTANT!		(Nlx_Odl is used internally by Neuralynx to generate additional compile information.)

// used only as a part of larger record
Nlx_Odl struct TetPoint {
	signed __int16		snADVal[TT_NUMELECTRODES];
};

// used only as a part of larger record
Nlx_Odl struct StereoPoint {
	signed __int16		snADVal[ST_NUMELECTRODES];
};

// used only as a part of larger record
Nlx_Odl struct SinglePoint {
	signed __int16		snADVal[SE_NUMELECTRODES];
};


// a tetrode (TTScAcqEnt record)
Nlx_Odl struct TTRec	{
		unsigned __int64	qwTimeStamp;			// TS
		unsigned __int32	dwScNumber;			// Channel number
		unsigned __int32	dwCellNumber;			// What cell was this calculated to be? filled in by online cluster analysis
		signed __int32		dnParams[MAX_PARAMS];		// Parameters calculated from snData by the clustering algoritm
		TetPoint		snData[SPIKE_NUMPOINTS];	// The A-D data samples

};

// a stereotrode
Nlx_Odl struct STRec	{
		unsigned __int64	qwTimeStamp;				// See TTRec, above
		unsigned __int32	dwScNumber;
		unsigned __int32	dwCellNumber;
		signed __int32		dnParams[MAX_PARAMS];
		StereoPoint		snData[SPIKE_NUMPOINTS];
};

// a single electrode
Nlx_Odl struct SERec	{
		unsigned __int64	qwTimeStamp;				// See TTRec, above
		unsigned __int32	dwScNumber;
		unsigned __int32	dwCellNumber;
		signed __int32		dnParams[MAX_PARAMS];
		SinglePoint		snData[SPIKE_NUMPOINTS];
};

//This is totally generic for passing around the appropriate information in function calls
//May be either a TT, ST, or SE rec.  //Waveforms can be acessed as (short*)(&(*SCRec[1]))
Nlx_Odl struct SCRec	{
		unsigned __int64	qwTimeStamp;
		unsigned __int32	dwScNumber;
		unsigned __int32	dwCellNumber;
		signed __int32	 dnParams[MAX_PARAMS];
};


Nlx_Odl struct SCMin {
		unsigned __int64	qwTimeStamp;
		unsigned __int32	dwCellNumber;		
};
////////////////////////////////
// Continuous Sampling Datatypes
////////////////////////////////

const int MAX_CSC_SAMPLES = 512;


// a CSCAcqEnt record for EEG data
Nlx_Odl struct CRRec	{
		unsigned __int64	qwTimeStamp;			// TS
		unsigned __int32	dwChannelNum;			// Channel number
		unsigned __int32	dwSampleFreq;			// freq in hertz of sampling rate
		unsigned __int32	dwNumValidSamples;		// number of snSamples containing useful data
		signed __int16	 snSamples[MAX_CSC_SAMPLES];		// the A-D data samples
};


//////////////////
// Event Datatypes
//////////////////

const int EVENT_NUM_EXTRAS = 8;
const int NLX_EventRecStringSize = 128;

//This enum is used to determine the source of the event and
//replaces all of the const ints declared in various locations.
//Since the event string for TTL records describes the source
//in more detail, this value is no longer necessary.  However,
//we are going to keep it around for backward compatability.

//The following event IDs are not used in Cheetah, and
//since this file was probably the only place these values were
//used, I'll keep them around in case some other program needs 
//to look up the old values for some reason.  These values are
//skipped in the enumeration in case they need to be added back in.
//CUBE_EVENTID = 0x69
//OPTIONAL_DIGITALINPUT_EVENTID = 10
//CONFIGFILE_DEFAULT_EVENTID = 5
//VIDEOEVENTGEN_EVENTID = 6

enum EventRecID {
	ERID_UNKNOWN = -1,
	ERID_DIGITAL_LYNX = 0,
	ERID_DT3010_BOARD1 = 1, //replaces DT3010_BOARD1_EVENTID
	ERID_DT3010_BOARD2 = 2, //replaces DT3010_BOARD2_EVENTID 
	ERID_RAW_DATA_FILE = 3,
	ERID_MANUAL_EVENT_ENTRY = 4, // any event generated by the user at the keyboard gets this id, replaces KEYBOARD_EVENTID
	ERID_DT3010_BOARD1_DAC_OUTPUT_STARTED = 7, // ttl value will be the wave buffer #, replaces DT3010_BOARD1_DAC_OUTPUT_STARTED
	ERID_DT3010_BOARD2_DAC_OUTPUT_STARTED = 8, // ttl value will be the wave buffer #, replaces DT3010_BOARD2_DAC_OUTPUT_STARTED
	ERID_LYNX_SX = 11,
	ERID_DT3010_BOARD1_ADC_ACQ_STARTED = 12, //acq started for dt3010 board 1, replaces DT3010_BOARD1_ADC_ACQ_STARTED
	ERID_DT3010_BOARD2_ADC_ACQ_STARTED = 14, //acq started for dt3010 board 2, replaces DT3010_BOARD2_ADC_ACQ_STARTED
	ERID_FALCON = 15,
	ERID_SPARROW = 16,
	ERID_ADDITIONAL_DIGITAL_IO = 17, //from addon DIO board, replaces ADDONDIGITALINPUT_EVENTID
	ERID_HAWK = 18,
	ERID_CHEETAH_160 = 119 //replaces DCDCB_EVENTID
};

// record for EventAcqEnt objects
Nlx_Odl struct EventRec	{
	short   	nstx;             			// always 800 from DCDCB
	short   	npkt_id;          			// DCDCB ID  1002, etc.
	short   	npkt_data_size;   			// always 2
	__int64	qwTimeStamp;
	short   	nevent_id;       			// just an id value
	short   	nttl;            			// TTL input value
	short   	ncrc;            			// from the DCDCB
	short   	ndummy1;          			// just a place holder
	short   	ndummy2;          			// just a place holder
	__int32	dnExtra[EVENT_NUM_EXTRAS];       			// extra "bit values"
	char		EventString[NLX_EventRecStringSize]; 	// char string for user input events
};


//////////////////////////
// Video Tracker Datatypes
//////////////////////////

const int NLX_VTREC_NUM_POINTS = 400;
const int NLX_VTREC_NUM_TARGETS = 50;
const unsigned __int16 NLX_VTREC_SWSTX(0x800);

// record for TrackerAcqEnt objects
Nlx_Odl struct VideoRec	{
	unsigned __int16	swstx;					// should be 0x800
	unsigned __int16	swid;					// the ID of the VT as assigned by Cheetah
	unsigned __int16	swdata_size;				// size of the VT record in bytes
	unsigned __int64	qwTimeStamp;				// TS
	unsigned __int32	dwPoints[NLX_VTREC_NUM_POINTS];		// the points with color bit values x&y - note: this is a bit field!
	signed __int16	sncrc;					// ignored, relic from Cheetah160VT
	signed __int32	dnextracted_x;				// calculated x coordinate from our extraction algorithm
	signed __int32	dnextracted_y;				// calculated y coordinate from our extraction algorithm
	signed __int32	dnextracted_angle;			// calculated head direction in degrees from the Y axis
	signed __int32	dntargets[NLX_VTREC_NUM_TARGETS];	// colored targets with same format as the points
};


// bit masks for isolating data from bitfields contained in 'dwPoints' and 'dntargets'

const unsigned int VREC_COLOR_MASK = 0x7000F000;  	// logical OR of all the colors
const unsigned int VREC_LU_MASK = 0x8000;		// luminance mask
const unsigned int VREC_RR_MASK = 0x4000;		// pure & raw RGB masks
const unsigned int VREC_RG_MASK = 0x2000;
const unsigned int VREC_RB_MASK = 0x1000;
const unsigned int VREC_PR_MASK = 0x40000000;
const unsigned int VREC_PG_MASK = 0x20000000;
const unsigned int VREC_PB_MASK = 0x10000000;
const unsigned int VREC_RS_MASK = 0x80000000;		// reserved bit mask
const unsigned int VREC_X_MASK = 0x00000FFF;		// x value mask
const unsigned int VREC_Y_MASK = 0x0FFF0000;		// y val


#pragma pack(pop, before_nlx_datatypes) // back to old packing scheme

#endif  //_NLX_DATATYPE_H

//////
// EOF
//////