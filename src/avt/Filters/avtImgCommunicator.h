#ifndef AVT_IMG_COMMUNICATOR_H
#define AVT_IMG_COMMUNICATOR_H

#include <filters_exports.h>
#include <pipeline_exports.h>
#include <avtSamplePointExtractor.h>
#include <algorithm>

#ifdef PARALLEL
#   include <mpi.h>
#endif

#define MSG_DATA 100
#define MSG_RESULT 101



class avtImgCommunicator
{
	int totalPatches;
	int *processorPatchesCount;
	imgMetaData *allRecvPatches;
	imgData *allRecvImgData;

    int 		num_procs;
    int 		my_id;

    imgMetaData setImg(int _inUse, int _procId, int _patchNumber, float dim_x, float dim_y, float screen_ll_x, float screen_ll_y, float screen_ur_x, float screen_ur_y, float _avg_z);

   

public:
	avtImgCommunicator();
	~avtImgCommunicator();

	void sendNumPatches(int destId, int numPatches);
	void sendPatchMetaData(int destId, imgMetaData tempImg);
	void sendPatchImgData(int destId, int arraySize, float *sendMsgBuffer);

	void masterRecvNumPatches();
	void masterRecvPatchMetaData();
	void masterRecvPatchImgData();

	void syncAllProcs();

	void composeImages(int width, int height, unsigned char *wholeImage);
	void setPatchImg(int procId, int patchNumber, int bufferSize, float buffer[]);

	void printPatches();

	int GetNumProcs(){ return num_procs;};
	int GetMyId(){ return my_id;};

	
#ifdef PARALLEL
	MPI_Status status;
	MPI_Datatype _img_mpi;

	MPI_Datatype createMetaDataType();
#endif
};

#endif