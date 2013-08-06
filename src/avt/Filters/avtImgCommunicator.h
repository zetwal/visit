#ifndef AVT_IMG_COMMUNICATOR_H
#define AVT_IMG_COMMUNICATOR_H

#include <filters_exports.h>
#include <pipeline_exports.h>
#include <avtSamplePointExtractor.h>
#include <algorithm>
#include <string>

#ifdef PARALLEL
#   include <mpi.h>
#endif

#define MSG_DATA 100
#define MSG_RESULT 101



class avtImgCommunicator
{
	int totalPatches;
	int numPatchesToCompose;
	int *processorPatchesCount;
	imgMetaData *allRecvPatches;
	imgData *allRecvImgData;

	
	iotaMeta *allRecvIotaMeta;
	std::vector<int> procToSend;
	


    int 		num_procs;
    int 		my_id;

    imgMetaData setImg(int _inUse, int _procId, int _patchNumber, float dim_x, float dim_y, float screen_ll_x, float screen_ll_y, float screen_ur_x, float screen_ur_y, float _avg_z);
    iotaMeta setIota(int _procId, int _patchNumber, float _imgArea, float _avg_z);
    int getDataPatchID(int procID, int patchID);
   

public:
	avtImgCommunicator();
	~avtImgCommunicator();

	void init();

	void gatherNumPatches(int numPatches);
	void gatherIotaMetaData(int arraySize, float *allIotaMetadata);

	void patchDecisionallocation();		// decides which processor should get which patches and tell each processor how many patches it will receive

	void receiveNumPatchesToCompose();
	void sendRecvandRecvInfo();
	void recvDataforDataToRecv(int &totalSendData, int *informationToSendArray, int &totalRecvData, int *informationToRecvArray);
	void sendnRecvPatchesMetanData();

	void sendPatchImgData(int destId, int arraySize, float *sendMsgBuffer);
	void masterRecvPatchImgData();

	void syncAllProcs();

	void composeImages(int width, int height, unsigned char *wholeImage, unsigned char background[3]);
	void setPatchImg(int procId, int patchNumber, int bufferSize, float buffer[]);

	void printPatches();
	

	int GetNumProcs(){ return num_procs;};
	int GetMyId(){ return my_id;};

	void gatherMetaData(int arraySize, float *allIotaMetadata);
	
#ifdef PARALLEL
	MPI_Status status;
	MPI_Datatype _img_mpi;

	MPI_Datatype createMetaDataType();
#endif
};



void createPpm(float array[], int dimx, int dimy, std::string filename);
std::string NumbToString (int Number);
#endif