#ifndef AVT_IMG_COMMUNICATOR_H
#define AVT_IMG_COMMUNICATOR_H

#include <filters_exports.h>
#include <pipeline_exports.h>
#include <avtSamplePointExtractor.h>
#include <algorithm>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#ifdef PARALLEL
#   include <mpi.h>
#endif

#define MSG_DATA 100
#define MSG_RESULT 101

struct imageBuffer{
	float *image;
	float depth;
};

struct code{
	float count;	// should be int but float makes it easier to send with MPI!
	float color[4];
};

class avtImgCommunicator
{
	int totalPatches;
	int numPatchesToCompose;
	int *processorPatchesCount;

	// each processor maintains a list of its neighbours that resides on the same node - local compositing takes advantage of that
	int numMyNeighboursOnNode;
	std::vector<int> neighboursId;

	float *imgBuffer;
	iotaMeta *allRecvIotaMeta;

	std::set<float> all_avgZ_proc0;
	std::vector<std::vector<float> > boundsPerBlockVec;

	int *patchesToSendArray, *patchesToRecvArray;
	int *numPatchesToSendArray, *numPatchesToRecvArray;
	int *recvDisplacementForProcs, *sendDisplacementForProcs;

	int *numPatchesToSendRecvArray;
	float *boundsPerBlockArray;
	int *blockDisplacementForProcs;
	int *numBlocksPerProc;

	int* patchesToCompositeLocallyArray;
	int* numPatchesToCompositeLocally;
	int* compositeDisplacementForProcs;

	int *compressedSizePerDiv;	//size of each division
	
	unsigned char background[3];

    int 		num_procs;
    int 		my_id;
    std::string hostname;

    imgMetaData setImg(int _inUse, int _procId, int _patchNumber, float dim_x, float dim_y, float screen_ll_x, float screen_ll_y, float screen_ur_x, float screen_ur_y, float _avg_z);
    iotaMeta setIota(int _procId, int _patchNumber, int dim_x, int dim_y, int screen_ll_x, int screen_ll_y, float _avg_z);
    int getDataPatchID(int procID, int patchID);

    int getHostname(hostname);

public:
	avtImgCommunicator();
	~avtImgCommunicator();

	void init();

	void gatherNumPatches(int numPatches);
	void gatherIotaMetaData(int arraySize, float *allIotaMetadata);

	void patchAllocationLogic();		// decides which processor should get which patches and tell each processor how many patches it will receive

	void scatterNumDataToCompose(int& totalSendData, int& totalRecvData, int& numDivisions, int& totalPatchesToCompositeLocally);
	void scatterDataToCompose(	int& totalSendData, int* informationToSendArray, 
								int& totalRecvData, int* informationToRecvArray, 
								int& numDivisions, float* blocksPerProc,
								int& totalPatchesToCompositeLocally, int* patchesToCompositeLocally);

	void sendPointToPoint(imgMetaData toSendMetaData, imgData toSendImgData, int tag);	// Send out the patches and receive them
	void recvPointToPoint(imgMetaData &recvMetaData, imgData &recvImgData);

	void recvPointToPointMetaData(imgMetaData &recvMetaData, int tag);
	void recvPointToPointImgData(imgMetaData recvMetaData, imgData &recvImgData, int tag);

	void gatherEncodingSizes(int *sizeEncoding, int numDivisions);
	void gatherAndAssembleEncodedImages(int sizex, int sizey, int sizeSending, float *image, int numDivisions);		// do the compositing of the subpatches

	void gatherAndAssembleImages(int sizex, int sizey, float *image, int numDivisions);		// do the compositing of the subpatches

	void getcompositedImage(int imgBufferWidth, int imgBufferHeight, unsigned char *wholeImage);	// get the final composited image



	// RLE Encoding
	int rleEncodeAll(int dimsX, int dimsY, int numDivs, float *imgArray,  float *& encoding, int *& sizeOfEncoding);
	void rleDecode(int encSize, float *encoding, int offset, float *& img);


	void syncAllProcs();
	int GetNumProcs(){ return num_procs;};
	int GetMyId(){ return my_id;};

	float clamp(float x);
	void setBackground(unsigned char _background[3]){ for (int i=0; i<3; i++) background[i] = _background[i]; }
	
#ifdef PARALLEL
	MPI_Status status;
	MPI_Datatype _img_mpi;

	MPI_Datatype createMetaDataType();
#endif
};



void createPpm(float array[], int dimx, int dimy, std::string filename);
void createPpmWithOffset(float array[], int dimx, int dimy, std::string filename, int offset);
std::string NumbToString (int Number);
#endif