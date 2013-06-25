#include <algorithm>
#include <cmath>

#include <avtParallel.h>
#include <ImproperUseException.h>

#include <stdio.h>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#ifdef PARALLEL
#   include <mpi.h>
#endif

#include <avtImgCommunicator.h>



/*
struct imgMetaData{
  int inUse;        // 1: in use/ 0:not in use
  int procId;       // to be removed!!!!!
  int patchNumber;
  int dims[2];      // height, width
  int screen_ll[2];
  int screen_ur[2];
  float avg_z;
};
*/

void displayStruct_imgRecv(int id, int from_id, imgMetaData temp){
	printf("\n || Recv Id: %d from: %d \n In  use: %d, Proc Id: %d , Patch Id: %d, \n Size: %d, %d,  \n Pos: %d, %d   %d, %d , \n Depth: %.2f \n\n",
			id, from_id,
			temp.inUse,
			temp.procId, temp.patchNumber, 
			temp.dims[0], temp.dims[1], 
			temp.screen_ll[0], temp.screen_ll[1], temp.screen_ur[0],temp.screen_ur[1], 
			temp.avg_z);
}


void displayStruct_imgSend(int id, int to_id, imgMetaData temp){
	printf("\n || Send Id: %d to: %d \n In  use: %d, Proc Id: %d , Patch Id: %d, \n Size: %d, %d,  \n Pos: %d, %d   %d, %d , \n Depth: %.2f \n\n",
			id, to_id,
			temp.inUse,
			temp.procId, temp.patchNumber, 
			temp.dims[0], temp.dims[1], 
			temp.screen_ll[0], temp.screen_ll[1], temp.screen_ur[0],temp.screen_ur[1], 
			temp.avg_z);
}


void displayStruct_meta(imgMetaData temp){
	printf("\n || In  use: %d, Proc Id: %d , Patch Id: %d, \n Size: %d, %d,  \n Pos: %d, %d   %d, %d , \n Depth: %.2f \n\n",
			temp.inUse,
			temp.procId, temp.patchNumber, 
			temp.dims[0], temp.dims[1], 
			temp.screen_ll[0], temp.screen_ll[1], temp.screen_ur[0],temp.screen_ur[1], 
			temp.avg_z);
}

bool sortImgByDepth(imgMetaData const& before, imgMetaData const& after){
	return before.avg_z < after.avg_z;
}


// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
avtImgCommunicator::avtImgCommunicator(){
	int ierr;

#ifdef PARALLEL
	ierr = MPI_Comm_size(VISIT_MPI_COMM, &num_procs);
	ierr = MPI_Comm_rank(VISIT_MPI_COMM, &my_id);

	// create a new mpi data type matching the one we need
	_img_mpi = createMetaDataType();
	MPI_Type_commit(&_img_mpi);
#else
    num_procs = 1;
    my_id = 0;
#endif
    
    std::cout << "!!! Id: " << my_id << "   Total: : " << num_procs << std::endl;

    if (my_id == 0){
		processorPatchesCount = new int[num_procs];
		totalPatches = 0;
	}
}



// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
avtImgCommunicator::~avtImgCommunicator(){
	if (my_id == 0){
		delete []processorPatchesCount;
		delete []allRecvPatches;
	}
}



// ****************************************************************************
//  Method: avtImgCommunicator::waitToSync
//
//  Purpose:
//			Wait for all processors to hearch here before continuing
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
void avtImgCommunicator::syncAllProcs(){

#ifdef PARALLEL
	MPI_Barrier(MPI_COMM_WORLD);
#endif

}


// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
void avtImgCommunicator::sendNumPatches(int destId, int numPatches){
	int patchesProc[2];
	patchesProc[0] = my_id;	patchesProc[1] = numPatches;

#ifdef PARALLEL
	MPI_Send(patchesProc, 2, MPI_INT, 0, MSG_DATA, MPI_COMM_WORLD);
#endif

}

/*
struct imgData{
  float patchNumber;
  float *imagePatch;
};
*/

// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
void avtImgCommunicator::sendPatchMetaData(int destId, imgMetaData tempImg){

#ifdef PARALLEL
	MPI_Send(&tempImg, 1, _img_mpi, 0, MSG_DATA, MPI_COMM_WORLD);

	//if (PAR_Rank() == 4)
	//	displayStruct_imgSend(PAR_Rank(),destId,tempImg);
#endif
}



// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
void avtImgCommunicator::sendPatchImgData(int destId, int arraySize, float *sendMsgBuffer){
#ifdef PARALLEL
	MPI_Send(sendMsgBuffer, arraySize, MPI_FLOAT, destId, MSG_DATA, MPI_COMM_WORLD);
#endif
}


// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
void avtImgCommunicator::masterRecvNumPatches(){
#ifdef PARALLEL
	if (my_id == 0){
		int tempRecvBuffer[2];
		for (int i=0; i<num_procs; i++){
			MPI_Recv(&tempRecvBuffer, 2, MPI_INT, MPI_ANY_SOURCE, MSG_DATA, MPI_COMM_WORLD, &status);	// receive
			processorPatchesCount[tempRecvBuffer[0]] = tempRecvBuffer[1];	// convert to the correct stucture
			totalPatches += processorPatchesCount[tempRecvBuffer[0]];		// count the number of patches
			
			std::cout << "!!! Recv: " << tempRecvBuffer[0] << "   Patch: " << processorPatchesCount[tempRecvBuffer[0]] << std::endl;
		}	
		
		allRecvPatches = new imgMetaData[totalPatches];
	}
#endif

}


// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
void avtImgCommunicator::masterRecvPatchMetaData(){
	#ifdef PARALLEL
	imgMetaData tempRecvImg;

	for (int i=0; i<totalPatches; i++){
		// Receive
		MPI_Recv (&tempRecvImg, 1, _img_mpi, MPI_ANY_SOURCE, MSG_DATA, MPI_COMM_WORLD, &status);
			
		int patchIndex = tempRecvImg.procId * processorPatchesCount[tempRecvImg.procId] + tempRecvImg.patchNumber;
		allRecvPatches[patchIndex] = setImg(tempRecvImg.inUse, tempRecvImg.procId,  tempRecvImg.patchNumber, tempRecvImg.dims[0], tempRecvImg.dims[1],  
							   tempRecvImg.screen_ll[0],tempRecvImg.screen_ll[1], tempRecvImg.screen_ur[0],tempRecvImg.screen_ur[1], tempRecvImg.avg_z);

		//if (tempRecvImg.procId == 4)
		//	displayStruct_imgRecv(0,tempRecvImg.procId, allRecvPatches[i]);
	}
	#endif
}



// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
void avtImgCommunicator::masterRecvPatchImgData(){
#ifdef PARALLEL
	int patchId = 0;
	allRecvImgData = new imgData[totalPatches];

	for (int i=0; i<num_procs; i++){ // start from 1 to ignore self
		for (int j=0; j<processorPatchesCount[i]; j++){
			//printf("\n proc: %d  			processor patch count: %d     patchID: %d    dims: %d, %d ",i, processorPatchesCount[i], patchId, allRecvPatches[patchId].dims[0], allRecvPatches[patchId].dims[1]);
			int imgSize = (allRecvPatches[patchId].dims[0] * allRecvPatches[patchId].dims[1]);

			float *recvMsgBuffer = new float[imgSize + 1];
			MPI_Recv(recvMsgBuffer,imgSize + 1, MPI_FLOAT, i, MSG_DATA, MPI_COMM_WORLD, &status);

			allRecvImgData[patchId].imagePatch = new float[imgSize];

			allRecvImgData[patchId].patchNumber = (int)recvMsgBuffer[0];
			for (int k=0; k<imgSize; k++){
				allRecvImgData[patchId].imagePatch[k] = recvMsgBuffer[k+1];
			}
			patchId++;

			//printf("\ntotal patches: %d    patchID: %d    dims: %d, %d ",totalPatches, patchId, allRecvPatches[patchId].dims[0], allRecvPatches[patchId].dims[1]);

			delete []recvMsgBuffer;
		}
	}
#endif
}


void avtImgCommunicator::printPatches(){
	for (int i=0; i<totalPatches; i++)
		displayStruct_meta(allRecvPatches[i]);

	for (int i=0; i<num_procs; i++)
		printf("\n Processor: %d   Patch count: %d \n",i,processorPatchesCount[i]);
}

// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
void avtImgCommunicator::composeImages(int imgBufferWidth, int imgBufferHeight, unsigned char *wholeImage){
	// sort the images by depth
	std::sort(allRecvPatches, allRecvPatches + totalPatches, &sortImgByDepth);

	// compose
	float *buffer = new float[imgBufferWidth * imgBufferHeight * 4];
	wholeImage = new unsigned char[imgBufferWidth * imgBufferHeight];

	for (int i=0; i<(imgBufferWidth * imgBufferHeight * 4); i++)
		buffer[i] = 0.0;

	for (int i=0; i<totalPatches; i++){
		int startingY = allRecvPatches[i].screen_ur[1];
		int startingX = allRecvPatches[i].screen_ll[0];

		int patchIndex = allRecvPatches[i].procId * processorPatchesCount[allRecvPatches[i].procId] + allRecvPatches[i].patchNumber;

		for (int j=0; j<allRecvPatches[i].dims[1]; j++)
			for (int k=0; k<allRecvPatches[i].dims[0]; k++){
				int subImgIndex = allRecvPatches[i].dims[0]*4*j +  k*4;
				int indexBuffer = (imgBufferWidth*4*startingY +  imgBufferWidth*4*j) + (startingX*4) + k*4;

				//Front to back compositing: composited = source * (1.0 - destination.a) + destination; 
				for (int color=0; color<4; color++)
					buffer[indexBuffer+color] = allRecvImgData[indexBuffer].imagePatch[subImgIndex+color] * (1.0 - buffer[indexBuffer+3]) + buffer[indexBuffer+color];
				

				// convert the image to unsigned char to send back
				for (int color=0; color<3; color++)
					wholeImage[indexBuffer+color] = (255 * buffer[indexBuffer+color]) *  buffer[indexBuffer+3];
				
			}
		
	}
	delete []buffer;
}


// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
imgMetaData avtImgCommunicator::setImg(int _inUse, int _procId, int _patchNumber, float dim_x, float dim_y, float screen_ll_x, float screen_ll_y, float screen_ur_x, float screen_ur_y, float _avg_z){
	imgMetaData temp;
	temp.inUse = _inUse;
	temp.procId = _procId;
	temp.patchNumber = _patchNumber;
	temp.dims[0] = dim_x;					temp.dims[1] = dim_y;
	temp.screen_ll[0] = screen_ll_x;		temp.screen_ll[1] = screen_ll_y;
	temp.screen_ur[0] = screen_ur_x;		temp.screen_ur[1] = screen_ur_y;
	temp.avg_z = _avg_z;
	
	return temp;
}

// ****************************************************************************
//  Method: avtImgCommunicator::
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
#ifdef PARALLEL
MPI_Datatype avtImgCommunicator::createMetaDataType(){
	MPI_Datatype _imgMeta_mpi;
	const int numItems = 7;
	int blockLengths[numItems] = {1, 1, 1, 2, 2, 2, 1};
	MPI_Datatype type[numItems] = { MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_FLOAT };
	MPI_Aint offsets[numItems] = {0, sizeof(int), sizeof(int)*2, sizeof(int)*3, sizeof(int)*5, sizeof(int)*7, sizeof(int)*9 };
	MPI_Type_struct(numItems, blockLengths,  offsets, type, &_imgMeta_mpi);
	
	return _imgMeta_mpi;
}

#endif


