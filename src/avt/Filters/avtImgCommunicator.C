
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
#include <fstream>



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


bool sortImgByDepth(imgMetaData const& before, imgMetaData const& after){ return before.avg_z > after.avg_z; }

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
    
    totalPatches = 0;

    processorPatchesCount = NULL;
	allRecvImgData= NULL;
	allRecvPatches= NULL;
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
		if (processorPatchesCount != NULL)
			delete []processorPatchesCount;

		if (allRecvImgData != NULL){
			for (int i=0; i<totalPatches; i++)
				delete []allRecvImgData[i].imagePatch;

			delete []allRecvImgData;
		}

		if (allRecvPatches != NULL)
			delete []allRecvPatches;
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
int avtImgCommunicator::getDataPatchID(int procID, int patchID){
	int sumPatches = 0;
	for (int i=0; i<procID; i++)
		sumPatches += processorPatchesCount[i];
	
	return (sumPatches+patchID);

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
void avtImgCommunicator::init(){
	if (my_id == 0)
		processorPatchesCount = new int[num_procs];
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
void avtImgCommunicator::gatherNumPatches(int numPatches){
	int patchesProc[2];
	patchesProc[0] = my_id;	patchesProc[1] = numPatches;

	#ifdef PARALLEL
		int *tempRecvBuffer = NULL;
		//std::cout << " " << my_id << "  num_procs: " << num_procs << "!!! Send: " << patchesProc[0] << " - " << patchesProc[1] << std::endl;

		// Only proc 0 receives data
		if (my_id == 0)		
			tempRecvBuffer = new int[num_procs*2];

		MPI_Gather(patchesProc, 2, MPI_INT,   tempRecvBuffer, 2,MPI_INT,         0, MPI_COMM_WORLD);		// all send to proc 0

		if (my_id == 0){		
			for (int i=0; i<num_procs; i++){
				processorPatchesCount[tempRecvBuffer[i*2]] = tempRecvBuffer[i*2 + 1];	// enter the number of patches for each processor
				totalPatches += processorPatchesCount[tempRecvBuffer[i*2]];				// count the number of patches

				std::cout << "!!! Recv: " << tempRecvBuffer[i*2] << "   Patch: " << processorPatchesCount[tempRecvBuffer[i*2]] << std::endl;
			}

			allRecvPatches = new imgMetaData[totalPatches];

			if (tempRecvBuffer != NULL)
				delete []tempRecvBuffer;
			tempRecvBuffer = NULL;
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
void avtImgCommunicator::gatherIotaMetaData(int arraySize, float *allIotaMetadata){
	int *recvSizePerProc = NULL;

	#ifdef PARALLEL
		float *tempRecvBuffer = NULL;
		int *offsetBuffer = NULL;

		if (my_id == 0){
			tempRecvBuffer = new float[totalPatches*4]; // x4: procId, patchNumber, imgArea, avg_z
			recvSizePerProc = new int[num_procs];	
			offsetBuffer = new int[num_procs];	
			for (int i=0; i<num_procs; i++){
				if (i == 0)
					offsetBuffer[i] = 0;
				else
					offsetBuffer[i] = offsetBuffer[i-1] + recvSizePerProc[i-1];

				recvSizePerProc[i] = processorPatchesCount[i]*4;
			}
		}

		syncAllProcs();
		MPI_Gatherv(allIotaMetadata, arraySize, MPI_FLOAT,   tempRecvBuffer, recvSizePerProc, offsetBuffer,MPI_FLOAT,         0, MPI_COMM_WORLD);		// all send to proc 0

		if (my_id == 0){
			iotaMeta tempPatch;
			std::vector<iotaMeta> vectorList;

			for (int i=0; i<totalPatches; i++){
				tempPatch.procId = 		(int) tempRecvBuffer[i*4 + 0];
				tempPatch.patchNumber = (int) tempRecvBuffer[i*4 + 1];
				tempPatch.imgArea = 	(int) tempRecvBuffer[i*4 + 2];
				tempPatch.avg_z = 			  tempRecvBuffer[i*4 + 3];

				vectorList.push_back(tempPatch);
			}

			//for (int i=0; i<vectorList.size(); i++){
			//	std::cout << "procId: " << vectorList[i].procId << "   patch: " << vectorList[i].patchNumber  << "   imgArea: " << vectorList[i].imgArea << "    " << vectorList[i].avg_z << std::endl;
			//}

			vectorList.clear();

			delete []recvSizePerProc;
			delete []offsetBuffer;
			delete []tempRecvBuffer;
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

void avtImgCommunicator::gatherMetaData(int arraySize, float *allIotaMetadata){
	int *recvSizePerProc = NULL;

	#ifdef PARALLEL
		float *tempRecvBuffer = NULL;
		int *offsetBuffer = NULL;

		if (my_id == 0){
			tempRecvBuffer = new float[totalPatches*10]; // x 10: procId, patchNumber, dims[0], dims[1], screen_ll[0], screen_ll[1], screen_ur[0], screen_ur[1], avg_z
			recvSizePerProc = new int[num_procs];	
			offsetBuffer = new int[num_procs];	
			for (int i=0; i<num_procs; i++){
				if (i == 0)
					offsetBuffer[i] = 0;
				else
					offsetBuffer[i] = offsetBuffer[i-1] + recvSizePerProc[i-1];

				recvSizePerProc[i] = processorPatchesCount[i]*10;
			}
		}

		syncAllProcs();
		MPI_Gatherv(allIotaMetadata, arraySize, MPI_FLOAT,   tempRecvBuffer, recvSizePerProc, offsetBuffer,MPI_FLOAT,         0, MPI_COMM_WORLD);		// all send to proc 0

		if (my_id == 0){
			imgMetaData tempRecvImg;

			for (int i=0; i<totalPatches; i++){
				tempRecvImg.procId = 		(int) tempRecvBuffer[i*10 + 0];
				tempRecvImg.patchNumber = 	(int) tempRecvBuffer[i*10 + 1];
				tempRecvImg.inUse = 		(int) tempRecvBuffer[i*10 + 2];
				tempRecvImg.dims[0] = 		(int) tempRecvBuffer[i*10 + 3];
				tempRecvImg.dims[1] = 		(int) tempRecvBuffer[i*10 + 4];
				tempRecvImg.screen_ll[0] = 	(int) tempRecvBuffer[i*10 + 5];
				tempRecvImg.screen_ll[1] = 	(int) tempRecvBuffer[i*10 + 6];
				tempRecvImg.screen_ur[0] = 	(int) tempRecvBuffer[i*10 + 7];
				tempRecvImg.screen_ur[1] = 	(int) tempRecvBuffer[i*10 + 8];
				tempRecvImg.avg_z = 			  tempRecvBuffer[i*10 + 9];

				int patchIndex = getDataPatchID(tempRecvImg.procId, tempRecvImg.patchNumber);
				allRecvPatches[patchIndex] = setImg(tempRecvImg.inUse, tempRecvImg.procId,  tempRecvImg.patchNumber, tempRecvImg.dims[0], tempRecvImg.dims[1],  
								   				tempRecvImg.screen_ll[0],tempRecvImg.screen_ll[1], tempRecvImg.screen_ur[0],tempRecvImg.screen_ur[1], tempRecvImg.avg_z);
			}

			//for (int i=0; i<vectorList.size(); i++){
			//	std::cout << "procId: " << vectorList[i].procId << "   patch: " << vectorList[i].patchNumber  << "   imgArea: " << vectorList[i].imgArea << "    " << vectorList[i].avg_z << std::endl;
			//}

			// Now that we have all the patches, create space to store the img data
			allRecvImgData = new imgData[totalPatches];

			delete []recvSizePerProc;
			delete []offsetBuffer;
			delete []tempRecvBuffer;
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
void avtImgCommunicator::sendPatchImgData(int destId, int arraySize, float *sendMsgBuffer){
#ifdef PARALLEL
	MPI_Send(&arraySize, 1, MPI_INT, destId, 1, MPI_COMM_WORLD);				// #tag 1
	MPI_Send(sendMsgBuffer, arraySize, MPI_FLOAT, destId, 2, MPI_COMM_WORLD);	// #tag 2
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
	int startingProcessor = 1;
	int patchId;

	for (int i=startingProcessor; i<num_procs; i++){ 		// start from 1 to ignore self
		for (int j=0; j<processorPatchesCount[i]; j++){		// for each patch

			int bufferSize;
			MPI_Recv(&bufferSize, 1, MPI_INT, i,1, MPI_COMM_WORLD, &status);				// #tag 1
			float *recvMsgBuffer = new float[bufferSize];
			MPI_Recv(recvMsgBuffer,bufferSize, MPI_FLOAT, i, 2, MPI_COMM_WORLD, &status);	// #tag 2

			patchId = getDataPatchID((int)recvMsgBuffer[0], (int)recvMsgBuffer[1]);

			allRecvImgData[patchId].procId = (int)recvMsgBuffer[0];
			allRecvImgData[patchId].patchNumber = (int)recvMsgBuffer[1];
			allRecvImgData[patchId].imagePatch = new float[bufferSize-2];

			for (int k=0; k<bufferSize-2; k++)
				allRecvImgData[patchId].imagePatch[k] = recvMsgBuffer[k+2];
			
			delete []recvMsgBuffer;
		}
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
void avtImgCommunicator::composeImages(int imgBufferWidth, int imgBufferHeight, unsigned char *wholeImage, unsigned char background[3]){
	// sort the metadata image array by depth
	std::sort(allRecvPatches, allRecvPatches + totalPatches, &sortImgByDepth);

	// Create a buffer to store the composited image
	float *buffer = new float[imgBufferWidth * imgBufferHeight * 4];
	for (int i=0; i<(imgBufferWidth * imgBufferHeight * 4); i+=4){
		buffer[i] = background[0]/255.0;	buffer[i+1] = background[1]/255.0; buffer[i+2] = background[2]/255.0; buffer[i+3] = 0.0;
	}

	int imgCount = 0;
	for (int i=0; i<totalPatches; i++){
		int startingX = allRecvPatches[i].screen_ll[0];
		int startingY = allRecvPatches[i].screen_ll[1]; 

		// get the patchId from the unsorted imageData array
		int patchId = getDataPatchID(allRecvPatches[i].procId, allRecvPatches[i].patchNumber);


		//////////////////////////////////
		////// DEBUG - save image

		// std::ofstream myfile;
		// if (imgCount == 0){
  // 			myfile.open ("/home/pascal/Desktop/example0.txt");
  // 			std::cout << " One " << allRecvPatches[i].dims[0] << " x " << allRecvPatches[i].dims[1] << " ~ " << allRecvPatches[i].screen_ll[0] << " , " << allRecvPatches[i].screen_ll[1] << "  ~  " << allRecvPatches[i].screen_ur[0] << " , " << allRecvPatches[i].screen_ur[1] << " ~ " << allRecvPatches[i].avg_z  << std::endl;
  		
  // 			createPpm(allRecvImgData[patchId].imagePatch, allRecvPatches[i].dims[0], allRecvPatches[i].dims[1], "/home/pascal/Desktop/example0");
  // 		}
		// else
		//  	if (imgCount == 1){
		//  		myfile.open ("/home/pascal/Desktop/example1.txt");
		//  		std::cout << " Two " << allRecvPatches[i].dims[0] << " x " << allRecvPatches[i].dims[1] << " ~ " << allRecvPatches[i].screen_ll[0] << " , " << allRecvPatches[i].screen_ll[1] << "  ~  " << allRecvPatches[i].screen_ur[0] << " , " << allRecvPatches[i].screen_ur[1] << " ~ " << allRecvPatches[i].avg_z  << std::endl;
  			
  // 				createPpm(allRecvImgData[patchId].imagePatch, allRecvPatches[i].dims[0], allRecvPatches[i].dims[1], "/home/pascal/Desktop/example1");
  // 			}
  			
  		//////////////////////////////////



		for (int j=0; j<allRecvPatches[i].dims[1]; j++){
			for (int k=0; k<allRecvPatches[i].dims[0]; k++){

				if ((startingX + k) > imgBufferWidth)
					continue;

				if ((startingY + j) > imgBufferHeight)
					continue;
				
				int subImgIndex = allRecvPatches[i].dims[0]*j*4 + k*4;										// index in the subimage 
				int bufferIndex = (startingY*imgBufferWidth*4 + j*imgBufferWidth*4) + (startingX*4 + k*4);	// index in the big buffer


				//Front to back compositing: 
				//composited = source * (1.0 - destination.a) + destination; 
				//buffer[bufferIndex+0] = allRecvImgData[patchId].imagePatch[subImgIndex+0] * (1.0 - buffer[bufferIndex+3]) + buffer[bufferIndex+0];
				//buffer[bufferIndex+1] = allRecvImgData[patchId].imagePatch[subImgIndex+1] * (1.0 - buffer[bufferIndex+3]) + buffer[bufferIndex+1];
				//buffer[bufferIndex+2] = allRecvImgData[patchId].imagePatch[subImgIndex+2] * (1.0 - buffer[bufferIndex+3]) + buffer[bufferIndex+2];
				//buffer[bufferIndex+3] = allRecvImgData[patchId].imagePatch[subImgIndex+3] * (1.0 - buffer[bufferIndex+3]) + buffer[bufferIndex+3];

				
				//////////////////////////////////
				////// DEBUG - save image
				
				//if (imgCount < 2){
				//	std::cout << j << ",  " << k << "   Buffer: " << buffer[bufferIndex+0] << " ,  " << buffer[bufferIndex+1] << " , " << buffer[bufferIndex+2] << " ,  " << buffer[bufferIndex+3];
				//	std::cout << "   incoming fragment: " << allRecvImgData[patchId].imagePatch[subImgIndex+0] << " ,  " << allRecvImgData[patchId].imagePatch[subImgIndex+1] << " ,  " << allRecvImgData[patchId].imagePatch[subImgIndex+2] << " ,  " << allRecvImgData[patchId].imagePatch[subImgIndex+3];
				//} 
				//if (imgCount < 2){
				//	myfile << allRecvImgData[patchId].imagePatch[subImgIndex+0] << " " << allRecvImgData[patchId].imagePatch[subImgIndex+1] << " " << allRecvImgData[patchId].imagePatch[subImgIndex+2] << " " << allRecvImgData[patchId].imagePatch[subImgIndex+3] << "\n";
				//}
				
				//////////////////////////////////

				
				//Back to Front compositing: 
				//composited_i = composited_i-1 * (1.0 - alpha_i) + incoming
				//alpha = alpha_i-1 * (1- alpha_i)
				buffer[bufferIndex+0] = (buffer[bufferIndex+0] * (1.0 - allRecvImgData[patchId].imagePatch[subImgIndex+3])) + allRecvImgData[patchId].imagePatch[subImgIndex+0];
				buffer[bufferIndex+1] = (buffer[bufferIndex+1] * (1.0 - allRecvImgData[patchId].imagePatch[subImgIndex+3])) + allRecvImgData[patchId].imagePatch[subImgIndex+1];
				buffer[bufferIndex+2] = (buffer[bufferIndex+2] * (1.0 - allRecvImgData[patchId].imagePatch[subImgIndex+3])) + allRecvImgData[patchId].imagePatch[subImgIndex+2];
				buffer[bufferIndex+3] = (buffer[bufferIndex+3] * (1.0 - allRecvImgData[patchId].imagePatch[subImgIndex+3]));
			}
			//if (imgCount < 2){
			//		myfile <<  "\n ";
			//}
		}

		//////////////////////////////////
		////// DEBUG - save image
		
		// if (imgCount == 0 || imgCount == 1){
  // 			myfile.close();
  // 		}

		// if (imgCount == 0){
  // 			createPpm(buffer, imgBufferWidth, imgBufferHeight, "/home/pascal/Desktop/one");
  // 		}

  // 		if (imgCount == 1){
  // 			createPpm(buffer, imgBufferWidth, imgBufferHeight, "/home/pascal/Desktop/two");
  // 		}

  		imgCount++;
  		
  		//////////////////////////////////
	}

	//Copy to main buffer
	for (int i=0; i< imgBufferHeight; i++)
		for (int j=0; j<imgBufferWidth; j++){
			int bufferIndex = (imgBufferWidth*4*i) + (j*4);
			int wholeImgIndex = (imgBufferWidth*3*i) + (j*3);

			wholeImage[wholeImgIndex+0] = (buffer[bufferIndex+0] ) * 255;
			wholeImage[wholeImgIndex+1] = (buffer[bufferIndex+1] ) * 255;
			wholeImage[wholeImgIndex+2] = (buffer[bufferIndex+2] ) * 255;
		}

	//////////////////////////////////
	////// DEBUG - save image

    //std::string imgFilenameFinal = "/home/pascal/Desktop/FinalBuffer.ppm";
    //createPpm(buffer, imgBufferWidth, imgBufferHeight, imgFilenameFinal);
  	
  	//////////////////////////////////

	delete []buffer;
	buffer = NULL;


	if (my_id == 0){
		if (processorPatchesCount != NULL)
			delete []processorPatchesCount;
		processorPatchesCount = NULL;

		if (allRecvImgData != NULL){
			for (int i=0; i<totalPatches; i++)
				delete []allRecvImgData[i].imagePatch;

			delete []allRecvImgData;
		}
		allRecvImgData = NULL;


		if (allRecvPatches != NULL)
			delete []allRecvPatches;
		allRecvPatches = NULL;
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
void avtImgCommunicator::setPatchImg(int procId, int patchNumber, int bufferSize, float buffer[]){

	int patchId = (int)buffer[1];

	allRecvImgData[patchId].procId = (int)buffer[0];
	allRecvImgData[patchId].patchNumber = (int)buffer[1];


	allRecvImgData[patchId].imagePatch = new float[bufferSize-2];
	for (int i=0; i<bufferSize-2; i++)
		allRecvImgData[patchId].imagePatch[i] = buffer[i+2];
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
void avtImgCommunicator::printPatches(){
	for (int i=0; i<totalPatches; i++)
		displayStruct_meta(allRecvPatches[i]);

	for (int i=0; i<num_procs; i++)
		printf("\n Processor: %d   Patch count: %d \n",i,processorPatchesCount[i]);
}



void createPpm(float array[], int dimx, int dimy, std::string filename){
    int i, j;
    std::cout << "createPpm2  dims: " << dimx << ", " << dimy << " -  " << filename.c_str() << std::endl;
    FILE *fp = fopen(filename.c_str(), "wb"); // b - binary mode 
    (void) fprintf(fp, "P6\n%d %d\n255\n", dimx, dimy);
    for (j = 0; j < dimy; ++j){
        for (i = 0; i < dimx; ++i){
            static unsigned char color[3];
            color[0] = array[j*(dimx*4) + i*4 + 0] * 255;  // red
            color[1] = array[j*(dimx*4) + i*4 + 1] * 255;  // green
            color[2] = array[j*(dimx*4) + i*4 + 2] * 255;  // blue 
            (void) fwrite(color, 1, 3, fp);
        }
    }
    (void) fclose(fp);
    std::cout << "End createPpm: " << std::endl;
}

std::string NumbToString (int Number)
{
     std::ostringstream ss;
     ss << Number;
     return ss.str();
}
