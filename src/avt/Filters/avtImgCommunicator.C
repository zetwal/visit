
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


bool sortImgByDepth(imgMetaData const& before, imgMetaData const& after){ return before.avg_z < after.avg_z; }

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

		if (allRecvImgData != NULL)
			delete []allRecvImgData;

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
void avtImgCommunicator::sendNumPatches(int destId, int numPatches){
	int patchesProc[2];
	patchesProc[0] = my_id;	patchesProc[1] = numPatches;

#ifdef PARALLEL
	MPI_Send(patchesProc, 2, MPI_INT, 0, MSG_DATA, MPI_COMM_WORLD);
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
			
			//std::cout << "!!! Recv: " << tempRecvBuffer[0] << "   Patch: " << processorPatchesCount[tempRecvBuffer[0]] << std::endl;
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

	// Now that we have all the patches, create space to store the img data
	allRecvImgData = new imgData[totalPatches];
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

	//MPI_Send(&curLength, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
	//		MPI_Send(pszFileBuffer+curStartNum, curLength, MPI_CHAR, i, 2, MPI_COMM_WORLD);
	MPI_Send(&arraySize, 1, MPI_INT, destId, 1, MPI_COMM_WORLD);
	MPI_Send(sendMsgBuffer, arraySize, MPI_FLOAT, destId, 2, MPI_COMM_WORLD);

	//if ((int)sendMsgBuffer[0] == 35){
	//	for (int i=1; i<arraySize; i+=4)
	//		printf("\n 35: %d <> %.6f %.6f %.6f %.6f \n",i/4, sendMsgBuffer[i],sendMsgBuffer[i+1],sendMsgBuffer[i+2],sendMsgBuffer[i+3]);
	//}
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
	//allRecvImgData = new imgData[totalPatches];

	for (int i=startingProcessor; i<num_procs; i++){ // start from 1 to ignore self
		for (int j=0; j<processorPatchesCount[i]; j++){

/*
			printf("\n proc: %d  			processor patch count: %d     patchID: %d    dims: %d, %d ",i, processorPatchesCount[i], patchId, allRecvPatches[patchId].dims[0], allRecvPatches[patchId].dims[1]);
			int imgSize = (allRecvPatches[patchId].dims[0] * allRecvPatches[patchId].dims[1] * 4);

			float *recvMsgBuffer = new float[imgSize + 2];

			int bufferSize;
			MPI_Recv(&bufferSize, 1, MPI_INT, 0,1, MPI_COMM_WORLD, &status);
			MPI_Recv(&bufferSize,imgSize + 1, MPI_FLOAT, i, 1, MPI_COMM_WORLD, &status);
			MPI_Recv(recvMsgBuffer,imgSize + 1, MPI_FLOAT, i, 2, MPI_COMM_WORLD, &status);

			allRecvImgData[patchId].imagePatch = new float[imgSize];

			allRecvImgData[patchId].procId = (int)recvMsgBuffer[0];
			allRecvImgData[patchId].patchNumber = (int)recvMsgBuffer[1];


			for (int k=0; k<imgSize; k++){
				allRecvImgData[patchId].imagePatch[k] = recvMsgBuffer[k+1];
			}
			patchId++;


*/

			//printf("\n proc: %d  			processor patch count: %d     patchID: %d    dims: %d, %d ",i, processorPatchesCount[i], patchId, allRecvPatches[patchId].dims[0], allRecvPatches[patchId].dims[1]);
			//int imgSize = (allRecvPatches[patchId].dims[0] * allRecvPatches[patchId].dims[1] * 4);

			//float *recvMsgBuffer = new float[imgSize + 2];

			int bufferSize;
			MPI_Recv(&bufferSize, 1, MPI_INT, i,1, MPI_COMM_WORLD, &status);
			float *recvMsgBuffer = new float[bufferSize];
			MPI_Recv(recvMsgBuffer,bufferSize, MPI_FLOAT, i, 2, MPI_COMM_WORLD, &status);

			patchId = 0;
			for (int k=0; k<((int)recvMsgBuffer[0]); k++)
				patchId += processorPatchesCount[k];
			patchId += ((int)recvMsgBuffer[1]);

			
			allRecvImgData[patchId].procId = (int)recvMsgBuffer[0];
			allRecvImgData[patchId].patchNumber = (int)recvMsgBuffer[1];
			allRecvImgData[patchId].imagePatch = new float[bufferSize-2];
			for (int k=0; k<bufferSize-2; k++){
				allRecvImgData[patchId].imagePatch[k] = recvMsgBuffer[k+2];
			}
			

			
			//
			// Debug
			// 
			/*
			if ((int)recvMsgBuffer[0] == 5 && (int)recvMsgBuffer[1] == 35){
				printf("\n");
				printf("\n Recv Proc: %d    Patch: %d    PatchId: %d",                 (int)recvMsgBuffer[0],               (int)recvMsgBuffer[1], patchId);
				printf("\n Img  Proc: %d    Patch: %d    PatchId: %d",        allRecvImgData[patchId].procId, allRecvImgData[patchId].patchNumber, patchId);
				printf("\n Meta Proc: %d    Patch: %d    PatchId: %d \n",     allRecvPatches[patchId].procId, allRecvPatches[patchId].patchNumber, patchId);
				for (int i=2; i<bufferSize-2; i+=4)
					printf("\n 35: %d ** %.6f %.6f %.6f %.6f \n",i/4, recvMsgBuffer[i],recvMsgBuffer[i+1],recvMsgBuffer[i+2],recvMsgBuffer[i+3]);

				//spPatchId = patchId;
			}
			*/

			//if (allRecvPatches[patchId].procId == 5 && allRecvPatches[patchId].patchNumber == 35){
            //    std::string imgFilename = "/home/pascal/Desktop/examplePtEx_in_avtImgComm.ppm";
            //    createPpm(allRecvImgData[patchId].imagePatch, allRecvPatches[patchId].dims[0], allRecvPatches[patchId].dims[1], imgFilename);
            //}

			/*
			//if (i==0 && j ==0){
			//if (i==5 && j ==35){
			if (allRecvPatches[i].procId == 5 && allRecvPatches[i].patchNumber == 35){
				std::cout << "\n\nPar_Rank(): " << allRecvPatches[i].procId << "   patch: " << allRecvPatches[i].patchNumber << std::endl;

                std::cout << "avtImgCommunicator height: " << allRecvPatches[i].dims[1] << "   width: " << allRecvPatches[i].dims[0] << std::endl;
                std::cout << "avtImgCommunicator screen_ll[0]: " << allRecvPatches[i].screen_ll[0] << "   screen_ll[1]: " << allRecvPatches[i].screen_ll[1] << std::endl;
                std::cout << "avtImgCommunicator screen_ur[0]: " << allRecvPatches[i].screen_ur[0] << "   screen_ur[1]: " << allRecvPatches[i].screen_ur[1] << std::endl;
                std::cout << "avtImgCommunicator avg_z: " << allRecvPatches[i].avg_z << std::endl;

				for (int l=0; l<allRecvPatches[0].dims[1]; l++){
	                for (int m=0; m<allRecvPatches[0].dims[0]; m++){

	                    int index = l*(4*allRecvPatches[0].dims[0]) + m*4;
	                    std::cout << index/4 << " ~ "  << allRecvImgData[0].imagePatch[index]<< ", " << allRecvImgData[0].imagePatch[index+1] << ", " << allRecvImgData[0].imagePatch[index+2] << ", " << allRecvImgData[0].imagePatch[index+3] << "  \n  ";
	                }
	                std::cout << "\n";
	            }
	            std::cout << "\n";
        	}
        	
			//printf("\ntotal patches: %d    patchID: %d    dims: %d, %d ",totalPatches, patchId, allRecvPatches[patchId].dims[0], allRecvPatches[patchId].dims[1]);
	*/
			delete []recvMsgBuffer;
			//}
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
void avtImgCommunicator::composeImages(int imgBufferWidth, int imgBufferHeight, unsigned char *wholeImage){
	// sort the images by depth
	std::sort(allRecvPatches, allRecvPatches + totalPatches, &sortImgByDepth);

	// Compose

	// Create a buffer to store the compositing
	float *buffer = new float[imgBufferWidth * imgBufferHeight * 4];
	for (int i=0; i<(imgBufferWidth * imgBufferHeight * 4); i++)
		buffer[i] = 0.0;


	//for (int i=0; i<1; i++){
	for (int i=0; i<totalPatches; i++){
		int startingY = allRecvPatches[i].screen_ll[1];  // need to invert
		int startingX = allRecvPatches[i].screen_ll[0];

		int patchId = 0;
		for (int k=0; k<allRecvPatches[i].procId; k++)
			patchId += processorPatchesCount[k];
		patchId += allRecvPatches[i].patchNumber;

		//printf("\n startingX %d - startingY %d \n", allRecvPatches[i].screen_ll[0], allRecvPatches[i].screen_ll[1]);

		//if (allRecvPatches[i].procId==5 && allRecvPatches[i].patchNumber==35){
        //    std::string imgFilename = "/home/pascal/Desktop/examplePtEx_in_avtImgComm_compose.ppm";
        //    createPpm(allRecvImgData[patchId].imagePatch, allRecvPatches[i].dims[0], allRecvPatches[i].dims[1], imgFilename);
        //}	


		for (int j=0; j<allRecvPatches[i].dims[1]; j++){
			for (int k=0; k<allRecvPatches[i].dims[0]; k++){
				int subImgIndex = allRecvPatches[i].dims[0]*j*4 + k*4;
				int bufferIndex = (startingY*imgBufferWidth*4 + j*imgBufferWidth*4) + (startingX*4 + k*4);

				//Front to back compositing: 
				//composited = source * (1.0 - destination.a) + destination; 
				buffer[bufferIndex+0] = allRecvImgData[patchId].imagePatch[subImgIndex+0] * (1.0 - buffer[bufferIndex+3]) + buffer[bufferIndex+0];
				buffer[bufferIndex+1] = allRecvImgData[patchId].imagePatch[subImgIndex+1] * (1.0 - buffer[bufferIndex+3]) + buffer[bufferIndex+1];
				buffer[bufferIndex+2] = allRecvImgData[patchId].imagePatch[subImgIndex+2] * (1.0 - buffer[bufferIndex+3]) + buffer[bufferIndex+2];
				buffer[bufferIndex+3] = allRecvImgData[patchId].imagePatch[subImgIndex+3] * (1.0 - buffer[bufferIndex+3]) + buffer[bufferIndex+3];

				//Back to Front compositing: 
				//composited_i = composited_i-1 * (1.0 - alpha_i) + incoming
				//alpha = alpha_i-1 * (1- alpha_i)
				//buffer[bufferIndex+0] = (buffer[bufferIndex+0] * (1.0 - allRecvImgData[patchId].imagePatch[subImgIndex+3])) + allRecvImgData[patchId].imagePatch[subImgIndex+0];
				//buffer[bufferIndex+0] = (buffer[bufferIndex+1] * (1.0 - allRecvImgData[patchId].imagePatch[subImgIndex+3])) + allRecvImgData[patchId].imagePatch[subImgIndex+1];
				//buffer[bufferIndex+0] = (buffer[bufferIndex+2] * (1.0 - allRecvImgData[patchId].imagePatch[subImgIndex+3])) + allRecvImgData[patchId].imagePatch[subImgIndex+2];
				//buffer[bufferIndex+3] = buffer[bufferIndex+3]  *(1.0 - allRecvImgData[patchId].imagePatch[subImgIndex+3]);
			
			}
		}
	}

	//copy to main buffer
	for (int i=0; i< imgBufferHeight; i++)
		for (int j=0; j<imgBufferWidth; j++){
			int bufferIndex = (imgBufferWidth*4*i) + (j*4);
			int wholeImgIndex = (imgBufferWidth*3*i) + (j*3);

			wholeImage[wholeImgIndex+0] = (buffer[bufferIndex+0]*255);
			wholeImage[wholeImgIndex+1] = (buffer[bufferIndex+1]*255);
			wholeImage[wholeImgIndex+2] = (buffer[bufferIndex+2]*255);

			//printf("\n i: %d, j: %d,   bufferIndex: %d   wholeImgIndex: %d  - rgb %.2f  %.2f  %.2f    -  rgb %d %d %d",i,j,bufferIndex,wholeImgIndex, 
			//	buffer[bufferIndex+0]*buffer[bufferIndex+3]*255, buffer[bufferIndex+1]*buffer[bufferIndex+3]*255, buffer[bufferIndex+2]*buffer[bufferIndex+3]*255,
			//	 					wholeImage[wholeImgIndex+0],					 wholeImage[wholeImgIndex+1],						wholeImage[wholeImgIndex+2]);
		}


    std::string imgFilenameFinal = "/home/pascal/Desktop/FinalBuffer.ppm";
    createPpm(buffer, imgBufferWidth, imgBufferHeight, imgFilenameFinal);
  

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