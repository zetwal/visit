
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

#ifdef PARALLEL
MPI_Datatype createImgDataType(const int dataSize){
	MPI_Datatype _img_mpi;
	const int numItems = 8;
	int blockLengths[numItems] = {1, 1, 1, 1, 2, 2, 2, 1};
	MPI_Datatype type[numItems] = { MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_FLOAT};
	MPI_Aint offsets[numItems] = {0, sizeof(int), sizeof(int)*2, sizeof(int)*3, sizeof(int)*4, sizeof(int)*6, sizeof(int)*8, sizeof(int)*10};
	MPI_Type_struct(numItems, blockLengths,  offsets, type, &_img_mpi);
	
	return _img_mpi;
}



MPI_Datatype createImgDataType(){	
	MPI_Datatype _img_mpi;
	const int numItems = 8;
	int blockLengths[numItems] = {1, 1,   1, 1,   2, 2, 2, 1};
	MPI_Datatype type[numItems] = { MPI_INT, MPI_INT,    MPI_INT, MPI_INT,   MPI_INT, MPI_INT, MPI_INT,   MPI_FLOAT};
	MPI_Aint offsets[numItems] = {0, sizeof(int), sizeof(int)*2, sizeof(int)*3, sizeof(int)*4, sizeof(int)*6, sizeof(int)*8, sizeof(int)*10};
	MPI_Type_struct(numItems, blockLengths,  offsets, type, &_img_mpi);
	
	return _img_mpi;
}

#endif

bool sortImgByDepth(imgMetaData const& before, imgMetaData const& after){ return before.avg_z > after.avg_z; }
bool sortImgByDepthIota(iotaMeta const& before, iotaMeta const& after){	return before.avg_z < after.avg_z; }
bool value_comparer(const std::pair<int,int> &before, const std::pair<int,int> &after){ return before.second < after.second; }
bool sortByVecSize(const std::vector<iotaMeta> &before, const std::vector<iotaMeta> &after){return before.size() > after.size();}












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
    numPatchesToCompose = 0;

    processorPatchesCount = NULL;
	imgBuffer = NULL;

	patchesToSendArray = NULL;
	patchesToRecvArray = NULL;
	numPatchesToSendArray = NULL;
	numPatchesToRecvArray = NULL;
	recvDisplacementForProcs = NULL;
	sendDisplacementForProcs = NULL;
	numPatchesToSendRecvArray = NULL;
	boundsPerBlockArray = NULL;
	blockDisplacementForProcs = NULL;
	numBlocksPerProc = NULL;

	all_avgZ_proc0.clear();
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

		if (allRecvIotaMeta != NULL)
			delete []allRecvIotaMeta;

		if (imgBuffer != NULL)
			delete []imgBuffer;
	}

	all_avgZ_proc0.clear();
	//boundsPerBlockVec.clear();
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
		// Only proc 0 receives data
		if (my_id == 0)		
			tempRecvBuffer = new int[num_procs*2];

		MPI_Gather(patchesProc, 2, MPI_INT,   tempRecvBuffer, 2,MPI_INT,         0, MPI_COMM_WORLD);		// all send to proc 0

		if (my_id == 0){		
			for (int i=0; i<num_procs; i++){
				processorPatchesCount[tempRecvBuffer[i*2]] = tempRecvBuffer[i*2 + 1];	// enter the number of patches for each processor
				totalPatches += processorPatchesCount[tempRecvBuffer[i*2]];				// count the number of patches

				//std::cout << "!!! Recv: " << tempRecvBuffer[i*2] << "   Patch: " << processorPatchesCount[tempRecvBuffer[i*2]] << std::endl;
			}

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

		MPI_Gatherv(allIotaMetadata, arraySize, MPI_FLOAT,   tempRecvBuffer, recvSizePerProc, offsetBuffer,MPI_FLOAT,         0, MPI_COMM_WORLD);		// all send to proc 0

		if (my_id == 0){
			// Allocate space to receive the many patches
			allRecvIotaMeta = new iotaMeta[totalPatches];

			iotaMeta tempPatch;

			for (int i=0; i<totalPatches; i++){
				tempPatch.procId = 		(int) tempRecvBuffer[i*4 + 0];
				tempPatch.patchNumber = (int) tempRecvBuffer[i*4 + 1];
				tempPatch.imgArea = 	(int) tempRecvBuffer[i*4 + 2];
				tempPatch.avg_z = 			  tempRecvBuffer[i*4 + 3];

				int patchIndex = getDataPatchID(tempPatch.procId, tempPatch.patchNumber);
				allRecvIotaMeta[patchIndex] = setIota(tempPatch.procId, tempPatch.patchNumber,tempPatch.imgArea, tempPatch.avg_z);
				all_avgZ_proc0.insert(tempPatch.avg_z); //insert avg_zs into the set to keep a count of the total number of avg_zs
			}

			if (recvSizePerProc != NULL)
				delete []recvSizePerProc;
			recvSizePerProc = NULL;

			if (offsetBuffer != NULL)
				delete []offsetBuffer;
			offsetBuffer = NULL;

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
int calculatePatchDivision (const std::vector<iotaMeta>& imgVector, std::vector<int>& procsAlreadyInList){
std::map<int,int> numPatchesPerProc;
	std::pair<std::map<int,int>::iterator, bool> isPresent;

	if (imgVector.size() == 0) return 0;

	for (int i = 0; i < imgVector.size(); ++i){
		const bool is_in = ((std::find(procsAlreadyInList.begin(), procsAlreadyInList.end(), imgVector[i].procId)) != (procsAlreadyInList.end()));
		if(!is_in){
			isPresent = numPatchesPerProc.insert (	std::pair<int,int>(imgVector[i].procId, imgVector[i].imgArea) );
			if(isPresent.second == false)			
				numPatchesPerProc[imgVector[i].procId] += imgVector[i].imgArea;
		}
	}

	if(numPatchesPerProc.size() == 0){
		std::cout << "---------->Proc Repeated! ";
		for (int i = 0; i < imgVector.size(); ++i){
			isPresent = numPatchesPerProc.insert (  std::pair<int,int>(imgVector[i].procId, imgVector[i].imgArea) );
			if(isPresent.second == false)			
				numPatchesPerProc[imgVector[i].procId] += imgVector[i].imgArea;
		}
	}

	return std::max_element(numPatchesPerProc.begin(), numPatchesPerProc.end(), value_comparer)->first;
}



void avtImgCommunicator::patchAllocationLogic(){

	// 1. do not send the numpatchespreproc vector!
	// 2. do away with numcompositedPatches test and replace with totalRecvPatches

	std::vector<std::vector<iotaMeta> > all_patches_sorted_avgZ_proc0(num_procs); 
	std::vector<int> procToSend			(num_procs);
	std::vector<int> numAvgZEachBlock	(num_procs);

	numBlocksPerProc = new int[num_procs]();
	boundsPerBlockVec.resize(num_procs);

	int num_avgZ_proc0 = all_avgZ_proc0.size();

	// calculate the range of avg_zs
	int num_avgZ_perBlock = (num_avgZ_proc0 / num_procs);
	int rem_avgZ = num_avgZ_proc0 % num_procs;

	printf("num_avg_z: %d num_procs: %d num_avgZ_perBlock: %d rem_avgZ: %d \n\n", num_avgZ_proc0, num_procs, num_avgZ_perBlock, rem_avgZ);
	
	//procToSend.resize 					(num_procs);
	//all_patches_sorted_avgZ_proc0.resize(num_procs);

	for(int i = 0; i < num_procs; i++){
		numAvgZEachBlock[i] = num_avgZ_perBlock + (rem_avgZ-- > 0 ? 1 : 0);
		std::cout << "i: " << i << " numAvgZEachBlock: " << numAvgZEachBlock[i] << std::endl;
	}

	// Sorting the patches
	std::sort(allRecvIotaMeta, allRecvIotaMeta + totalPatches, &sortImgByDepthIota);

	// Populate the all_patches_sorted_avgZ_proc0 vector with data from allPatchesMeta
	int current_avgZ_Block = 0;
	float prevAvgZ = totalPatches > 0 ? allRecvIotaMeta[0].avg_z : 0;
	int num_avgZ_thisBlock = 1;

	for (int i = 0; i < totalPatches; ++i){
		float currentAvgZ = allRecvIotaMeta[i].avg_z;

		if (currentAvgZ != prevAvgZ && num_avgZ_thisBlock == numAvgZEachBlock[current_avgZ_Block]){
			num_avgZ_thisBlock = 1;
			current_avgZ_Block++;
			prevAvgZ = currentAvgZ;
		} else if (currentAvgZ != prevAvgZ){
			prevAvgZ = currentAvgZ;
			num_avgZ_thisBlock++;
		}

		all_patches_sorted_avgZ_proc0[current_avgZ_Block].push_back(allRecvIotaMeta[i]);
	}

	//sort the avg_z vector by size
	std::sort (all_patches_sorted_avgZ_proc0.begin(), all_patches_sorted_avgZ_proc0.end(), sortByVecSize);

	// Calculate which block of avg_z to send to which processor
	// ith block is sent to the processor in procToSend[i]

	//printf("num_avg_z: %d num_procs: %d\n\n", num_avgZ_proc0, num_procs);
	for (int i = 0; i < num_procs; ++i){
		procToSend[i] = calculatePatchDivision(all_patches_sorted_avgZ_proc0[i], procToSend);

		numBlocksPerProc[procToSend[i]] += 2; // 2 for lower and upper bound

		int size = all_patches_sorted_avgZ_proc0[i].size();
		if(size >= 1)	boundsPerBlockVec[procToSend[i]].push_back(all_patches_sorted_avgZ_proc0[i][0].avg_z);
		if(size >= 2) 	boundsPerBlockVec[procToSend[i]].push_back(all_patches_sorted_avgZ_proc0[i][size-1].avg_z);
		printf("division: %d, procToSend: %d \n", i, procToSend[i]);
	}

	// Printing for error check
	for(int i = 0; i < num_procs; ++i){
		printf("Block %d\n", i );
		for(std::vector<iotaMeta>::iterator it = all_patches_sorted_avgZ_proc0[i].begin(); it != all_patches_sorted_avgZ_proc0[i].end(); ++it){
			printf("\t %.5f\n", it->avg_z);
		}
	}


	std::vector< std::vector<int> >	 patchesToSendVec (num_procs);
	std::vector< std::map<int,int> > patchesToRecvMap (num_procs); 
	std::pair< std::map<int,int>::iterator, bool > isPresent;

	for (int procId = 0; procId < num_procs; procId++)	patchesToRecvMap[procId].clear();
	
	// Populate recvMap and sendVec
	for (int i = 0; i < all_patches_sorted_avgZ_proc0.size(); i++){

		int destProcId = procToSend[i];
	 	for (int j = 0; j < all_patches_sorted_avgZ_proc0[i].size(); j++){

	 		int originProcId = all_patches_sorted_avgZ_proc0[i][j].procId;
	 		if(originProcId != destProcId){

		 		// Array of the form [0 1 2 3 ...]
		 		// 					 even numbers(0,2..): patchNumber 
		 		//					 odd numbers(1,3...): destProcId
	 		    patchesToSendVec[originProcId].push_back( all_patches_sorted_avgZ_proc0[i][j].patchNumber );
	 			patchesToSendVec[originProcId].push_back( destProcId );

	 			// Array of the form [0 1 2 3 ...]
		 		// 					 even numbers(0,2..): procId 
		 		//					 odd numbers(1,3...): numPatches
					isPresent = patchesToRecvMap[destProcId].insert( std::pair<int,int>(originProcId, 1) );
					if(isPresent.second == false)		++patchesToRecvMap[destProcId][originProcId];
					
					//printf("Inserted to dest: %d the origin: %d avg_z: %.2f\n", destProcId, originProcId, allPatchesMeta[patchIndex].avg_z );
	 		}
	 	}

	 	//std::cout << "------division:" << i << " procToSend:" << procToSend[i] << std::endl;
	}


	recvDisplacementForProcs 	= new int[num_procs]();
	sendDisplacementForProcs 	= new int[num_procs]();
	numPatchesToSendArray 		= new int[num_procs]();
	numPatchesToRecvArray 		= new int[num_procs]();
	blockDisplacementForProcs 	= new int[num_procs]();
	numPatchesToSendRecvArray 	= new int[3*num_procs]();

	int totalRecvSize = 0;
	int totalSendSize = 0;
	int totalBlockSize = 0;

	for (int currentProcId = 0; currentProcId < num_procs; currentProcId++){

		numPatchesToSendArray[currentProcId] = patchesToSendVec[currentProcId].size();
		numPatchesToRecvArray[currentProcId] = 2*patchesToRecvMap[currentProcId].size();

		numPatchesToSendRecvArray[currentProcId*3] 		= numPatchesToSendArray[currentProcId];
		numPatchesToSendRecvArray[currentProcId*3 + 1] 	= numPatchesToRecvArray[currentProcId];
		numPatchesToSendRecvArray[currentProcId*3 + 2]	= numBlocksPerProc[currentProcId];

		totalRecvSize += numPatchesToRecvArray[currentProcId];
		totalSendSize += numPatchesToSendArray[currentProcId];
		totalBlockSize+= numBlocksPerProc[currentProcId];

		if(currentProcId!=0) {
			recvDisplacementForProcs[currentProcId] = recvDisplacementForProcs[currentProcId-1] + numPatchesToRecvArray[currentProcId-1];
			sendDisplacementForProcs[currentProcId] = sendDisplacementForProcs[currentProcId-1] + numPatchesToSendArray[currentProcId-1];
			blockDisplacementForProcs[currentProcId] = blockDisplacementForProcs[currentProcId-1] + numBlocksPerProc[currentProcId-1];
		}

	}

	patchesToRecvArray 	= new int[totalRecvSize];
	patchesToSendArray 	= new int[totalSendSize];
	boundsPerBlockArray = new float[totalBlockSize];

	for (int currentProcId = 0; currentProcId < num_procs; currentProcId++){

		//printf("\n\n### dest: %d size: %ld\n ", currentProcId, patchesToRecvMap[currentProcId].size());
		int count = 0;
		for (std::map<int,int>::iterator it = patchesToRecvMap[currentProcId].begin(); it!=patchesToRecvMap[currentProcId].end(); ++it){
				patchesToRecvArray[recvDisplacementForProcs[currentProcId] + (count++)] = it->first;
				patchesToRecvArray[recvDisplacementForProcs[currentProcId] + (count++)] = it->second;
		}

		count = 0;
		for (std::vector<int>::iterator it = patchesToSendVec[currentProcId].begin(); it!=patchesToSendVec[currentProcId].end(); ++it)
			patchesToSendArray[sendDisplacementForProcs[currentProcId] + (count++)] = *it;

		count = 0;
		for (std::vector<float>::iterator it = boundsPerBlockVec[currentProcId].begin(); it!=boundsPerBlockVec[currentProcId].end(); ++it)
			boundsPerBlockArray[blockDisplacementForProcs[currentProcId] + (count++)] = *it;

	}

	//printf("\n ..................................................\n");

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

void avtImgCommunicator::scatterNumDataToCompose(int &totalSendData, int &totalRecvData, int &numDivisions){
	int totalSendRecvData[3];
	#ifdef PARALLEL
		MPI_Scatter(numPatchesToSendRecvArray, 3, MPI_INT, totalSendRecvData, 3, MPI_INT, 0, MPI_COMM_WORLD);
	#endif

	totalSendData = totalSendRecvData[0];
	totalRecvData = totalSendRecvData[1];
	numDivisions = totalSendRecvData[2]; // 2 times num divisions

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
void avtImgCommunicator::scatterDataToCompose(int &totalSendData, int *informationToSendArray, int &totalRecvData, int *informationToRecvArray, int &numDivisions, float *blocksPerProc){
	#ifdef PARALLEL
        MPI_Scatterv(patchesToSendArray, numPatchesToSendArray, sendDisplacementForProcs, MPI_INT, informationToSendArray, totalSendData, MPI_INT, 0,MPI_COMM_WORLD);
        MPI_Scatterv(patchesToRecvArray, numPatchesToRecvArray, recvDisplacementForProcs, MPI_INT, informationToRecvArray, totalRecvData, MPI_INT, 0,MPI_COMM_WORLD);
        MPI_Scatterv(boundsPerBlockArray, numBlocksPerProc, blockDisplacementForProcs, MPI_FLOAT, blocksPerProc, numDivisions, MPI_FLOAT, 0,MPI_COMM_WORLD);
	#endif

    if(my_id == 0){
        if(patchesToSendArray!=NULL) delete[] patchesToSendArray;	patchesToSendArray = NULL;
        if(patchesToRecvArray!=NULL) delete[] patchesToRecvArray;
        if(numPatchesToSendArray!=NULL) delete[] numPatchesToSendArray;
        if(numPatchesToRecvArray!=NULL) delete[] numPatchesToRecvArray;
        if(recvDisplacementForProcs!=NULL) delete[] recvDisplacementForProcs;
        if(sendDisplacementForProcs!=NULL) delete[] sendDisplacementForProcs;
        if(blockDisplacementForProcs!=NULL) delete[] blockDisplacementForProcs;
        if(numPatchesToSendRecvArray!=NULL) delete[] numPatchesToSendRecvArray;
        if(boundsPerBlockArray!=NULL) delete[] boundsPerBlockArray;
        if(numBlocksPerProc!=NULL) delete[] numBlocksPerProc;
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
void avtImgCommunicator::sendPointToPoint(imgMetaData toSendMetaData, imgData toSendImgData, int tag){
	#ifdef PARALLEL
		// Commit the datatype
		MPI_Datatype _TestImg_mpi;
		_TestImg_mpi = createImgDataType();
		MPI_Type_commit(&_TestImg_mpi);

	    MPI_Send(&toSendMetaData, 1, _TestImg_mpi, toSendMetaData.destProcId, tag, MPI_COMM_WORLD); // 2
		MPI_Send(toSendImgData.imagePatch, toSendMetaData.dims[0]*toSendMetaData.dims[1]*4, MPI_FLOAT, toSendMetaData.destProcId, tag-1, MPI_COMM_WORLD); //send the image data // 1
   
   		MPI_Type_free(&_TestImg_mpi);
    #endif
}



void avtImgCommunicator::recvPointToPoint(imgMetaData &recvMetaData, imgData &recvImgData){
	#ifdef PARALLEL
		// Commit the datatype
		MPI_Datatype _TestImg_mpi;
		_TestImg_mpi = createImgDataType();
		MPI_Type_commit(&_TestImg_mpi);

        MPI_Recv (&recvMetaData, 1, _TestImg_mpi, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, &status);

        recvImgData.procId = recvMetaData.procId;
        recvImgData.patchNumber = recvMetaData.patchNumber;
        recvImgData.imagePatch = new float[recvMetaData.dims[0]*recvMetaData.dims[1] * 4];

        MPI_Recv(recvImgData.imagePatch, recvMetaData.dims[0]*recvMetaData.dims[1]*4, MPI_FLOAT, status.MPI_SOURCE, 1, MPI_COMM_WORLD, &status);

        MPI_Type_free(&_TestImg_mpi);
    #endif
}

void avtImgCommunicator::recvPointToPointMetaData(imgMetaData &recvMetaData, int tag){
	#ifdef PARALLEL
		// Commit the datatype
		MPI_Datatype _TestImg_mpi;
		_TestImg_mpi = createImgDataType();
		MPI_Type_commit(&_TestImg_mpi);

        MPI_Recv (&recvMetaData, 1, _TestImg_mpi, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status); // 2

        MPI_Type_free(&_TestImg_mpi);
    #endif
}


void avtImgCommunicator::recvPointToPointImgData(imgMetaData recvMetaData, imgData &recvImgData, int tag){
	#ifdef PARALLEL
        MPI_Recv(recvImgData.imagePatch, recvMetaData.dims[0]*recvMetaData.dims[1]*4, MPI_FLOAT, status.MPI_SOURCE, tag-1, MPI_COMM_WORLD, &status); // 1
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
void avtImgCommunicator::gatherAndAssembleImages(int sizex, int sizey, float *images, int numDivisions){
	#ifdef PARALLEL
		float *tempRecvBuffer = NULL;
		int *recvSizePerProc = NULL;
		int *offsetBuffer = NULL;
		int totalDivisions = 0;
		std::map<float,int> depthPartitions;	//float: z-value, int: division index

		// Only proc 0 receives data
		if (my_id == 0){
			int divIndex = 0;
			recvSizePerProc = new int[num_procs];	
			offsetBuffer = new int[num_procs];	

			for (int i=0; i<num_procs; i++){
				int numBoundsPerBlock = boundsPerBlockVec[i].size()/2;
				totalDivisions += numBoundsPerBlock;
				std::cout << i << " ~ "<< "boundsPerBlockVec[i].size(): " << numBoundsPerBlock << std::endl;
				recvSizePerProc[i] = numBoundsPerBlock*(sizex*sizey*4);

				for (int j=0; j<boundsPerBlockVec[i].size(); j+=2){
					depthPartitions.insert( std::pair< float,int > ( (boundsPerBlockVec[i][j] + boundsPerBlockVec[i][j+1]) * 0.5, divIndex ) );					
					divIndex++;
				}

				if (i == 0)
					offsetBuffer[i] = 0;
				else
					offsetBuffer[i] = offsetBuffer[i-1] + recvSizePerProc[i-1];
			}

			totalDivisions = totalDivisions;
			tempRecvBuffer = new float[ ( (sizex*sizey*4)*totalDivisions) ];

			//std::cout << "divIndex: " << divIndex << "   totalDivisions: " << totalDivisions << std::endl;	
		}

		//std::cout << my_id << " ~ numDivisions: " << numDivisions << std::endl;

        //  send   recv  others
		MPI_Gatherv(images, ((sizex*sizey*4)*numDivisions), MPI_FLOAT,    tempRecvBuffer, recvSizePerProc, offsetBuffer,MPI_FLOAT,        0, MPI_COMM_WORLD);		// all send to proc 0

		if (my_id == 0){
			// Create a buffer to store the composited image
			imgBuffer = new float[sizex * sizey * 4];
			
			// Front to back
			for (int i=0; i<(sizex * sizey * 4); i+=4){
				imgBuffer[i+0] = background[0]/255.0;	
				imgBuffer[i+1] = background[1]/255.0; 
				imgBuffer[i+2] = background[2]/255.0; 
				imgBuffer[i+3] = 1.0;
			}

			std::map< float,int >::iterator it;
			it=depthPartitions.end();
			int count = 0;
			do{
				--it;
				std::cout << it->first << " => " << it->second << '\n';
				imageBuffer temp;
				temp.image = new float[sizex*sizey*4];
				memcpy(temp.image, &tempRecvBuffer[(count*(sizex*sizey*4))], (sizex*sizey*4)*sizeof(float) );

				//std::string imgFilenameFinal = "/home/pascal/Desktop/Gathered_on_0_" + NumbToString(it->second) + "_Buffer.ppm";
	        	//createPpm(temp.image, sizex, sizey, imgFilenameFinal);

	        	for (int j=0; j<sizey; j++){
					for (int k=0; k<sizex; k++){
						int imgIndex = sizex*4*j + k*4;										// index in the image 

						// Back to front compositing
						imgBuffer[imgIndex+0] = (imgBuffer[imgIndex+0] * (1.0 - temp.image[imgIndex+3])) + temp.image[imgIndex+0];
						imgBuffer[imgIndex+1] = (imgBuffer[imgIndex+1] * (1.0 - temp.image[imgIndex+3])) + temp.image[imgIndex+1];
						imgBuffer[imgIndex+2] = (imgBuffer[imgIndex+2] * (1.0 - temp.image[imgIndex+3])) + temp.image[imgIndex+2];
					}
				}

				//std::string imgFilename_Final = "/home/pascal/Desktop/_"+ NumbToString(count) + "_Numfinal.ppm";
	        	//createPpm(imgBuffer, sizex, sizey, imgFilename_Final);

	        	delete []temp.image;
	        	count++;
			}while( it!=depthPartitions.begin());
			
			//std::string imgFilename_Final = "/home/pascal/Desktop/_final.ppm";
	        //createPpm(imgBuffer, sizex, sizey, imgFilename_Final);

			if (tempRecvBuffer != NULL)
				delete []tempRecvBuffer;
			tempRecvBuffer = NULL;
		}

		syncAllProcs();

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
void avtImgCommunicator::getcompositedImage(int imgBufferWidth, int imgBufferHeight, unsigned char *wholeImage){

	for (int i=0; i< imgBufferHeight; i++)
		for (int j=0; j<imgBufferWidth; j++){
			int bufferIndex = (imgBufferWidth*4*i) + (j*4);
			int wholeImgIndex = (imgBufferWidth*3*i) + (j*3);

			wholeImage[wholeImgIndex+0] = (imgBuffer[bufferIndex+0] ) * 255;
			wholeImage[wholeImgIndex+1] = (imgBuffer[bufferIndex+1] ) * 255;
			wholeImage[wholeImgIndex+2] = (imgBuffer[bufferIndex+2] ) * 255;
		}

	if (imgBuffer != NULL)
		delete []imgBuffer;
	imgBuffer = NULL;
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
	temp.destProcId = _procId;
	temp.patchNumber = _patchNumber;
	temp.dims[0] = dim_x;					temp.dims[1] = dim_y;
	temp.screen_ll[0] = screen_ll_x;		temp.screen_ll[1] = screen_ll_y;
	temp.screen_ur[0] = screen_ur_x;		temp.screen_ur[1] = screen_ur_y;
	temp.avg_z = _avg_z;
	
	return temp;
}


iotaMeta avtImgCommunicator::setIota(int _procId, int _patchNumber, float _imgArea, float _avg_z){
	iotaMeta temp;
	temp.procId = _procId;
	temp.patchNumber = _patchNumber;
	temp.imgArea = _imgArea;
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

void createPpmWithOffset(float array[], int dimx, int dimy, std::string filename, int offset){
    int i, j;
    std::cout << "createPpm2  dims: " << dimx << ", " << dimy << " -  " << filename.c_str() << std::endl;
    FILE *fp = fopen(filename.c_str(), "wb"); // b - binary mode 
    (void) fprintf(fp, "P6\n%d %d\n255\n", dimx, dimy);
    for (j = 0; j < dimy; ++j){
        for (i = 0; i < dimx; ++i){
            static unsigned char color[3];
            color[0] = array[offset+ j*(dimx*4) + i*4 + 0] * 255;  // red
            color[1] = array[offset+ j*(dimx*4) + i*4 + 1] * 255;  // green
            color[2] = array[offset+ j*(dimx*4) + i*4 + 2] * 255;  // blue 
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
