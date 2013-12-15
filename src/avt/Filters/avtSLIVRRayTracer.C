/*****************************************************************************
*
* Copyright (c) 2000 - 2013, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* LLNL-CODE-442911
* All rights reserved.
*
* This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
* full copyright notice is contained in the file COPYRIGHT located at the root
* of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
*
* Redistribution  and  use  in  source  and  binary  forms,  with  or  without
* modification, are permitted provided that the following conditions are met:
*
*  - Redistributions of  source code must  retain the above  copyright notice,
*    this list of conditions and the disclaimer below.
*  - Redistributions in binary form must reproduce the above copyright notice,
*    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
*    documentation and/or other materials provided with the distribution.
*  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
*    be used to endorse or promote products derived from this software without
*    specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
* AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
* IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
* ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
* LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
* DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
* SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
* CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
* LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
* OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
* DAMAGE.
*
*****************************************************************************/


// ************************************************************************* //
//                           avtSLIVRRayTracer.C                              //
// ************************************************************************* //

#include <avtSLIVRRayTracer.h>


#include <avtImage.h>
#include <avtSamplePoints.h>
#include <avtVolume.h>

#include <ImproperUseException.h>

avtSLIVRRayTracer::avtSLIVRRayTracer(){
	std::cout << "avtSLIVRRayTracer::avtSLIVRRayTracer()" << std::endl;

    proc = 0;
    numProcs = 1;
}

avtSLIVRRayTracer::~avtSLIVRRayTracer(){
	std::cout << "avtSLIVRRayTracer::~avtSLIVRRayTracer()" << std::endl;
}


void
avtSLIVRRayTracer::Execute(void){
	std::cout << "avtSLIVRRayTracer::Execute()... " << std::endl;
	int height, width, fullWidth, fullHeight;
	height = width = 100;
	fullWidth = fullHeight = 500;
	vtkImageData *image = avtImageRepresentation::NewImage(width, height);

    std::cout << "@#@numPatches: " << numPatches << std::endl;

    //
    // Populate an initial image, either with the background or with an
    // opaque image that is to be inserted into the middle of the rendering.
    //
    unsigned char *data = (unsigned char *)image->GetScalarPointer(0, 0, 0);
    int nPixels = width*height;
    double *zbuffer = new double[nPixels];
    for (int i = 0 ; i < nPixels ; i++)
    {
        zbuffer[i] = 1.;
    }

    //
    // Draw the initial background into the image.
    //
    //int fullHeight = volume->GetVolumeHeight();
    //int fullWidth  = volume->GetVolumeWidth();

    vtkImageData *fullImage = avtImageRepresentation::NewImage(fullWidth, fullHeight);
    unsigned char *fulldata = (unsigned char *) fullImage->GetScalarPointer(0, 0, 0);

    //
    // Tell our output what its new image is.
    //
    SetOutputImage(image);

    //
    // Clean up memory.
    //
    image->Delete();
    fullImage->Delete();
    delete [] zbuffer;

	std::cout << "... avtSLIVRRayTracer::Execute()  End!!!" << std::endl;
}

/*
void avtSLIVRRayTracer::compositing(void){
    bool parallelOn = false;
    if (numProcs > 1)
        parallelOn = true;
       
    //
    // Single Processor
    //
    if (parallelOn == false){
        
        //
        // Get the metadata
        //
        std::vector<imgMetaData> allImgMetaData;          // contains the metadata to composite the image
        int numPatches = extractor.getImgPatchSize();     // get the number of patches
    
        for (int i=0; i<numPatches; i++){
            imgMetaData temp;
            temp = extractor.getImgMetaPatch(i);
            allImgMetaData.push_back(temp);
        }
        debug5 << "Number of patches: " << numPatches << std::endl;

        //
        // Sort with the largest z first
        //
        std::sort(allImgMetaData.begin(), allImgMetaData.end(), &sortImgMetaDataByDepth);

        //
        // Composite the images
        //
        
        // Creates a buffer to store the composited image
        float *composedData = new float[screen[0] * screen[1] * 4];

        for (int i=0; i<(screen[0] * screen[1] * 4); i+=4){
            composedData[i+0] = background[0]/255.0;   
            composedData[i+1] = background[1]/255.0; 
            composedData[i+2] = background[2]/255.0; 
            composedData[i+3] = 1.0;
        }

        for (int i=0; i<numPatches; i++){
            imgMetaData currentPatch = allImgMetaData[i];

            imgData tempImgData;
            tempImgData.imagePatch = new float[currentPatch.dims[0] * currentPatch.dims[1] * 4];
            extractor.getnDelImgData(currentPatch.patchNumber, tempImgData);

            for (int j=0; j<currentPatch.dims[1]; j++){
                for (int k=0; k<currentPatch.dims[0]; k++){

                    int startingX = currentPatch.screen_ll[0];
                    int startingY = currentPatch.screen_ll[1]; 

                    if ((startingX + k) > screen[0])
                        continue;

                    if ((startingY + j) > screen[1])
                        continue;
                    
                    int subImgIndex = currentPatch.dims[0]*j*4 + k*4;                                   // index in the subimage 
                    int bufferIndex = (startingY*screen[0]*4 + j*screen[0]*4) + (startingX*4 + k*4);    // index in the big buffer

                    // back to Front compositing: composited_i = composited_i-1 * (1.0 - alpha_i) + incoming; alpha = alpha_i-1 * (1- alpha_i)
                    composedData[bufferIndex+0] = imgComm.clamp((composedData[bufferIndex+0] * (1.0 - tempImgData.imagePatch[subImgIndex+3])) + tempImgData.imagePatch[subImgIndex+0]);
                    composedData[bufferIndex+1] = imgComm.clamp((composedData[bufferIndex+1] * (1.0 - tempImgData.imagePatch[subImgIndex+3])) + tempImgData.imagePatch[subImgIndex+1]);
                    composedData[bufferIndex+2] = imgComm.clamp((composedData[bufferIndex+2] * (1.0 - tempImgData.imagePatch[subImgIndex+3])) + tempImgData.imagePatch[subImgIndex+2]);
                    composedData[bufferIndex+3] = imgComm.clamp((composedData[bufferIndex+3] * (1.0 - tempImgData.imagePatch[subImgIndex+3])) + tempImgData.imagePatch[subImgIndex+3]);
                }
            }

            if (tempImgData.imagePatch != NULL)
                delete []tempImgData.imagePatch;
            tempImgData.imagePatch = NULL;
        }
        allImgMetaData.clear();


        // Creates an image structure to hold the image
        avtImage_p whole_image;
        whole_image = new avtImage(this);

        float *zbuffer = new float[screen[0] * screen[1]];
        unsigned char *imgTest = NULL;

        vtkImageData *img = avtImageRepresentation::NewImage(screen[0], screen[1]);
        whole_image->GetImage() = img;

        imgTest = new unsigned char[screen[0] * screen[1] * 3];
        imgTest = whole_image->GetImage().GetRGBBuffer();

        zbuffer = new float[screen[0] * screen[1]]();
        for (int s=0; s<screen[0] * screen[1]; s++)
            zbuffer[s] = 20.0;
        zbuffer = whole_image->GetImage().GetZBuffer();


        // Get the composited image
        for (int i=0; i< screen[1]; i++)
            for (int j=0; j<screen[0]; j++){
                int bufferIndex = (screen[0]*4*i) + (j*4);
                int wholeImgIndex = (screen[0]*3*i) + (j*3);

                imgTest[wholeImgIndex+0] = (composedData[bufferIndex+0] ) * 255;
                imgTest[wholeImgIndex+1] = (composedData[bufferIndex+1] ) * 255;
                imgTest[wholeImgIndex+2] = (composedData[bufferIndex+2] ) * 255;
            }

        img->Delete();

        if (zbuffer != NULL)
            delete []zbuffer;

        SetOutput(whole_image);
        return;
    }


        //
        // Parallel
        //

        //
        // -- -- -- Timing -- 

        // imgComm.syncAllProcs();     // only required for time testing purposes

        visitTimer->StopTimer(timingVolToImg, "VolToImg");
        visitTimer->DumpTimings();
        
        int  timingComm = visitTimer->StartTimer();
        int  timingCommMeta = visitTimer->StartTimer();

       
        //
        // Getting the patches & send/receive the number of patches that each has to 0
        //
        int numPatches = extractor.getImgPatchSize();     // get the number of patches
        imgComm.gatherNumPatches(numPatches);

        debug5 << PAR_Rank() << "   avtRayTracer::Execute  - Getting the patches -    numPatches: " << numPatches << "   total assigned: " << extractor.getTotalAssignedPatches() << endl;


        //
        // Send/Receive the patches iota metadata to proc 0
        //

        float *tempSendBuffer = NULL;
        tempSendBuffer = new float[numPatches*7];

        std::multimap<int, imgMetaData> imgMetaDataMultiMap;
        imgMetaDataMultiMap.clear();
        for (int i=0; i<numPatches; i++){
            imgMetaData temp;
            temp = extractor.getImgMetaPatch(i);
            imgMetaDataMultiMap.insert(  std::pair<int, imgMetaData>   (temp.patchNumber, temp));

            tempSendBuffer[i*7 + 0] = temp.procId;
            tempSendBuffer[i*7 + 1] = temp.patchNumber;
            tempSendBuffer[i*7 + 2] = temp.dims[0];
            tempSendBuffer[i*7 + 3] = temp.dims[1];
            tempSendBuffer[i*7 + 4] = temp.screen_ll[0];
            tempSendBuffer[i*7 + 5] = temp.screen_ll[1];
            tempSendBuffer[i*7 + 6] = temp.avg_z;

            debug5 << PAR_Rank() << "   ~ has patch #: " << temp.patchNumber << "   avg_z: " << temp.avg_z << endl;
        }

        imgComm.gatherIotaMetaData(numPatches*7, tempSendBuffer); 

        delete []tempSendBuffer;
        tempSendBuffer = NULL;


        //
        // -- -- Timing --
        visitTimer->StopTimer(timingCommMeta, "Communicating metadata");
        visitTimer->DumpTimings();
        debug5 << PAR_Rank() << "   avtRayTracer::Execute  - Send/Receive the patches iota " << endl;

        int  timingCommDecision = visitTimer->StartTimer();


        //
        // Patch allocation Logic & inform other procs about logic
        // 
        if (PAR_Rank() == 0)
            imgComm.patchAllocationLogic();    

        debug5 << PAR_Rank() << "   avtRayTracer::Execute  - imgComm.patchAllocationLogic() " << endl;

        // informationToRecvArray:   (procId, numPatches)           (procId, numPatches)         (procId, numPatches)  ...
        // informationToSendArray:   (patchNumber, destProcId)     (patchNumber, destProcId)    (patchNumber, destProcId)  ...
        int totalSendData, totalRecvData, numZDivisions, totalPatchesToCompositeLocally;
        int *informationToRecvArray, *informationToSendArray, *patchesToCompositeLocally;
        float *divisionsArray = NULL;
        informationToRecvArray = informationToSendArray = patchesToCompositeLocally = NULL;
        totalSendData = totalRecvData = numZDivisions = totalPatchesToCompositeLocally = 0;

        //
        // Send info about which patch to receive and which patches to send & receive
        imgComm.scatterNumDataToCompose(totalSendData, totalRecvData, numZDivisions, totalPatchesToCompositeLocally);

        debug5 << PAR_Rank() << " ~  num patches to send: " << totalSendData/2 << " num processors to recv from: " << totalRecvData/2 << "    numZDivisions: " << numZDivisions << "   totalPatchesToCompositeLocally: " <<   totalPatchesToCompositeLocally << endl;


        //
        //receive the information about how many patches to receive from other processors
        informationToRecvArray = new int[totalRecvData];
        informationToSendArray = new int[totalSendData];
        divisionsArray = new float[numZDivisions];


        if (totalPatchesToCompositeLocally > 0)
            patchesToCompositeLocally = new int[totalPatchesToCompositeLocally];

        imgComm.scatterDataToCompose(totalSendData, informationToSendArray, totalRecvData, informationToRecvArray, numZDivisions, divisionsArray, totalPatchesToCompositeLocally, patchesToCompositeLocally);
        numZDivisions /= 2;   // halve it as it's the size of the array that contains start n end

        
        //
        // --- Timing -- //
        
        visitTimer->StopTimer(timingCommDecision, "Decision making and informing procs");
        visitTimer->DumpTimings();
        
        int  timingLocalCompositing = visitTimer->StartTimer();
        
        
        //
        // Local compositing if any is required
        // 
        int index = 1;
        int start, end;
        start = end = index;

        std::vector<imgData> compositedDataVec;
        debug5 << PAR_Rank() << " ~  totalPatchesToCompositeLocally: " << totalPatchesToCompositeLocally << endl;

        while (index <= totalPatchesToCompositeLocally){
            if (patchesToCompositeLocally[index] == -1 || index == totalPatchesToCompositeLocally){
                end = index-1;

                int dimsX, dimsY;

                imgMetaData composedPatch = imgMetaDataMultiMap.find(patchesToCompositeLocally[start])->second;
                imgMetaData lastPatch = imgMetaDataMultiMap.find(patchesToCompositeLocally[end])->second;

                imgData composedData;
                composedData.patchNumber    = composedPatch.patchNumber;
                composedData.procId         = composedPatch.procId;

                int imgBufferWidth = abs(lastPatch.screen_ll[0] + lastPatch.dims[0] - composedPatch.screen_ll[0]);
                int imgBufferHeight = abs(lastPatch.screen_ll[1] + lastPatch.dims[1] - composedPatch.screen_ll[1]);

                composedData.imagePatch = new float[((imgBufferWidth * imgBufferHeight * 4)) ]();  //size: imgBufferWidth * imgBufferHeight * 4, initialized to 0

                composedPatch.dims[0]   = imgBufferWidth;
                composedPatch.dims[1]   = imgBufferHeight;
                composedPatch.inUse     = false;


                for (int patchIndex=start; patchIndex<=end; patchIndex++){
                    imgMetaData currentPatch = imgMetaDataMultiMap.find(patchesToCompositeLocally[patchIndex])->second;

                    int startingX = currentPatch.screen_ll[0] - composedPatch.screen_ll[0];
                    int startingY = currentPatch.screen_ll[1] - composedPatch.screen_ll[1]; 

                    imgData tempImgData;
                    tempImgData.imagePatch = new float[currentPatch.dims[0] * currentPatch.dims[1] * 4];
                    extractor.getnDelImgData(currentPatch.patchNumber, tempImgData);

                    // Assemble the idivisions into 1 layer
                    for (int j=0; j<currentPatch.dims[1]; j++){
                        for (int k=0; k<currentPatch.dims[0]; k++){


                            if ( ((startingX + k) < 0) || ((startingX + k) > imgBufferWidth) )
                                continue;

                            if ( ((startingY + j) < 0) || ((startingY + j) > imgBufferHeight) )
                                continue;
                            
                            int subImgIndex = currentPatch.dims[0]*j*4 + k*4;                             // index in the subimage 
                            int bufferIndex = (startingY*imgBufferWidth*4 + j*imgBufferWidth*4) + (startingX*4 + k*4);  // index in the big buffer

                            // back to Front compositing: composited_i = composited_i-1 * (1.0 - alpha_i) + incoming; alpha = alpha_i-1 * (1- alpha_i)
                            composedData.imagePatch[bufferIndex+0] = (composedData.imagePatch[bufferIndex+0] * (1.0 - tempImgData.imagePatch[subImgIndex+3])) + tempImgData.imagePatch[subImgIndex+0];
                            composedData.imagePatch[bufferIndex+1] = (composedData.imagePatch[bufferIndex+1] * (1.0 - tempImgData.imagePatch[subImgIndex+3])) + tempImgData.imagePatch[subImgIndex+1];
                            composedData.imagePatch[bufferIndex+2] = (composedData.imagePatch[bufferIndex+2] * (1.0 - tempImgData.imagePatch[subImgIndex+3])) + tempImgData.imagePatch[subImgIndex+2];
                            composedData.imagePatch[bufferIndex+3] = (composedData.imagePatch[bufferIndex+3] * (1.0 - tempImgData.imagePatch[subImgIndex+3])) + tempImgData.imagePatch[subImgIndex+3];
                        }
                    }

                    if (tempImgData.imagePatch != NULL)
                        delete []tempImgData.imagePatch;
                    tempImgData.imagePatch = NULL;

                    imgMetaDataMultiMap.erase(patchesToCompositeLocally[patchIndex]);
                }
                start = index+1;

                imgMetaDataMultiMap.insert(std::pair<int,imgMetaData>(composedPatch.patchNumber, composedPatch));
                compositedDataVec.push_back(composedData);
            }
            index++;
        }

        //
        // --- Timing -- //

        visitTimer->StopTimer(timingLocalCompositing, "Local compositing");
        visitTimer->DumpTimings();
        
        int  timingLocalWork = visitTimer->StartTimer();
        

        //
        // Set the destination for each patch
        //

        for (int k = 0; k < totalSendData; k+=2){
            std::multimap<int,imgMetaData>::iterator it = imgMetaDataMultiMap.find(informationToSendArray[k]);
            (it->second).destProcId = informationToSendArray[k+1];
        }

        //
        // storage initialization
        //        
        int remainingPatches = 0;
        std::vector<imgMetaData> allImgMetaData;                        // to contain the metadata to composite
        std::multimap< std::pair<int,int>, imgData> imgDataToCompose;   // to contain the data to composite

        //std::multimap<int, imgMetaData> imgMetaDataMultiMap;          // contained all the patches it produced

        allImgMetaData.clear();
        imgDataToCompose.clear();


        //
        // Copying the patches that it will need
        //
        std::multimap<int,imgMetaData>::iterator it;
        for (it = imgMetaDataMultiMap.begin(); it != imgMetaDataMultiMap.end(); ++it ){

            imgMetaData tempImgMetaData = it->second;

            if (tempImgMetaData.destProcId == tempImgMetaData.procId ){
                imgData tempImgData;
                tempImgData.imagePatch = new float[it->second.dims[0] * it->second.dims[1] * 4];
                extractor.getnDelImgData(tempImgMetaData.patchNumber, tempImgData);

                allImgMetaData.push_back(tempImgMetaData);
                imgDataToCompose.insert( std::pair< std::pair<int,int>, imgData> (  std::pair<int,int>(tempImgMetaData.procId, tempImgMetaData.patchNumber), tempImgData)   );

                remainingPatches++;
            }    
        }


        //
        // --- Timing -- //
        visitTimer->StopTimer(timingLocalWork, "Local compositing patches to send");
        visitTimer->DumpTimings();
        
        int  timingSendReceive = visitTimer->StartTimer();
        

        //
        // Sending and receiving from other patches (does a kind of binary swap - half send, half receive and each list gets subdivided)
        // 
        int startingProc = 0;
        int endingProc = PAR_Size() - 1;

        while (startingProc != endingProc){
            int newStartingProc, newEndingProc, numInOtherHalf;
            int *procsInOtherList = NULL;

            int numProcessors = endingProc - startingProc + 1;
            int half = numProcessors/2;
            int middleProc = startingProc + (half-1);
            int doFirst = SEND;     //1: send 2: receive
            
            if (PAR_Rank() <= middleProc){
                numInOtherHalf = numProcessors - half;
                procsInOtherList = new int[numInOtherHalf];
                
                for (int i=0; i<numInOtherHalf; i++)
                    procsInOtherList[i] = startingProc+half+i;
                
                doFirst = SEND;
                newStartingProc = startingProc;
                newEndingProc = middleProc; 
            }
            else{
                numInOtherHalf = half;
                procsInOtherList = new int[numInOtherHalf];
                
                for (int i=0; i<half; i++)
                    procsInOtherList[i] = startingProc+i;
                
                doFirst = RECEIVE;
                newStartingProc = middleProc+1;
                newEndingProc = endingProc;
            }

            if (doFirst == SEND){
                //
                // Send
                std::set<int> senderListSet;     senderListSet.clear();

                for (int i=0; i<numInOtherHalf; i++)
                    for(int j = 0; j < totalSendData; j+=2) 
                        if (informationToSendArray[j+1] == procsInOtherList[i])
                            senderListSet.insert(informationToSendArray[j+1]);
                    
                std::multimap<int,imgMetaData>::iterator it;
                for (it = imgMetaDataMultiMap.begin(); it != imgMetaDataMultiMap.end(); ++it ){
                    imgMetaData tempImgMetaData = it->second;

                    const bool is_in = ((std::find(senderListSet.begin(), senderListSet.end(), it->second.destProcId)) != (senderListSet.end()));
                    if (is_in){
                        imgMetaData tempImgMetaData = it->second;
                        imgData tempImgData;
                        tempImgData.patchNumber = tempImgMetaData.patchNumber;
                        tempImgData.procId = tempImgMetaData.procId;
                        tempImgData.imagePatch = new float[tempImgMetaData.dims[0] * tempImgMetaData.dims[1] * 4];
                        

                        if(tempImgMetaData.inUse)
                            extractor.getnDelImgData(tempImgMetaData.patchNumber, tempImgData);
                        else{
                            const bool is_inC = (std::find(compositedDataVec.begin(), compositedDataVec.end(), tempImgData)) != compositedDataVec.end();  
                            if(is_inC) 
                                tempImgData = *(std::find(compositedDataVec.begin(), compositedDataVec.end(), tempImgData));
                            else 
                                debug5 << PAR_Rank() << " Ray casting: SLIVR uuuuuuuh it didn't find the patch" << endl;
                        }

                        imgComm.sendPointToPoint(tempImgMetaData,tempImgData, numProcessors);

                        if (tempImgData.imagePatch != NULL)
                            delete []tempImgData.imagePatch;
                        tempImgData.imagePatch = NULL;
                    }
                }

                senderListSet.clear();

                //
                // Receive

                //
                // counting how many to receive
                int numToReceive = 0;
                for (int i=0; i<totalRecvData/2; i++)
                    for (int j=0; j<numInOtherHalf; j++)
                        if ((informationToRecvArray[i*2] == procsInOtherList[j]) && (informationToRecvArray[i*2 + 1] > 0))
                            numToReceive+=informationToRecvArray[i*2 + 1];

                
                // does the receive
                for (int i=0; i<numToReceive; i++){
                    imgData tempImgData;
                    imgMetaData tempImgMetaData;

                    imgComm.recvPointToPointMetaData(tempImgMetaData, numProcessors);

                    tempImgData.procId = tempImgMetaData.procId;
                    tempImgData.patchNumber = tempImgMetaData.patchNumber;
                    tempImgData.imagePatch = new float[tempImgMetaData.dims[0]*tempImgMetaData.dims[1] * 4];
                    imgComm.recvPointToPointImgData(tempImgMetaData, tempImgData, numProcessors);

                    allImgMetaData.push_back(tempImgMetaData);
                    imgDataToCompose.insert( std::pair< std::pair<int,int>, imgData> (std::pair<int,int>(tempImgMetaData.procId, tempImgMetaData.patchNumber), tempImgData));
                }
            }else{
                //
                // Receive
        
                // counting how many to receive
                int numToReceive = 0;
                for (int i=0; i<totalRecvData/2; i++)
                  for (int j=0; j<numInOtherHalf; j++)
                        if ((informationToRecvArray[i*2] == procsInOtherList[j]) && (informationToRecvArray[i*2 + 1] > 0))
                            numToReceive+=informationToRecvArray[i*2 + 1];

                // does the receive
                for (int i=0; i<numToReceive; i++){
                    imgData tempImgData;
                    imgMetaData tempImgMetaData;

                    //imgComm.recvPointToPoint(tempImgMetaData, tempImgData);
                    imgComm.recvPointToPointMetaData(tempImgMetaData, numProcessors);

                    tempImgData.procId = tempImgMetaData.procId;
                    tempImgData.patchNumber = tempImgMetaData.patchNumber;
                    tempImgData.imagePatch = new float[tempImgMetaData.dims[0]*tempImgMetaData.dims[1] * 4];
                    imgComm.recvPointToPointImgData(tempImgMetaData, tempImgData, numProcessors);

                    allImgMetaData.push_back(tempImgMetaData);
                    imgDataToCompose.insert( std::pair< std::pair<int,int>, imgData> (std::pair<int,int>(tempImgMetaData.procId, tempImgMetaData.patchNumber), tempImgData));
                }


                //
                // Send
                std::set<int> senderListSet;    senderListSet.clear();

                for (int i=0; i<numInOtherHalf; i++)
                    for(int j = 0; j < totalSendData; j+=2) 
                        if (informationToSendArray[j+1] == procsInOtherList[i])
                            senderListSet.insert(informationToSendArray[j+1]);
                    
                
                std::multimap<int,imgMetaData>::iterator it;
                for (it = imgMetaDataMultiMap.begin(); it != imgMetaDataMultiMap.end(); ++it ){
                    imgMetaData tempImgMetaData = it->second;

                    const bool is_in = ((std::find(senderListSet.begin(), senderListSet.end(), it->second.destProcId)) != (senderListSet.end()));
                    if (is_in){
                        imgMetaData tempImgMetaData = it->second;
                        imgData tempImgData;
                        tempImgData.patchNumber = tempImgMetaData.patchNumber;
                        tempImgData.procId = tempImgMetaData.procId;
                        tempImgData.imagePatch = new float[it->second.dims[0] * it->second.dims[1] * 4];

                        if(tempImgMetaData.inUse)
                            extractor.getnDelImgData(tempImgMetaData.patchNumber, tempImgData);
                        else{
                            const bool is_inC = (std::find(compositedDataVec.begin(), compositedDataVec.end(), tempImgData)) != compositedDataVec.end();  
                            if(is_inC) 
                                tempImgData = *(std::find(compositedDataVec.begin(), compositedDataVec.end(), tempImgData));
                            else 
                                debug5 << PAR_Rank() << "Ray casting: SLIVR uuuuuuuh it didn't find the patch" << endl;
                        }
                        

                        imgComm.sendPointToPoint(tempImgMetaData,tempImgData, numProcessors);

                        if (tempImgData.imagePatch != NULL)
                            delete []tempImgData.imagePatch;
                        tempImgData.imagePatch = NULL;
                    }
                }

                senderListSet.clear();
            }
  
            if (procsInOtherList != NULL)
                delete []procsInOtherList;
            procsInOtherList = NULL;
    
            startingProc = newStartingProc;
            endingProc = newEndingProc;
        }

        if (informationToRecvArray != NULL)
            delete []informationToRecvArray;
        informationToRecvArray = NULL;

        if (informationToSendArray != NULL)
            delete []informationToSendArray;
        informationToSendArray = NULL;



        //
        // --- Timing -- 
        
        visitTimer->StopTimer(timingSendReceive, "Send Receive");
        visitTimer->DumpTimings();
        
        int  timingLocalCompositing2 = visitTimer->StartTimer();
        


        //
        // Each proc does local compositing and then sends
        //
        int imgBufferWidth = screen[0];
        int imgBufferHeight = screen[1];
        float *buffer = new float[((imgBufferWidth * imgBufferHeight * 4)) * numZDivisions]();  //size: imgBufferWidth * imgBufferHeight * 4, initialized to 0

        std::sort(allImgMetaData.begin(), allImgMetaData.end(), &sortImgMetaDataByDepth);

        std::multimap< std::pair<int,int>, imgData>::iterator itImgData;
        int bufferDivisionIndex = 0;
        int divIndex = 0;
        int totalSize = allImgMetaData.size();

        debug5 << PAR_Rank() << "   ~   totalSize to compose: " << totalSize << "    numZDivisions: " << numZDivisions << endl;

        for (int k=0; k<numZDivisions; k++){
            debug5 << PAR_Rank() << "   ~   division boundaries: " << divisionsArray[k*2 + 0] << " to " << divisionsArray[k*2 + 1] << endl;
        }

        for (int patchIndex=0; patchIndex<totalSize; patchIndex++){
            if (allImgMetaData[patchIndex].avg_z >= divisionsArray[divIndex*2] && allImgMetaData[patchIndex].avg_z <= divisionsArray[divIndex*2+1]){  //new index
            }else{
                for (int z=0; z<numZDivisions; z++){
                    if (allImgMetaData[patchIndex].avg_z >= divisionsArray[z*2] && allImgMetaData[patchIndex].avg_z <= divisionsArray[z*2+1]) {
                        divIndex = z;
                        break;
                    }  
                }
            }
            bufferDivisionIndex = (imgBufferWidth * imgBufferHeight * 4) * divIndex; 

            int startingX = allImgMetaData[patchIndex].screen_ll[0];
            int startingY = allImgMetaData[patchIndex].screen_ll[1]; 

            debug5 << PAR_Rank() << "   ~   divIndex: " << divIndex << "    composing patch #: " << allImgMetaData[patchIndex].procId << " ,  " << allImgMetaData[patchIndex].patchNumber << "   size: " << allImgMetaData[patchIndex].dims[0] << " x  " << allImgMetaData[patchIndex].dims[1] << "   avg_z: " << allImgMetaData[patchIndex].avg_z << endl;


            itImgData = imgDataToCompose.find( std::pair<int,int>(allImgMetaData[patchIndex].procId, allImgMetaData[patchIndex].patchNumber) );
            if (itImgData == imgDataToCompose.end()){
                debug5 << "Error in local compositing - shouldn't happen!!!  " << allImgMetaData[patchIndex].patchNumber << endl;
                continue;
            }
            
            // Assemble the idivisions into 1 layer
            for (int j=0; j<allImgMetaData[patchIndex].dims[1]; j++){
                for (int k=0; k<allImgMetaData[patchIndex].dims[0]; k++){
                    if ( ((startingX + k) < 0) || ((startingX + k) > imgBufferWidth) )
                        continue;

                    if ( ((startingY + j) < 0) || ((startingY + j) > imgBufferHeight) )
                        continue;
                    
                    int subImgIndex = allImgMetaData[patchIndex].dims[0]*j*4 + k*4;                             // index in the subimage 
                    int bufferIndex = (startingY*imgBufferWidth*4 + j*imgBufferWidth*4) + (startingX*4 + k*4);  // index in the big buffer

                    // back to Front compositing: composited_i = composited_i-1 * (1.0 - alpha_i) + incoming; alpha = alpha_i-1 * (1- alpha_i)
                    buffer[bufferDivisionIndex + bufferIndex+0] = (buffer[bufferDivisionIndex + bufferIndex+0] * (1.0 - itImgData->second.imagePatch[subImgIndex+3])) + itImgData->second.imagePatch[subImgIndex+0];
                    buffer[bufferDivisionIndex + bufferIndex+1] = (buffer[bufferDivisionIndex + bufferIndex+1] * (1.0 - itImgData->second.imagePatch[subImgIndex+3])) + itImgData->second.imagePatch[subImgIndex+1];
                    buffer[bufferDivisionIndex + bufferIndex+2] = (buffer[bufferDivisionIndex + bufferIndex+2] * (1.0 - itImgData->second.imagePatch[subImgIndex+3])) + itImgData->second.imagePatch[subImgIndex+2];
                    buffer[bufferDivisionIndex + bufferIndex+3] = (buffer[bufferDivisionIndex + bufferIndex+3] * (1.0 - itImgData->second.imagePatch[subImgIndex+3])) + itImgData->second.imagePatch[subImgIndex+3];
                }
            }

            if (itImgData->second.imagePatch != NULL)
                delete []itImgData->second.imagePatch;
            itImgData->second.imagePatch = NULL;
        }
}
*/