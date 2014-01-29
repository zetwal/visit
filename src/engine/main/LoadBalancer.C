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
//                                LoadBalancer.C                             //
// ************************************************************************* //

#include <LoadBalancer.h>
#include <avtCallback.h>

#include <stdlib.h>

#include <deque>
#include <set>

#include <avtDatabase.h>
#include <avtDatabaseMetaData.h>
#include <avtTerminatingSink.h>
#include <avtOriginatingSource.h>
#include <avtIOInformation.h>
#include <avtSILRestrictionTraverser.h>
#include <avtStreamingGhostGenerator.h>
#include <avtTypes.h>
#include <VisItException.h>
#include <DebugStream.h>
#include <AbortException.h>
#ifdef PARALLEL
    #include <mpi.h>
    #include <avtParallel.h>
#endif

#include <algorithm>
#include <map>
#include <deque>
#include <iterator>
#include <vector>
#include <sstream>



#include <avtDatabaseMetaData.h>

using     std::string;
using     std::vector;
using     std::deque;
using     std::set;


//
// Function Prototypes.
//

static avtDataRequest_p    ReduceCallback(void *,
                                                   avtContract_p);
static bool                      StreamingCheckerCallback(void *,
                                                   avtContract_p);
static bool                      ContinueCallback(void *, int);

//
// Static class data
//
ParAbortCallback        LoadBalancer::abortCallback = NULL;
void                   *LoadBalancer::abortCallbackArgs = NULL;
ProgressCallback        LoadBalancer::progressCallback = NULL;
void                   *LoadBalancer::progressCallbackArgs = NULL;

bool                    LoadBalancer::allowDynamic = false;
LoadBalanceScheme       LoadBalancer::scheme       =
                                       LOAD_BALANCE_CONTIGUOUS_BLOCKS_TOGETHER;

// ****************************************************************************
//  Method:  LoadBalancer::AllowDynamic
//
//  Purpose:
//    Allow dynamic load balancing when possible.
//
//  Programmer:  Jeremy Meredith
//  Creation:    September 20, 2001
//
//  Modifications:
//
//    Hank Childs, Sun Mar  6 08:42:50 PST 2005
//    Renamed method to AllowDynamic from ForceDynamic.
//
// ****************************************************************************

void
LoadBalancer::AllowDynamic()
{
    LoadBalancer::allowDynamic = true;
}

bool
LoadBalancer::GetAllowDynamic()
{
    return LoadBalancer::allowDynamic;
}

// ****************************************************************************
//  Method:  LoadBalancer::SetScheme
//
//  Purpose:
//      Sets the load balancing scheme used when static load balancing.
//
//  Programmer:  Hank Childs
//  Creation:    May 12, 2003
//
// ****************************************************************************

void
LoadBalancer::SetScheme(LoadBalanceScheme s)
{
    LoadBalancer::scheme = s;
}

// ****************************************************************************
// Method: LoadBalancer::GetScheme
//
// Purpose: 
//   Return the static load balancing scheme.
//
// Programmer: Brad Whitlock
// Creation:   Mon Oct 10 11:50:01 PDT 2011
//
// Modifications:
//   
// ****************************************************************************

LoadBalanceScheme
LoadBalancer::GetScheme()
{
    return LoadBalancer::scheme;
}

std::string
LoadBalancer::GetSchemeAsString()
{
    std::string str;
    switch(LoadBalancer::scheme)
    {
    case LOAD_BALANCE_CONTIGUOUS_BLOCKS_TOGETHER:
        str = "Contiguous Blocks Together";
        break;
    case LOAD_BALANCE_STRIDE_ACROSS_BLOCKS:
        str = "Stride Across Blocks";
        break;
    case LOAD_BALANCE_RANDOM_ASSIGNMENT:
        str = "Random Assignment";
        break;
    case LOAD_BALANCE_DBPLUGIN_DYNAMIC:
        str = "Database Plugin Dynamic";
        break;
    case LOAD_BALANCE_RESTRICTED:
        str = "Restricted";
        break;
    case LOAD_BALANCE_ABSOLUTE:
        str = "Absolute";
        break;
    case LOAD_BALANCE_STREAM:
        str = "Stream";
        break;
    case LOAD_BALANCE_KDTREE:
        str = "K-d tree based Assignment";
        break;
    default:
        break;
    }
    return str;
}



// ****************************************************************************
//  Method: LoadBalancer::chopPartition
//
//  Purpose:
//      Break a partition into 2 k-d tree style
//
//  Arguments:
//      parent      The parent partition.
//      childOne    The first child
//      childTwo    The second child
//      axisOrder   Order of the axes to cycle through
//
//  Programmer: Pascal Grosset  
//  Creation:   December 9, 2013
//
// ****************************************************************************
int          
LoadBalancer::kdtreeBuilding(int numDivisions, int logicalBounds[3], double minSpatialExtents[3], double maxSpatialExtents[3], 
                                               std::vector<patchMetaData> patches, std::vector<int> &list, bool amr){

    std::cout << rank << " ~~~ " << " LoadBalancer::kdtreeBuilding " << numDivisions << std::endl;

    /////////////////////////////////////////////////////////////////
    // Sort to find the order of axes from largest to smallest
    int axisOrder[3];
    
    std::multimap<int,int> myMultimap;
    for (int i=0; i<3; i++)
        myMultimap.insert(std::pair<int,int>(logicalBounds[i],i));
    
    int order=2;
    for (std::multimap<int,int>::iterator it=myMultimap.begin(); it!=myMultimap.end(); ++it){
        axisOrder[order]=(*it).second;
        order--;
    }
    std::cout << rank << " ~~~ " << "axisOrder: " << axisOrder[0] << ", " << axisOrder[1] << ", " << axisOrder[2] << std::endl << std::endl;


    /////////////////////////////////////////////////////////////////
    // Do the region splitting according to a k-d tree
    std::deque<partitionExtents> myPartitions;
    myPartitions.clear();

    partitionExtents parent, one, two;
    parent.axisIndex = 2;   // set it to the last one so that on the next iteration we get the first one! :)
    parent.head = true;
    parent.dims[0] = logicalBounds[0];  parent.dims[1] = logicalBounds[1];  parent.dims[2] = logicalBounds[2];
    parent.minExtents[0] = minSpatialExtents[0];    parent.minExtents[1] = minSpatialExtents[1];    parent.minExtents[2] = minSpatialExtents[2];
    parent.maxExtents[0] = maxSpatialExtents[0];    parent.maxExtents[1] = maxSpatialExtents[1];    parent.maxExtents[2] = maxSpatialExtents[2];

    myPartitions.push_back(parent);
    int startIndex = 0;

    while (myPartitions.size() != numDivisions){
        parent = myPartitions.front();   myPartitions.pop_front();
        //if (rank == 0)
        //     std::cout << rank << " ~~~ " << "Parent: "<< parent.axisIndex << "   - Extents (min-max):  " << parent.minExtents[0]<< ", " << parent.minExtents[1]<< ", " << parent.minExtents[2] << "   -   " << parent.maxExtents[0]<< ", " << parent.maxExtents[1]<< ", " << parent.maxExtents[2] << "  dims: " << parent.dims[0]<< ", " << parent.dims[1]<< ", " << parent.dims[2] << std::endl;
        
        chopPartition(parent,one,two,axisOrder);
        myPartitions.push_back(one);
        myPartitions.push_back(two);
        if (one.head == true)
        	startIndex = myPartitions.size()-2;
        else
        	startIndex = startIndex-1;

        // if (rank == 0){
        //     std::cout << rank << " ~~ " <<"One: "<<  one.axisIndex << "   - Extents (min-max):  " << one.minExtents[0]<< ", " << one.minExtents[1]<< ", " << one.minExtents[2] << "   -   " << one.maxExtents[0]<< ", " << one.maxExtents[1]<< ", " << one.maxExtents[2] << "  dims: " << one.dims[0]<< ", " << one.dims[1]<< ", " << one.dims[2] << std::endl;
        //     std::cout << rank << " ~~ " <<"two: "<<  two.axisIndex << "   - Extents (min-max):  " << two.minExtents[0]<< ", " << two.minExtents[1]<< ", " << two.minExtents[2] << "   -   " << two.maxExtents[0]<< ", " << two.maxExtents[1]<< ", " << two.maxExtents[2] << "  dims: " << two.dims[0]<< ", " << two.dims[1]<< ", " << two.dims[2] << std::endl;
        //     std::cout << std::endl;
        // }
    }
    
    //
    //Push head to the start of the list: this ensures more contiguity
    bool headfound = false;
    do{
    	partitionExtents tempExt = myPartitions.front();
    	if (tempExt.head == true)
    		headfound = true;
    	else{
    		// remove the front and put it at the back
    		myPartitions.pop_front();
    		myPartitions.push_back(tempExt);
    	}
    }while(headfound == true);


    /////////////////////////////////////////////////////////////////
    // Determine which patch is in which region

    // insert patches into a multimap (patchesTemp) sorted on id of patch
    std::multimap<int,patchMetaData> patchesTemp;
    for (std::vector<int>::iterator it=list.begin(); it!=list.end(); it++)
        patchesTemp.insert(std::pair<int,patchMetaData>((int)(*it),patches[(int)(*it)]));

    std::vector<int> patchesInsideList;  //list of patches that belong to that partition
    std::vector<int> patchesOverlapList; //list of division that overlap that partition

    double minX = myPartitions[rank].minExtents[0];  double maxX = myPartitions[rank].maxExtents[0];
    double minY = myPartitions[rank].minExtents[1];  double maxY = myPartitions[rank].maxExtents[1];
    double minZ = myPartitions[rank].minExtents[2];  double maxZ = myPartitions[rank].maxExtents[2];
    
    for (std::multimap<int,patchMetaData>::iterator it=patchesTemp.begin(); it!=patchesTemp.end(); it++){
        int key = (*it).first;

        bool done = false;
        double centroid[3];
        centroid[0] = ((*it).second.minSpatialExtents[0] + (*it).second.maxSpatialExtents[0])/2.0;
        centroid[1] = ((*it).second.minSpatialExtents[1] + (*it).second.maxSpatialExtents[1])/2.0;
        centroid[2] = ((*it).second.minSpatialExtents[2] + (*it).second.maxSpatialExtents[2])/2.0;

        double offset = 0.0000001;
        // Check if that patche belongs to this partition
        if (centroid[0]+offset >= minX   &&   centroid[0]+offset < maxX)
            if (centroid[1]+offset >= minY   &&   centroid[1]+offset < maxY)
                if (centroid[2]+offset >= minZ   &&   centroid[2]+offset < maxZ){
                    patchesInsideList.push_back(key);
                    done = true;
                }

        // else check if that patch overlaps with this partition
        if (done == false){
            if (patchOverlaps((*it).second.minSpatialExtents[0],(*it).second.maxSpatialExtents[0],
                              (*it).second.minSpatialExtents[1],(*it).second.maxSpatialExtents[1],
                              (*it).second.minSpatialExtents[2],(*it).second.maxSpatialExtents[2],
                              minX,maxX,minY,maxY,minZ,maxZ) == true)
                patchesOverlapList.push_back(key);
        }
    }

    /////////////////////////////////////////////////////////////////
    // Assign the patches to the list
    list.clear();
    std::sort(patchesInsideList.begin(),patchesInsideList.end()); // not specifically required

    // Adding those that belong
    for (int j=0; j<patchesInsideList.size(); j++)
        list.push_back(patchesInsideList[j]);


    // Adding those that are not specifically in but only overlap!
    if (amr)
        for (int j=0; j<patchesOverlapList.size(); j++)
            list.push_back(patchesOverlapList[j]);
    
    std::cout << rank << " ~ " << "  patchesInsideList.size(): " << patchesInsideList.size() << 
                                  "  patchesOverlapList.size(): " << patchesOverlapList.size() << 
                                  "  list.size(): " << list.size() << std::endl;

    return list.size();
}

// ****************************************************************************
//  Method: LoadBalancer::chopPartition
//
//  Purpose:
//      Break a partition into 2 k-d tree style
//
//  Arguments:
//      parent      The parent partition.
//      childOne    The first child
//      childTwo    The second child
//      axisOrder   Order of the axes to cycle through
//
//  Programmer: Pascal Grosset  
//  Creation:   December 9, 2013
//
// ****************************************************************************
int
LoadBalancer::chopPartition(partitionExtents parent, partitionExtents & childOne, partitionExtents & childTwo, int axisOrder[3]){
    int axisIndex = (parent.axisIndex+1)%3;
    int axis = axisOrder[axisIndex];
    int count = 0;
    
    while (true){
        if (parent.dims[axis] >= 1)
            break;
        else{
            axisIndex = (axisIndex+1)%3;
            axis = axisOrder[axisIndex];
        }
        
        count++;
        if (count == 3)
            break;
    }
    if (count == 3)
        return -1;  // We are going to be cycling forever here! So stop!
        
    childOne.axisIndex = axisIndex;
    childTwo.axisIndex = axisIndex;
    //std::cout << "axis: " << axis << std::endl;
        
    if (axis == 0){             // x-axis
        childOne.dims[1] = childTwo.dims[1] = parent.dims[1];
        childOne.dims[2] = childTwo.dims[2] = parent.dims[2];
        
        childOne.minExtents[1] = childTwo.minExtents[1] = parent.minExtents[1];
        childOne.maxExtents[1] = childTwo.maxExtents[1] = parent.maxExtents[1];
        
        childOne.minExtents[2] = childTwo.minExtents[2] = parent.minExtents[2];
        childOne.maxExtents[2] = childTwo.maxExtents[2] = parent.maxExtents[2];
        
        
        childOne.dims[0] = parent.dims[0]/2;
        childTwo.dims[0] = parent.dims[0]-childOne.dims[0];
        
        childOne.minExtents[0] = parent.minExtents[0];
        childOne.maxExtents[0] = parent.minExtents[0] + (parent.maxExtents[0]-parent.minExtents[0])/2.0;
        
        childTwo.minExtents[0] = childOne.maxExtents[0];
        childTwo.maxExtents[0] = parent.maxExtents[0];
    }
    else
        if (axis == 1){         // y-axis
            childOne.dims[0] = childTwo.dims[0] = parent.dims[0];
            childOne.dims[2] = childTwo.dims[2] = parent.dims[2];
            
            childOne.minExtents[0] = childTwo.minExtents[0] = parent.minExtents[0];
            childOne.maxExtents[0] = childTwo.maxExtents[0] = parent.maxExtents[0];
            
            childOne.minExtents[2] = childTwo.minExtents[2] = parent.minExtents[2];
            childOne.maxExtents[2] = childTwo.maxExtents[2] = parent.maxExtents[2];
            
            childOne.dims[1] = parent.dims[1]/2;
            childTwo.dims[1] = parent.dims[1]-childOne.dims[1];
            
            childOne.minExtents[1] = parent.minExtents[1];
            childOne.maxExtents[1] = parent.minExtents[1] + (parent.maxExtents[1]-parent.minExtents[1])/2.0;
            
            childTwo.minExtents[1] = childOne.maxExtents[1];
            childTwo.maxExtents[1] = parent.maxExtents[1];
        }
        else
            if (axis == 2){     // z-axis
                childOne.dims[0] = childTwo.dims[0] = parent.dims[0];
                childOne.dims[1] = childTwo.dims[1] = parent.dims[1];
                
                childOne.minExtents[0] = childTwo.minExtents[0] = parent.minExtents[0];
                childOne.maxExtents[0] = childTwo.maxExtents[0] = parent.maxExtents[0];
                
                childOne.minExtents[1] = childTwo.minExtents[1] = parent.minExtents[1];
                childOne.maxExtents[1] = childTwo.maxExtents[1] = parent.maxExtents[1];
                
                
                childOne.dims[2] = parent.dims[2]/2;
                childTwo.dims[2] = parent.dims[2]-childOne.dims[2];
                
                childOne.minExtents[2] = parent.minExtents[2];
                childOne.maxExtents[2] = parent.minExtents[2] + (parent.maxExtents[2]-parent.minExtents[2])/2.0;
                
                childTwo.minExtents[2] = childOne.maxExtents[2];
                childTwo.maxExtents[2] = parent.maxExtents[2];
            }
    if (parent.head == true)
    	childOne.head = true;
    else
    	childOne.head = false;
    childTwo.head = false;
    
    return 0;
}



// ****************************************************************************
//  Method: LoadBalancer::patchOverlaps
//
//  Purpose:
//      Determine if a patch overlaps in a partition
//
//  Arguments:
//      patch      The dimensions of the partition.
//      partition  The dimensions of the patch
//
//  Programmer: Pascal Grosset  
//  Creation:   December 22, 2013
//
// ****************************************************************************
bool 
LoadBalancer::patchOverlaps(float patchMinX, float patchMaxX, float patchMinY, float patchMaxY, float patchMinZ, float patchMaxZ,
                            float partitionMinX, float partitionMaxX, float partitionMinY, float partitionMaxY, float partitionMinZ, float partitionMaxZ){
    // Check if not outside; if it is not outside it has to be somewhere inside 
    if (patchMaxX < partitionMinX) return false;
    if (patchMinX > partitionMaxX) return false;
    if (patchMaxY < partitionMinY) return false;
    if (patchMinY > partitionMaxY) return false;
    if (patchMaxZ < partitionMinZ) return false;
    if (patchMinZ > partitionMaxZ) return false;
        return true;
}


// ****************************************************************************
//  Method: LoadBalancer::RegisterAbortCallback
//
//  Purpose:
//      Registers the AbortCallback.
//
//  Arguments:
//      pc      The abort callback.
//      args    The arguments to the abort callback.
//
//  Programmer: Jeremy Meredith
//  Creation:   September 20, 2001
//
// ****************************************************************************

void
LoadBalancer::RegisterAbortCallback(ParAbortCallback pc, void *args)
{
    abortCallback     = pc;
    abortCallbackArgs = args;
}

// ****************************************************************************
//  Method: LoadBalancer::RegisterProgressCallback
//
//  Purpose:
//      Registers the ProgressCallback.
//
//  Arguments:
//      pc      The progress callback.
//      args    The arguments to the progress callback.
//
//  Programmer: Jeremy Meredith
//  Creation:   September 19, 2001
//
// ****************************************************************************

void
LoadBalancer::RegisterProgressCallback(ProgressCallback pc, void *args)
{
    progressCallback     = pc;
    progressCallbackArgs = args;
}

// ****************************************************************************
//  Method:  LoadBalancer::CheckAbort
//
//  Purpose:
//    Check for an abort message
//
//  Programmer:  Jeremy Meredith
//  Creation:    September 20, 2001
//
// ****************************************************************************

bool
LoadBalancer::CheckAbort(bool informSlaves)
{
    return abortCallback(abortCallbackArgs, informSlaves);
}

// ****************************************************************************
//  Method:  LoadBalancer::UpdateProgress
//
//  Purpose:
//    Send a progress update using the registered callback.
//
//  Arguments:
//    current    number of parts processed
//    total      total number of parts
//
//  Programmer:  Jeremy Meredith
//  Creation:    September 19, 2001
//
// ****************************************************************************

void
LoadBalancer::UpdateProgress(int current, int total)
{
    if (progressCallback != NULL)
        progressCallback(progressCallbackArgs, "Calculating", "Calculating",
                         current, total);
}

// ****************************************************************************
//  Method: LoadBalancer constructor
//
//  Arguments:
//      np      The number of processors.
//      r       The index of this processor (0-origin).
//
//  Programmer: Hank Childs
//  Creation:   June 17, 2001
//
//  Modifications:
//    Jeremy Meredith, Thu Jul 26 12:29:29 PDT 2001
//    Use a data structure for pipeline info.
//
//    Jeremy Meredith, Thu Sep 20 00:49:26 PDT 2001
//    Added the registration of the DynamicChecker callback.
//
//    Hank Childs, Tue Feb 19 19:45:43 PST 2008
//    Rename "dynamic" to "streaming", since we really care about whether we
//    are streaming, not about whether we are doing dynamic load balancing.
//    And the two are no longer synonymous.
//
// ****************************************************************************

LoadBalancer::LoadBalancer(int np, int r)
{
    //
    // Pipeline index 0 is reserved for meta-data pipelines.
    //
    pipelineInfo.push_back(LBInfo("dummy_pipeline"));

    rank   = r;
    nProcs = np;

    amrOn = false;
    amrLevels.clear();

    //
    // Register callbacks with the avt pipeline.
    //
    avtOriginatingSource::SetLoadBalancer(ReduceCallback, this);
    avtOriginatingSource::SetStreamingChecker(StreamingCheckerCallback, this);
    avtTerminatingSink::SetGuideFunction(ContinueCallback, this);
}


// ****************************************************************************
//  Method: LoadBalancer::CheckDynamicLoadBalancing
//
//  Purpose:
//      Takes in the pipeline specification for the entire pipeline and
//      determines whether or not it will be dynamically load balanced
//
//  Arguments:
//      spec    A pipeline specification.
//
//  Returns:    true if dynamic; false if static or none
//
//  Programmer: Jeremy Meredith
//  Creation:   September 19, 2001
//
//  Modifications:
//    Jeremy Meredith, Fri Sep 21 14:39:58 PDT 2001
//    Added checks for forcing static/dynamic load balancing.
//
//    Hank Childs, Sat Feb 19 14:27:03 PST 2005
//    Allow for dynamic load balancing on serial engines.  Also allow for
//    dynamic load balancing to take place if the command line has 
//    "-allowdynamic".  Finally, cache our answer so that we don't make the
//    same calculation repeatedly in DLB mode.
//
//    Hank Childs, Thu Feb 14 15:54:36 PST 2008
//    This function is misnamed.  It really has to do with whether or not we
//    are *streaming* ... not necessarily doing dynamic load balancing 
//    (DLB implies streaming, but streaming does not imply DLB).
//    Return true if we are streaming.  Renaming will come later.
//
// ****************************************************************************

bool
LoadBalancer::CheckDynamicLoadBalancing(avtContract_p input)
{
    //
    // See if we have already have decided.  If so, just return our cached
    // decision.
    //
    int index = input->GetPipelineIndex();
    LBInfo &lbinfo = pipelineInfo[index];
    if (lbinfo.haveInitializedDLB)
        return lbinfo.doDLB;

    //
    // If the user has not explicitly asked for DLB, then don't do it.
    //
    if (!allowDynamic)
    {
        // Almost always false.
        lbinfo.doDLB = false || (scheme == LOAD_BALANCE_STREAM);
        lbinfo.haveInitializedDLB = true;
        return lbinfo.doDLB;
    }

    //
    // Some hard and fast rules:
    //
    // Pipeline index 0 is reserved for meta-data and inlined pipelines.  So
    // no DLB for those.
    //
    // Cannot dynamic load balance some pipelines because of the filters
    // in the pipeline.
    //
    // We cannot do dynamic load balancing if the database does not believe
    // we can do dynamic load balancing (for example because we need ghost
    // data communicated or materials reconstructed).
    //
    avtDataRequest_p data = input->GetDataRequest();
    std::string dbname = lbinfo.db;
    avtDatabase *db = dbMap[dbname];
    if (input->GetPipelineIndex() == 0 ||
        input->ShouldUseStreaming() == false || 
        db->CanDoStreaming(data) == false)
    {
        lbinfo.doDLB = false;
        lbinfo.haveInitializedDLB = true;
        return false;
    }

    //
    // Don't do DLB if we have 2 or 3 procs.  It's not worth it.
    //
    if (nProcs == 2 || nProcs == 3)
    {
        lbinfo.doDLB = false;
        lbinfo.haveInitializedDLB = true;
        return false;
    }

    //
    // The user has asked for DLB.  And nothing in the pipeline is prevent it.
    // Do it!
    //
    lbinfo.doDLB = true;
    lbinfo.haveInitializedDLB = true;
    return true;
}


// ****************************************************************************
//  Method: LoadBalancer::DetermineAppropriateScheme
//
//  Purpose: Permits a mesh to override the default load balance scheme
//
//  Programmer: Mark C. Miller
//  Creation:   October 4, 2005
//
//  Modifications:
//
//    Mark C. Miller, Thu Nov 17 11:46:43 PST 2005
//    Test for non-NULL mmd before dereferencing
//
//    Hank Childs, Fri Nov 18 16:15:12 PST 2005
//    Better accomodate CMFE ... don't assume that the database stored with
//    this pipeline index is the one we're referencing.
//
//    Hank Childs, Tue Mar  7 10:43:53 PST 2006
//    Check here to see if we should do DBPLUGIN_DYNAMIC, rather than
//    setting it as a global.
//    
//    Mark C. Miller, Wed Jun 17 14:25:59 PDT 2009
//    Replaced CATCH(...) with CATCHALL.
// ****************************************************************************

LoadBalanceScheme
LoadBalancer::DetermineAppropriateScheme(avtContract_p input)
{

    //
    // See if we have already have decided.  If so, just return our cached
    // decision.
    //
    int index = input->GetPipelineIndex();
    const LBInfo &lbinfo = pipelineInfo[index];
    std::string dbname = lbinfo.db;
    avtDatabase *db = dbMap[dbname];

    avtDataRequest_p data = input->GetDataRequest();
    avtDatabaseMetaData *md = db->GetMetaData(db->GetMostRecentTimestep());
    string meshName;

    TRY
    {
        meshName = md->MeshForVar(data->GetVariable());
    }
    CATCHALL
    {
        // Probably a CMFE.
        return scheme;
    }
    ENDTRY;

    if (md->GetFormatCanDoDomainDecomposition())
        return LOAD_BALANCE_DBPLUGIN_DYNAMIC;

    const avtMeshMetaData *mmd = md->GetMesh(meshName);

    // for (int p=0; p<mmd->numBlocks; p++){
    //     std::cout << "  ~ 2D Load balance  Parent: " << p << "   size: " << mmd->patch_parent[p].size() << std::endl;
    //     for (int j=0; j<mmd->patch_parent[p].size(); j++)
    //         std::cout << "  ~ 2D Vec Parent: " <<  p << "   child: " << mmd->patch_parent[p][j] << std::endl;
    // }

    if (mmd && mmd->loadBalanceScheme != LOAD_BALANCE_UNKNOWN)
    {
        debug1 << "Default load balance scheme \""
               << LoadBalanceSchemeToString(scheme).c_str() << "\""
               << " being overridden in favor of \""
               << LoadBalanceSchemeToString(mmd->loadBalanceScheme).c_str() << "\""
               << " for mesh \"" << meshName.c_str() << "\"" << endl;
        return mmd->loadBalanceScheme;
    }

    return scheme;
}

// ****************************************************************************
//  Method: LoadBalancer::Reduce
//
//  Purpose:
//      Takes in the pipeline specification for the entire pipeline and
//      determines which portion this processor should do on this pass.
//
//  Arguments:
//      spec    A pipeline specification.
//
//  Returns:    A data specification that is a subset of the data specification
//              in the pipeline specification.
//
//  Programmer: Hank Childs
//  Creation:   June 17, 2001
//
//  Modifications:
//    Jeremy Meredith, Thu Jul 26 12:21:57 PDT 2001
//    Added dynamic load balancing.
//
//    Hank Childs, Wed Aug  8 16:48:45 PDT 2001
//    Added intersection with original restriction when in parallel (to
//    preserve restrictions that don't involve domains).
//
//    Jeremy Meredith, Wed Aug 29 12:07:23 PDT 2001
//    Fixed a bug in static load balancing -- there was overlap in 
//    assignment of domains, preventing all domains from being assigned
//    to processors.
//
//    Jeremy Meredith, Mon Sep 17 23:07:43 PDT 2001
//    Made the master be aware of the domains already cached on each 
//    processor, and to choose those domains for a processor when possible.
//
//    Jeremy Meredith, Thu Sep 20 00:52:37 PDT 2001
//    Made it use the CheckDynamicLoadBalancing function.
//    Made it send progress updates in dynamic mode.
//    Added logic to make full use of IO hints.
//
//    Jeremy Meredith, Fri Sep 21 14:41:08 PDT 2001
//    Added support for aborting dynamic loadbalanced execution.
//    Changed dynamic loadbalancing to delay sending the "complete"
//    signal to slave processes until it is sure no abort has happened.
//
//    Hank Childs, Mon Dec  2 14:46:03 PST 2002
//    Use a SIL restriction traverser to find the domain list.
//
//    Hank Childs, Wed Dec  4 17:23:12 PST 2002
//    Made use of a more efficient call to intersect SIL restrictions.
//
//    Hank Childs, Mon May 12 19:34:54 PDT 2003
//    Account for different load balancing schemes.
//
//    Mark C. Miller, Wed Jun  9 21:50:12 PDT 2004
//    Eliminated use of MPI_ANY_TAG and modified to use GetUniqueMessageTags
//
//    Mark C. Miller, Tue Sep 28 19:57:42 PDT 2004
//    Added the very trivial load balance scheme LOAD_BALANCE_DBPLUGIN_DYNAMIC
//    where all processors are assigned the one and only block
//
//    Hank Childs, Sat Feb 19 14:27:03 PST 2005
//    Allow for dynamic load balancing with serial engines.
//
//    Hank Childs, Sat Mar  5 18:53:28 PST 2005
//    Take more care in setting up unique message tags.  These methods may be
//    called a different number of times.  So use a static.
//
//    Jeremy Meredith, Wed May 11 09:12:40 PDT 2005
//    Added "restricted" load balancing mode.  This is intended for
//    non-global filesystems and simulation-mode engines.  It occurs when
//    each processor can only access a limited subset of the domains.
//
//    Mark C. Miller, Thu Sep 15 11:30:18 PDT 2005
//    Added "absolute" load balancing mode where domains are assigned based
//    upon their absolute domain number modulo the number of processors.
//    This guarentees that domains are never re-read on other processors.
//    However, it obviously can negatively effect balance based on the
//    current SIL restriction. Nonetheless, assuming the user has chosen
//    a number of processors to achieve adaquate performance in the
//    worst-case plotting scanerio, this scheme will certainly do no worse
//    and might do better due to guarenteeing that no re-reading is done.
//    
//    Mark C. Miller, Wed Nov 16 10:46:36 PST 2005
//    Added call to DetermineAppropriateLoadBalanceScheme 
//    Changed DBPLUGIN_DYNAMIC scheme so that every processors gets the
//    complete list of domains
//
//    Mark C. Miller, Mon Jan 22 22:09:01 PST 2007
//    Changed MPI_COMM_WORLD to VISIT_MPI_COMM
//
//    Hank Childs, Sun Feb 10 09:43:42 PST 2008
//    Use the streaming ghost generator to direct the load balancing,
//    if there is a streaming ghost generator.
//
//    Dave Pugmire, Tue May 25 10:15:35 EDT 2010
//    Add domain single domain replication to all processors.
//
//    Dave Pugmire, Mon Jun 14 14:16:57 EDT 2010
//    Single domain replication needs to mark pipeline update complete.
//
// ****************************************************************************

avtDataRequest_p
LoadBalancer::Reduce(avtContract_p input)
{
    avtDataRequest_p data = input->GetDataRequest();

    //
    // It is difficult for the load balancer to communicate with the originating
    // source because it is done through callbacks.
    // So we do it by setting a Boolean in the contract.  Since there is only
    // one path that involves actually doing data replication, and many that don't,
    // we will unset the Boolean now and reset it in the case we actually do
    // data replication.
    //
    bool dataReplicationRequested = input->ReplicateSingleDomainOnAllProcessors();
    input->SetReplicateSingleDomainOnAllProcessors(false);

    //
    // Pipeline index 0 is reserved for meta-data.  It should already be
    // load balanced.
    //
    if (input->GetPipelineIndex() == 0)
    {
        return data;
    }

    //
    // Assess load balancing specially for serial engines.
    //
    if (nProcs <= 1)
    {
        bool doDynLB = CheckDynamicLoadBalancing(input);
        if (!doDynLB && scheme != LOAD_BALANCE_STREAM)
        {
            pipelineInfo[input->GetPipelineIndex()].complete = true;
            return data;
        }
        else
        {
            avtDataObjectSource::RegisterProgressCallback(NULL,NULL);
            avtSILRestriction_p orig_silr   = data->GetRestriction();
            avtSILRestriction_p silr        = new avtSILRestriction(orig_silr);
            avtDataRequest_p new_data = new avtDataRequest(data, silr);
            avtSILRestrictionTraverser trav(silr);

            vector<int> list;
            trav.GetDomainList(list);
    
            if (pipelineInfo[input->GetPipelineIndex()].current < 0)
                pipelineInfo[input->GetPipelineIndex()].current  = 0;
            int domain = list[pipelineInfo[input->GetPipelineIndex()].current];
            int sggDomain = avtStreamingGhostGenerator::LBGetNextDomain();
            if (sggDomain >= 0)
                domain = sggDomain;
            vector<int> domainList(1, domain);
            new_data->GetRestriction()
                                   ->RestrictDomainsForLoadBalance(domainList);
            UpdateProgress(pipelineInfo[input->GetPipelineIndex()].current,
                           (int)list.size());
            pipelineInfo[input->GetPipelineIndex()].current++;
            if (pipelineInfo[input->GetPipelineIndex()].current == list.size())
                pipelineInfo[input->GetPipelineIndex()].complete = true;
            return new_data;
        }
    }

#ifdef PARALLEL

    avtSILRestriction_p orig_silr   = data->GetRestriction();
    avtSILRestriction_p silr        = new avtSILRestriction(orig_silr);
    avtDataRequest_p new_data = new avtDataRequest(data, silr);
    avtSILRestrictionTraverser trav(silr);

    // set up MPI message tags
    static int lastDomDoneMsg = GetUniqueMessageTag();
    static int newDomToDoMsg = GetUniqueMessageTag();

    if (scheme == LOAD_BALANCE_STREAM)
    {
        if (pipelineInfo[input->GetPipelineIndex()].current < 0)
        {
            pipelineInfo[input->GetPipelineIndex()].current = 0;

            //
            // We probably want to do something more sophisticated in the future 
            // (like walking through a SIL).  For now, just use the "chunks"
            // mechanism set up with convenience methods.
            //
            vector<int> list;
            trav.GetDomainList(list);

            int amountPer = list.size() / nProcs;
            int oneExtraUntil = list.size() % nProcs;
            int lastDomain = 0;
            for (int i = 0 ; i < nProcs ; i++)
            {
                if (i == rank)
                {
                    int amount = amountPer + (i < oneExtraUntil ? 1 : 0);
                    for (int j = 0 ; j < amount ; j++)
                    {
                        domainListForStreaming.push_back(list[j+lastDomain]);
                    }
                }
                lastDomain += amountPer + (i < oneExtraUntil ? 1 : 0);
            }
        }

        int domain = domainListForStreaming[pipelineInfo[input->GetPipelineIndex()].current];
        int sggDomain = avtStreamingGhostGenerator::LBGetNextDomain();
        if (sggDomain >= 0)
            domain = sggDomain;
        vector<int> domainList(1, domain);
        new_data->GetRestriction()
                               ->RestrictDomainsForLoadBalance(domainList);
        UpdateProgress(pipelineInfo[input->GetPipelineIndex()].current,
                       domainListForStreaming.size());
        pipelineInfo[input->GetPipelineIndex()].current++;
        if (pipelineInfo[input->GetPipelineIndex()].current == domainListForStreaming.size())
        {
            pipelineInfo[input->GetPipelineIndex()].complete = true;
            domainListForStreaming.clear();
        }
    }
    // Can we do dynamic load balancing?
    else if (! CheckDynamicLoadBalancing(input))
    {
        //
        // We probably want to do something more sophisticated in the future 
        // (like walking through a SIL).  For now, just use the "chunks"
        // mechanism set up with convenience methods.
        //
        vector<int> list;
        vector<int> mylist;
        trav.GetDomainList(list);

        if (dataReplicationRequested && list.size() == 1)
        {
            silr->RestrictDomainsForLoadBalance(list);
            pipelineInfo[input->GetPipelineIndex()].complete = true;

            // Communicate back to the pipeline that we are replicating.
            input->SetReplicateSingleDomainOnAllProcessors(true);

            return data;
        }

        //
        // For variables (including meshes) that require specific types of
        // load balancing, we override the scheme here
        //
        LoadBalanceScheme theScheme = DetermineAppropriateScheme(input);
        std::cout << "theScheme: " << theScheme << std::endl;

        int index = input->GetPipelineIndex();
        const LBInfo &lbinfo = pipelineInfo[index];
        std::string dbname = lbinfo.db;
        avtDatabase *db = dbMap[dbname];

        avtDataRequest_p data = input->GetDataRequest();
        avtDatabaseMetaData *md = db->GetMetaData(db->GetMostRecentTimestep());
        string meshname = md->MeshForVar(new_data->GetVariable());
        std::cout << "!!!! !!!! meshname: " << meshname << std::endl;

        const avtMeshMetaData *mmd = md->GetMesh(meshname);
        patchesInfo = mmd->patches;

        if (mmd->levels.size() > 1){
            amrOn = true;
            amrLevels = mmd->levels;
            
            avtCallback::SetAMR(true);
        }else{
            amrOn = false;
            amrLevels.clear();

            avtCallback::SetAMR(false);
        }

        int logicalBounds[3];
        double minSpatialExtents[3], maxSpatialExtents[3];
        for (int i=0; i<3; i++){
            logicalBounds[i]=mmd->logicalBounds[i];
            minSpatialExtents[i]=mmd->minSpatialExtents[i];
            maxSpatialExtents[i]=mmd->maxSpatialExtents[i];
        }
        std::vector<int> templist;
        std::vector<int> numPatchesPerProc;
        templist = list;

        // if (rank ==0){
        //     std::cout << "List size: " << templist.size() << std::endl;
        //     for (int i=0; i<templist.size(); i++){
        //         std::cout << templist[i] << "  ";
        //     }
        
        //     std::cout << "amrLevels levels: " << amrLevels.size() << std::endl;
        //     for (int i=0; i<amrLevels.size(); i++){
        //         std::cout << amrLevels[i] << "  ";
        //     }
        // }

        // std::cout << rank << " ~ loadBalancer.C: Full dimensions: " << mmd->minSpatialExtents[0] << " , " << mmd->maxSpatialExtents[0] << 
        //                          "      " << mmd->minSpatialExtents[1] << " , " << mmd->maxSpatialExtents[1] << 
        //                          "      " << mmd->minSpatialExtents[2] << " , " << mmd->maxSpatialExtents[2] << std::endl;

        if (theScheme == LOAD_BALANCE_KDTREE)
        {
            std::cout << "||| K-d tree load balancing ||| " << patchesInfo.size() << std::endl;

            int numPatches = kdtreeBuilding(nProcs, logicalBounds, minSpatialExtents, maxSpatialExtents, patchesInfo,templist,amrOn);
            std::cout << rank << "  numPatches: " << numPatches << std::endl;
            
            mylist.clear();
            mylist.reserve(numPatches);
            for (int j=0; j<numPatches; j++)
                mylist.push_back(templist[j]);
        }
        else if (theScheme == LOAD_BALANCE_CONTIGUOUS_BLOCKS_TOGETHER)
        {
            int amountPer = list.size() / nProcs;
            int oneExtraUntil = list.size() % nProcs;
            int lastDomain = 0;
            for (int i = 0 ; i < nProcs ; i++)
            {
                if (i == rank)
                {
                    int amount = amountPer + (i < oneExtraUntil ? 1 : 0);
                    for (int j = 0 ; j < amount ; j++)
                    {
                        mylist.push_back(list[j+lastDomain]);
                    }
                }
                lastDomain += amountPer + (i < oneExtraUntil ? 1 : 0);
            }
        }
        else if (theScheme == LOAD_BALANCE_STRIDE_ACROSS_BLOCKS)
        {
            for (int j = 0 ; j < list.size() ; j++)
            {
                if (j % nProcs == rank)
                    mylist.push_back(list[j]);
            }
        }
        else if (theScheme == LOAD_BALANCE_ABSOLUTE)
        {
            for (int j = 0 ; j < list.size() ; j++)
            {
                if (list[j] % nProcs == rank)
                    mylist.push_back(list[j]);
            }
        }
        else if (theScheme == LOAD_BALANCE_RESTRICTED)
        {
            LBInfo &lbInfo(pipelineInfo[input->GetPipelineIndex()]);
            IOInfo &ioInfo(ioMap[lbInfo.db]);
            const HintList &hints(ioInfo.ioInfo.GetHints());

            for (int j = 0 ; j < list.size() ; j++)
            {
                if (hints.size() >= rank)
                {
                    const vector<int> &doms = hints[rank];
                    int ndoms = doms.size();
                    for (int h=0; h<ndoms; h++)
                    {
                        if (doms[h] == list[j])
                        {
                            mylist.push_back(list[j]);
                            break;
                        }
                    }
                }
            }
        }
        else if (theScheme ==  LOAD_BALANCE_RANDOM_ASSIGNMENT)
        {
            // all procs randomly jumble the list of domain ids
            // all procs compute same jumbled list due to same seed
            // [ which won't be true on a heterogeneous platform ]
            int j;
            vector<int> jumbledList = list;
            srand(0xDeadBeef);
            for (j = 0 ; j < list.size() * 5; j++)
            {
               int i1 = rand() % list.size();
               int i2 = rand() % list.size();
               int tmp = jumbledList[i1];
               jumbledList[i1] = jumbledList[i2];
               jumbledList[i2] = tmp;
            }
            // now, do round-robin assignment from the jumbled list
            for (j = 0 ; j < list.size() ; j++)
            {
                if (j % nProcs == rank)
                    mylist.push_back(jumbledList[j]);
            }
        }
        else if (theScheme == LOAD_BALANCE_DBPLUGIN_DYNAMIC)
        {
            // Every processor gets the complete list
            mylist = list;
        }

        //Display the list:
        std::stringstream ss;
        ss << rank << " ~ size: " << mylist.size() << "  patches: ";
        for (int i=0; i<mylist.size(); i++)
            ss << mylist[i] << ", ";
        //std::cout << ss.str() << std::endl;

        silr->RestrictDomainsForLoadBalance(mylist);
        pipelineInfo[input->GetPipelineIndex()].complete = true;
    }
    else
    {
        // disable progress updates from the filters this time around
        avtDataObjectSource::RegisterProgressCallback(NULL,NULL);

        LBInfo &lbInfo(pipelineInfo[input->GetPipelineIndex()]);
        IOInfo &ioInfo(ioMap[lbInfo.db]);
        if (rank == 0)
        {
            // -------------------------------------
            //     MASTER LOADBALANCER PROCESSES
            // -------------------------------------

            // Allocate enough space to hold the completed domains
            ioInfo.domains.resize(nProcs);
            ioInfo.files.resize(nProcs);
            bool validFileMap = (ioInfo.fileMap.size() != 0);

            // Get the list of domains to process
            vector<int> domainList;
            trav.GetDomainList(domainList);

            // Make a work list and a completed list
            int         totaldomains = domainList.size();
            deque<int>  incomplete(domainList.begin(), domainList.end());
            vector<int> complete;

            debug5 << "LoadBalancer Master -- starting with " 
                   << incomplete.size() << " domains\n";

            // pull from the incomplete list and push onto the complete list
            // until all domains are complete
            bool abort = false;
            int domain;
            UpdateProgress(0,0);
            while (complete.size() < totaldomains)
            {
                // check for an abort
                if (!abort &&
                    CheckAbort(false))
                {
                    abort = true;
                    totaldomains -= incomplete.size();
                    incomplete.clear();
                }

                // update the progress
                UpdateProgress(complete.size() + (domainList.size() - incomplete.size()),
                               domainList.size()*2);


                // get the completed domain number
                MPI_Status stat;
                MPI_Recv(&domain, 1, MPI_INT, MPI_ANY_SOURCE,
                         lastDomDoneMsg, VISIT_MPI_COMM, &stat);
                int processor = stat.MPI_SOURCE;

                // -1 means the first pass by the slave; nothing completed yet
                if (domain != -1)
                {
                    // add it to the complete list
                    complete.push_back(domain);
                }

                // figure out what to tell this processor to do
                if (incomplete.empty())
                    continue;

                // find a cached domain for next processor
                deque<int>::iterator i;
                for (i = incomplete.begin(); i != incomplete.end(); i++)
                {
                    if (ioInfo.domains[processor].find(*i) != 
                        ioInfo.domains[processor].end())
                        break;
                }
                // if no match, try to find one that is in a file
                // already opened by this processor
                if (i == incomplete.end())
                {
                    for (i = incomplete.begin(); i != incomplete.end(); i++)
                    {
                        int fileno = 0;
                        if (validFileMap)
                            fileno = ioInfo.fileMap[*i];
                        if (ioInfo.files[processor].count(fileno) > 0)
                            break;
                    }
                }
                // if still no match, find one that is in a file
                // opened by the fewest number of processors
                if (i == incomplete.end())
                {
                    int mindomain = -1;
                    int minopen   = 999999999;
                    for (i = incomplete.begin(); i != incomplete.end(); i++)
                    {
                        int fileno = 0;
                        if (validFileMap)
                            fileno = ioInfo.fileMap[*i];
                        // count the number of processors which have
                        // this file opened
                        int nopen = 0;
                        for (int j=0; j<ioInfo.files.size(); j++)
                            if (ioInfo.files[j].count(fileno) > 0)
                                nopen++;
                        if (nopen < minopen)
                        {
                            mindomain = *i;
                            minopen   = nopen;
                        }
                    }
                    for (i = incomplete.begin(); i != incomplete.end(); i++)
                    {
                        if (*i == mindomain)
                            break;
                    }
                }                    

                // if no match, just take the next one in line
                if (i == incomplete.end())
                    i=incomplete.begin();

                domain = *i;
                incomplete.erase(i);

                ioInfo.domains[processor].insert(domain);
                if (validFileMap)
                    ioInfo.files[processor].insert(ioInfo.fileMap[domain]);
                else
                    ioInfo.files[processor].insert(0);

                // send the new domain number to that processor
                debug5 << "LoadBalancer Master: sending domain " 
                       << domain << " to processor "<<processor<<"\n";
                MPI_Send(&domain, 1, MPI_INT, processor, newDomToDoMsg, VISIT_MPI_COMM);
            }

            // we're all done -- -2 means to abort, -1 means to send results
            int status = abort ? -2 : -1;
            for (int i=1; i<nProcs; i++)
                MPI_Send(&status, 1, MPI_INT, i, newDomToDoMsg,VISIT_MPI_COMM);

            if (abort)
                EXCEPTION0(AbortException);

            // all work is done
            UpdateProgress(1,0);
            lbInfo.complete = true;
            new_data->GetRestriction()->TurnOffAll();
            MPI_Barrier(VISIT_MPI_COMM);


        }
        else
        {
            // -------------------------------------
            //            SLAVE PROCESSES
            // -------------------------------------

            // send our last completed domain to the master
            int domain = lbInfo.current;
            MPI_Send(&domain, 1, MPI_INT, 0, lastDomDoneMsg, VISIT_MPI_COMM);

            // get our new work unit
            MPI_Status stat;
            MPI_Recv(&domain, 1, MPI_INT, 0, newDomToDoMsg, VISIT_MPI_COMM, &stat);
            lbInfo.current = domain;

            if (domain == -2)
            {
                EXCEPTION0(AbortException);
            }
            else if (domain == -1)
            {
                //  -1 is a tag for "no work" -- we are all done
                lbInfo.complete = true;
                new_data->GetRestriction()->TurnOffAll();
                MPI_Barrier(VISIT_MPI_COMM);
            }
            else
            {
                vector<int> domainList(1, domain);
                new_data->GetRestriction()
                                   ->RestrictDomainsForLoadBalance(domainList);
            }
        }
    }

    // By intersecting with the original restriction, we will ensure that
    // we are catching restrictions beyond domains, like materials, etc.
    // See comments in SIL restriction code regarding 'FastIntersect'.
    new_data->GetRestriction()->FastIntersect(orig_silr);

    return new_data;
#else
    EXCEPTION1(VisItException, "nprocs was > 1 in a non-parallel code");
#endif

}


// ****************************************************************************
//  Method: LoadBalancer::AddDatabase
//
//  Purpose:
//      Add a database to the load balancer.  Also adds the I/O information 
//      that should be used when balancing a load.
//
//  Arguments:
//      name     The name of the database.
//      ioinfo   Information about which domains should be grouped together on
//               the same processor.
//
//  Notes:  This will need to expand to support IO restrictions for
//          clustered (non-global) file systems.
//
//  Programmer:  Jeremy Meredith
//  Creation:    July 26, 2001
//
//  Modifications:
//    Jeremy Meredith, Thu Sep 20 00:54:58 PDT 2001
//    Added setting of fileMap from the iohints.
//
//    Mark C. Miller, Tue Sep 28 19:57:42 PDT 2004
//    Added code to set the load balance scheme if the metadata indicates
//    plugin can do its own decomposition
//
//    Hank Childs, Sun Feb 27 11:12:44 PST 2005
//    Added avtDatabase argument.
//
//    Mark C. Miller, Wed Nov 16 10:46:36 PST 2005
//    Removed avtIOInformation and avtDatabaseMetaData args because the
//    are obtainable from the database_ptr
//
//    Hank Childs, Tue Mar  7 10:43:53 PST 2006
//    Make the decision to do DBPLUGIN_DYNAMIC load balancing on a per
//    input basis.
//
// ****************************************************************************

void
LoadBalancer::AddDatabase(const string &db, avtDatabase *db_ptr, int time)
{
    const avtIOInformation& io = db_ptr->GetIOInformation(time);

    dbMap[db] = db_ptr;
    ioMap[db].ioInfo = io;
    ioMap[db].fileMap.resize(io.GetNDomains());

    debug4 << "LoadBalancer::AddDatabase - db=" << db.c_str() << endl;
    debug4 << "    iohints=[";
    const HintList &hints = io.GetHints();
    for (int i=0; i<hints.size(); i++)
    {
        debug4 << " {";
        for (int j=0; j<hints[i].size(); j++)
        {
            ioMap[db].fileMap[hints[i][j]] = i;
            debug4 << hints[i][j];
            if (j<hints[i].size()-1) debug4 << ",";
        }
        debug4 << "}";
        if (i<hints.size()-1)
            debug4 << "\n             ";
    }
    debug4 << "]  " << endl;
}

// ****************************************************************************
//  Method: LoadBalancer::AddPipeline
//
//  Purpose:
//      Creates a unique pipeline index for an avt pipeline.  Also tells the
//      load balancer what database is associated with that pipeline.
//
//  Arguments:
//      name     The name of the database.
//
//  Returns:     A unique index for a pipeline.
//
//  Programmer:  Jeremy Meredith
//  Creation:    July 26, 2001
//
// ****************************************************************************

int
LoadBalancer::AddPipeline(const string &db)
{
    int index = pipelineInfo.size();
    pipelineInfo.push_back(LBInfo(db));
    return index;
}

 
// ****************************************************************************
//  Method: LoadBalancer::ResetPipeline
//
//  Purpose:
//      Resets the status of a pipeline so it can re-execute.
//
//  Arguments:
//      index   The pipeline index.
//
//  Programmer: Hank Childs
//  Creation:   November 21, 2001
//
// ****************************************************************************

void
LoadBalancer::ResetPipeline(int index)
{
    if (index < 0 || index >= pipelineInfo.size())
    {
        debug1 << "Given an invalid pipeline index to reset (" << index << ")."
               << endl;
        return;
    }

    pipelineInfo[index].complete = false;
    pipelineInfo[index].current  = -1;
}

   
// ****************************************************************************
//  Method: LoadBalancer::ContinueExecute
//
//  Purpose:
//      Determines if all of the data has been read in for a pipeline and if
//      it can stop executing.  Data object sinks call this method.
//
//  Arguments:
//      index    A pipeline index.
//
//  Returns:     true if that pipeline needs to continue executing.
//
//  Programmer:  Jeremy Meredith
//  Creation:    July 26, 2001
//
// ****************************************************************************

bool
LoadBalancer::ContinueExecute(int index)
{
    return (! pipelineInfo[index].complete);
}


// ****************************************************************************
//  Function: ReduceCallback
//
//  Purpose:
//      A C function that can be registered with avt pipelines and can serve
//      as a callback when a pipeline needs to be balanced.  This should be
//      registered by a LoadBalancer when that object is contructed.
//
//  Arguments:
//      ptr      A pointer to the load balancer.
//      spec     The pipeline specification to balance.
//
//  Returns:     The data specification.
//
//  Programmer:  Hank Childs
//  Creation:    June 17, 2001
//
// ****************************************************************************

static avtDataRequest_p
ReduceCallback(void *ptr, avtContract_p spec)
{
    LoadBalancer *lb = (LoadBalancer *) ptr;
    return lb->Reduce(spec);
}


// ****************************************************************************
//  Function: DynamicCheckerCallback
//
//  Purpose:
//      A C function that can be registered with avt pipelines and can serve
//      as a callback when one needs to determine if the loadbalancer will
//      perform dynamic load balancing.  This should be registered by a 
//      LoadBalancer when that object is contructed.
//
//  Arguments:
//      ptr      A pointer to the load balancer.
//      spec     The pipeline specification to balance.
//
//  Returns:     The data specification.
//
//  Programmer:  Jeremy Meredith
//  Creation:    September 19, 2001
//
//  Modifications:
//
//    Hank Childs, Tue Feb 19 19:45:43 PST 2008
//    Rename "dynamic" to "streaming", since we really care about whether we
//    are streaming, not about whether we are doing dynamic load balancing.
//    And the two are no longer synonymous.
//
// ****************************************************************************

static bool
StreamingCheckerCallback(void *ptr, avtContract_p spec)
{
    LoadBalancer *lb = (LoadBalancer *) ptr;
    bool rv = lb->CheckDynamicLoadBalancing(spec);
    return rv;
}


// ****************************************************************************
//  Function: ContinueCallback
//
//  Purpose:
//      A C function that can be registered with avt pipelines and can serve as
//      a callback when a pipeline is not sure if it should continue executing.
//      This should be registered by a LoadBalancer when that object is
//      constructed.
//
//  Arguments:
//      ptr       A pointer to the load balancer.
//      index     A pipeline index.
//
//  Returns:      true if the object should continue executing.
//
//  Programmer:   Hank Childs
//  Creation:     June 17, 2001
//
// ****************************************************************************

static bool
ContinueCallback(void *ptr, int index)
{
    LoadBalancer *lb = (LoadBalancer *) ptr;
    return lb->ContinueExecute(index);
}



