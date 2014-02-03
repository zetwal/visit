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
//                          avtRayExtractor.C                        //
// ************************************************************************* //

#include <avtRayExtractor.h>

#include <float.h>

#include <vtkCellData.h>
#include <vtkDataSet.h>
#include <vtkHexahedron.h>
#include <vtkPixel.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPyramid.h>
#include <vtkQuad.h>
#include <vtkQuadraticHexahedron.h>
#include <vtkTetra.h>
#include <vtkTriangle.h>
#include <vtkUnsignedCharArray.h>
#include <vtkVoxel.h>
#include <vtkWedge.h>
#include <vtkImageData.h>

#include <avtCellList.h>
#include <avtDatasetExaminer.h>
#include <avtHexahedronExtractor.h>
#include <avtHexahedron20Extractor.h>
#include <avtMassVoxelExtractor.h>
#include <avtParallel.h>
#include <avtPointExtractor.h>
#include <avtPyramidExtractor.h>
#include <avtRayFunction.h>
#include <avtRelativeValueSamplePointArbitrator.h>
#include <avtSamplePoints.h>
#include <avtTetrahedronExtractor.h>
#include <avtVolume.h>
#include <avtWedgeExtractor.h>
#include <avtCallback.h>

#include <DebugStream.h>
#include <InvalidCellTypeException.h>
#include <TimingsManager.h>

#include <Utility.h>
#include <DebugStream.h>

#include <avtImgCommunicator.h>
#include <avtImageSource.h>
#include <avtDatasetToImageFilter.h>

bool sortImgMetaDataByDepthCopy(imgMetaData const& before, imgMetaData const& after){ return before.avg_z > after.avg_z; }
bool sortImgMetaDataByDepthLargest(imgMetaData const& before, imgMetaData const& after){ return before.avg_z < after.avg_z; }

bool sortImgByCoordinatesX(imgMetaData const& before, imgMetaData const& after){ return (before.screen_ll[0] < after.screen_ll[0]); }
bool sortImgByCoordinatesLastX(imgMetaData const& before, imgMetaData const& after){ return (before.screen_ur[0] > after.screen_ur[0]); }
bool sortImgByCoordinatesY(imgMetaData const& before, imgMetaData const& after){ return (before.screen_ll[1] < after.screen_ll[1]); }
bool sortImgByCoordinatesLastY(imgMetaData const& before, imgMetaData const& after){ return (before.screen_ur[1] > after.screen_ur[1]); }

// ****************************************************************************
//  Method: avtRayExtractor constructor
//
//  Arguments:
//      w       The width.
//      h       The height.
//      d       The depth.
//
//  Programmer: Hank Childs
//  Creation:   December 5, 2000
//     
//  Modifications:
//
//    Hank Childs, Thu Nov 15 15:39:48 PST 2001
//    Moved construction of cell list to Execute to account new limitations of
//    sample points involving multiple variables.
//
//    Hank Childs, Tue Jan  1 10:01:20 PST 2002
//    Initialized sendCells.
//
//    Hank Childs, Sun Dec 14 11:07:56 PST 2003
//    Initialized massVoxelExtractor.
//
//    Hank Childs, Fri Nov 19 13:57:02 PST 2004
//    Initialized rectilinearGridsAreInWorldSpace.
//
//    Hank Childs, Fri Dec 10 09:59:57 PST 2004
//    Initialized shouldDoTiling.
//
//    Hank Childs, Wed Feb  2 08:56:00 PST 2005
//    Initialize modeIs3D.
//
//    Hank Childs, Sun Dec  4 19:12:42 PST 2005
//    Initialize kernelBasedSampling.
//
//    Hank Childs, Tue Jan 24 16:42:40 PST 2006
//    Added point extractor.
//
//    Timo Bremer, Thu Sep 13 14:02:40 PDT 2007
//    Added hex20 extractor.
//
//    Hank Childs, Tue Jan 15 14:26:06 PST 2008
//    Initialize members for sample point arbitration.
//
//    Hank Childs, Fri Jan  9 14:10:25 PST 2009
//    Initialize jitter.
//
// ****************************************************************************

avtRayExtractor::avtRayExtractor(int w, int h, int d)
{
    width  = w;
    height = h;
    depth  = d;

    currentNode = 0;
    totalNodes  = 0;

    massVoxelExtractor  = NULL;

    jitter           = false;

    rectilinearGridsAreInWorldSpace = false;
    aspect = 1.;

    shouldDoTiling = false;
    modeIs3D = true;
    rootGathersAll = false;

    shouldSetUpArbitrator    = false;
    arbitratorPrefersMinimum = false;
    arbitrator               = NULL;

    patchCount = 0;
    totalAssignedPatches = 0;

    imgComm.setDoneVolumeRendering(true);


    rayCastingSLIVR = false;
    trilinearInterpolation = false;
    lighting = false;
    lightPosition[0] = lightPosition[1] = lightPosition[2] = 0.0;   lightPosition[3] = 1.0;
    materialProperties[0] = 0.4; materialProperties[1] = 0.75; materialProperties[3] = 0.0; materialProperties[3] = 15.0;
}


// ****************************************************************************
//  Method: avtRayExtractor destructor
//
//  Programmer: Hank Childs
//  Creation:   December 8, 2000
//      
//  Modifications:
//
//    Hank Childs, Sun Dec 14 11:07:56 PST 2003
//    Deleted massVoxelExtractor.
//
//    Hank Childs, Tue Jan 24 16:42:40 PST 2006
//    Deleted pointExtractor.
//
//    Timo Bremer, Thu Sep 13 14:02:40 PDT 2007
//    Deleted hex20Extractor.
//
//    Hank Childs, Tue Jan 15 21:25:01 PST 2008
//    Delete arbitrator.
//
// ****************************************************************************

avtRayExtractor::~avtRayExtractor()
{

    if (massVoxelExtractor != NULL)
    {
        delete massVoxelExtractor;
        massVoxelExtractor = NULL;
    }

    if (arbitrator != NULL)
    {
        delete arbitrator;
        arbitrator = NULL;
    }

    delImgPatches();
}




// ****************************************************************************
//  Method: avtRayExtractor::SetRectilinearGridsAreInWorldSpace
//
//  Purpose:
//      Tells this object that any input rectilinear grids are in world space,
//      not image space.
//
//  Programmer: Hank Childs
//  Creation:   November 19, 2004
//
//  Modifications:
//    Jeremy Meredith, Thu Feb 15 12:14:04 EST 2007
//    Set rectilinearGridsAreInWorldSpace based on the passed in value
//    (which might be false), not to true.
//
// ****************************************************************************

void
avtRayExtractor::SetRectilinearGridsAreInWorldSpace(bool val,
                 const avtViewInfo &v, double a)
{
    rectilinearGridsAreInWorldSpace = val;
    viewInfo = v;
    aspect = a;
}


// ****************************************************************************
//  Method: avtRayExtractor::RestrictToTile
//
//  Purpose:
//      Tells the extractor whether or not it should only sample within a tile.
//
//  Programmer: Hank Childs
//  Creation:   November 19, 2004
//
// ****************************************************************************

void
avtRayExtractor::RestrictToTile(int wmin, int wmax, int hmin, int hmax)
{
    shouldDoTiling = true;
    width_min  = wmin;
    width_max  = wmax;
    height_min = hmin;
    height_max = hmax;
    modified = true;
}


// ****************************************************************************
//  Method: avtRayExtractor::Execute
//
//  Purpose:
//      This is the real execute method that gets the sample points out of a
//      dataset.
//
//  Programmer: Hank Childs
//  Creation:   December 5, 2000
//
//  Modifications:
//
//    Kathleen Bonnell, Sat Apr 21 13:06:27 PDT 2001 
//    Moved major portion of code to recursive Execute method that walks 
//    down the data tree.
//
//    Hank Childs, Wed Jun  6 10:31:04 PDT 2001
//    Removed domain list argument.  Blew away outdated comments.
//
//    Eric Brugger, Mon Nov  5 13:46:04 PST 2001
//    Modified to always compile the timing code.
//
//    Hank Childs, Thu Nov 15 15:39:48 PST 2001
//    Set up the extractors here (instead of constructor), since they must
//    know how many variables they are working on.
//
//    Hank Childs, Mon Nov 19 14:56:40 PST 2001
//    Gave progress while resampling.
//
// ****************************************************************************

void
avtRayExtractor::Execute(void)
{
    //std::cout << PAR_Rank() << " ~  avtRayExtractor::Execute" << std::endl;
    debug5 << " ~  avtRayExtractor::Execute" << std::endl;

    int timingsIndex = visitTimer->StartTimer();

    SetUpExtractors();

    if (rayCastingSLIVR == true){
        massVoxelExtractor->SetLighting(lighting);
        massVoxelExtractor->SetLightDirection(lightDirection);
        massVoxelExtractor->SetMatProperties(materialProperties);
        massVoxelExtractor->SetModelViewMatrix(modelViewMatrix);
        massVoxelExtractor->SetTransferFn(transferFn1D);
        massVoxelExtractor->SetViewDirection(view_direction);
        massVoxelExtractor->SetViewUp(view_up);
        massVoxelExtractor->SetMeshDims(meshMin,meshMax);
    }

    //
    // Determine partition extents
    //
    //std::cout << PAR_Rank() << "   Partition extents: " << currentPartitionExtents[0] << ", " << currentPartitionExtents[1] << ", " << currentPartitionExtents[2] << ", "
    //                                                    << currentPartitionExtents[3] << ", " << currentPartitionExtents[4] << ", " << currentPartitionExtents[5] << std::endl; 

    avtDataTree_p tree = GetInputDataTree();
    totalNodes = tree->GetNumberOfLeaves();
    currentNode = 0;
    partitionExtentsComputationDone = false;
    ExecuteTree(tree);

    visitTimer->StopTimer(timingsIndex, "Ray point extraction");

    // avtImage_p image = ExecuteRayTracer();
    // SetOutputImage(image);
}


// ****************************************************************************
//  Method: avtRayExtractor::SetUpExtractors
//
//  Purpose:
//      Sets up the extractors and tell them which volume to extract into.
//
//  Programmer: Hank Childs
//  Creation:   November 15, 2001
//
//  Modifications:
//
//    Hank Childs, Tue Jan  1 10:01:20 PST 2002
//    Tell the extractors whether they should extract from large cells.
//
//    Hank Childs, Sun Dec 14 11:07:56 PST 2003
//    Set up massVoxelExtractor.
//
//    Hank Childs, Fri Dec 10 09:59:57 PST 2004
//    Do the sampling in tiles if necessary.
//
//    Hank Childs, Sun Dec  4 19:12:42 PST 2005
//    Add support for kernel based sampling.
//
//    Timo Bremer, Thu Sep 13 14:02:40 PDT 2007
//    Added hex20 extractor.
//
//    Hank Childs, Fri Jan  9 14:11:24 PST 2009
//    Tell extractors whether or not to jitter.  Also remove call to 
//    massVoxelExtractor regarding "sendCellsMode", as it does not participate
//    in that mode ... so the call was worthless.
//
// ****************************************************************************

void
avtRayExtractor::SetUpExtractors(void)
{
    avtImage_p output = GetTypedOutput();
    
    if (massVoxelExtractor != NULL)
        delete massVoxelExtractor;

    //
    // Set up the extractors and tell them which cell list to use.
    //
    avtCellList *cl;    // not needed but avtMassVoxelExtractor requires it
    avtVolume *volume;  // not needed but avtMassVoxelExtractor requires it
    massVoxelExtractor = new avtMassVoxelExtractor(width, height, depth, volume,cl);
    massVoxelExtractor->SetTrilinear(trilinearInterpolation);
    massVoxelExtractor->SetRayCastingSLIVR(rayCastingSLIVR);
    massVoxelExtractor->SetJittering(jitter);

    if (shouldDoTiling)
        massVoxelExtractor->Restrict(width_min, width_max-1, height_min, height_max-1);
}


// ****************************************************************************
//  Method: avtRayExtractor::SetUpArbitrator
//
//  Purpose:
//      Tells this module that it should set up an arbitrator.
//
//  Programmer: Hank Childs
//  Creation:   January 15, 2008
//
// ****************************************************************************

void
avtRayExtractor::SetUpArbitrator(std::string &name, bool pm)
{
    arbitratorVarName        = name;
    arbitratorPrefersMinimum = pm;
    shouldSetUpArbitrator    = true;
}


// ****************************************************************************
//  Method: avtRayExtractor::PreExecute
//
//  Purpose:
//      Determines how many points we have if we have a point mesh.  This will
//      allow us to choose an appropriate radius.
//
//  Programmer: Hank Childs
//  Creation:   February 28, 2006
//
//  Modifications:
//
//    Hank Childs, Tue Jan 15 21:23:49 PST 2008
//    Set up the sample point arbitrator.
//
//    Hank Childs, Sat Nov 21 13:29:21 PST 2009
//    Add support for long longs.
//
// ****************************************************************************

void
avtRayExtractor::PreExecute(void)
{
    avtDatasetToImgFilter::PreExecute();

    if (GetInput()->GetInfo().GetAttributes().GetTopologicalDimension() == 0)
    {
        avtDataset_p ds = GetTypedInput();
        VISIT_LONG_LONG nzones = avtDatasetExaminer::GetNumberOfZones(ds);
        VISIT_LONG_LONG total_nzones;
        SumLongLongArrayAcrossAllProcessors(&nzones, &total_nzones, 1);
        
        if (total_nzones == 0)
        {
            point_radius = 0.1;
            return;
        }

        // In image space, the total volume will be 4 (-1->+1 in X,-1->+1 in Y,
        // 0->+1 in Z).  But: we want to treat all dimensions evenly.  So
        // use 8 (doubling Z) and then correct for it later (when we use the
        // number).
        int dim = GetInput()->GetInfo().GetAttributes().GetSpatialDimension();
        double start_vol = (dim == 3 ? 8. : 4.);
        double vol_per_point = start_vol / total_nzones;
        double exp = (dim == 3 ? 0.333333 : 0.5);
        double side_length = pow(vol_per_point, exp) / 2;
        point_radius = side_length * 1.1; // a little extra
    }
}


// ****************************************************************************
//  Method: avtRayExtractor::PostExecute
//
//  Purpose:
//      Unregisters the sample point arbitrator
//
//  Programmer: Hank Childs
//  Creation:   January 15, 2008
//
// ****************************************************************************

void
avtRayExtractor::PostExecute(void)
{
    avtDatasetToImgFilter::PostExecute();

    if (shouldSetUpArbitrator)
    {
        avtRay::SetArbitrator(NULL);
        if (arbitrator != NULL)
        {
            delete arbitrator;
            arbitrator = NULL;
        }
    }
}


// ****************************************************************************
//  Method: avtRayExtractor::ExecuteTree
//
//  Purpose:
//      This is the recursive execute method that gets the sample points 
//      out of a dataset.  
//
//  Arguments:
//      dt      The data tree that should be processed.
//
//  Programmer: Kathleen Bonnell 
//  Creation:   April 21, 2001. 
//
//  Modifications:
//
//    Hank Childs, Wed Jun  6 10:22:48 PDT 2001
//    Renamed ExecuteTree.
//
//    Hank Childs, Tue Jun 19 19:24:39 PDT 2001
//    Put in logic to handle bad data trees.
//
//    Hank Childs, Tue Nov 13 15:22:07 PST 2001
//    Add support for multiple variables.
//
//    Hank Childs, Mon Nov 19 14:56:40 PST 2001
//    Gave progress while resampling.
//
//    Hank Childs, Mon Apr 15 15:34:43 PDT 2002
//    Give clearer error messages.
//
//    Hank Childs, Fri Jul 18 11:42:10 PDT 2003
//    Do not sample ghost zones.  This gives slightly better performance.
//    And ghost zones occasionally have the wrong value (due to problems with
//    the code that produced it).
//
//    Hank Childs, Sun Dec 14 11:07:56 PST 2003
//    Make use of massVoxelExtractor.
//
//    Hank Childs, Fri Aug 27 16:47:45 PDT 2004
//    Rename ghost data arrays.
//
//    Hank Childs, Fri Nov 19 13:57:02 PST 2004
//    If the rectilinear grids are in world space, then let the mass voxel
//    extractor know about it.
//
//    Hank Childs, Sat Jan 29 13:32:54 PST 2005
//    Added 2D extractors.
//
//    Hank Childs, Sun Jan  1 10:53:09 PST 2006
//    Moved raster based sampling into its own method.  Added support for
//    kernel based sampling.
//
//    Manasa Prasad, 
//    Converted the recursive function to iteration
//
// ****************************************************************************

void
avtRayExtractor::ExecuteTree(avtDataTree_p dt)
{
    if (*dt == NULL){
        return;
    }

    if (dt->GetNChildren() <= 0 && (!(dt->HasData()))){
        return;
    }

    unsigned long m_size, m_rss;
    GetMemorySize(m_size, m_rss);

    debug5 << PAR_Rank() << " ~ avtRayExtractor::ExecuteTree  .. .  " 
           << "    Memory use before: " << m_size << "  rss (MB): " << m_rss/(1024*1024) << endl;

    totalAssignedPatches = dt->GetNChildren();

    patchCount = 0;
    imageMetaPatchVector.clear();
    imgDataHashMap.clear();
    
    if ((totalAssignedPatches != 0) && (dt->ChildIsPresent(0) && !( *(dt->GetChild(0)) == NULL))){
    }else
        totalAssignedPatches = 0;

    for (int i = 0; i < totalAssignedPatches; i++) {
        if (dt->ChildIsPresent(i) && !( *(dt->GetChild(i)) == NULL))
        {
            avtDataTree_p child = dt->GetChild(i);

            vtkDataSet *ds = child->GetDataRepresentation().GetDataVTK();
            RasterBasedSample(ds, i);

            UpdateProgress(10*currentNode+9, 10*totalNodes);
            currentNode++;
        }
    }

    debug5 << PAR_Rank() << " ~ avtRayExtractor::ExecuteTree  .. .  " << "  patchCount:" << patchCount << endl;

    GetMemorySize(m_size, m_rss);
    debug5 << PAR_Rank() << " ~ Memory use after: " << m_size << "  rss (MB): " << m_rss/(1024*1024)
           <<  "   ... avtRayExtractor::ExecuteTree done@!!!" << endl;
    
}

//
// Previous recursive equivalent
//
// void
// avtRayExtractor::ExecuteTree(avtDataTree_p dt)
// {
//     if (*dt == NULL)
//     {
//         return;
//     }
//     if (dt->GetNChildren() <= 0 && (!(dt->HasData())))
//     {
//         return;
//     }

//     if (dt->GetNChildren() != 0)
//     {
//         for (int i = 0; i < dt->GetNChildren(); i++)
//         {
//             if (dt->ChildIsPresent(i))
//                 ExecuteTree(dt->GetChild(i));
//         }

//         return;
//     }

//     //
//     // Get the dataset for this leaf in the tree.
//     //
//     vtkDataSet *ds = dt->GetDataRepresentation().GetDataVTK();

//     //
//     // Iterate over all cells in the mesh and call the appropriate 
//     // extractor for each cell to get the sample points.
//     //
//     if (kernelBasedSampling)
//         KernelBasedSample(ds);
//     else
//         RasterBasedSample(ds);

//     UpdateProgress(10*currentNode+9, 10*totalNodes);
//     currentNode++;
// }



// ****************************************************************************
//  Method: avtRayExtractor::delImgPatches
//
//  Purpose:
//      allocates space to the pointer address and copy the image generated to it
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
void
avtRayExtractor::delImgPatches(){
    imageMetaPatchVector.clear();

    for (iter_t it=imgDataHashMap.begin(); it!=imgDataHashMap.end(); it++){
        if ((*it).second.imagePatch != NULL)
            delete [](*it).second.imagePatch;

        (*it).second.imagePatch = NULL;
    }
    imgDataHashMap.clear();
}



// ****************************************************************************
//  Method: avtRayExtractor::getImgData
//
//  Purpose:
//      copies a patchover
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
void 
avtRayExtractor::getnDelImgData(int patchId, imgData &tempImgData){
    iter_t it = imgDataHashMap.find(patchId);
    if (it == imgDataHashMap.end())
        std::cout << "########################### ERROR!!!! ################################################" << std::endl;

    tempImgData.procId = it->second.procId;
    tempImgData.patchNumber = it->second.patchNumber;
    memcpy(tempImgData.imagePatch,it->second.imagePatch,imageMetaPatchVector[patchId].dims[0] * 4 * imageMetaPatchVector[patchId].dims[1] * sizeof(float));

    if (it->second.imagePatch != NULL)
        delete [](*it).second.imagePatch;
    it->second.imagePatch = NULL;
}


void 
avtRayExtractor::getImgData(int patchId, imgData &tempImgData){
    iter_t it = imgDataHashMap.find(patchId);

     if (it == imgDataHashMap.end())
        std::cout << "########################### ERROR!!!! ################################################" << std::endl;

    tempImgData.procId = it->second.procId;
    tempImgData.patchNumber = it->second.patchNumber;
    memcpy(tempImgData.imagePatch,it->second.imagePatch,imageMetaPatchVector[patchId].dims[0] * 4 * imageMetaPatchVector[patchId].dims[1] * sizeof(float));
}

void 
avtRayExtractor::delImgData(int patchId){
    iter_t it = imgDataHashMap.find(patchId);

     if (it == imgDataHashMap.end())
        std::cout << "########################### ERROR!!!! ################################################" << std::endl;

    if (it->second.imagePatch != NULL)
        delete [](*it).second.imagePatch;
    it->second.imagePatch = NULL;
}

// ****************************************************************************
//  Method: avtRayExtractor::
//
//  Purpose:
//      allocates space to the pointer address and copy the image generated to it
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
imgMetaData
avtRayExtractor::initMetaPatch(int id){
    imgMetaData temp;
    temp.inUse = 0;
    temp.procId = PAR_Rank();
    temp.destProcId = PAR_Rank();
    temp.patchNumber = id;
    temp.dims[0] = temp.dims[1] = -1;
    temp.screen_ll[0] = temp.screen_ll[1] = -1;
    temp.screen_ur[0] = temp.screen_ur[1] = -1;
    temp.avg_z = -1.0;
    
    return temp;
}


// ****************************************************************************
//  Method: avtRayExtractor::RasterBasedSample
//
//  Purpose:
//      Does raster based sampling.
//
//  Programmer: Hank Childs
//  Creation:   January 1, 2006
//
//  Modifications:
//    Jeremy Meredith, Thu Feb 15 11:44:28 EST 2007
//    Added support for rectilinear grids with an inherent transform.
//
//    Hank Childs, Fri Jun  1 12:50:45 PDT 2007
//    Added support for non-scalars.
//
//    Timo Bremer, Thu Sep 13 14:02:40 PDT 2007
//    Added support for hex-20s.
//
//    Hank Childs, Mon Oct 29 20:29:55 PST 2007
//    Ignore surface primitives in 3D.
//
// ****************************************************************************

void
avtRayExtractor::RasterBasedSample(vtkDataSet *ds, int num)
{
    if (ds->GetDataObjectType() == VTK_RECTILINEAR_GRID){

        bool completelyInside = true;

        //
        // Check if that patch is completely inside or not
        double patchExtents[6], centroid[3];
        ds->GetBounds(patchExtents);
        centroid[0] = (patchExtents[0] + patchExtents[1])/2.0;
        centroid[1] = (patchExtents[2] + patchExtents[3])/2.0;
        centroid[2] = (patchExtents[4] + patchExtents[5])/2.0;

        if (centroid[0] >= currentPartitionExtents[0] && centroid[0]<currentPartitionExtents[1])
            if (centroid[1] >= currentPartitionExtents[2] && centroid[0]<currentPartitionExtents[3])
                if (centroid[2] >= currentPartitionExtents[4] && centroid[0]<currentPartitionExtents[5])
                    completelyInside = false;

        //
        // if AMR Check which level patch it is (does it have parents) 

        //std::cout << PAR_Rank() << " ~ " << avtCallback::UseTestVal() << std::endl;
        //std::cout << PAR_Rank() << " ~ " << num << "   Patch extents: " << patchExtents[0] << ", " << patchExtents[1] << "   " << patchExtents[2] << ", " << patchExtents[3] << "   " << patchExtents[4] << ", " << patchExtents[5] << std::endl; 

        avtDataAttributes &atts = GetInput()->GetInfo().GetAttributes();
        const double *xform = NULL;
        if (atts.GetRectilinearGridHasTransform())
            xform = atts.GetRectilinearGridTransform();
        massVoxelExtractor->SetGridsAreInWorldSpace(
           rectilinearGridsAreInWorldSpace, viewInfo, aspect, xform);

        std::vector<std::string> varnames;
        std::vector<int>         varsizes;
        varnames.push_back(varName);
        varsizes.push_back(1);

        double scRange[2];
        ds->GetScalarRange(scRange);

        double minUsedScalar = transferFn1D->GetMinUsedScalar();
        double maxUsedScalar = transferFn1D->GetMaxUsedScalar();

        debug5 << PAR_Rank() << " ~ " << num << "  |  Scalar range: " << scRange[0] << ", " << scRange[1] << "   used scalar: " << minUsedScalar << ", " << maxUsedScalar << " || Tests: " << (scRange[1] < minUsedScalar && scRange[0] < minUsedScalar) << " , " << (scRange[0] > maxUsedScalar && scRange[1] > maxUsedScalar) << endl;
    
        if (scRange[1] < minUsedScalar && scRange[0] < minUsedScalar)
            return;

        if (scRange[0] > maxUsedScalar && scRange[1] > maxUsedScalar)
            return;
        
         if (partitionExtentsComputationDone == false){
            int minPos[2], maxPos[2];
            float minDepth, maxDepth;
            double minExtents[4], maxExtents[4];
            minExtents[0]=currentPartitionExtents[0];
            minExtents[1]=currentPartitionExtents[1];
            minExtents[2]=currentPartitionExtents[2];
            minExtents[3]=1.0;

            maxExtents[0]=currentPartitionExtents[3];
            maxExtents[1]=currentPartitionExtents[4];
            maxExtents[2]=currentPartitionExtents[5];
            maxExtents[3]=1.0;

            massVoxelExtractor->world_to_screen(minExtents, width,height, minPos,minDepth);
            massVoxelExtractor->world_to_screen(maxExtents, width,height, maxPos,maxDepth);

            //std::cout << PAR_Rank() << " ~ " << minExtents[0] << ", " << minExtents[1] << ", " << minExtents[2] << "  -  " << maxExtents[0] << ", " << maxExtents[1] << ", " << maxExtents[2] << "     screen width, height: " << width << ", " << height  << "   Partitions  min: " << minPos[0] << ", " << minPos[1] << ", " << minDepth << "     max: " << maxPos[0] << ", " << maxPos[1] << ", " << maxDepth << std::endl;
            partitionExtentsComputationDone = true;
        }
        massVoxelExtractor->SetPartitionExtents(currentPartitionExtents);  // minX, minY, minZ,    maxX, maxY,maxZ
        massVoxelExtractor->SetLogicalBounds(logicalBounds[0],logicalBounds[1],logicalBounds[2]);
        massVoxelExtractor->setProcIdPatchID(PAR_Rank(),num);
        massVoxelExtractor->SetAMR(avtCallback::UseAMR());
        massVoxelExtractor->Extract((vtkRectilinearGrid *) ds, varnames, varsizes);

        imgMetaData tmpImageMetaPatch;
        tmpImageMetaPatch = initMetaPatch(patchCount);



        massVoxelExtractor->getImageDimensions(tmpImageMetaPatch.inUse, tmpImageMetaPatch.dims, tmpImageMetaPatch.screen_ll, tmpImageMetaPatch.screen_ur, tmpImageMetaPatch.avg_z, tmpImageMetaPatch.fullyInside);
        if (tmpImageMetaPatch.inUse == 1){

            tmpImageMetaPatch.destProcId = tmpImageMetaPatch.procId;
            for (int i=0; i<6; i++)
                tmpImageMetaPatch.extents[i] = patchExtents[i];
            imageMetaPatchVector.push_back(tmpImageMetaPatch);

            imgData tmpImageDataHash;
            tmpImageDataHash.procId = tmpImageMetaPatch.procId;           tmpImageDataHash.patchNumber = tmpImageMetaPatch.patchNumber;         tmpImageDataHash.imagePatch = NULL;
            tmpImageDataHash.imagePatch = new float[(tmpImageMetaPatch.dims[0]*4)*tmpImageMetaPatch.dims[1]];

            massVoxelExtractor->getComputedImage(tmpImageDataHash.imagePatch);
            imgDataHashMap.insert( std::pair<int, imgData> (tmpImageDataHash.patchNumber , tmpImageDataHash) );

            patchCount++;

            //std::string imgFilename_Final = "/home/pascal/Desktop/imgTests/_"+ NumbToString(tmpImageDataHash.procId) + "_" + NumbToString(tmpImageDataHash.patchNumber) + ".ppm";
            //createPpm(tmpImageDataHash.imagePatch, tmpImageMetaPatch.dims[0], tmpImageMetaPatch.dims[1], imgFilename_Final);
        }
        
    }
    return;
}


// ****************************************************************************
//  Method: avtRayExtractor::ExecuteRayTracer
//
//  Purpose:
//      Does SLIVR ray tracing.
//
//  Programmer: Pascal Grosset and Manasa Prasad
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
avtImage_p
avtRayExtractor::ExecuteRayTracer(){
    int  timingVolToImg;

    if (parallelOn)
        timingVolToImg = visitTimer->StartTimer();

    //int screen[2];
    // screen[0] = screen0;
    // screen[1] = screen1;

    // unsigned char background[3];
    // background[0] = bg0;
    // background[1] = bg1;
    // background[2] = bg2;

    debug5 << "parallelOn: " << parallelOn << std::endl;
    // std::cout << "screen: " << screen[0] << " " << screen[1] << std::endl;
    // std::cout << "screen: " << background[0] << " " << background[1] << " " << background[2] << std::endl;
    // std::cout << "patch size: " << getImgPatchSize() << std::endl;

    if (avtCallback::UseKdTreeLoadBalancer() == true)
        debug5 << "Used kd tree load balancer" << std::endl;
    else
        debug5 << "NOT used kd tree load balancer" << std::endl;

    if (avtCallback::UseAMR() == true)
        debug5 << "AMR" << std::endl;
    else
        debug5 << "NOT AMR" << std::endl;

    //
    // Single Processor
    //
    if (parallelOn == false){
        
        debug5 << "In serial mode" << std::endl;
        //
        // Get the metadata
        //
        std::vector<imgMetaData> allImgMetaData;    // contains the metadata to composite the image
        int numPatches = getImgPatchSize();         // get the number of patches
    
        for (int i=0; i<numPatches; i++){
            imgMetaData temp;
            temp = getImgMetaPatch(i);
            allImgMetaData.push_back(temp);
        }
        debug5 << "Number of patches: " << numPatches << std::endl;


        //
        // Sort with the largest z first
        //
        std::sort(allImgMetaData.begin(), allImgMetaData.end(), &sortImgMetaDataByDepthCopy);


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
            getnDelImgData(currentPatch.patchNumber, tempImgData);

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

        //SetOutput(whole_image);
        return whole_image;
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
    int numPatches = getImgPatchSize();     // get the number of patches
    imgComm.gatherNumPatches(numPatches);

    debug5 << PAR_Rank() << "   avtRayTracer::Execute  - Getting the patches -    numPatches: " << numPatches << "   total assigned: " << getTotalAssignedPatches() << endl;
    //std::cout << PAR_Rank() << "   avtRayTracer::Execute  - Getting the patches -    numPatches: " << numPatches << "   total assigned: " << getTotalAssignedPatches() << endl;


    //
    // Send/Receive the patches iota metadata to proc 0
    //

    float *tempSendBuffer = NULL;
    tempSendBuffer = new float[numPatches*7];

    std::multimap<int, imgMetaData> imgMetaDataMultiMap;
    imgMetaDataMultiMap.clear();
    for (int i=0; i<numPatches; i++){
        imgMetaData temp;
        temp = getImgMetaPatch(i);
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

    if (tempSendBuffer != NULL)
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
    //std::cout << PAR_Rank() << "   avtRayTracer::Execute  - imgComm.patchAllocationLogic() " << endl;

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
    //std::cout << PAR_Rank() << " ~  num patches to send: " << totalSendData/2 << " num processors to recv from: " << totalRecvData/2 << "    numZDivisions: " << numZDivisions << "   totalPatchesToCompositeLocally: " <<   totalPatchesToCompositeLocally << endl;


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

            //std::cout << PAR_Rank() << " ~ imgBufferWidth: " << imgBufferWidth << "   imgBufferHeight: " << imgBufferHeight << std::endl;
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
                getnDelImgData(currentPatch.patchNumber, tempImgData);

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
            getnDelImgData(tempImgMetaData.patchNumber, tempImgData);

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
                        getnDelImgData(tempImgMetaData.patchNumber, tempImgData);
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
                        getnDelImgData(tempImgMetaData.patchNumber, tempImgData);
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
    float *buffer = NULL;
    buffer = new float[((imgBufferWidth * imgBufferHeight * 4)) * numZDivisions]();  //size: imgBufferWidth * imgBufferHeight * 4, initialized to 0

    std::sort(allImgMetaData.begin(), allImgMetaData.end(), &sortImgMetaDataByDepthCopy);

    std::multimap< std::pair<int,int>, imgData>::iterator itImgData;
    int bufferDivisionIndex = 0;
    int divIndex = 0;
    int totalSize = allImgMetaData.size();

    debug5 << PAR_Rank() << "  ~  totalSize to compose: " << totalSize << "    numZDivisions: " << numZDivisions << endl;

    for (int k=0; k<numZDivisions; k++){
        debug5 << PAR_Rank() << "  ~  division boundaries: " << divisionsArray[k*2 + 0] << " to " << divisionsArray[k*2 + 1] << endl;
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

        debug5 << PAR_Rank() << "  ~  divIndex: " << divIndex << "    composing patch #: " << allImgMetaData[patchIndex].procId << " ,  " << allImgMetaData[patchIndex].patchNumber << "   size: " << allImgMetaData[patchIndex].dims[0] << " x  " << allImgMetaData[patchIndex].dims[1] << "   avg_z: " << allImgMetaData[patchIndex].avg_z << endl;


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

    allImgMetaData.clear();
    imgDataToCompose.clear();

    debug5 << PAR_Rank() << "  ~ composing patch done: " << endl;


    //
    // --- Timing -- 
    visitTimer->StopTimer(timingLocalCompositing2, "Local Comnpositing of received patches");
    visitTimer->DumpTimings();
    
    int  timingRLE = visitTimer->StartTimer();
    


    //
    // RLE Encoding
    //

    float *encoding = NULL;
    int *sizeEncoding = NULL;

    int totalEncodingSize = imgComm.rleEncodeAll(imgBufferWidth,imgBufferHeight, numZDivisions,buffer,  encoding,sizeEncoding);

    if (buffer != NULL)
        delete []buffer;
    buffer = NULL;
    
    debug5 << PAR_Rank() << "   ~ encoding done!  "<< endl;

    //
    // --- Timing -- 
    visitTimer->StopTimer(timingRLE, "RLE ");
    visitTimer->DumpTimings();
    
    int  finalSend = visitTimer->StartTimer();
    

    //
    // Proc 0 recieves and does the final assmebly
    //
    if (PAR_Rank() == 0)
        imgComm.setBackground(background);

    // Gather all the images
    imgComm.gatherEncodingSizes(sizeEncoding, numZDivisions);                                                       // size of images
    imgComm.gatherAndAssembleEncodedImages(screen[0], screen[1], totalEncodingSize*5, encoding, numZDivisions);     // data from each processor

    debug5 << PAR_Rank() << "   ~ gatherEncodingSizes " << endl;


    if (encoding != NULL)
        delete []encoding;
    encoding = NULL;

    if (sizeEncoding != NULL)
        delete []sizeEncoding;
    sizeEncoding = NULL;

    if (divisionsArray != NULL)
        delete []divisionsArray;
    divisionsArray = NULL;
   

    //
    // --- Timing -- 
    visitTimer->StopTimer(finalSend, "Final send ");
    visitTimer->DumpTimings();
    
    visitTimer->StopTimer(timingComm, "Communicating");
    visitTimer->DumpTimings();

    int  timingCompositinig = visitTimer->StartTimer();


    //
    // Compositing
    //

    // create images structures to hold these
    avtImage_p whole_image, tempImage;

    tempImage = new avtImage(this);     // for processors other than proc 0 ; a dummy
    unsigned char *imgTest = NULL;
    
    // Processor 0 does a special compositing
    if (PAR_Rank() == 0)
    {
        whole_image = new avtImage(this);

        float *zbuffer = new float[screen[0] * screen[1]];
        

        // creates input for the
        vtkImageData *img = avtImageRepresentation ::NewImage(screen[0], screen[1]);
        whole_image->GetImage() = img;


        imgTest = new unsigned char[screen[0] * screen[1] * 3];
        imgTest = whole_image->GetImage().GetRGBBuffer();

        zbuffer = new float[screen[0] * screen[1]]();
        for (int s=0; s<screen[0] * screen[1]; s++)
            zbuffer[s] = 20.0;
        zbuffer = whole_image->GetImage().GetZBuffer();

        
        // Get the composited image
        imgComm.getcompositedImage(screen[0], screen[1], imgTest); 
        img->Delete();

        if (zbuffer != NULL)
            delete []zbuffer;
    }
    imgComm.syncAllProcs();

    if (PAR_Rank() == 0)
        tempImage->Copy(*whole_image);

    visitTimer->StopTimer(timingCompositinig, "Compositing");
    visitTimer->DumpTimings();

    return tempImage;
}



// ****************************************************************************
//  Method: LoadBalancer::chopPartitionRT
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
//  Creation:   December 23, 2013
//
// ****************************************************************************
int
avtRayExtractor::chopPartitionRT(partitionExtents parent, partitionExtents & childOne, partitionExtents & childTwo, int axisOrder[3]){
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
    if (count == 3){
    	debug5 << "LoadBalancer::chopPartition Error in kdtree" << std::endl;
        return -1;  // We are going to be cycling forever here! So stop!
    }
    
    childOne.axisIndex = childTwo.axisIndex = axisIndex;
        
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
//  Method: LoadBalancer::getPartitionExtents
//
//  Purpose:
//      Do a k-d tree partitioning
//
//  Arguments:
//      numDivisions      
//      logicalBounds    
//      minSpatialExtents    
//      maxSpatialExtents
//      extents
//
//  Programmer: Pascal Grosset  
//  Creation:   December 9, 2013
//
// ****************************************************************************
void          
avtRayExtractor::getPartitionExtents(int id, int numDivisions, int logicalBounds[3], double minSpatialExtents[3], double maxSpatialExtents[3], double extents[6]){
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

    while (myPartitions.size() != numDivisions){
        parent = myPartitions.front();   myPartitions.pop_front();

        chopPartitionRT(parent,one,two,axisOrder);
        myPartitions.push_back(one);
        myPartitions.push_back(two);
    }


    //
    //Push head to the start of the list: this ensures more contiguity
    bool headfound = false;
    int tempCount = 0;
    do{
        partitionExtents tempExt = myPartitions.front();
        if (tempExt.head == true){
            headfound = true;
            break;
        }
        else{
            // remove the front and put it at the back
            myPartitions.pop_front();
            myPartitions.push_back(tempExt);
        }

        tempCount++;
        if (tempCount > myPartitions.size()){    // should not have to resort to this!!!!
            debug5 << "avtRayExtractor::getPartitionExtents kdtree has had to abort an infinite loop!!!! shouldn't be happening" << std::endl;
            break;
        }
    }while(headfound == false);

    /////////////////////////////////////////////////////////////////
    // Determine which patch is in which region
    extents[0] = myPartitions[id].minExtents[0];  extents[3] = myPartitions[id].maxExtents[0];
    extents[1] = myPartitions[id].minExtents[1];  extents[4] = myPartitions[id].maxExtents[1];
    extents[2] = myPartitions[id].minExtents[2];  extents[5] = myPartitions[id].maxExtents[2];
}


bool 
avtRayExtractor::patchOverlap(float patchMinX, float patchMaxX, float patchMinY, float patchMaxY, float patchMinZ, float patchMaxZ,
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


void createPpmCopyComp(unsigned char array[], int dimx, int dimy, std::string filename){
    int i, j;
    //std::cout << "createPpm2  dims: " << dimx << ", " << dimy << " -  " << filename.c_str() << std::endl;
    FILE *fp = fopen(filename.c_str(), "wb"); // b - binary mode 
    (void) fprintf(fp, "P6\n%d %d\n255\n", dimx, dimy);
    for (j = 0; j < dimy; ++j){
        for (i = 0; i < dimx; ++i){
            static unsigned char color[3];
            color[0] = array[j*(dimx*3) + i*3 + 0] ;  // red
            color[1] = array[j*(dimx*3) + i*3 + 1] ;  // green
            color[2] = array[j*(dimx*3) + i*3 + 2] ;  // blue 
            (void) fwrite(color, 1, 3, fp);
        }
    }
    (void) fclose(fp);
    //std::cout << "End createPpm: " << std::endl;
}


void createPpmRE_RGBA(unsigned char array[], int dimx, int dimy, std::string filename){
    int i, j;
    //std::cout << "createPpm2  dims: " << dimx << ", " << dimy << " -  " << filename.c_str() << std::endl;
    FILE *fp = fopen(filename.c_str(), "wb"); // b - binary mode 
    (void) fprintf(fp, "P6\n%d %d\n255\n", dimx, dimy);
    for (j = 0; j < dimy; ++j){
        for (i = 0; i < dimx; ++i){
            static unsigned char color[3];
            float alpha = array[j*(dimx*4) + i*4 + 3] / 255.0;
            color[0] = array[j*(dimx*4) + i*4 + 0]  * alpha;  // red
            color[1] = array[j*(dimx*4) + i*4 + 1] * alpha;  // green
            color[2] = array[j*(dimx*4) + i*4 + 2] * alpha;  // blue 
            (void) fwrite(color, 1, 3, fp);
        }
    }
    (void) fclose(fp);
    //std::cout << "End createPpm: " << std::endl;
}

avtImage_p
avtRayExtractor::ExecuteRayTracerLB(){
    int timingVolToImg;

    //
    // Get the metadata
    //
    std::vector<imgMetaData> allImgMetaData;    // contains the metadata to composite the image
    int numPatches = getImgPatchSize();         // get the number of patches

    for (int i=0; i<numPatches; i++){
        imgMetaData temp;
        temp = getImgMetaPatch(i);
        allImgMetaData.push_back(temp);
    }
    debug5 << PAR_Rank() << " ~ avtRayTracer::ExecuteRayTracerLB  - Getting the patches - num patches used: " << numPatches << "   total assigned: " << getTotalAssignedPatches() << endl;
    //std::cout << PAR_Rank() << " ~ avtRayTracer::ExecuteRayTracerLB  - Getting the patches - num patches used : " << numPatches << "   total assigned: " << getTotalAssignedPatches() << endl;

    int localNodeCompositingTiming;
    localNodeCompositingTiming = visitTimer->StartTimer();
    
    int imgBufferWidth, imgBufferHeight;
    int startX, startY, endX, endY;
    float avg_z = 0;
    imgBufferWidth = imgBufferHeight = 0;
    if (numPatches > 0){ 
        //
        // Sort to find extents of patches
        std::sort(allImgMetaData.begin(), allImgMetaData.end(), &sortImgByCoordinatesX);
        startX = allImgMetaData[0].screen_ll[0];

        std::sort(allImgMetaData.begin(), allImgMetaData.end(), &sortImgByCoordinatesLastX);
        endX = allImgMetaData[0].screen_ur[0];

        std::sort(allImgMetaData.begin(), allImgMetaData.end(), &sortImgByCoordinatesY);
        startY = allImgMetaData[0].screen_ll[1];

        std::sort(allImgMetaData.begin(), allImgMetaData.end(), &sortImgByCoordinatesLastY);
        endY = allImgMetaData[0].screen_ur[1];

        //
        // Sort with the largest z first
        std::sort(allImgMetaData.begin(), allImgMetaData.end(), &sortImgMetaDataByDepthLargest);

        imgBufferWidth = endX - startX + 1;
        imgBufferHeight = endY - startY + 1;

        imgComm.setHasImageToComposite(true);

        debug5 << PAR_Rank() << " ~ done sorting and extents: screen extents X: " << startX << ", " << endX <<  "  ~~~~  screen extents Y: " << startY << ", " << endY << std::endl;
        //std::cout << PAR_Rank() << " ~ done sorting and extents: screen extents X: " << startX << ", " << endX <<  "  ~~~~  screen extents Y: " << startY << ", " << endY << std::endl;
    }else{
        imgComm.setHasImageToComposite(false);
    }
    imgComm.setDoneVolumeRendering(true);

    //
    // Local Compositing within one processor
    //
    
    // Creates a buffer to store the composited image
    float *localBuffer = NULL;
    localBuffer = new float[imgBufferWidth * imgBufferHeight * 4]();

    int localProcCompsitingTiming;
    localProcCompsitingTiming = visitTimer->StartTimer();

    for (int i=0; i<numPatches; i++){
        imgMetaData currentPatch = allImgMetaData[i];
        imgData tempImgData;
        tempImgData.imagePatch = NULL;
        tempImgData.imagePatch = new float[currentPatch.dims[0] * currentPatch.dims[1] * 4]();

        getImgData(currentPatch.patchNumber, tempImgData);

        int startingX = currentPatch.screen_ll[0] - startX;
        int startingY = currentPatch.screen_ll[1] - startY; 

        if (i == 0)
            avg_z = currentPatch.avg_z;

        for (int j=0; j<currentPatch.dims[1]; j++){
            for (int k=0; k<currentPatch.dims[0]; k++){
                
                if ((startingX + k) > imgBufferWidth) continue;
                if ((startingY + j) > imgBufferHeight) continue;
                
                int subImgIndex = currentPatch.dims[0]*j*4 + k*4;                                           // index in the subimage 
                int bufferIndex = (startingY*imgBufferWidth*4 + j*imgBufferWidth*4) + (startingX*4 + k*4);  // index in the big buffer

                if (localBuffer[bufferIndex+3] > 1.0) continue;
                if (tempImgData.imagePatch[subImgIndex+3] <= 0.0) continue;

                // Front to Back
                localBuffer[bufferIndex+0] = imgComm.clamp( (tempImgData.imagePatch[subImgIndex+0] * (1.0 - localBuffer[bufferIndex+3])) + localBuffer[bufferIndex+0] );
                localBuffer[bufferIndex+1] = imgComm.clamp( (tempImgData.imagePatch[subImgIndex+1] * (1.0 - localBuffer[bufferIndex+3])) + localBuffer[bufferIndex+1] );
                localBuffer[bufferIndex+2] = imgComm.clamp( (tempImgData.imagePatch[subImgIndex+2] * (1.0 - localBuffer[bufferIndex+3])) + localBuffer[bufferIndex+2] );
                localBuffer[bufferIndex+3] = imgComm.clamp( (tempImgData.imagePatch[subImgIndex+3] * (1.0 - localBuffer[bufferIndex+3])) + localBuffer[bufferIndex+3] ); 
            }
        }

        //std::string imgFilename_comp = "/home/pascal/Desktop/imgTests/_composed_ " + NumbToString(PAR_Rank()) + "_"+ NumbToString(i) +"_"+ NumbToString(currentPatch.patchNumber) + "_.ppm";
        //createPpm(localBuffer, imgBufferWidth, imgBufferHeight, imgFilename_comp);

        if (tempImgData.imagePatch != NULL)
            delete []tempImgData.imagePatch;
        tempImgData.imagePatch = NULL;

        delImgData(currentPatch.patchNumber);
    }

    //
    // No longer need patches at this point, so doing some clean up and memory release
    delImgPatches();
    allImgMetaData.clear();
    imageMetaPatchVector.clear();
    imgDataHashMap.clear();

    visitTimer->StopTimer(localProcCompsitingTiming, "Local Proc Compositing Timing");
    visitTimer->DumpTimings();

    //imgComm.syncAllProcs();

    debug5 << "Local compositing done :  num patches: " << numPatches << "   size: " << imgBufferWidth << " x " << imgBufferHeight  << "  avg_z: " <<  avg_z << std::endl;
    //std::cout << PAR_Rank()  << " ~ Local compositing done :  num patches: " << numPatches << "   size: " << imgBufferWidth << " x " << imgBufferHeight  << "  avg_z: " <<  avg_z << std::endl;

    //std::string imgFilename_comp = "/home/pascal/Desktop/imgTests/_proc_ " + NumbToString(PAR_Rank()) + "_.ppm";
    //createPpm(localBuffer, imgBufferWidth, imgBufferHeight, imgFilename_comp);

    //
    // Compositing
    //

	//
    // create images structures to hold these
    avtImage_p whole_image, tempImage;
    whole_image = new avtImage(this);
    
    tempImage = new avtImage(this);     	// for processors other than proc 0 ; a dummy
    float *zbuffer = new float[screen[0] * screen[1]];

    if (PAR_Rank() == 0)
        imgComm.setBackground(background);

    int sendingTags[2] = {15,16};
    if (avtCallback::UseusingIcet() == false){
        if (PAR_Size() > 1){
            if (avtCallback::UseTree() == true){
                //std::cout << PAR_Rank() << " ~ Tree compositing "  << endl;
                debug5 << PAR_Rank() << " ~ Tree compositing "  << endl;
            	debug5 << PAR_Rank() << " ~ Do compositing on one node ...................... " << numPatches << std::endl << std::endl << std::endl;
                //
                // Compositing among contiguous processors on one node
                //
                int compositingNodeTiming;
                compositingNodeTiming = visitTimer->StartTimer();

                std::vector<int>collocatedProcs;
                collocatedProcs.clear();
                for (std::list<int>::iterator it=contiguousMergingProcs.begin(); it != contiguousMergingProcs.end(); ++it)
                    collocatedProcs.push_back((int)*it);
                
                int internalTags[3]={42,43,44};
                imgComm.doNodeCompositing(collocatedProcs, startX, startY, imgBufferWidth, imgBufferHeight, avg_z, localBuffer, internalTags);
                
                visitTimer->StopTimer(compositingNodeTiming, "Compositing Node Timing");
                visitTimer->DumpTimings();
				debug5 << PAR_Rank() << " ~ Done with compositing on one node!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << numPatches << std::endl << std::endl;

                //
                // Compositing across nodes
                //
                debug5 << PAR_Rank() << " ~ Do compositing across nodes ........................ " << numPatches << std::endl << std::endl << std::endl;
                
                int compositingAcrossNodesTiming;
                compositingAcrossNodesTiming = visitTimer->StartTimer();

                    int externalTags[3]={52,53,54};
                    imgComm.doNodeCompositing(processorCompositingOrder, startX, startY, imgBufferWidth, imgBufferHeight, avg_z, localBuffer, externalTags);

                visitTimer->StopTimer(compositingAcrossNodesTiming, "Compositing Across Nodes");
                visitTimer->DumpTimings();

                debug5 << PAR_Rank() << " ~ Done with compositing across nodes!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << numPatches << std::endl << std::endl;
                
                

                //
                // Root Node timing
                //
                debug5 << PAR_Rank() << " ~ Do gather on 0 .............................. " << numPatches << std::endl << std::endl << std::endl;
                int rootNodeFinalTiming;
                rootNodeFinalTiming = visitTimer->StartTimer();
                    imgComm.finalAssemblyOnRoot(screen[0], screen[1], startX, startY, imgBufferWidth, imgBufferHeight, localBuffer, sendingTags);
                visitTimer->StopTimer(rootNodeFinalTiming, "Root Node final timing");
                visitTimer->DumpTimings();

                debug5 << PAR_Rank() << " ~ Done gather on 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!" << numPatches << std::endl << std::endl;
                imgComm.syncAllProcs();
                debug5 << PAR_Rank() << " ~ Tree Compositing - All done. Drawing the image now!" << std::endl << std::endl;
            }
            else{
                debug5 << PAR_Rank() << "  ~ Direct Send compositing "  << endl;
                //
                // RLE Encoding
                //
                float *encoding = NULL;
                int *sizeEncoding = NULL;
                
                int totalEncodingSize = imgComm.rleEncodeAll(imgBufferWidth,imgBufferHeight, 1,localBuffer,  encoding,sizeEncoding);

                debug5 << PAR_Rank() << "  ~ encoding done!  initial size: " << imgBufferWidth * imgBufferHeight * 4 << "    now: " << totalEncodingSize << endl;
                //std::cout << PAR_Rank() << "  ~ encoding done!  initial size: " << imgBufferWidth * imgBufferHeight * 4 << "    now: " << totalEncodingSize << endl;


                //
                // --- Timing -- 
                int  finalSend = visitTimer->StartTimer();
                
                //
                // Proc 0 recieves and does the final assmebly
                //
                imgComm.gatherEncodingSizesLB(sizeEncoding, 1);   

                int dataToSend[4];                  // size of images
                dataToSend[0] = imgBufferWidth;
                dataToSend[1] = imgBufferHeight;
                dataToSend[2] = startX;
                dataToSend[3] = startY;
                imgComm.gatherAndAssembleEncodedImagesLB(screen[0], screen[1], dataToSend, totalEncodingSize*5, encoding, 1, avg_z);     // data from each processor

                debug5 << PAR_Rank() << " ~ Direct Send - All done. Drawing the image now!" << std::endl << std::endl;

                if (encoding != NULL)
                    delete []encoding;
                encoding = NULL;

                if (sizeEncoding != NULL)
                    delete []sizeEncoding;
                sizeEncoding = NULL;

                //
                // --- Timing -- 
                visitTimer->StopTimer(finalSend, "Final send ");
                visitTimer->DumpTimings();
            }
        }else{
            //std::cout << " PAR_Size() 1" << std::endl;
            imgComm.finalAssemblyOnRoot(screen[0], screen[1], startX, startY, imgBufferWidth, imgBufferHeight, localBuffer, sendingTags);
            //std::cout << "Done final assembly " << std::endl;
        }

        int  timingCompositinig = visitTimer->StartTimer();
        unsigned char *imgTest = NULL;
        
        // Processor 0 does a special compositing
        if (PAR_Rank() == 0)
        {
            //std::cout << " Done final assembly 2" << std::endl;
            whole_image = new avtImage(this);

            // creates input for the
            vtkImageData *img = avtImageRepresentation ::NewImage(screen[0], screen[1]);
            whole_image->GetImage() = img;

            imgTest = new unsigned char[screen[0] * screen[1] * 3];
            imgTest = whole_image->GetImage().GetRGBBuffer();

            zbuffer = new float[screen[0] * screen[1]]();
            for (int s=0; s<screen[0] * screen[1]; s++)
                zbuffer[s] = 20.0;

            //std::cout << " Done final assembly 2.5" << std::endl;
            zbuffer = whole_image->GetImage().GetZBuffer();
            //std::cout << " Done final assembly 3" << std::endl;
            // Get the composited image
            imgComm.getcompositedImage(screen[0], screen[1], imgTest); 
            img->Delete();
            //std::cout << " Done final assembly 4" << std::endl;
            debug5 << PAR_Rank() << "   ~ final: " << endl;

            if (zbuffer != NULL)
                delete []zbuffer;
        }
        imgComm.syncAllProcs();

        if (PAR_Rank() == 0)
            tempImage->Copy(*whole_image);

        visitTimer->StopTimer(timingCompositinig, "Compositing");
        visitTimer->DumpTimings();
    } 
    else   // Using iceT
    {
        for (int i = 0 ; i < screen[0] * screen[1] ; i++)
            zbuffer[i] = avg_z;

        vtkImageData *vtk = toVTKImage(localBuffer, imgBufferWidth, imgBufferHeight, 0, 0, -1);
        avtImageRepresentation *vtk_image = new avtImageRepresentation(vtk, zbuffer);

        if(imgBufferWidth == 0 || imgBufferHeight == 0){
            vtk_image->SetOrigin(0, 0);
            vtk_image->SetBoundingSize(1, 1);
        }else{
            vtk_image->SetOrigin(startX, startY);
            vtk_image->SetBoundingSize(imgBufferWidth, imgBufferHeight);
        }

        imgComm.syncAllProcs();

        whole_image->GetImage() = *vtk_image;
        tempImage->Copy(*whole_image);
    }

    if (localBuffer != NULL)
        delete []localBuffer;
    localBuffer = NULL;

    if (zbuffer != NULL)
        delete []zbuffer;
    zbuffer = NULL;

    return tempImage;
}



vtkImageData* 
avtRayExtractor::toVTKImage(float* buffer, int width, int height, int startX, int startY, float avg_z ){
    //
    // Create an image that we can place each pixel into.
    //

    int fullHeight;
    int fullWidth;

    if(avg_z == -1){
      fullWidth = width;//screen[0];
      fullHeight = height;//screen[1];
    }else{
      fullWidth = screen[0];
      fullHeight = screen[1];
    }

    if(fullWidth == 0 || fullHeight == 0){
        fullWidth = 1;
        fullHeight = 1;

        int nPixels = fullWidth*fullHeight;

        vtkImageData *image = avtImageRepresentation::NewImage(fullWidth, fullHeight);
        image->AllocateScalars(VTK_UNSIGNED_CHAR, 4);
        unsigned char *data = (unsigned char *)image->GetScalarPointer(0, 0, 0);

        for (int i = 0 ; i < 4*nPixels ; i++)
            data[i] = 0;

        return image;
    }

    int nPixels = fullWidth*fullHeight;

    debug5 << PAR_Rank() << "\t\t fullWidth: " << fullWidth << "    fullHeight: " << fullHeight << std::endl; 

    vtkImageData *image = avtImageRepresentation::NewImage(fullWidth, fullHeight);
    image->AllocateScalars(VTK_UNSIGNED_CHAR, 4);
    //image->SetNumberOfScalarComponents(4, image->GetInformation());

    //
    // Populate an initial image, either with the background or with an
    // opaque image that is to be inserted into the middle of the rendering.
    //
    unsigned char *data = (unsigned char *)image->GetScalarPointer(0, 0, 0);

    //
    // Now that we have the background in the full image, copy it into what
    // we need for this image.
    //
    for (int i = 0 ; i < 4*nPixels ; i++)
        data[i] = 0;

    //  std::cout << PAR_Rank() <<  "   before writing buffer to data" << std::endl;

    for (int j = 0 ; j < height ; j++){
        for (int k = 0 ; k < width ; k++){
            int index = width*4*j + k*4;

            if ((startX + k) > fullWidth)
                continue;

            if ((startY + j) > fullHeight)
            continue;

            int bufferIndex = (startY*fullWidth*4 + j*fullWidth*4) + (startX*4+ k*4); 

            if(bufferIndex > fullWidth * fullHeight * 4) {
                continue;
            }
            if(index > width * height * 4) continue;

            data[bufferIndex    ] = buffer[index ] * 255;
            data[bufferIndex + 1] = buffer[index +1] * 255;
            data[bufferIndex + 2] = buffer[index +2] * 255;
            data[bufferIndex + 3] = buffer[index +3] * 255;
        }
    }

    avtImageRepresentation *temp = new avtImageRepresentation(image);

    return image;
}


void    
avtRayExtractor::InsertOpaqueImage(avtImage_p img)
{
    opaqueImage = img;
    //std::cout << PAR_Rank() << "        inserting opaque image" << std::endl;
}

void
avtRayExtractor::GetContiguousNodeList()
{
    contiguousMergingProcs.clear();
	int id = PAR_Rank();
    int position = -1;
    int myPos, myId;
    
    std::vector<int>::iterator it;
    it = find(processorCompositingOrder.begin(), processorCompositingOrder.end(), id);
    if (it != processorCompositingOrder.end())
        position = it-processorCompositingOrder.begin();
        
    contiguousMergingProcs.push_back(id);
    //
	// Check if the ones around me are on my node: two directions up and down
	
	// Up
	myPos = position;
	myId = id;
    bool found = false;
	do{
		int nextUp = myPos-1;
		if (nextUp < 0)
			break;
		int nodeId = processorCompositingOrder[nextUp];

		if (imgComm.checkIfProcessorIsOnMyNode(nodeId)){
			contiguousMergingProcs.push_front(nodeId);
			found = true;
			myPos = nextUp;
			myId = nodeId;
		}else
            found = false;

	}while(found == true);
        
    // Down
	myPos = position;
	myId = id;
    found = false;
	do{
		int nextDown = myPos+1;
		if (nextDown >= processorCompositingOrder.size())
			break;
			
		int nodeId = processorCompositingOrder[nextDown];

		if (imgComm.checkIfProcessorIsOnMyNode(nodeId)){
			contiguousMergingProcs.push_back(nodeId);
			found = true;
			myPos = nextDown;
			myId = nodeId;
		}else
            found = false;

	}while(found == true);
	
    std::stringstream ss;
    ss << PAR_Rank() << " ~ Contiguous procs size: " << contiguousMergingProcs.size() << "  patches:  \n";
    for (std::list<int>::iterator it=contiguousMergingProcs.begin(); it != contiguousMergingProcs.end(); ++it)
        ss << *it << "\n";
    //std::cout << ss.str() << std::endl;
    debug5 << ss.str() << std::endl;
}
