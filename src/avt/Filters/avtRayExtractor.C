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


#ifdef PARALLEL
    sendCells        = true;
#else
    sendCells        = false;
#endif
    jitter           = false;
    rayfoo           = NULL;

    rectilinearGridsAreInWorldSpace = false;
    aspect = 1.;

    shouldDoTiling = false;



    SetTrilinear(false);

    shouldSetUpArbitrator    = false;
    arbitratorPrefersMinimum = false;
    arbitrator               = NULL;

    patchCount = 0;
    totalAssignedPatches = 0;

    rayCastingSLIVR = false;
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
    std::cout << "          In Execute" << std::endl;

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

    avtDataTree_p tree = GetInputDataTree();
    totalNodes = tree->GetNumberOfLeaves();
    currentNode = 0;
    ExecuteTree(tree);

    visitTimer->StopTimer(timingsIndex, "Ray point extraction");
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
    avtSamplePoints_p output = GetTypedOutput();
    // if (kernelBasedSampling)
    //     output->SetUseWeightingScheme(true);

    //
    // This will always be NULL the first time through.  For subsequent tiles
    // (provided we are doing tiling) will not have this issue.
    //
    if (output->GetVolume() == NULL)
        output->SetVolume(width, height, depth);
    else
        output->GetVolume()->ResetSamples();
    output->ResetCellList();
    avtVolume *volume = output->GetVolume();
    if (shouldDoTiling)
        volume->Restrict(width_min, width_max-1, height_min, height_max-1);


    if (massVoxelExtractor != NULL)
    {
        delete massVoxelExtractor;
    }


    //
    // Set up the extractors and tell them which cell list to use.
    //
    avtCellList *cl = output->GetCellList();
   
    massVoxelExtractor = new avtMassVoxelExtractor(width, height, depth, volume,cl);
    massVoxelExtractor->SetTrilinear(trilinearInterpolation);
    massVoxelExtractor->SetRayCastingSLIVR(rayCastingSLIVR);
    massVoxelExtractor->SetJittering(jitter);


    if (shouldDoTiling)
    {
       
        massVoxelExtractor->Restrict(width_min, width_max-1,
                                     height_min, height_max-1);
    }
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
    avtDatasetToSamplePointsFilter::PreExecute();

    if (shouldSetUpArbitrator)
    {
        avtSamplePoints_p samples = GetTypedOutput();

        std::cout << "avtRayExtractor::GetTypedOutput : " << samples->GetNumberOfRealVariables() << std::endl;


        int nvars = samples->GetNumberOfRealVariables();
        int theMatch = -1;
        int tmpIndex = 0;
        for (int i = 0 ; i < nvars ; i++)
        {
            bool foundMatch = false;
            if (samples->GetVariableName(i) == arbitratorVarName)
                foundMatch = true;

            if (foundMatch)
            {
                theMatch = tmpIndex;
                break;
            }
            else
                tmpIndex += samples->GetVariableSize(i);
        }

        if (theMatch != -1)
        {
            arbitrator = new avtRelativeValueSamplePointArbitrator(
                                      arbitratorPrefersMinimum, tmpIndex);
            avtRay::SetArbitrator(arbitrator);
        }
    }

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
    avtDatasetToSamplePointsFilter::PostExecute();

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
    if (*dt == NULL)
    {
        return;
    }
    if (dt->GetNChildren() <= 0 && (!(dt->HasData())))
    {
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
    
    if (rayCastingSLIVR == true)
        if ((totalAssignedPatches != 0) && (dt->ChildIsPresent(0) && !( *(dt->GetChild(0)) == NULL))){
        }else
            totalAssignedPatches = 0;

    for (int i = 0; i < totalAssignedPatches; i++) {
        if (dt->ChildIsPresent(i) && !( *(dt->GetChild(i)) == NULL))
        {
            avtDataTree_p child = dt->GetChild(i);

            //
            // Get the dataset for this leaf in the tree.
            //
            vtkDataSet *ds = child->GetDataRepresentation().GetDataVTK();

            //
            // Iterate over all cells in the mesh and call the appropriate 
            // extractor for each cell to get the sample points.
            //
            // if (kernelBasedSampling)
            //     KernelBasedSample(ds);
            // else
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

    tempImgData.procId = it->second.procId;
    tempImgData.patchNumber = it->second.patchNumber;
    memcpy(tempImgData.imagePatch,it->second.imagePatch,imageMetaPatchVector[patchId].dims[0] * 4 * imageMetaPatchVector[patchId].dims[1] * sizeof(float));

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
    // if (modeIs3D && ds->GetDataObjectType() == VTK_RECTILINEAR_GRID)
    // {
        avtDataAttributes &atts = GetInput()->GetInfo().GetAttributes();
        const double *xform = NULL;
        if (atts.GetRectilinearGridHasTransform())
            xform = atts.GetRectilinearGridTransform();
        massVoxelExtractor->SetGridsAreInWorldSpace(
           rectilinearGridsAreInWorldSpace, viewInfo, aspect, xform);
        avtSamplePoints_p samples = GetTypedOutput();
        int numVars = samples->GetNumberOfRealVariables(); 
        std::vector<std::string> varnames;
        std::vector<int>         varsizes;
        for (int i = 0 ; i < numVars ; i++)
        {
            varnames.push_back(samples->GetVariableName(i));
            varsizes.push_back(samples->GetVariableSize(i));
        }


        if (rayCastingSLIVR == true){
            double scRange[2];
            ds->GetScalarRange(scRange);

            double minUsedScalar = transferFn1D->GetMinUsedScalar();
            double maxUsedScalar = transferFn1D->GetMaxUsedScalar();

            debug5 << PAR_Rank() << " ~ " << num << "  |  Scalar range: " << scRange[0] << ", " << scRange[1] << "   used scalar: " << minUsedScalar << ", " << maxUsedScalar << " || Tests: " << (scRange[1] < minUsedScalar && scRange[0] < minUsedScalar) << " , " << (scRange[0] > maxUsedScalar && scRange[1] > maxUsedScalar) << endl;
        
            if (scRange[1] < minUsedScalar && scRange[0] < minUsedScalar)
                return;

            if (scRange[0] > maxUsedScalar && scRange[1] > maxUsedScalar)
                return;
        }

        massVoxelExtractor->setProcIdPatchID(PAR_Rank(),num);
        massVoxelExtractor->Extract((vtkRectilinearGrid *) ds,
                                    varnames, varsizes);
        if (rayCastingSLIVR == true){
            imgMetaData tmpImageMetaPatch;
            tmpImageMetaPatch = initMetaPatch(patchCount);

            massVoxelExtractor->getImageDimensions(tmpImageMetaPatch.inUse, tmpImageMetaPatch.dims, tmpImageMetaPatch.screen_ll, tmpImageMetaPatch.screen_ur, tmpImageMetaPatch.avg_z);
            if (tmpImageMetaPatch.inUse == 1){
                double bounds[6];
                ds->GetBounds(bounds);
                
                tmpImageMetaPatch.destProcId = tmpImageMetaPatch.procId;
                for (int i=0; i<6; i++)
                    tmpImageMetaPatch.extents[i] = bounds[i];
                imageMetaPatchVector.push_back(tmpImageMetaPatch);

                //std::cout << tmpImageMetaPatch.procId << " ~ " << tmpImageMetaPatch.patchNumber << "  |   bounds: " << bounds[0] << ", " << bounds[1] << "   ;   " << bounds[2] << ", " << bounds[3] << "   ;   " << bounds[4] << ", " << bounds[5] << std::endl;

                
                imgData tmpImageDataHash;
                tmpImageDataHash.procId = tmpImageMetaPatch.procId;           tmpImageDataHash.patchNumber = tmpImageMetaPatch.patchNumber;         tmpImageDataHash.imagePatch = NULL;
                tmpImageDataHash.imagePatch = new float[(tmpImageMetaPatch.dims[0]*4)*tmpImageMetaPatch.dims[1]];

                massVoxelExtractor->getComputedImage(tmpImageDataHash.imagePatch);
                imgDataHashMap.insert( std::pair<int, imgData> (tmpImageDataHash.patchNumber , tmpImageDataHash) );

                patchCount++;
            }
        }

        // if (rayCastingSLIVR == true){
        //     imgMetaData tmpImageMetaPatch;
        //     tmpImageMetaPatch = initMetaPatch(patchCount);

        //     massVoxelExtractor->getImageDimensions(tmpImageMetaPatch.inUse, tmpImageMetaPatch.dims, tmpImageMetaPatch.screen_ll, tmpImageMetaPatch.screen_ur, tmpImageMetaPatch.avg_z);
        //     if (tmpImageMetaPatch.inUse == 1){
        //         tmpImageMetaPatch.destProcId = tmpImageMetaPatch.procId;
        //         imageMetaPatchVector.push_back(tmpImageMetaPatch);
                
        //         imgData tmpImageDataHash;
        //         tmpImageDataHash.procId = tmpImageMetaPatch.procId;           tmpImageDataHash.patchNumber = tmpImageMetaPatch.patchNumber;         tmpImageDataHash.imagePatch = NULL;
        //         tmpImageDataHash.imagePatch = new float[(tmpImageMetaPatch.dims[0]*4)*tmpImageMetaPatch.dims[1]];

        //         massVoxelExtractor->getComputedImage(tmpImageDataHash.imagePatch);
        //         imgDataHashMap.insert( std::pair<int, imgData> (tmpImageDataHash.patchNumber , tmpImageDataHash) );

        //         patchCount++;
        //     }
        // }

        //debug5 << PAR_Rank() << " ~ avtRayExtractor::RasterBasedSample  .. .  " << "  patchCount:" << patchCount << endl;

        return;
    //}
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

void
avtRayExtractor::ExecuteRayTracer(bool parallelOn, int screen0, int screen1, unsigned char bg0, unsigned char bg1, unsigned char bg2){

}

//when you wake up, run with everything commented out, and check if it works
// if not, then try removing arguments, one by one

avtImage_p
avtRayExtractor::ExecuteRayTracer(bool parallelOn, int screen0, int screen1, 
                                  unsigned char bg0, unsigned char bg1, unsigned char bg2,
                                  int timingVolToImg, int timingIndex){

    int screen[2];
    screen[0] = screen0;
    screen[1] = screen1;

    unsigned char background[3];
    background[0] = bg0;
    background[1] = bg1;
    background[2] = bg2;

    debug5 << "parallelOn: " << parallelOn << std::endl;
    // std::cout << "screen: " << screen[0] << " " << screen[1] << std::endl;
    // std::cout << "screen: " << background[0] << " " << background[1] << " " << background[2] << std::endl;
    // std::cout << "patch size: " << getImgPatchSize() << std::endl;

    if (avtCallback::UseKdTreeLoadBalancer() == true)
        std::cout << "Used kd tree load balancer" << std::endl;
    else
        std::cout << "NOT used kd tree load balancer" << std::endl;

    if (avtCallback::UseAMR() == true)
        std::cout << "AMR" << std::endl;
    else
        std::cout << "NOT AMR" << std::endl;

    //
    // Single Processor
    //
    if (parallelOn == false){
        
        std::cout << "In serial mode" << std::endl;
        //
        // Get the metadata
        //
        std::vector<imgMetaData> allImgMetaData;          // contains the metadata to composite the image
        int numPatches = getImgPatchSize();     // get the number of patches
    
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

        //return;

        //return imgTest;
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
    std::cout << PAR_Rank() << "   avtRayTracer::Execute  - imgComm.patchAllocationLogic() " << endl;

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

    std::cout << PAR_Rank() << " ~  num patches to send: " << totalSendData/2 << " num processors to recv from: " << totalRecvData/2 << "    numZDivisions: " << numZDivisions << "   totalPatchesToCompositeLocally: " <<   totalPatchesToCompositeLocally << endl;


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
    float *buffer = new float[((imgBufferWidth * imgBufferHeight * 4)) * numZDivisions]();  //size: imgBufferWidth * imgBufferHeight * 4, initialized to 0

    std::sort(allImgMetaData.begin(), allImgMetaData.end(), &sortImgMetaDataByDepthCopy);

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

    allImgMetaData.clear();
    imgDataToCompose.clear();

    debug5 << PAR_Rank() << "   ~ composing patch done: " << endl;


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

        debug5 << PAR_Rank() << "   ~ final: " << endl;
        std::cout << PAR_Rank() << "   ~ final: " << std::endl;

        if (zbuffer != NULL)
            delete []zbuffer;
    }
    imgComm.syncAllProcs();

    if (PAR_Rank() == 0)
        tempImage->Copy(*whole_image);
    //SetOutput(tempImage);
    


    //**************************************





    visitTimer->StopTimer(timingCompositinig, "Compositing");
    visitTimer->DumpTimings();

    visitTimer->StopTimer(timingIndex, "Ray Tracing");
    visitTimer->DumpTimings();

    // if (imgTest != NULL)
    //   delete []imgTest;

    //delImgPatches();

    return tempImage;

    //return imgTest;
    //return;

}

