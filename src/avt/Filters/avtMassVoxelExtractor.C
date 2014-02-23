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
//                            avtMassVoxelExtractor.C                        //
// ************************************************************************* //

#include <avtMassVoxelExtractor.h>
#include <avtImgCommunicator.h>

#include <float.h>

#include <avtAccessor.h>
#include <avtCellList.h>
#include <avtVolume.h>

#include <vtkDataArray.h>
#include <vtkCamera.h>
#include <vtkCellData.h>
#include <vtkMatrix4x4.h>
#include <vtkPointData.h>
#include <vtkRectilinearGrid.h>
#include <vtkTemplateAliasMacro.h>
#include <vtkUnsignedCharArray.h>

#include <DebugStream.h>
#include <avtCallback.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm> 


void normalize(float vec[3]){
    float inverse_sqrt_sum_squared = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    if (inverse_sqrt_sum_squared != 0)
        inverse_sqrt_sum_squared = 1.0/inverse_sqrt_sum_squared;

    for (int i=0;i<3; i++)
        vec[i] = vec[i]*inverse_sqrt_sum_squared;
}


float dot(float vecA[3], float vecB[3]){
    return ((vecA[0]*vecB[0]) + (vecA[1]*vecB[1]) + (vecA[2]*vecB[2]));
}


// ****************************************************************************
//  Method: avtMassVoxelExtractor constructor
//
//  Arguments:
//     w     The number of sample points in the x direction (width).
//     h     The number of sample points in the y direction (height).
//     d     The number of sample points in the z direction (depth).
//     vol   The volume to put samples into.
//     cl    The cell list to put cells whose sampling was deferred.
//
//  Programmer: Hank Childs
//  Creation:   December 14, 2003
//
//  Modifications:
//
//    Hank Childs, Fri Nov 19 14:50:58 PST 2004
//    Initialize gridsAreInWorldSpace.
//
//    Jeremy Meredith, Thu Feb 15 13:11:34 EST 2007
//    Added an ability to extract voxels using the world-space version
//    even when they're really in image space.
//
//    Hank Childs, Wed Aug 27 11:11:28 PDT 2008
//    Initialize spatial coordinates array.
//
//    Hank Childs, Wed Dec 24 11:22:43 PST 2008
//    Remove reference to ProportionSpaceToZBufferSpace data member.
//
//    Kathleen Biagas, Fri Jul 13 09:23:55 PDT 2012
//    Use double instead of float.
//
// ****************************************************************************

avtMassVoxelExtractor::avtMassVoxelExtractor(int w, int h, int d,
                                             avtVolume *vol, avtCellList *cl)
    : avtExtractor(w, h, d, vol, cl)
{
    fullImgWidth = w;
    fullImgHeight = h;
    gridsAreInWorldSpace = false;
    pretendGridsAreInWorldSpace = false;
    trilinearInterpolation = false;
    rayCastingSLIVR = false;
    isAMR = false;

    aspect = 1;
    view_to_world_transform = vtkMatrix4x4::New();
    world_to_view_transform = vtkMatrix4x4::New();
    X = NULL;
    Y = NULL;
    Z = NULL;
    divisors_X = NULL;
    divisors_Y = NULL;
    divisors_Z = NULL;

    lighting = false;
    fullyInside = false;
    lightPosition[0] = lightPosition[1] = lightPosition[2] = 0.0;   lightPosition[3] = 1.0;
    materialProperties[0] = 0.4; materialProperties[1] = 0.75; materialProperties[3] = 0.0; materialProperties[3] = 15.0;

    proc = patch = 0;
    patchDrawn = 0;
    imgDims[0] = imgDims[1] = 0;                // size of the patch
    imgLowerLeft[0] = imgLowerLeft[1] = 0;      // coordinates in the whole image
    imgUpperRight[0] = imgUpperRight[1] = 0;    //
    imgDepth = -1;                              // from the depth buffer
    imgArray = NULL;                            // the image data


    // threads
    enableThreads = false;
    numThreads = 0;
    threadHandles = NULL;
    threadArgument = NULL;
    allPatchesProcessed = false;
}


// ****************************************************************************
//  Method: avtMassVoxelExtractor destructor
//
//  Purpose:
//      Defines the destructor.  Note: this should not be inlined in the header
//      because it causes problems for certain compilers.
//
//  Programmer: Hank Childs
//  Creation:   February 5, 2004
//
//  Modifications:
//
//    Hank Childs, Sun Nov 21 10:35:40 PST 2004
//    Delete the view to world transform.
//
//    Hank Childs, Wed Aug 27 11:10:51 PDT 2008
//    Delete the spatial coordinate arrays.
//
//    Hank Childs, Wed Dec 24 11:22:43 PST 2008
//    Remove reference to ProportionSpaceToZBufferSpace data member.
//
// ****************************************************************************

avtMassVoxelExtractor::~avtMassVoxelExtractor()
{
    view_to_world_transform->Delete();
    world_to_view_transform->Delete();

    if (X != NULL)
        delete [] X;
    if (Y != NULL)
        delete [] Y;
    if (Z != NULL)
        delete [] Z;
    if (divisors_X != NULL)
        delete [] divisors_X;
    if (divisors_Y != NULL)
        delete [] divisors_Y;
    if (divisors_Z != NULL)
        delete [] divisors_Z;
    if (imgArray != NULL)
        delete []imgArray;

    imgArray = NULL;
}


// ****************************************************************************
//  Function:  AssignEight
//
//  Purpose:
//      A templated function to assign 8 values to a floating point array.
//
//  Programmer: Hank Childs
//  Creation:   August 26, 2008
//
//  Modifications:
//    Kathleen Biagas, Fri Jul 13 09:23:55 PDT 2012
//    Use double instead of float.
//
// ****************************************************************************

template <class T> static void
AssignEight(double *vals, int *index, int s, int m, T *array)
{
    for (int i = 0 ; i < 8 ; i++)
        vals[i] = (double) array[s*index[i]+m];
}

static void
AssignEight(int vartype, double *vals, int *index, int s, int m, void *array)
{
    switch (vartype)
    {
      case VTK_CHAR:
        AssignEight(vals, index, s, m, (char *) array);
        break;
      case VTK_UNSIGNED_CHAR:
        AssignEight(vals, index, s, m, (unsigned char *) array);
        break;
      case VTK_SHORT:
        AssignEight(vals, index, s, m, (short *) array);
        break;
      case VTK_UNSIGNED_SHORT:
        AssignEight(vals, index, s, m, (unsigned short *) array);
        break;
      case VTK_INT:
        AssignEight(vals, index, s, m, (int *) array);
        break;
      case VTK_UNSIGNED_INT:
        AssignEight(vals, index, s, m, (unsigned int *) array);
        break;
      case VTK_UNSIGNED_LONG:
        AssignEight(vals, index, s, m, (long *) array);
        break;
      case VTK_FLOAT:
        AssignEight(vals, index, s, m, (float *) array);
        break;
      case VTK_DOUBLE:
        AssignEight(vals, index, s, m, (double *) array);
        break;
      case VTK_ID_TYPE:
        AssignEight(vals, index, s, m, (vtkIdType *) array);
        break;
    }
}

// ****************************************************************************e
//  Function:  ConvertToDouble
//
//  Purpose:
//     A function that performs a cast and conversion to a double.
//
//  Programmer: Hank Childs
//  Creation:   August 26, 2008
//
//  Modifications:
//    Kathleen Biagas, Fri Jul 13 09:23:55 PDT 2012
//    Use double instead of float.
//
// ****************************************************************************

static double
ConvertToDouble(int vartype, int index, int s, int m, void *array)
{
    switch (vartype)
    {
      case VTK_CHAR:
        return (double) ((char*)array)[s*index+m];
      case VTK_UNSIGNED_CHAR:
        return (double) ((unsigned char*)array)[s*index+m];
      case VTK_SHORT:
        return (double) ((short*)array)[s*index+m];
      case VTK_UNSIGNED_SHORT:
        return (double) ((unsigned short*)array)[s*index+m];
      case VTK_INT:
        return (double) ((int*)array)[s*index+m];
      case VTK_UNSIGNED_INT:
        return (double) ((unsigned int*)array)[s*index+m];
      case VTK_UNSIGNED_LONG:
        return (double) ((unsigned long*)array)[s*index+m];
      case VTK_FLOAT:
        return (double) ((float*)array)[s*index+m];
      case VTK_DOUBLE:
        return ((double*)array)[s*index+m];
      case VTK_ID_TYPE:
        return (double) ((vtkIdType*)array)[s*index+m];
    }

    return 0.;
}


// ****************************************************************************
//  Function: FindMatch
//
//  Purpose:
//      Traverses an ordered array in logarithmic time.
//
//  Programmer: Hank Childs
//  Creation:   November 22, 2004
//
//  Modifications:
//    Kathleen Biagas, Fri Jul 13 09:23:55 PDT 2012
//    Use double instead of float.
//
// ****************************************************************************

static inline int FindMatch(const double *A, const double &a, const int &nA)
{
    if ((a < A[0]) || (a > A[nA-1]))
        return -1;

    int low = 0;
    int hi  = nA-1;
    while ((hi - low) > 1)
    {
        int guess = (hi+low)/2;
        if (A[guess] == a)
            return guess;
        if (a < A[guess])
            hi = guess;
        else
            low = guess;
    }


    return low;
}



// ****************************************************************************
//  Method: avtMassVoxelExtractor::RegisterGrid
//
//  Purpose:
//      Registers a rectilinear grid.  This is called in combination with
//      SampleAlongSegment, which assumes that data members in this class have
//      been set up.
//
//  Programmer: Hank Childs
//  Creation:   November 20, 2004
//
//  Modifications:
//
//    Hank Childs, Fri Jun  1 15:37:33 PDT 2007
//    Add support for non-scalars.
//
//    Hank Childs, Wed Aug 27 11:06:27 PDT 2008
//    Add support for non-floats.
//
//    Hank Childs, Thu Aug 28 10:52:32 PDT 2008
//    Make sure we only sample the variables that were requested.
//
//    Kathleen Biagas, Fri Jul 13 09:23:55 PDT 2012
//    Use double instead of float.
//
// ****************************************************************************

void
avtMassVoxelExtractor::RegisterGrid(vtkRectilinearGrid *rgrid,
                                    std::vector<std::string> &varorder,
                                    std::vector<int> &varsize)
{
    int  i, j, k;

    rgrid->GetDimensions(dims);
    if (X != NULL)
        delete [] X;
    if (Y != NULL)
        delete [] Y;
    if (Z != NULL)
        delete [] Z;

    // dims is the size of each of the small 3D patches e.g. 52x16x16 (or grid)
    // X, Y & Z store the "real" coordinates each point in the grid (the above grid) e.g. 0.61075, 0.19536, 0.01936 for 0,0,0
    X = new double[dims[0]];
    for (i = 0 ; i < dims[0] ; i++)
        X[i] = rgrid->GetXCoordinates()->GetTuple1(i);
    Y = new double[dims[1]];
    for (i = 0 ; i < dims[1] ; i++)
        Y[i] = rgrid->GetYCoordinates()->GetTuple1(i);
    Z = new double[dims[2]];
    for (i = 0 ; i < dims[2] ; i++)
        Z[i] = rgrid->GetZCoordinates()->GetTuple1(i);

    vtkDataArray *arr = rgrid->GetCellData()->GetArray("avtGhostZones");
    if (arr != NULL)
        ghosts = (unsigned char *) arr->GetVoidPointer(0);  // here
    else
        ghosts = NULL;

    ncell_arrays = 0;
    for (i = 0 ; i < rgrid->GetCellData()->GetNumberOfArrays() ; i++)
    {
        vtkDataArray *arr = rgrid->GetCellData()->GetArray(i);
        const char *name = arr->GetName();
        int idx = -1;
        for (j = 0 ; j < varorder.size() ; j++)
        {
            if (varorder[j] == name)
            {
                idx = 0;
                for (k = 0 ; k < j ; k++)
                    idx += varsize[k];
                break;
            }
        }
        if (idx < 0)
            continue;
        cell_index[ncell_arrays] = idx;
        cell_vartypes[ncell_arrays] = arr->GetDataType();
        cell_size[ncell_arrays] = arr->GetNumberOfComponents();
        cell_arrays[ncell_arrays++] = arr->GetVoidPointer(0);
    }

    npt_arrays = 0;
    for (i = 0 ; i < rgrid->GetPointData()->GetNumberOfArrays() ; i++)
    {
        vtkDataArray *arr = rgrid->GetPointData()->GetArray(i);
        const char *name = arr->GetName();
        int idx = -1;
        for (j = 0 ; j < varorder.size() ; j++)
        {
            if (varorder[j] == name)
            {
                idx = 0;
                for (k = 0 ; k < j ; k++)
                    idx += varsize[k];
                break;
            }
        }
        if (idx < 0)
            continue;
        pt_index[npt_arrays] = idx;
        pt_vartypes[npt_arrays] = arr->GetDataType();
        pt_size[npt_arrays] = arr->GetNumberOfComponents();
        pt_arrays[npt_arrays++] = arr->GetVoidPointer(0);
    }

    if (divisors_X != NULL)
        delete [] divisors_X;
    if (divisors_Y != NULL)
        delete [] divisors_Y;
    if (divisors_Z != NULL)
        delete [] divisors_Z;

    //
    // We end up dividing by the term A[i+1]/A[i] a whole bunch.  So store
    // out its inverse so that we can do cheap multiplication.  This gives us
    // a 5% performance boost.
    //
    divisors_X = new double[dims[0]-1];
    for (i = 0 ; i < dims[0] - 1 ; i++)
        divisors_X[i] = (X[i+1] == X[i] ? 1. : 1./(X[i+1]-X[i]));
    divisors_Y = new double[dims[1]-1];
    for (i = 0 ; i < dims[1] - 1 ; i++)
        divisors_Y[i] = (Y[i+1] == Y[i] ? 1. : 1./(Y[i+1]-Y[i]));
    divisors_Z = new double[dims[2]-1];
    for (i = 0 ; i < dims[2] - 1 ; i++)
        divisors_Z[i] = (Z[i+1] == Z[i] ? 1. : 1./(Z[i+1]-Z[i]));
}



// ****************************************************************************
//  Method: avtMassVoxelExtractor::Extract
//
//  Purpose:
//      Extracts the grid into the sample points.
//
//  Programmer: Hank Childs
//  Creation:   November 19, 2004
//
//  Modifications:
//    Jeremy Meredith, Thu Feb 15 13:11:34 EST 2007
//    Added an ability to extract voxels using the world-space version
//    even when they're really in image space.
//
//    Hank Childs, Fri Jun  1 16:40:10 PDT 2007
//    Added support for non-scalars.
//
// ****************************************************************************

void
avtMassVoxelExtractor::Extract(vtkRectilinearGrid *rgrid,
                std::vector<std::string> &varnames, std::vector<int> &varsizes)
{
    if (gridsAreInWorldSpace || pretendGridsAreInWorldSpace)
        if (rayCastingSLIVR)
            simpleExtractWorldSpaceGrid(rgrid, varnames, varsizes);
        else
          ExtractWorldSpaceGrid(rgrid, varnames, varsizes);
    else
        ExtractImageSpaceGrid(rgrid, varnames, varsizes);
}



// ****************************************************************************
//  Method: avtMassVoxelExtractor::ExtractWorldSpaceGrid
//
//  Purpose:
//      Extracts a grid that has already been put into world space.  This case
//      typically corresponds to ray casting.
//
//  Programmer: Hank Childs
//  Creation:   November 19, 2004
//
//  Modifications:
//
//    Hank Childs, Fri Dec 10 10:36:40 PST 2004
//    Use the "restricted" screen space.  This will normally be the whole
//    screen space, but if we are tiling, then it will correspond to only
//    the tile.
//
//    Hank Childs, Fri Jun  1 15:45:58 PDT 2007
//    Add support for non-scalars.
//
//    Kathleen Biagas, Fri Jul 13 09:23:55 PDT 2012
//    Use double instead of float.
//
// ****************************************************************************

void
avtMassVoxelExtractor::ExtractWorldSpaceGrid(vtkRectilinearGrid *rgrid,
                 std::vector<std::string> &varnames, std::vector<int> &varsize)
{
    //
    // Some of our sampling routines need a chance to pre-process the data.
    // Register the grid here so we can do that.
    //
    RegisterGrid(rgrid, varnames, varsize);

    //
    // Set up a list of ranges to look at.
    //
    const int max_ranges = 100; // this should be bigger than log(max(W,H))
    int width_min[max_ranges];
    int width_max[max_ranges];
    int height_min[max_ranges];
    int height_max[max_ranges];
    int curRange = 0;
    width_min[curRange] = restrictedMinWidth;
    width_max[curRange] = restrictedMaxWidth+1;
    height_min[curRange] = restrictedMinHeight;
    height_max[curRange] = restrictedMaxHeight+1;
    curRange++;

    while (curRange > 0)
    {
        //
        // Determine what range we are dealing with on this iteration.
        //
        curRange--;
        int w_min = width_min[curRange];
        int w_max = width_max[curRange];
        int h_min = height_min[curRange];
        int h_max = height_max[curRange];

        //
        // Let's find out if this range can even intersect the dataset.
        // If not, just skip it.
        //
        if (!FrustumIntersectsGrid(w_min, w_max, h_min, h_max))
        {
            continue;
        }

        int num_rays = (w_max-w_min)*(h_max-h_min);
        if (num_rays > 5)
        {
            //
            // Break the range down into smaller subchunks and recurse.
            //
            int w_mid = (w_min+w_max) / 2;
            int h_mid = (h_min+h_max) / 2;
            width_min[curRange] = w_min; width_max[curRange] = w_mid;
            height_min[curRange] = h_min; height_max[curRange] = h_mid;
            curRange++;
            width_min[curRange] = w_mid; width_max[curRange] = w_max;
            height_min[curRange] = h_min; height_max[curRange] = h_mid;
            curRange++;
            width_min[curRange] = w_min; width_max[curRange] = w_mid;
            height_min[curRange] = h_mid; height_max[curRange] = h_max;
            curRange++;
            width_min[curRange] = w_mid; width_max[curRange] = w_max;
            height_min[curRange] = h_mid; height_max[curRange] = h_max;
            curRange++;
        }
        else
        {
            //
            // We have a small amount of rays, so just evaluate them.
            //
            for (int i = w_min ; i < w_max ; i++)
                for (int j = h_min ; j < h_max ; j++)
                {
                    double origin[4];
                    double terminus[4];
                    GetSegment(i, j, origin, terminus);
                    SampleAlongSegment(origin, terminus, i, j);
                }
        }
    }
}


// ****************************************************************************
//  Method: avtMassVoxelExtractor::simpleExtractWorldSpaceGrid
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************

void
avtMassVoxelExtractor::simpleExtractWorldSpaceGrid(vtkRectilinearGrid *rgrid,
                 std::vector<std::string> &varnames, std::vector<int> &varsize)
{
    //std::cout << " avtMassVoxelExtractor::simpleExtractWorldSpaceGrid "  << std::endl;
    patchDrawn = 0;

    //
    // Some of our sampling routines need a chance to pre-process the data.
    // Register the grid here so we can do that.
    //
    RegisterGrid(rgrid, varnames, varsize);   // stores the values in a structure so that it can be used

    //
    // Determine what range we are dealing with on this iteration.
    //
    int w_min = restrictedMinWidth;
    int w_max = restrictedMaxWidth+1;
    int h_min = restrictedMinHeight;
    int h_max = restrictedMaxHeight+1;

    imgWidth = imgHeight = 0;
    //std::cout << "w_min: " << w_max << "   w_max: " << w_max << "   h_min: " << h_min << "   h_max: " << h_max << std::endl;

    //
    // Let's find out if this range can even intersect the dataset.
    // If not, just skip it.
    //
    if (!FrustumIntersectsGrid(w_min, w_max, h_min, h_max))
       return;
    
    double scRange[2], bounds[6];
    int dDims[3];
    rgrid->GetScalarRange(scRange);
    rgrid->GetBounds(bounds);
    rgrid->GetDimensions(dDims);


    //
    // AMR: set where to clip off in world space
    //
    bool minFound, maxFound;
   
    minIndex[0] = minIndex[1] = minIndex[2] = 0;
    minIndex[0] = dims[0]-1;
    minIndex[1] = dims[1]-1;
    minIndex[2] = dims[2]-1;

    //
    // Check whether fully inside or not
    if (((X[0] >= currentPartitionExtents[0]         && X[0] < currentPartitionExtents[3]) && 
         (X[dims[0]-1] >= currentPartitionExtents[0] && X[dims[0]-1] < currentPartitionExtents[3]))  &&
         
        ((Y[0] >= currentPartitionExtents[1]         && Y[0] < currentPartitionExtents[4]) && 
         (Y[dims[1]-1] >= currentPartitionExtents[1] && Y[dims[1]-1] < currentPartitionExtents[3]))  &&
         
        ((Z[0] >= currentPartitionExtents[2]         && Z[0] < currentPartitionExtents[5]) && 
         (Z[dims[2]-1] >= currentPartitionExtents[2] && Z[dims[2]-1] < currentPartitionExtents[3])) )
         fullyInside = true;

    if (isAMR == true && avtCallback::UseaMRDuplication() == true){
        // 
        // x
        minFound = maxFound = false;

        // fully in
        if (X[0] >= currentPartitionExtents[0] && X[0] < currentPartitionExtents[3]){ 
            minIndex[0] = 0;
            minFound = true;
        }

        if (X[dims[0]-1] >= currentPartitionExtents[0] && X[dims[0]-1] < currentPartitionExtents[3]){ 
            maxIndex[0] = dims[0]-1;
            maxFound = true;
        }

        // other cases
        if (minFound == false)
            minIndex[0] = 0;
            for (int i=0; i<dims[0]-1; i++){
                if (X[i] >= currentPartitionExtents[0]){
                    minFound = true;
                    minIndex[0] = i;
                    break;
                }else{
                    minIndex[0] = i;
                }
            }

        if (maxFound == false)
            maxIndex[0] = dims[0]-1;
            for (int i=dims[0]-1; i>=0; i--){
                if (X[i] < currentPartitionExtents[3]){
                    maxFound = true;
                    break;
                }else{
                    maxIndex[0] = i;
                }
            }



        //
        // y
        minFound = maxFound = false;

        // fully in on one side
        if (Y[0] >= currentPartitionExtents[1] && Y[0] < currentPartitionExtents[4]){ 
            minIndex[1] = 0;
            minFound = true;
        }

        if (Y[dims[1]-1] >= currentPartitionExtents[1] && Y[dims[1]-1] < currentPartitionExtents[4]){ 
            maxIndex[1] = dims[1]-1;
            maxFound = true;
        }

        // other cases
        if (minFound == false)
            minIndex[1] = 0;
            for (int i=0; i<dims[1]-1; i++){
                if (Y[i] >= currentPartitionExtents[1]){
                    minFound = true;
                    minIndex[1] = i;
                    break;
                }else{
                    minIndex[1] = i;
                }
            }

        if (maxFound == false)
            maxIndex[1] = dims[1]-1;
            for (int i=dims[1]-1; i>=0; i--){
                if (Y[i] < currentPartitionExtents[4]){
                    maxFound = true;
                    break;
                }else{
                    maxIndex[1] = i;
                }
            }


        //
        // z
        minFound = maxFound = false;

        // fully in on one side
        if (Z[0] >= currentPartitionExtents[2] && Z[0] < currentPartitionExtents[5]){ // fully in on one side
            minIndex[2] = 0;
            minFound = true;
        }

        if (Z[dims[2]-1] >= currentPartitionExtents[2] && Z[dims[2]-1] < currentPartitionExtents[5]){ // fully in on the other side
            maxIndex[2] = dims[2]-1;
            maxFound = true;
        }

        // other cases
        if (minFound == false)
            minIndex[2] = 0;
            for (int i=0; i<dims[2]-1; i++){
                if (Z[i] >= currentPartitionExtents[2]){
                    minFound = true;
                    minIndex[2] = i;
                    break;
                }else{
                    minIndex[2] = i;
                }
            }

        if (maxFound == false)
            maxIndex[2] = dims[2]-1;
            for (int i=dims[2]-1; i>=0; i--){
                if (Z[i] < currentPartitionExtents[5]){
                    maxFound = true;
                    break;
                }else{
                    maxIndex[2] = i;
                }
            }

        debug5 << proc << " ~ " << patch << 
             "  orig ind: 0 , 0, 0  to " << dims[0]-1 << ", " << dims[1]-1 << ", " << dims[2]-1 << 
             "  new ind: " << minIndex[0] << ", " << minIndex[1] << ", " << minIndex[2] << "  to " << maxIndex[0] << ", " << maxIndex[1] << ", " << maxIndex[2] <<
             "  patch extents: " << X[0] << ", " << Y[0] << ", " << Z[0] << "  to  " << X[dims[0]-1] << ", " << Y[dims[1]-1] << ", " << Z[dims[2]-1] <<
             "  new patch extents: " << X[minIndex[0]] << ", " << Y[minIndex[1]] << ", " << Z[minIndex[2]] << "  to  " << X[maxIndex[0]] << ", " << Y[maxIndex[1]] << ", " << Z[maxIndex[2]] <<
             "  partition extents: " << currentPartitionExtents[0] << ", " << currentPartitionExtents[1] << ", " << currentPartitionExtents[2] << "  to  " << currentPartitionExtents[3] << ", " << currentPartitionExtents[4] << ", " << currentPartitionExtents[5] << std::endl;
    }

    debug5 << proc << " ~ " << patch << 
             "  orig ind: 0 , 0, 0  to " << dims[0]-1 << ", " << dims[1]-1 << ", " << dims[2]-1 << 
             "  patch extents: " << X[0] << ", " << Y[0] << ", " << Z[0] << "  to  " << X[dims[0]-1] << ", " << Y[dims[1]-1] << ", " << Z[dims[2]-1] <<
             "  partition extents: " << currentPartitionExtents[0] << ", " << currentPartitionExtents[1] << ", " << currentPartitionExtents[2] << "  to  " << currentPartitionExtents[3] << ", " << currentPartitionExtents[4] << ", " << currentPartitionExtents[5] << std::endl;


    //
    // Determine the screen size of the patch being processed
    //
    xMin = yMin = 1000000;
    xMax = yMax = -1000000;
    
    float coordinates[8][3];
    coordinates[0][0] = X[0];           coordinates[0][1] = Y[0];           coordinates[0][2] = Z[0];
    coordinates[1][0] = X[dims[0]-1];   coordinates[1][1] = Y[0];           coordinates[1][2] = Z[0];
    coordinates[2][0] = X[dims[0]-1];   coordinates[2][1] = Y[dims[1]-1];   coordinates[2][2] = Z[0];
    coordinates[3][0] = X[0];           coordinates[3][1] = Y[dims[1]-1];   coordinates[3][2] = Z[0];

    coordinates[4][0] = X[0];           coordinates[4][1] = Y[0];           coordinates[4][2] = Z[dims[2]-1];
    coordinates[5][0] = X[dims[0]-1];   coordinates[5][1] = Y[0];           coordinates[5][2] = Z[dims[2]-1];
    coordinates[6][0] = X[dims[0]-1];   coordinates[6][1] = Y[dims[1]-1];   coordinates[6][2] = Z[dims[2]-1];
    coordinates[7][0] = X[0];           coordinates[7][1] = Y[dims[1]-1];   coordinates[7][2] = Z[dims[2]-1];

    //std::cout << proc << " ~ " << patch << "  X[0]: " << X[0] << "   X[dims[0]-1]: " << X[dims[0]-1] << "   meshMin: " << meshMin[0] << "   meshMax: " << meshMax[0] << std::endl;
    double _world[4], _view[4];
    float offset, offset_0, error_correction;
    offset = offset_0 = error_correction = 0.0f;
    _world[3] = 1.0;
    imgDepth = 0;

    for (int i=0; i<8; i++){
        _world[0] = coordinates[i][0]; 
        _world[1] = coordinates[i][1]; 
        _world[2] = coordinates[i][2];
        
        int screenPos[2]; 
        float tempImgDepth;
        world_to_screen(_world, fullImgWidth, fullImgHeight, screenPos, tempImgDepth);

        if (xMin > screenPos[0]) xMin = screenPos[0];
        if (xMax < screenPos[0]) xMax = screenPos[0];
        if (yMin > screenPos[1]) yMin = screenPos[1];
        if (yMax < screenPos[1]) yMax = screenPos[1];

        if (i == 0)
            imgDepth = tempImgDepth;
        else
            if (imgDepth < tempImgDepth)
                imgDepth = tempImgDepth;
        //std::cout << "origin: " << coordinates[i][0] << ", " << coordinates[i][1] << ", " << coordinates[i][2] << "   screenPos: " << screenPos[0] << ", " << screenPos[1] << std::endl;
    }

    xMin = xMin - error_correction;
    yMin = yMin - error_correction;
    xMax = xMax + error_correction;
    yMax = yMax + error_correction;

    xMin = xMin - offset;
    yMin = yMin - offset;
    xMax = xMax + offset;
    yMax = yMax + offset;

    if (xMin < 0) xMin = 0;
    if (yMin < 0) yMin = 0;
    if (xMax < 0) xMax = 0;
    if (yMax < 0) yMax = 0;

    if (xMin >= w_max) xMin = w_max-1;
    if (yMin >= h_max) yMin = h_max-1;
    if (xMax >= w_max) xMax = w_max-1;
    if (yMax >= h_max) yMax = h_max-1;

    int intOffset = 0;
    xMin = xMin - intOffset;    xMax = xMax + intOffset;    imgWidth =  xMax - xMin + 1;
    yMin = yMin - intOffset;    yMax = yMax + intOffset;    imgHeight = yMax - yMin + 1;

    imgArray = new float[((imgWidth)*4) * imgHeight](); // declare and initialize to 0.0 ()

    imgDims[0] = imgWidth;       imgDims[1] = imgHeight;
    imgLowerLeft[0] = xMin;      imgLowerLeft[1] = yMin;
    imgUpperRight[0] = xMax;     imgUpperRight[1] = yMax;

    if (pretendGridsAreInWorldSpace==true)
        debug5 << "pretendGridsAreInWorldSpace: true" << std::endl;
    else
        debug5 << "pretendGridsAreInWorldSpace: false" << std::endl;
    debug5 << "distance between points: " <<  screen_to_WorldDistance(xMin, yMin, xMin+1,yMin) << std::endl;
    

    task tempTask;
    if (enableThreads){
        int sizePartition = (xMax - xMin)/numThreads;
        threadArguments *threadArgs = NULL;
        threadArgs = new threadArguments[numThreads];

        for (int i=0; i<numThreads; i++){
            threadArgs[i].thisPtr = this;
            threadArgs[i].arg0 = i;
            threadArgs[i].arg1 = xMin + sizePartition*i;

            if (i != numThreads-1)
                threadArgs[i].arg2 = threadArgs[i].arg1 + sizePartition;
            else
                threadArgs[i].arg2 = xMax;

            threadArgs[i].arg3 = yMin;
            threadArgs[i].arg4 = yMax;

            if ( pthread_create(&threadHandles[i], NULL, runThread, (void *)&threadArgs[i]) )
                std::cout << "Could NOT create thread " << i << std::endl;


            tempTask.xMin = xMin + sizePartition*i;   
            tempTask.xMax = tempTask.xMin + sizePartition;
            tempTask.yMin = yMin;   tempTask.yMax = yMax;
            taskList.push_back(tempTask);
        }

        for (int i=0; i<numThreads; i++)
            pthread_join(threadHandles[i],NULL);
    
        if (threadArgs != NULL)
            delete []threadArgs;
        threadArgs = NULL;

        //closeThreads();

    }else{
        for (int i = xMin ; i < xMax ; i++)
            for (int j = yMin ; j < yMax ; j++)
            {
                double origin[4];                               // starting point where we start sampling
                double terminus[4];                             // ending point where we stop sampling
                GetSegment(i, j, origin, terminus);             // find the starting point & ending point of the ray
                SampleAlongSegment(origin, terminus, i, j);     // Go get the segments along this ray and store them in 
            }
    }

    if (patchDrawn == 0){
        if (imgArray != NULL)
            delete []imgArray;

        imgArray = NULL;
    }
}



// ****************************************************************************
//  Method: avtMassVoxelExtractor::runThread
//
//  Purpose:
//      Total hack! Since pthreads does not like C++ and the only thing it
//      accepts is a static function, we create one as an entry point into the
//      class!
//      So this function receives the arguments that the thread was called with
//      and dispatches it to the class.
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
void * avtMassVoxelExtractor::runThread(void *arg){
    threadArguments *temp = (threadArguments *)arg;

    (temp->thisPtr)->sampleImage(temp->arg0,temp->arg1,temp->arg2,temp->arg3,temp->arg4);
}


// ****************************************************************************
//  Method: avtMassVoxelExtractor::sampleImage
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
void avtMassVoxelExtractor::sampleImage(int threadId, int x_Min, int x_Max, int y_Min, int y_Max){
    debug5 << "Running with threads " << threadId << " of " << avtCallback::UseNumThreads() <<
              "Min/max: " << x_Min << ", " << x_Max << "  " << y_Min << ", " << y_Max << std::endl;

    for (int i = x_Min ; i < x_Max ; i++)
        for (int j = y_Min ; j < y_Max ; j++)
        {
            double origin[4];                               // starting point where we start sampling
            double terminus[4];                             // ending point where we stop sampling
            GetSegment(i, j, origin, terminus);             // find the starting point & ending point of the ray
            SampleAlongSegment(origin, terminus, i, j);     // Go get the segments along this ray and store them in 
        }
}

// ****************************************************************************
//  Method: avtMassVoxelExtractor::GetSegment
//
//  Purpose:
//      Gets a line segment based on a pixel location.
//
//  Arguments:
//      w           x; of the width of the render window 
//      h           y; of the height of the render window
//      origin      
//      terminus
//
//  Programmer: Hank Childs
//  Creation:   November 21, 2004
//
//  Modifications:
// 
//    Jeremy Meredith, Fri Feb  9 14:00:51 EST 2007
//    Flip back across the x axis if pretendGridsAreInWorldSpace is set.
//
//    Hank Childs, Fri Jan  9 14:07:49 PST 2009
//    Add support for jittering.
//
//    Kathleen Biagas, Fri Jul 13 09:23:55 PDT 2012
//    Use double instead of float.
//
// ****************************************************************************

void
avtMassVoxelExtractor::GetSegment(int w, int h, double *origin, double *terminus)
    const
{
    double view[4];

    //
    // The image is being reflected across a center vertical line.  This is the
    // easiest place to fix it.  I haven't tracked down where it is, but it
    // would appear to a problem with the interface from our images to VTK's.
    //
    // Note: If we're only pretending that grids are in world space -- i.e.
    // we're resampling, not raycasting an image -- then flipping across
    // the x axis like this is wrong.  Take this into account here.
    //
    view[0] = (w - width/2.)/(width/2.);
    if (pretendGridsAreInWorldSpace)
        view[0] *= -1;
    view[1] = (h - height/2.)/(height/2.);
    view[2] = cur_clip_range[0];
    view[3] = 1.;
    view_to_world_transform->MultiplyPoint(view, origin);
    
    //debug5 << "View: " << view[0] << ", " << view[1] << ", " << view[2] << ", " << view[3] << std::endl;
    //debug5 << "cur_clip_range: " << cur_clip_range[0] << ", " << cur_clip_range[1] << std::endl;
    
    if (origin[3] != 0.)
    {
        origin[0] /= origin[3];
        origin[1] /= origin[3];
        origin[2] /= origin[3];
    }

    view[0] = (w - width/2.)/(width/2.);
    if (pretendGridsAreInWorldSpace)
        view[0] *= -1;
    view[1] = (h - height/2.)/(height/2.);
    view[2] = cur_clip_range[1];
    view[3] = 1.;
    view_to_world_transform->MultiplyPoint(view, terminus);
    if (terminus[3] != 0.)
    {
        terminus[0] /= terminus[3];
        terminus[1] /= terminus[3];
        terminus[2] /= terminus[3];
    }

    if (jitter)
    {
        int reliable_random_number = (13*w*h + 14*w*w + 79*h*h + 247*w + 779*h)%513;
        double jitter = (1.0/depth) * ((reliable_random_number-256) / (256.0));
        double dir[3];
        dir[0] = (terminus[0] - origin[0])*jitter;
        dir[1] = (terminus[1] - origin[1])*jitter;
        dir[2] = (terminus[2] - origin[2])*jitter;
        origin[0] += dir[0];
        origin[1] += dir[1];
        origin[2] += dir[2];
        terminus[0] += dir[0];
        terminus[1] += dir[1];
        terminus[2] += dir[2];
    }
}


// ****************************************************************************
//  Method: avtMassVoxelExtractor::SampleAlongSegment
//
//  Purpose:
//      Samples the grid along a line segment.
//
//  Programmer: Hank Childs
//  Creation:   November 20, 2004
//
//  Modifications:
//
//    Hank Childs, Tue Jan  3 17:26:11 PST 2006
//    Fix bug that ultimately led to UMR where sampling occurred along 
//    invalid values.
//
//    Hank Childs, Wed Dec 24 11:22:43 PST 2008
//    No longer use the ProportionSpaceToZBufferSpace data member, as we now 
//    do our sampling in even intervals (wbuffer).
//
//    Kathleen Biagas, Fri Jul 13 09:23:55 PDT 2012
//    Use double instead of float.
//
// ****************************************************************************

void
avtMassVoxelExtractor::SampleAlongSegment(const double *origin, 
                                          const double *terminus, int w, int h)
{
    int first = 0;
    int last = 0;
    bool hasIntersections = FindSegmentIntersections(origin, terminus,
                                                     first, last);
    if (!hasIntersections)
        return;
    //std::cout <<  proc << " *# hasIntersection #* " <<patch  << std::endl;

    bool foundHit = false;
    int curX = -1;
    int curY = -1;
    int curZ = -1;
    bool xGoingUp = (terminus[0] > origin[0]);
    bool yGoingUp = (terminus[1] > origin[1]);
    bool zGoingUp = (terminus[2] > origin[2]);

    double x_dist = (terminus[0]-origin[0]);
    double y_dist = (terminus[1]-origin[1]);
    double z_dist = (terminus[2]-origin[2]);

    double          *prop_buffer = NULL;
    int             *ind_buffer = NULL;
    bool            *valid_sample = NULL;

    prop_buffer   = new double[3*depth];
    ind_buffer    = new int[3*depth];
    valid_sample  = new bool[depth];

    double pt[3];
    bool hasSamples = false;
 
    for (int i = first ; i < last ; i++)
    {
        int *ind = ind_buffer + 3*i;
        double *dProp = prop_buffer + 3*i;
        valid_sample[i] = false;

        double proportion = ((double)i)/((double)depth);
        pt[0] = origin[0] + proportion*x_dist;
        pt[1] = origin[1] + proportion*y_dist;
        pt[2] = origin[2] + proportion*z_dist;

        ind[0] = -1;
        ind[1] = -1;
        ind[2] = -1;

        if (!foundHit)
        {
            //
            // We haven't found any hits previously.  Exhaustively search
            // through arrays and try to find a hit.
            //
            ind[0] = FindMatch(X, pt[0], dims[0]);
            if (ind[0] >= 0)
                dProp[0] = (pt[0] - X[ind[0]]) * divisors_X[ind[0]];
            ind[1] = FindMatch(Y, pt[1], dims[1]);
            if (ind[1] >= 0)
                dProp[1] = (pt[1] - Y[ind[1]]) * divisors_Y[ind[1]];
            ind[2] = FindMatch(Z, pt[2], dims[2]);
            if (ind[2] >= 0)
                dProp[2] = (pt[2] - Z[ind[2]]) * divisors_Z[ind[2]];
        }
        else
        {
            //
            // We have found a hit before.  Try to locate the next sample 
            // based on what we already found.
            //
            if (xGoingUp)
            {
                for ( ; curX < dims[0]-1 ; curX++)
                {
                    if (pt[0] >= X[curX] && pt[0] <= X[curX+1])
                    {
                        dProp[0] = (pt[0] - X[curX]) * divisors_X[curX];
                        ind[0] = curX;
                        break;
                    }
                }
            }
            else
            {
                for ( ; curX >= 0 ; curX--)
                {
                    if (pt[0] >= X[curX] && pt[0] <= X[curX+1])
                    {
                        dProp[0] = (pt[0] - X[curX]) * divisors_X[curX];
                        ind[0] = curX;
                        break;
                    }
                }
            }
            if (yGoingUp)
            {
                for ( ; curY < dims[1]-1 ; curY++)
                {
                    if (pt[1] >= Y[curY] && pt[1] <= Y[curY+1])
                    {
                        dProp[1] = (pt[1] - Y[curY]) * divisors_Y[curY];
                        ind[1] = curY;
                        break;
                    }
                }
            }
            else
            {
                for ( ; curY >= 0 ; curY--)
                {
                    if (pt[1] >= Y[curY] && pt[1] <= Y[curY+1])
                    {
                        dProp[1] = (pt[1] - Y[curY]) * divisors_Y[curY];
                        ind[1] = curY;
                        break;
                    }
                }
            }
            if (zGoingUp)
            {
                for ( ; curZ < dims[2]-1 ; curZ++)
                {
                    if (pt[2] >= Z[curZ] && pt[2] <= Z[curZ+1])
                    {
                        dProp[2] = (pt[2] - Z[curZ]) * divisors_Z[curZ];
                        ind[2] = curZ;
                        break;
                    }
                }
            }
            else
            {
                for ( ; curZ >= 0 ; curZ--)
                {
                    if (pt[2] >= Z[curZ] && pt[2] <= Z[curZ+1])
                    {
                        dProp[2] = (pt[2] - Z[curZ]) * divisors_Z[curZ];
                        ind[2] = curZ;
                        break;
                    }
                }
            }
        }

        bool intersectedDataset = !(ind[0] < 0 || ind[1] < 0 || ind[2] < 0);
        if (!intersectedDataset)
        {
            if (!foundHit) 
            {
                // We still haven't found the start.  Keep looking.
                continue;
            }
            else
            {
                // This is the true terminus.
                last = i;
                break;
            }
        }
        else  // Did intersect data set.
        {
            if (!foundHit)
            {
                // This is the first true sample.  "The true start"
                first = i;
            }
        }

        valid_sample[i] = true;
        foundHit = true;
        hasSamples = true;

        curX = ind[0];
        curY = ind[1];
        curZ = ind[2];
    }

    if (hasSamples){
        SampleVariable(first, last, w, h, prop_buffer, ind_buffer, valid_sample);
    }

    if (prop_buffer != NULL)
        delete [] prop_buffer;
    if (ind_buffer != NULL)
        delete [] ind_buffer;
    if (valid_sample != NULL)
        delete [] valid_sample;
}


// ****************************************************************************
//  Method: avtMassVoxelExtractor::SampleVariable
//
//  Purpose:
//      Actually samples the variable into our temporaray structure.
//
//  Programmer: Hank Childs
//  Creation:   November 20, 2004
//
//  Modifications:
//
//    Hank Childs, Fri Jun  1 15:45:58 PDT 2007
//    Add support for non-scalars.
//
//    Kathleen Biagas, Fri Jul 13 09:23:55 PDT 2012
//    Use double instead of float.
//
// ****************************************************************************
void
avtMassVoxelExtractor::SampleVariable(int first, int last, int w, int h, double * prop_buffer, int * ind_buffer, bool * valid_sample)
{
    // if (proc == 5 && patch == 24)
    //     std::cout << "proc: " << proc << "   patch: " << patch << "  first: " << first << "   last: " << last << "   w: " << w << "  h:" << h << "   dims: " << dims[0] << ", " << dims[1] << ", " << dims[2] << std::endl;

    bool inrun = false;
    int  count = 0;

    avtRay *ray;
    if (rayCastingSLIVR == false)
        ray = volume->GetRay(w, h);

    int myInd[3];
    bool calc_cell_index = ((ncell_arrays > 0) || (ghosts != NULL));

    float gradientOffset = 0.5;
    bool alphaLimitReached = false;

    double dest_rgb[4] = {0.0,0.0,0.0, 0.0};     // to store the computed color
    for (int i = first ; i < last ; i++)
    {
        if (alphaLimitReached)
            break;

        // index computation
        const int *ind = ind_buffer + 3*i;       // x,y,z indices into the cell
        const double *prop = prop_buffer + 3*i;  // x: distance from left of cell |  y: distance from bottom of cell |  z: distance from the front of cell
        int index = 0;
        if (calc_cell_index)
            index = ind[2]*((dims[0]-1)*(dims[1]-1)) + ind[1]*(dims[0]-1) + ind[0]; // overall index

        //if (proc == 5 && patch == 24)
        //    std::cout << "\tTrying indices: " << ind[0] << ", " << ind[1] << ", " << ind[2] << "   prop: " << prop[0] << ", " << prop[1] << ", " << prop[2] << "  index: " << index << std::endl;

        if (rayCastingSLIVR && (isAMR && avtCallback::UseaMRDuplication() == true)){
            if (ind[0] < minIndex[0] || ind[0] > maxIndex[0])
                valid_sample[i] = false;

            if (ind[1] < minIndex[1] || ind[1] > maxIndex[1])
                valid_sample[i] = false;

            if (ind[2] < minIndex[2] || ind[2] > maxIndex[2])
                valid_sample[i] = false;
        }

        if (ghosts != NULL)
            if (ghosts[index] != 0)
                valid_sample[i] = false;
           
        if (!valid_sample[i])
            continue;

        if (!valid_sample[i] && inrun)
            if (rayCastingSLIVR == false){
                ray->SetSamples(i-count, i-1, tmpSampleList);
                inrun = false;
                count = 0;
            }

        int index_left, index_right,            index_top, index_bottom,         index_front, index_back;
        float dist_from_left, dist_from_right,  dist_from_top,dist_from_bottom,  dist_from_front, dist_from_back;
        int indices[6]; //0 left, 1 right  |  2: bottom, 3 top  | 4: front, 5 back

        int offsetLow[3], offsetHigh[3];
        offsetLow[0] = offsetLow[1] = offsetLow[2] = 1;
        offsetHigh[0] = offsetHigh[1] = offsetHigh[2] = 2;  

        float x_left = prop[0];         float x_right = 1. - x_left;    
        float y_bottom = prop[1];       float y_top = 1. - y_bottom;   
        float z_front = prop[2];        float z_back = 1. - z_front;   

        // get the index and distance from the center of the neighbouring cells
        if (ncell_arrays > 0){  // cell centered
            getIndexandDistFromCenter(prop[0], ind[0], index_left, index_right,   dist_from_left, dist_from_right);
            getIndexandDistFromCenter(prop[1],   ind[1], index_bottom,index_top,    dist_from_bottom,dist_from_top);
            getIndexandDistFromCenter(prop[2],  ind[2], index_front, index_back,   dist_from_front, dist_from_back);
        }else{                  // node centered
            index_left = ind[0];        index_right = index_left+1;
            index_bottom = ind[1];      index_top = index_bottom+1;
            index_front = ind[2];       index_back = index_front+1;

            dist_from_left = prop[0];       dist_from_right = 1.0 - prop[0];
            dist_from_bottom = prop[1];     dist_from_top = 1.0 - prop[1];
            dist_from_front = prop[2];      dist_from_back = 1.0 - prop[2];
        }
        indices[0] = index_left;        indices[1] = index_right;
        indices[2] = index_bottom;      indices[3] = index_top;
        indices[4] = index_front;       indices[5] = index_back;


        if (logicalBounds[0] == 1)
            indices[0] = indices[1] = 0;
        
        if (logicalBounds[1] == 1)
            indices[2] = indices[3] = 0;
        
        if (logicalBounds[2] == 1)
            indices[4] = indices[5] = 0;
        

        if (trilinearInterpolation){
            if (logicalBounds[0] > 1)
                if (indices[0] < 0 || indices[0]>dims[0]-2)
                    valid_sample[i] = false;

            if (logicalBounds[0] > 1)
                if (indices[1] < 0 || indices[1]>dims[0]-2)
                    valid_sample[i] = false;


            if (logicalBounds[1] > 1)
                if (indices[2] < 0 || indices[2]>dims[1]-2)
                    valid_sample[i] = false;

            if (logicalBounds[1] > 1)
                if (indices[3] < 0 || indices[3]>dims[1]-2)
                    valid_sample[i] = false;


            if (logicalBounds[2] > 1)
                if (indices[4] < 0 || indices[4]>dims[2]-2)
                    valid_sample[i] = false;

            if (logicalBounds[2] > 1)
                if (indices[5] < 0 || indices[5]>dims[2]-2)
                    valid_sample[i] = false;
        }

        if (rayCastingSLIVR && (isAMR && avtCallback::UseaMRDuplication() == true)){
            if (indices[0] < minIndex[0] || indices[0] > maxIndex[0])
                valid_sample[i] = false;

            if (indices[1] < minIndex[0] || indices[1] > maxIndex[0])
                valid_sample[i] = false;


            if (indices[2] < minIndex[1] || indices[2] > maxIndex[1])
                valid_sample[i] = false;

            if (indices[3] < minIndex[1] || indices[3] > maxIndex[1])
                valid_sample[i] = false;


            if (indices[4] < minIndex[2] || indices[4] > maxIndex[2])
                valid_sample[i] = false;

            if (indices[5] < minIndex[2] || indices[5] > maxIndex[2])
                valid_sample[i] = false;
        }


        if (rayCastingSLIVR){
            if (logicalBounds[0] > 1)
                if (!(ind[0] >= offsetLow[0] && ind[0] <= (dims[0]-1) - offsetHigh[0]))
                    valid_sample[i] = false;

            if (logicalBounds[1] > 1)
                if (!(ind[1] >= offsetLow[1] && ind[1] <= (dims[1]-1) - offsetHigh[1]))
                    valid_sample[i] = false;

            if (logicalBounds[2] > 1)
                if (!(ind[2] >= offsetLow[2] && ind[2] <= (dims[2]-1) - offsetHigh[2]))
                    valid_sample[i] = false;
        }

        if (!valid_sample[i])
            continue;

        // if (proc == 5 && patch == 24)
        //     std::cout << "\t\tIn indices: " << ind[0] << ", " << ind[1] << ", " << ind[2] << "   prop: " << prop[0] << ", " << prop[1] << ", " << prop[2] << std::endl;

        if (trilinearInterpolation || rayCastingSLIVR){
            //
            // Cell centered data
            //
            if (ncell_arrays > 0){
                int indexT[8];
                computeIndices(dims, indices, indexT);      // get the index of each of the 8 corners around the point in the data

                for (int l = 0 ; l < ncell_arrays ; l++)    // ncell_arrays: usually 1
                {
                    void  *cell_array = cell_arrays[l];
                    double values[8];

                    for (int m = 0 ; m < cell_size[l] ; m++){       // cell_size[l] usually 1
                        AssignEight(cell_vartypes[l], values, indexT, cell_size[l], m, cell_array);         
                        double scalarValue = trilinearInterpolate(values, dist_from_left, dist_from_bottom, dist_from_front);

                        // if (proc == 0)
                        //     std::cout << proc << " ~ " << patch << "  |   Values: " << values[0] << ", " << values[1] << ",  " << values[2] << ", " << values[3] << ", " << values[4] << ", " << values[5] << ", " << values[6] << ", " << values[7] << "   _ _  "
                        //                                         << "  |   dists: " << dist_from_left << ", " << dist_from_bottom << ", " << dist_from_front << "  |  scalar value: " << scalarValue
                        //                                         << "  |   pos: " << x_right << ", " << y_top << ", " << z_back << "  |  w, h: " << w << ", " << h 
                        //                                         << "  |   dims: " << dims[0] << ", " << dims[1] << ", " << dims[2]<< "  |   ind: " << ind[0] << ", " << ind[1] << ", " << ind[2]
                        //                                         << "  |   indices: " << indices[0] << ", " << indices[1] << ", " << indices[2] << ", " << indices[3] << ", " << indices[4] << ", " << indices[5]
                        //                                         << "  |   indexT: " << indexT[0] << ", " << indexT[1] << ", " << indexT[2] << ", " << indexT[3] << ", " << indexT[4] << ", " << indexT[5] << ", " << indexT[6] << ", " << indexT[7] << std::endl << std::endl << std::endl;
                        // // }
                        
                        if (rayCastingSLIVR){
                            double source_rgb[4];
                            float gradient[3];
                            int retVal = transferFn1D->QueryTF(scalarValue,source_rgb);

                            if ( ((retVal == 0)||(source_rgb[3]==0)) || (source_rgb[0]==0 && source_rgb[1]==0 && source_rgb[2]==0) ){
                                // no need to do anything more if there will be no color
                            }
                            else{
                                //
                                // Compute Lighting (if needed)
                                //
                                if (lighting == true){
                                    double vals[6];

                                    // h = offset = 1/2 the distance between grids
                                    // grad = 1/2*h * ( f(x+h,y,z)-f(x-h,y,z)    f(x,y+h,z)-f(x,y-h,z)   f(x,y,z-h)-f(x,y,z-h)  )
                                    
                                    float distFromRight, distFromLeft, distFromTop, distFromBottom, distFromFront, distFromBack;
                                    int indexLeft, indexRight, indexTop, indexBottom, indexFront, indexBack;
                                    

                                    double gradVals[8];
                                    int indexGrad[8], gradInd[3], gradIndices[6];

                                    void  *cell_array = cell_arrays[0];

                                    //
                                    // X
                                    //
                                    for (int i=0; i<6; i++)
                                        gradIndices[i] = indices[i];

                                    //
                                    // find x-h
                                    float tempX = prop[0];
                                    gradInd[0] = ind[0];
                                    if (prop[0] - gradientOffset < 0.0){
                                        tempX = (prop[0] - gradientOffset)+1.0;
                                        gradInd[0] = ind[0]-1;
                                    }
                                    
                                    getIndexandDistFromCenter(tempX, gradInd[0],  indexLeft, indexRight,  distFromLeft, distFromRight);
                                    gradIndices[0] = indexLeft;    gradIndices[1] = indexRight;
                                    computeIndices(dims, gradIndices, indexGrad);
                                    AssignEight(cell_vartypes[0], gradVals, indexGrad, 1, 0, cell_array);
                                    vals[0] = trilinearInterpolate(gradVals, distFromLeft, dist_from_bottom, dist_from_front);

                                    //
                                    // find x+h
                                    tempX = prop[0];
                                        gradInd[0] = ind[0];
                                    if (prop[0] + gradientOffset > 1.0){
                                        tempX = (prop[0] + gradientOffset)-1.0;
                                        gradInd[0] = ind[0]+1;
                                    }

                                    getIndexandDistFromCenter(tempX, gradInd[0],  indexLeft, indexRight,  distFromLeft, distFromRight);
                                    gradIndices[0] = indexLeft;    gradIndices[1] = indexRight;
                                    computeIndices(dims, gradIndices, indexGrad);
                                    AssignEight(cell_vartypes[0], gradVals, indexGrad, 1, 0, cell_array);
                                    vals[1] = trilinearInterpolate(gradVals, distFromLeft, dist_from_bottom, dist_from_front);

                                    

                                    //
                                    // Y
                                    //
                                    for (int i=0; i<6; i++)
                                        gradIndices[i] = indices[i];

                                    //
                                    // find y-h
                                    float tempY = prop[1];
                                    gradInd[1] = ind[1];
                                    if (prop[1] - gradientOffset < 0.0){
                                        tempY = (prop[1] - gradientOffset)+1.0;
                                        gradInd[1] = ind[1]-1;
                                    }
                                    
                                    getIndexandDistFromCenter(tempY, gradInd[1],  indexBottom, indexTop,  distFromBottom, distFromTop);
                                    gradIndices[2] = indexBottom ;    gradIndices[3] = indexTop;
                                    computeIndices(dims, gradIndices, indexGrad);
                                    AssignEight(cell_vartypes[0], gradVals, indexGrad, 1, 0, cell_array);
                                    vals[2] = trilinearInterpolate(gradVals, dist_from_left, distFromBottom, dist_from_front);

                                    //
                                    // find y+h
                                    tempY = prop[1];
                                    gradInd[1] = ind[1];
                                    if (prop[1] + gradientOffset > 1.0){
                                        tempY = (prop[1] + gradientOffset)-1.0;
                                        gradInd[1] = ind[1]+1;
                                    }

                                    getIndexandDistFromCenter(tempY, gradInd[1],  indexBottom, indexTop,  distFromBottom, distFromTop);
                                    gradIndices[2] = indexBottom;    gradIndices[3] = indexTop;
                                    computeIndices(dims, gradIndices, indexGrad);
                                    AssignEight(cell_vartypes[0], gradVals, indexGrad, 1, 0, cell_array);
                                    vals[3] = trilinearInterpolate(gradVals, dist_from_left, distFromBottom, dist_from_front);

                                    
                                    //
                                    // Z
                                    //
                                    for (int i=0; i<6; i++)
                                        gradIndices[i] = indices[i];

                                    //
                                    // find z-h
                                    float tempZ = prop[2];
                                    gradInd[2] = ind[2];
                                    if (prop[2] - gradientOffset < 0.0){
                                        tempZ = (prop[2] - gradientOffset)+1.0;
                                        gradInd[2] = ind[2]-1;
                                    }
                                    
                                    getIndexandDistFromCenter(tempZ, gradInd[2],  indexFront, indexBack,  distFromFront, distFromBack);
                                    gradIndices[4] = indexFront;    gradIndices[5] = indexBack;
                                    computeIndices(dims, gradIndices, indexGrad);
                                    AssignEight(cell_vartypes[0], gradVals, indexGrad, 1, 0, cell_array);
                                    vals[4] = trilinearInterpolate(gradVals, dist_from_left, dist_from_bottom, distFromFront);

                                    //
                                    // find z+h
                                    tempZ = prop[2];
                                    gradInd[2] = ind[2];
                                    if (prop[2] + gradientOffset > 1.0){
                                        tempZ = (prop[2] + gradientOffset)-1.0;
                                        gradInd[2] = ind[2]+1;
                                    }

                                    getIndexandDistFromCenter(tempZ, gradInd[2],  indexFront, indexBack,  distFromFront, distFromBack);
                                    gradIndices[4] = indexFront;    gradIndices[5] = indexBack;
                                    computeIndices(dims, gradIndices, indexGrad);
                                    AssignEight(cell_vartypes[0], gradVals, indexGrad, 1, 0, cell_array);
                                    vals[5] = trilinearInterpolate(gradVals, dist_from_left, dist_from_bottom, distFromFront);

                                    gradient[0] = (1.0/(2.0*gradientOffset)) * (vals[1] - vals[0]);
                                    gradient[1] = (1.0/(2.0*gradientOffset)) * (vals[3] - vals[2]);
                                    gradient[2] = (1.0/(2.0*gradientOffset)) * (vals[5] - vals[4]);

                                    normalize(gradient);
                                }

                                computePixelColor(source_rgb, dest_rgb, gradient);
                                if (dest_rgb[3] >= 1.0)
                                    alphaLimitReached = true;
                            }
                        }
                        else
                            tmpSampleList[count][cell_index[l]+m] = scalarValue; 
                    }
                }
            }

            //
            // Node centered data
            //
            if (npt_arrays > 0)
            {
                int indexT[8];
                computeIndicesVert(dims, indices, indexT);

                for (int l = 0 ; l < npt_arrays ; l++)
                {
                    void  *pt_array = pt_arrays[l];
                    double values[8];

                    for (int m = 0 ; m < pt_size[l] ; m++)
                    {
                        AssignEight(pt_vartypes[l], values, indexT, pt_size[l], m, pt_array);
                        double scalarValue = trilinearInterpolate(values, dist_from_left, dist_from_bottom, dist_from_front);

                        if (rayCastingSLIVR){
                            double source_rgb[4];
                            float gradient[3];
                            int retVal = transferFn1D->QueryTF(scalarValue,source_rgb);
                            if ( ((retVal == 0)||(source_rgb[3]==0)) || (source_rgb[0]==0 && source_rgb[1]==0 && source_rgb[2]==0) ){
                                // no need to do anything more if there will be no color
                            }
                            else{
                                //
                                // Compute Lighting (if needed)
                                //
                                if (lighting == true){
                                    double vals[6];

                                    // h = offset = 1/2 the distance between grids
                                    // grad = 1/2*h * ( f(x+h,y,z)-f(x-h,y,z)    f(x,y+h,z)-f(x,y-h,z)   f(x,y,z-h)-f(x,y,z-h)  )
                                    
                                    float distFromLeft, distFromBottom, distFromFront;
                                    int indexLeft, indexBottom, indexFront;
                                    double gradVals[8];
                                    int indexGrad[8], gradInd[3], gradIndices[6];

                                    //
                                    // X
                                    //
                                    for (int i=0; i<6; i++)
                                        gradIndices[i] = indices[i];

                                    //
                                    // find x-h
                                    distFromLeft = prop[0];
                                    indexLeft = ind[0];
                                    if (prop[0] - gradientOffset < 0.0){
                                        distFromLeft = (prop[0] - gradientOffset)+1.0;
                                        indexLeft = ind[0]-1;
                                    }

                                    gradIndices[0] = indexLeft;    gradIndices[1] = indexLeft+1;
                                    computeIndicesVert(dims, gradIndices, indexGrad);
                                    AssignEight(pt_vartypes[0], gradVals, indexGrad, 1, 0, pt_array);
                                    vals[0] = trilinearInterpolate(gradVals, distFromLeft, dist_from_bottom, dist_from_front);

                                    //
                                    // find x+h
                                    distFromLeft = prop[0];
                                    indexLeft = ind[0];
                                    if (prop[0] + gradientOffset > 1.0){
                                        distFromLeft = (prop[0] + gradientOffset)-1.0;
                                        indexLeft = ind[0]+1;
                                    }

                                    gradIndices[0] = indexLeft;    gradIndices[1] = indexLeft+1;
                                    computeIndicesVert(dims, gradIndices, indexGrad);
                                    AssignEight(pt_vartypes[0], gradVals, indexGrad, 1, 0, pt_array);
                                    vals[1] = trilinearInterpolate(gradVals, distFromLeft, dist_from_bottom, dist_from_front);


                                    
                                    //
                                    // Y
                                    //
                                    for (int i=0; i<6; i++)
                                        gradIndices[i] = indices[i];

                                    //
                                    // find y-h
                                    distFromBottom = prop[1];
                                    indexBottom = ind[1];
                                    if (prop[1] - gradientOffset < 0.0){
                                        distFromBottom = (prop[1] - gradientOffset)+1.0;
                                        indexBottom = ind[1]-1;
                                    }
                                    
                                    gradIndices[2] = indexBottom ;    gradIndices[3] = indexBottom+1;
                                    computeIndices(dims, gradIndices, indexGrad);
                                    AssignEight(pt_vartypes[0], gradVals, indexGrad, 1, 0, pt_array);
                                    vals[2] = trilinearInterpolate(gradVals, dist_from_left, distFromBottom, dist_from_front);

                                    //
                                    // find y+h
                                    distFromBottom = prop[1];
                                    indexBottom = ind[1];
                                    if (y_top + gradientOffset > 1.0){
                                        distFromBottom = (prop[1] + gradientOffset)-1.0;
                                        indexBottom = ind[1]+1;
                                    }

                                    gradIndices[2] = indexBottom;    gradIndices[3] = indexBottom+1;
                                    computeIndices(dims, gradIndices, indexGrad);
                                    AssignEight(pt_vartypes[0], gradVals, indexGrad, 1, 0, pt_array);
                                    vals[3] = trilinearInterpolate(gradVals, dist_from_left, distFromBottom, dist_from_front);

                                    

                                    //
                                    // Z
                                    //
                                    for (int i=0; i<6; i++)
                                        gradIndices[i] = indices[i];

                                    //
                                    // z-h
                                    distFromFront = prop[2];
                                    indexFront = ind[2];
                                    if (prop[2] + gradientOffset > 1.0){
                                        distFromFront = (prop[2]  - gradientOffset)+1.0;
                                        indexFront = ind[2]-1;
                                    }
    
                                    gradIndices[4] = indexFront;    gradIndices[5] = indexFront+1;
                                    computeIndices(dims, gradIndices, indexGrad);
                                    AssignEight(pt_vartypes[0], gradVals, indexGrad, 1, 0, pt_array);
                                    vals[4] = trilinearInterpolate(gradVals, dist_from_left, dist_from_bottom, distFromFront);

                                    //
                                    // z+h
                                    distFromFront = prop[2];
                                    indexFront = ind[2];
                                    if (prop[2] + gradientOffset > 1.0){
                                        distFromFront = (prop[2] + gradientOffset)-1.0;
                                        gradInd[2] = ind[2]+1;
                                    }

                                    gradIndices[4] = indexFront;    gradIndices[5] = indexFront+1;
                                    computeIndices(dims, gradIndices, indexGrad);
                                    AssignEight(pt_vartypes[0], gradVals, indexGrad, 1, 0, pt_array);
                                    vals[5] = trilinearInterpolate(gradVals, dist_from_left, dist_from_bottom, distFromFront);

                                    gradient[0] = (1.0/(2.0*gradientOffset)) * (vals[1] - vals[0]);
                                    gradient[1] = (1.0/(2.0*gradientOffset)) * (vals[3] - vals[2]);
                                    gradient[2] = (1.0/(2.0*gradientOffset)) * (vals[5] - vals[4]);

                                    normalize(gradient);
                                }

                                //
                                // Compute the color
                                //
                                computePixelColor(source_rgb, dest_rgb, gradient);
                                if (dest_rgb[3] >= 1.0)
                                    alphaLimitReached = true;
                            }
                        }
                        else
                            tmpSampleList[count][pt_index[l]+m] = scalarValue;
                    }
                }
            }
        }
        else{
            if (ncell_arrays > 0)
                for (int l = 0 ; l < ncell_arrays ; l++)
                    for (int m = 0 ; m < cell_size[l] ; m++)
                        tmpSampleList[count][cell_index[l]+m] = ConvertToDouble(cell_vartypes[l], index, cell_size[l], m, cell_arrays[l]);

            if (npt_arrays > 0)
            {
                int index[8];
                index[0] = (ind[2])*dims[0]*dims[1] + (ind[1])*dims[0] + (ind[0]);
                index[1] = (ind[2])*dims[0]*dims[1] + (ind[1])*dims[0] + (ind[0]+1);
                index[2] = (ind[2])*dims[0]*dims[1] + (ind[1]+1)*dims[0] + (ind[0]);
                index[3] = (ind[2])*dims[0]*dims[1] + (ind[1]+1)*dims[0] + (ind[0]+1);
                index[4] = (ind[2]+1)*dims[0]*dims[1] + (ind[1])*dims[0] + (ind[0]);
                index[5] = (ind[2]+1)*dims[0]*dims[1] + (ind[1])*dims[0] + (ind[0]+1);
                index[6] = (ind[2]+1)*dims[0]*dims[1] + (ind[1]+1)*dims[0] + (ind[0]);
                index[7] = (ind[2]+1)*dims[0]*dims[1] + (ind[1]+1)*dims[0] + (ind[0]+1);
                double x_right = prop[0];
                double x_left = 1. - prop[0];
                double y_top = prop[1];
                double y_bottom = 1. - prop[1];
                double z_back = prop[2];
                double z_front = 1. - prop[2];
                for (int l = 0 ; l < npt_arrays ; l++)
                {
                    void  *pt_array = pt_arrays[l];
                    int    s = pt_size[l];
                    for (int m = 0 ; m < s ; m++)
                    {
                        double vals[8];
                        AssignEight(pt_vartypes[l], vals, index, s, m, pt_array);
                        double val = 
                            x_left*y_bottom*z_front*vals[0] +
                            x_right*y_bottom*z_front*vals[1] +
                            x_left*y_top*z_front*vals[2] +
                            x_right*y_top*z_front*vals[3] +
                            x_left*y_bottom*z_back*vals[4] +
                            x_right*y_bottom*z_back*vals[5] +
                            x_left*y_top*z_back*vals[6] +
                            x_right*y_top*z_back*vals[7];
                        tmpSampleList[count][pt_index[l]+m] = val;
                    }
                }    
            }
        }

        inrun = true;
        count++;
    }

    //
    // Make sure we get runs at the end.
    //
    if (rayCastingSLIVR){
       imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 0] = std::min(std::max(dest_rgb[0],0.0),1.0);
       imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 1] = std::min(std::max(dest_rgb[1],0.0),1.0);
       imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 2] = std::min(std::max(dest_rgb[2],0.0),1.0);
       imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 3] = std::min(std::max(dest_rgb[3],0.0),1.0);



        // if (proc == 0){  // red
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 0] = 1;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 1] = 0;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 2] = 0;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 3] = 1;
        // }

        // if (proc == 1){ // green
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 0] = 0;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 1] = 1;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 2] = 0;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 3] = 1;
        // }

        // if (proc == 2){ // blue
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 0] = 0;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 1] = 0;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 2] = 1;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 3] = 1;
        // }

        // if (proc == 3){  //yellow
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 0] = 1;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 1] = 1;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 2] = 0;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 3] = 1;
        // }

        // if (proc == 4){ // cyan (blue-vert)
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 0] = 0;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 1] = 1;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 2] = 1;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 3] = 1;
        // }

        // if (proc == 5){ // magenta
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 0] = 1;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 1] = 0;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 2] = 1;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 3] = 1;
        // }

        // if (proc == 6){ // grey
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 0] = 0.5;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 1] = 0.5;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 2] = 0.5;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 3] = 1;
        // }

        // if (proc == 7){ // black
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 0] = 0;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 1] = 0;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 2] = 0;
        //     imgArray[(h-yMin)*(imgWidth*4) + (w-xMin)*4 + 3] = 1;
        // }
    }
    else
        if (inrun)
            ray->SetSamples(last-count, last-1, tmpSampleList);
}


// ****************************************************************************
//  Method: avtMassVoxelExtractor::world_to_screen
//
//  Purpose:
//          http://www.songho.ca/opengl/gl_transform.html
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************

void 
avtMassVoxelExtractor::world_to_screen(double _world[4], int imgWidth, int imgHeight, int screenPos[2], float &depth)
{
    double _view[4] = {0,0,0,0};

    // from world to view
    world_to_view_transform->MultiplyPoint(_world, _view);

    // perspective divide -> normalized screen space: -1,-1,-1 to 1,1,1
    if (_view[3] != 0.){
        _view[0] /= _view[3];
        _view[1] /= _view[3];
        _view[2] /= _view[3];
    }

    // viewport transform: screen coordinates
    screenPos[0] = int(_view[0]*(imgWidth/2.)  + (imgWidth /2.) + 0.5);     // 0.5 added so that proper rounding is done
    screenPos[1] = int(_view[1]*(imgHeight/2.) + (imgHeight/2.) + 0.5);

    depth = _view[2];
}


// ****************************************************************************
//  Method: avtMassVoxelExtractor::screen_to_WorldDistance
//
//  Purpose:
//
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
float 
avtMassVoxelExtractor::screen_to_WorldDistance(int x1, int y1, int x2, int y2){
    double origin1[4], origin2[4];
    double view[4];
    int w = x1;
    int h = y1;

    view[0] = (w - width/2.)/(width/2.);
    if (pretendGridsAreInWorldSpace)
        view[0] *= -1;
    view[1] = (h - height/2.)/(height/2.);
    view[2] = cur_clip_range[0];
    view[3] = 1.;
    view_to_world_transform->MultiplyPoint(view, origin1);

    w = x2;
    h = y2;
    view[0] = (w - width/2.)/(width/2.);
    if (pretendGridsAreInWorldSpace)
        view[0] *= -1;
    view[1] = (h - height/2.)/(height/2.);
    view[2] = cur_clip_range[0];
    view[3] = 1.;
    view_to_world_transform->MultiplyPoint(view, origin2);

    return fabs(origin1[0]-origin2[0]);
}



// ****************************************************************************
//  Method: avtMassVoxelExtractor::SetGridsAreInWorldSpace
//
//  Purpose:
//      Tells the MVE whether or not it is extracting a world space or image
//      space grid.
//
//  Programmer: Hank Childs
//  Creation:   November 20, 2004
//
//  Modifications:
//    Jeremy Meredith, Thu Feb 15 13:11:34 EST 2007
//    Added an ability to extract voxels using the world-space version
//    even when they're really in image space.
//
//    Hank Childs, Wed Dec 24 11:22:43 PST 2008
//    No longer calculate deprecated data member ProportionSpaceToZBufferSpace.
//
//    Hank Childs, Fri Nov 18 08:32:45 PST 2011
//    Fix ray cast of rectilinear with panning.
//
// ****************************************************************************

void
avtMassVoxelExtractor::SetGridsAreInWorldSpace(bool val, const avtViewInfo &v,
                                               double asp, const double *xform)
{
    gridsAreInWorldSpace = val;

    if (!gridsAreInWorldSpace)
    {
        if (xform)
        {
            // We're essentially doing resampling, but we cannot
            // use the faster version because there is another
            // transform to take into consideration.  In this case,
            // just revert to the world space algorithm and
            // fake the necessary parameters.
            pretendGridsAreInWorldSpace = true;
        }
        else
        {
            // We can use the faster version that doesn't depend on the
            // rest of the parameters set in this function, so return;
            return;
        }
    }

    view = v;
    aspect = asp;

    if (pretendGridsAreInWorldSpace)
    {
        view = avtViewInfo();
        view.setScale = true;
        view.parallelScale = 1;
        view.nearPlane = 1;
        view.farPlane = 2;

        aspect = 1.0;
    }

    //
    // Set up a VTK camera.  This will allow us to get the direction of
    // each ray and also the origin of each ray (this is simply the
    // position of the camera for perspective projection).
    //
    vtkCamera *cam = vtkCamera::New();
    // We have *= -1.0 throughout the code.  Here is "yet another".
    view.imagePan[0] *= -1.0;
    view.imagePan[1] *= -1.0;
    view.SetCameraFromView(cam);
    cam->GetClippingRange(cur_clip_range);
    vtkMatrix4x4 *mat = cam->GetCompositeProjectionTransformMatrix(aspect,
                                         cur_clip_range[0], cur_clip_range[1]);

    if (xform)
    {
        vtkMatrix4x4 *rectTrans = vtkMatrix4x4::New();
        rectTrans->DeepCopy(xform);
        vtkMatrix4x4::Multiply4x4(mat, rectTrans, view_to_world_transform);
        world_to_view_transform->DeepCopy(view_to_world_transform);
        view_to_world_transform->Invert();
        rectTrans->Delete();
    }
    else
    {
        // being executed for now for raycasting slivr
        vtkMatrix4x4::Invert(mat, view_to_world_transform);
        world_to_view_transform->DeepCopy(mat);
    }
    cam->Delete();
}


// ****************************************************************************
//  Method: avtMassVoxelExtractor::FrustumIntersectsGrid
//
//  Purpose:
//      Determines if a frustum intersects the grid.
//
//  Programmer: Hank Childs
//  Creation:   November 21, 2004
//
//  Modifications:
//    Kathleen Biagas, Fri Jul 13 09:23:55 PDT 2012
//    Use double instead of float.
//
// ****************************************************************************

bool
avtMassVoxelExtractor::FrustumIntersectsGrid(int w_min, int w_max, int h_min,
                                             int h_max) const
{
    //
    // Start off by getting the segments corresponding to the bottom left (bl),
    // upper left (ul), bottom right (br), and upper right (ur) rays.
    //
    double bl_start[4];
    double bl_end[4];
    GetSegment(w_min, h_min, bl_start, bl_end);

    double ul_start[4];
    double ul_end[4];
    GetSegment(w_min, h_max, ul_start, ul_end);

    double br_start[4];
    double br_end[4];
    GetSegment(w_max, h_min, br_start, br_end);

    double ur_start[4];
    double ur_end[4];
    GetSegment(w_max, h_max, ur_start, ur_end);

    //
    // Now use those segments to construct bounding planes.  If the grid is
    // not on the plus side of the bounding planes, then none of the frustum
    // will intersect the grid.
    //
    // Note: the plus side of the plane is dependent on the order that these
    // points are sent into the routine "FindPlaneNormal".  There are some
    // subtleties with putting the arguments in the right order.
    //
    double normal[3];
    FindPlaneNormal(bl_start, bl_end, ul_start, normal);
    if (!GridOnPlusSideOfPlane(bl_start, normal))
        return false;
    FindPlaneNormal(bl_start, br_start, br_end, normal);
    if (!GridOnPlusSideOfPlane(bl_start, normal))
        return false;
    FindPlaneNormal(ur_start, ul_start, ur_end, normal);
    if (!GridOnPlusSideOfPlane(ur_start, normal))
        return false;
    FindPlaneNormal(ur_start, ur_end, br_start, normal);
    if (!GridOnPlusSideOfPlane(ur_start, normal))
        return false;

    return true;
}


// ****************************************************************************
//  Method: avtMassVoxelExtractor::FindPlaneNormal
//
//  Purpose:
//      Finds the normal to a plane using three points.
//
//  Programmer: Hank Childs
//  Creation:   November 21, 2004
//
//  Modifications:
//    Kathleen Biagas, Fri Jul 13 09:23:55 PDT 2012
//    Use double instead of float.
//
// ****************************************************************************

void
avtMassVoxelExtractor::FindPlaneNormal(const double *pt1, const double *pt2,
                                       const double *pt3, double *normal)
{
    //
    // Set up vectors P1P2 and P1P3.
    //
    double v1[3];
    double v2[3];

    v1[0] = pt2[0] - pt1[0];
    v1[1] = pt2[1] - pt1[1];
    v1[2] = pt2[2] - pt1[2];
    v2[0] = pt3[0] - pt1[0];
    v2[1] = pt3[1] - pt1[1];
    v2[2] = pt3[2] - pt1[2];

    //
    // The normal is the cross-product of these two vectors.
    //
    normal[0] = v1[1]*v2[2] - v1[2]*v2[1];
    normal[1] = v1[2]*v2[0] - v1[0]*v2[2];
    normal[2] = v1[0]*v2[1] - v1[1]*v2[0];
}


// ****************************************************************************
//  Method: avtMassVoxelExtractor::GridOnPlusSideOfPlane
//
//  Purpose:
//      Determines if a grid is on the plus side of a plane.
//
//  Programmer: Hank Childs
//  Creation:   November 21, 2004
//
//  Modifications:
//    Jeremy Meredith, Thu Feb 15 13:11:34 EST 2007
//    Added an ability to extract voxels using the world-space version
//    even when they're really in image space.
//
//    Kathleen Biagas, Fri Jul 13 09:23:55 PDT 2012
//    Use double instead of float.
//
// ****************************************************************************

bool
avtMassVoxelExtractor::GridOnPlusSideOfPlane(const double *origin, 
                                             const double *normal) const
{
    double x_min = X[0];
    double x_max = X[dims[0]-1];
    double y_min = Y[0];
    double y_max = Y[dims[1]-1];
    double z_min = Z[0];
    double z_max = Z[dims[2]-1];

    for (int i = 0 ; i < 8 ; i++)
    {
        double pt[3];
        pt[0] = (i & 1 ? x_max : x_min);
        pt[1] = (i & 2 ? y_max : y_min);
        pt[2] = (i & 4 ? z_max : z_min);

        //
        // The plane is of the form Ax + By + Cz - D = 0.
        //
        // Using the origin, we can calculate D:
        // D = A*origin[0] + B*origin[1] + C*origin[2]
        //
        // We want to know if 'pt' gives:
        // A*pt[0] + B*pt[1] + C*pt[2] - D >=? 0.
        //
        // We can substitute in D to get
        // A*(pt[0]-origin[0]) + B*(pt[1]-origin[1]) + C*(pt[2-origin[2]) ?>= 0
        //
        double val  = normal[0]*(pt[0] - origin[0])
                   + normal[1]*(pt[1] - origin[1])
                   + normal[2]*(pt[2] - origin[2]);

        // Note: If we're only pretending that grids are in world space -- i.e.
        // we're resampling, not raycasting an image -- then flipping across
        // the x axis also negates the proper normal values here.
        if (pretendGridsAreInWorldSpace)
            val *= -1;

        if (val >= 0)
            return true;
    }

    return false;
}


// ****************************************************************************
//  Method: avtMassVoxelExtractor::FindSegmentIntersections
//
//  Purpose:
//      Finds the intersection points of a line segment and a rectilinear grid.
//  
//  Programmer: Hank Childs
//  Creation:   November 21, 2004
//
//  Modifications:
//
//    Hank Childs, Tue Feb  5 15:44:53 PST 2008
//    Fix bugs with the origin or the terminus of the segment being inside 
//    the volume.
//
//    Hank Childs, Wed Dec 24 11:21:57 PST 2008
//    Change the logic for perspective projections to account for w-buffering
//    (that is, even sampling in space).
//
//    Hank Childs, Wed Dec 31 09:08:50 PST 2008
//    For the case where the segment intersects the volume one time, return
//    that the entire segment should be examined, as it is probably more
//    likely that we have floating point error than we have intersected a
//    corner of the data set.
//
//    Kathleen Biagas, Fri Jul 13 09:23:55 PDT 2012
//    Use double instead of float.
//
// ****************************************************************************

bool
avtMassVoxelExtractor::FindSegmentIntersections(const double *origin, 
                                  const double *terminus, int &start, int &end)
{
    double  t, x, y, z;

    int num_hits = 0;
    double hits[6]; // Should always be 2 or 0.

    double x_min = X[0];
    double x_max = X[dims[0]-1];
    double y_min = Y[0];
    double y_max = Y[dims[1]-1];
    double z_min = Z[0];
    double z_max = Z[dims[2]-1];

    if (x_min <= origin[0] && origin[0] <= x_max &&
        y_min <= origin[1] && origin[1] <= y_max &&
        z_min <= origin[2] && origin[2] <= z_max)
    {
        hits[num_hits++] = 0.0;
    }
    if (x_min <= terminus[0] && terminus[0] <= x_max &&
        y_min <= terminus[1] && terminus[1] <= y_max &&
        z_min <= terminus[2] && terminus[2] <= z_max)
    {
        hits[num_hits++] = 1.0;
    }

    //
    // If the terminus and the origin have the same X, then we will find
    // the intersection at another face.
    //
    if (terminus[0] != origin[0])
    {
        //
        // See if we hit the X-min face.
        //
        t = (x_min - origin[0]) / (terminus[0] - origin[0]);
        y = origin[1] + t*(terminus[1] - origin[1]);
        z = origin[2] + t*(terminus[2] - origin[2]);
        if (y_min <= y && y <= y_max && z_min <= z && z <= z_max &&
            t > 0. && t < 1.)
        {
            hits[num_hits++] = t;
        }

        //
        // See if we hit the X-max face.
        //
        t = (x_max - origin[0]) / (terminus[0] - origin[0]);
        y = origin[1] + t*(terminus[1] - origin[1]);
        z = origin[2] + t*(terminus[2] - origin[2]);
        if (y_min <= y && y <= y_max && z_min <= z && z <= z_max &&
            t > 0. && t < 1.)
        {
            hits[num_hits++] = t;
        }
    }

    //
    // If the terminus and the origin have the same Y, then we will find
    // the intersection at another face.
    //
    if (terminus[1] != origin[1])
    {
        //
        // See if we hit the Y-min face.
        //
        t = (y_min - origin[1]) / (terminus[1] - origin[1]);
        x = origin[0] + t*(terminus[0] - origin[0]);
        z = origin[2] + t*(terminus[2] - origin[2]);
        if (x_min <= x && x <= x_max && z_min <= z && z <= z_max &&
            t > 0. && t < 1.)
        {
            hits[num_hits++] = t;
        }

        //
        // See if we hit the Y-max face.
        //
        t = (y_max - origin[1]) / (terminus[1] - origin[1]);
        x = origin[0] + t*(terminus[0] - origin[0]);
        z = origin[2] + t*(terminus[2] - origin[2]);
        if (x_min <= x && x <= x_max && z_min <= z && z <= z_max &&
            t > 0. && t < 1.)
        {
            hits[num_hits++] = t;
        }
    }

    //
    // If the terminus and the origin have the same Z, then we will find
    // the intersection at another face.
    //
    if (terminus[2] != origin[2])
    {
        //
        // See if we hit the Z-min face.
        //
        t = (z_min - origin[2]) / (terminus[2] - origin[2]);
        x = origin[0] + t*(terminus[0] - origin[0]);
        y = origin[1] + t*(terminus[1] - origin[1]);
        if (x_min <= x && x <= x_max && y_min <= y && y <= y_max &&
            t > 0. && t < 1.)
        {
            hits[num_hits++] = t;
        }

        //
        // See if we hit the Z-max face.
        //
        t = (z_max - origin[2]) / (terminus[2] - origin[2]);
        x = origin[0] + t*(terminus[0] - origin[0]);
        y = origin[1] + t*(terminus[1] - origin[1]);
        if (x_min <= x && x <= x_max && y_min <= y && y <= y_max &&
            t > 0. && t < 1.)
        {
            hits[num_hits++] = t;
        }
    }

    if (num_hits == 0)
        return false;

    //
    // We are expecting exactly two hits.  If we don't get that, then
    // we probably hit an edge of the dataset.  Give up on this optimization.
    //
    if (num_hits != 2)
    {
        start = 0;
        end = depth-1;
        return true;
    }

    if (hits[0] > hits[1])
    {
        double t = hits[0];
        hits[0] = hits[1];
        hits[1] = t;
    }

    if (hits[0] < 0 && hits[1] < 0)
        // Dataset on back side of camera -- no intersection.
        return false;
    if (hits[0] > 1. && hits[1] > 1.)
        // Dataset past far clipping plane -- no intersection.
        return false;

    // This is the correct calculation whether we are using orthographic or
    // perspective projection ... because with perspective we are using
    // w-buffering.
    start = (int) floor(depth*hits[0]);
    end   = (int) ceil(depth*hits[1]);
    if (start < 0)
        start = 0;
    if (end > depth)
        end = depth;

    return true;
}


// ****************************************************************************
//  Method: avtImgCommunicator::getIndexandDistFromCenter
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************

void 
avtMassVoxelExtractor::getIndexandDistFromCenter(float dist, int index,    int &index_before, int &index_after,    float &dist_before, float &dist_after){
    float center = 0.5;
    if (dist < center){
        index_before = index-1;
        index_after = index;
        dist_before = dist + center;
        dist_after = 1.0 - dist_before;
    }else{
        index_before = index;
        index_after = index+1;
        dist_before = dist - center;
        dist_after = 1.0 - dist_before;
    }
}


// ****************************************************************************
//  Method: avtImgCommunicator::computeIndices
//
//  Purpose:
//
//  Parameters:
//              indices: 0 left, 1 right  |  2: bottom, 3 top  | 4: front, 5 back
//              dims: dimensions x, y & z
//
//                                 z       y      x
//              returnIndices: 0:front, bottom, left
//                             1:front, bottom, right
//                             2:front, top,    left
//                             3:front, top,    right
//                             4:back,  bottom, left
//                             5:back,  bottom, right
//                             6:back,  top,    left
//                             7:back,  top,    right
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
void 
avtMassVoxelExtractor::computeIndices(int dims[3], int indices[6], int returnIndices[8]){
    returnIndices[0] = (indices[4])*((dims[0]-1)*(dims[1]-1)) + (indices[2])*(dims[0]-1) + (indices[0]);
    returnIndices[1] = (indices[4])*((dims[0]-1)*(dims[1]-1)) + (indices[2])*(dims[0]-1) + (indices[1]);

    returnIndices[2] = (indices[4])*((dims[0]-1)*(dims[1]-1)) + (indices[3])*(dims[0]-1) + (indices[0]);
    returnIndices[3] = (indices[4])*((dims[0]-1)*(dims[1]-1)) + (indices[3])*(dims[0]-1) + (indices[1]);

    returnIndices[4] = (indices[5])*((dims[0]-1)*(dims[1]-1)) + (indices[2])*(dims[0]-1) + (indices[0]);
    returnIndices[5] = (indices[5])*((dims[0]-1)*(dims[1]-1)) + (indices[2])*(dims[0]-1) + (indices[1]);

    returnIndices[6] = (indices[5])*((dims[0]-1)*(dims[1]-1)) + (indices[3])*(dims[0]-1) + (indices[0]);
    returnIndices[7] = (indices[5])*((dims[0]-1)*(dims[1]-1)) + (indices[3])*(dims[0]-1) + (indices[1]);
}


void 
avtMassVoxelExtractor::computeIndicesVert(int dims[3], int indices[6], int returnIndices[8]){
    returnIndices[0] = (indices[4])*((dims[0])*(dims[1])) + (indices[2])*(dims[0]) + (indices[0]);
    returnIndices[1] = (indices[4])*((dims[0])*(dims[1])) + (indices[2])*(dims[0]) + (indices[1]);

    returnIndices[2] = (indices[4])*((dims[0])*(dims[1])) + (indices[3])*(dims[0]) + (indices[0]);
    returnIndices[3] = (indices[4])*((dims[0])*(dims[1])) + (indices[3])*(dims[0]) + (indices[1]);

    returnIndices[4] = (indices[5])*((dims[0])*(dims[1])) + (indices[2])*(dims[0]) + (indices[0]);
    returnIndices[5] = (indices[5])*((dims[0])*(dims[1])) + (indices[2])*(dims[0]) + (indices[1]);

    returnIndices[6] = (indices[5])*((dims[0])*(dims[1])) + (indices[3])*(dims[0]) + (indices[0]);
    returnIndices[7] = (indices[5])*((dims[0])*(dims[1])) + (indices[3])*(dims[0]) + (indices[1]);
}


// ****************************************************************************
//  Method: avtImgCommunicator::trilinearInterpolate
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
//  Reference: http://www.scratchapixel.com/lessons/3d-advanced-lessons/interpolation/trilinear-interpolation/
//
// ****************************************************************************

double 
avtMassVoxelExtractor::trilinearInterpolate(double vals[8], float dist_from_left, float dist_from_bottom, float dist_from_front){
    float dist_from_right = 1.0 - dist_from_left;
    float dist_from_top = 1.0 - dist_from_bottom;
    float dist_from_back = 1.0 - dist_from_front;
    
    double val =    dist_from_right     * dist_from_top         * dist_from_back * vals[0] + 
                    dist_from_left      * dist_from_top         * dist_from_back * vals[1] +
                    dist_from_right     * dist_from_bottom      * dist_from_back * vals[2] +
                    dist_from_left      * dist_from_bottom      * dist_from_back * vals[3] +

                    dist_from_right     * dist_from_top         * dist_from_front * vals[4] +
                    dist_from_left      * dist_from_top         * dist_from_front * vals[5] +
                    dist_from_right     * dist_from_bottom      * dist_from_front * vals[6] +
                    dist_from_left      * dist_from_bottom      * dist_from_front * vals[7];
    return val;
}


// ****************************************************************************
//  Method: avtMassVoxelExtractor::computePixelColor
//
//  Purpose:
//      Computes color
//
//  Programmer: Pascal Grosset
//  Creation:   June 10, 2013
//
//  Modifications:
//      Need to take into account lighting
//      Need to take into accoujnt multiple light sources
//
// ****************************************************************************

void
avtMassVoxelExtractor::computePixelColor(double source_rgb[4], double dest_rgb[4], float gradient[3]){
    // Phong Shading
    if (lighting == true){
        float dir[3];           // The view "right" vector.
        double view_right[3];   // view_direction cross view_up
                                
        dir[0] = -view_direction[0];
        dir[1] = -view_direction[1];
        dir[2] = -view_direction[2];
        
        double gradientDouble[3], transformedGradient[3];
        for (int i=0; i<3; i++)
            transformedGradient[i] = gradientDouble[i] = gradient[i];

        // cos(angle) = a.b;  angle between normal and light
        float normal_dot_light = dot(gradient,dir);   // angle between light and normal;
        normal_dot_light = std::max(0.0, std::min((double)fabs(normal_dot_light),1.0) );

        // Calculate color using phong shading
        // I = (I  * ka) + [ (I_i  * kd * (L.N)) + (Ia_i * ks * (R.V)^ns) ]_for each light source i
        // I = (I  * ka) +   (I  * kd*abs(cos(angle))) + (Ia * ks*abs(cos(angle))^ns)
        for (int i=0; i<3; i++)
            source_rgb[i] = source_rgb[i] * materialProperties[0] + 
                            source_rgb[i] * materialProperties[1] * normal_dot_light;

        for (int i=0; i<3; i++)
            source_rgb[i] += materialProperties[2] * pow((double)normal_dot_light,materialProperties[3]) * source_rgb[3];   // I  * kd*abs(cos(angle))
    }

    // front to back compositing
    for (int i=0; i<4; i++){
        if (source_rgb[i] > 1.0)
            source_rgb[i] = 1.0;
        
        dest_rgb[i] = source_rgb[i] * (1.0 - dest_rgb[3]) + dest_rgb[i];
    }

    // back to front
    //    dest_rgb[i] = dest_rgb[i] * (1.0 - source_rgb[3]) + source_rgb[i];

    patchDrawn = 1;
}







// ****************************************************************************
//  Method: avtMassVoxelExtractor::getImageDimensions
//
//  Purpose:
//      transfers the metadata of the patch
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
void
avtMassVoxelExtractor::getImageDimensions(int &inUse, int dims[2], int screen_ll[2], int screen_ur[2], float &avg_z, bool &_fullyInside)
{
    inUse = patchDrawn;

    dims[0] = imgDims[0];    dims[1] = imgDims[1];

    screen_ll[0] = imgLowerLeft[0];     screen_ll[1] = imgLowerLeft[1];
    screen_ur[0] = imgUpperRight[0];    screen_ur[1] = imgUpperRight[1];

    avg_z = imgDepth;
    _fullyInside = fullyInside;
}


// ****************************************************************************
//  Method: avtMassVoxelExtractor::getComputedImage
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
avtMassVoxelExtractor::getComputedImage(float *image)
{
    memcpy(image,imgArray, imgDims[0]*4*imgDims[1]*sizeof(float));
   
    if (imgArray != NULL)
        delete []imgArray;
    imgArray = NULL;
}



// ****************************************************************************
//  Method: avtMassVoxelExtractor::FindIndex
//
//  Purpose:
//      Finds the index that corresponds to a point.
//
//  Programmer:   Hank Childs
//  Creation:     December 14, 2003
//
//  Modifications:
//    Kathleen Biagas, Fri Jul 13 07:44:46 PDT 2012
//    Templatized to handle coordinates of various types.
//
// ****************************************************************************

template <class T> inline int
FindIndex(const double &pt, const int &last_hit, const int &n,
          T *vals)
{
    for (int i = last_hit ; i < n-1 ; ++i)
    {
        if (pt >= (double)vals[i] && (pt <= (double)vals[i+1]))
            return i;
    }

    for (int i = 0 ; i < last_hit ; ++i)
    {
        if (pt >= (double)vals[i] && (pt <= (double)vals[i+1]))
            return i;
    }

    return -1;
}

inline int
FindIndex(vtkDataArray *coordArray,const double &pt, const int &last_hit, 
          const int &n)
{
    switch(coordArray->GetDataType())
    {
        vtkTemplateAliasMacro(return FindIndex(pt, last_hit, n,
            static_cast<VTK_TT *>(coordArray->GetVoidPointer(0))));
        default:    return -1;
    }
}


// ****************************************************************************
//  Function:  FindRange
//
//  Purpose:
//      A templated function to find a range.
//
//  Programmer: Kathleen Biagas 
//  Creation:   July 12, 2012
//
// ****************************************************************************

template <class T> inline void
FindRange(int ind, double c, double &min, double &max, T *coord)
{
    double range = (double)coord[ind+1] - (double)coord[ind];
    min = 1. - (c - (double)coord[ind])/range;
    max = 1. - min;
}

inline void
FindRange(vtkDataArray *coordArray, int ind, double c, double &min, double &max)
{
    switch(coordArray->GetDataType())
    {
        vtkTemplateAliasMacro(FindRange(ind, c, min, max,
            static_cast<VTK_TT *>(coordArray->GetVoidPointer(0))));
        default:
            EXCEPTION1(VisItException, "Unknown Coordinate type");
    }
}

// ****************************************************************************
//  Method: avtMassVoxelExtractor::ExtractImageSpaceGrid
//
//  Purpose:
//      Extracts a grid that has already been put into image space.  This case
//      typically corresponds to resampling.
//
//  Programmer: Hank Childs
//  Creation:   December 14, 2003
//
//  Modifications:
//
//    Hank Childs, Fri Aug 27 16:00:57 PDT 2004
//    Rename ghost data array.
//
//    Hank Childs, Fri Nov 19 14:50:58 PST 2004
//    Renamed from Extract.
//
//    Hank Childs, Mon Jul 11 14:01:28 PDT 2005
//    Fix indexing issue with ghost zones ['5712].
//
//    Hank Childs, Mon Feb  6 12:36:51 PST 2006
//    Fix another issue with ghost zones that only comes up with AMR grids
//    ['6940].
//
//    Hank Childs, Fri Jun  1 15:45:58 PDT 2007
//    Add support for non-scalars.
//
//    Hank Childs, Wed Aug 27 11:07:04 PDT 2008
//    Add support for non-floats.
//
//    Hank Childs, Thu Aug 28 10:52:32 PDT 2008
//    Make sure we only sample the variables that were requested.
//
//    Hank Childs, Sun Nov  1 14:32:46 CST 2009
//    Fix bug where ghost data could cause an extra sample to be put in the
//    avtVolume, with that sample's data be uninitialized memory.
//
//    Kathleen Biagas, Fri Jul 13 09:30:53 PDT 2012
//    Handle Coordinates as their native type. Use double internally.
//
// ****************************************************************************

void
avtMassVoxelExtractor::ExtractImageSpaceGrid(vtkRectilinearGrid *rgrid,
                std::vector<std::string> &varnames, std::vector<int> &varsizes)
{
    int  i, j, k, l, m;

    int dims[3];
    rgrid->GetDimensions(dims);
    const int nX = dims[0];
    const int nY = dims[1];
    const int nZ = dims[2];


    int last_x_hit = 0;
    int last_y_hit = 0;
    int last_z_hit = 0;

    vtkUnsignedCharArray *ghosts = (vtkUnsignedCharArray *)rgrid->GetCellData()->GetArray("avtGhostZones");  //here
    std::vector<void *>  cell_arrays;
    std::vector<int>     cell_vartypes;
    std::vector<int>     cell_size;
    std::vector<int>     cell_index;
    for (i = 0 ; i < rgrid->GetCellData()->GetNumberOfArrays() ; i++)
    {
        vtkDataArray *arr = rgrid->GetCellData()->GetArray(i);
        const char *name = arr->GetName();
        int idx = -1;
        for (j = 0 ; j < varnames.size() ; j++)
        {
            if (varnames[j] == name)
            {
                idx = 0;
                for (k = 0 ; k < j ; k++)
                    idx += varsizes[k];
                break;
            }
        }
        if (idx < 0)
            continue;
        cell_index.push_back(idx);
        cell_size.push_back(arr->GetNumberOfComponents());
        cell_vartypes.push_back(arr->GetDataType());
        cell_arrays.push_back(arr->GetVoidPointer(0));
    }

    std::vector<void *>  pt_arrays;
    std::vector<int>     pt_vartypes;
    std::vector<int>     pt_size;
    std::vector<int>     pt_index;
    for (i = 0 ; i < rgrid->GetPointData()->GetNumberOfArrays() ; i++)
    {
        vtkDataArray *arr = rgrid->GetPointData()->GetArray(i);
        const char *name = arr->GetName();
        int idx = -1;
        for (j = 0 ; j < varnames.size() ; j++)
        {
            if (varnames[j] == name)
            {
                idx = 0;
                for (k = 0 ; k < j ; k++)
                    idx += varsizes[k];
                break;
            }
        }
        if (idx < 0)
            continue;
        pt_index.push_back(idx);
        pt_size.push_back(arr->GetNumberOfComponents());
        pt_vartypes.push_back(arr->GetDataType());
        pt_arrays.push_back(arr->GetVoidPointer(0));
    }


    vtkDataArray *xarray = rgrid->GetXCoordinates();
    vtkDataArray *yarray = rgrid->GetYCoordinates();
    vtkDataArray *zarray = rgrid->GetZCoordinates();

    int startX = SnapXLeft(xarray->GetTuple1(0));
    int stopX  = SnapXRight(xarray->GetTuple1(nX-1));
    int startY = SnapYBottom(yarray->GetTuple1(0));
    int stopY  = SnapYTop(yarray->GetTuple1(nY-1));
    int startZ = SnapZFront(zarray->GetTuple1(0));
    int stopZ  = SnapZBack(zarray->GetTuple1(nZ-1));

    for (j = startY ; j <= stopY ; j++)
    {
        double yc = YFromIndex(j);
        int yind = FindIndex(yarray, yc, last_y_hit, nY);
        if (yind == -1)
            continue;
        last_y_hit = yind;

        double y_bottom  = 0.;
        double y_top = 1.;
        if (pt_arrays.size() > 0) 
        {
            FindRange(yarray, yind, yc, y_bottom, y_top);
        }
        for (i = startX ; i <= stopX ; i++)
        {
            double xc = XFromIndex(i);
            int xind = FindIndex(xarray, xc, last_x_hit, nX);
            if (xind == -1)
                continue;
            last_x_hit = xind;

            double x_left  = 0.;
            double x_right = 1.;
            if (pt_arrays.size() > 0)
            {
                FindRange(xarray, xind, xc, x_left, x_right);
            }

            last_z_hit = 0;
            int count = 0;
            int firstZ = -1;
            int lastZ  = stopZ;
            for (k = startZ ; k <= stopZ ; k++)
            {
                double zc = ZFromIndex(k);
                int zind = FindIndex(zarray, zc, last_z_hit, nZ);
                if (zind == -1)
                {
                    if (firstZ == -1)
                        continue;
                    else
                    {
                        lastZ = k-1;
                        break;
                    }
                }
                if ((count == 0) && (firstZ == -1))
                    firstZ = k;
                last_z_hit = zind;

                //
                // Don't sample from ghost zones.
                //
                if (ghosts != NULL)
                {
                    int index = zind*((nX-1)*(nY-1)) + yind*(nX-1) + xind;
                    if (ghosts->GetValue(index) != 0)
                    {
                        if (count > 0)
                        {
                            avtRay *ray = volume->GetRay(i, j);
                            ray->SetSamples(firstZ, k-1, tmpSampleList);
                        }
                        firstZ = -1;
                        count = 0;
                        continue;
                    }
                }

                double z_front  = 0.;
                double z_back = 1.;
                if (pt_arrays.size() > 0)
                {
                    FindRange(zarray, zind, zc, z_front, z_back);
                }

                int var_index = 0;
                for (l = 0 ; l < cell_arrays.size() ; l++)
                {
                    int index = zind*((nX-1)*(nY-1)) + yind*(nX-1) + xind;
                    for (m = 0 ; m < cell_size[l] ; m++)
                        tmpSampleList[count][cell_index[l]+m] = 
                                  ConvertToDouble(cell_vartypes[l],index,
                                              cell_size[l], m, cell_arrays[l]);
                }
                if (pt_arrays.size() > 0)
                {
                    int index[8];
                    index[0] = (zind)*nX*nY + (yind)*nX + (xind);
                    index[1] = (zind)*nX*nY + (yind)*nX + (xind+1);
                    index[2] = (zind)*nX*nY + (yind+1)*nX + (xind);
                    index[3] = (zind)*nX*nY + (yind+1)*nX + (xind+1);
                    index[4] = (zind+1)*nX*nY + (yind)*nX + (xind);
                    index[5] = (zind+1)*nX*nY + (yind)*nX + (xind+1);
                    index[6] = (zind+1)*nX*nY + (yind+1)*nX + (xind);
                    index[7] = (zind+1)*nX*nY + (yind+1)*nX + (xind+1);
                    for (l = 0 ; l < pt_arrays.size() ; l++)
                    {
                        void  *pt_array = pt_arrays[l];
                        int    s        = pt_size[l];
                        for (m = 0 ; m < s ; m++)
                        {
                            double vals[8];
                            AssignEight(pt_vartypes[l], vals, index, 
                                        s, m, pt_array);
                            double val = 
                                  x_left*y_bottom*z_front*vals[0] +
                                  x_right*y_bottom*z_front*vals[1] +
                                  x_left*y_top*z_front*vals[2] +
                                  x_right*y_top*z_front*vals[3] +
                                  x_left*y_bottom*z_back*vals[4] +
                                  x_right*y_bottom*z_back*vals[5] +
                                  x_left*y_top*z_back*vals[6] +
                                  x_right*y_top*z_back*vals[7];
                            tmpSampleList[count][pt_index[l]+m] = val;
                        }
                    }    
                }
                count++;
            }

            if (count > 0)
            {
                avtRay *ray = volume->GetRay(i, j);
                ray->SetSamples(firstZ, lastZ, tmpSampleList);
            }
        }
    }
}




// ****************************************************************************
//  Method: avtMassVoxelExtractor::initThreads
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
void avtMassVoxelExtractor::initThreads(){
    if (numThreads > 0)
        enableThreads = true;

    threadHandles = new pthread_t[numThreads];
    threadArgument = new threadArgument[numThreads];

    pthread_mutex_init(&mutexPatchAvailable, NULL);
    pthread_cond_init(&condPatchAvailable, NULL);

    for (int i=0; i<numThreads; i++){
        threadArgument[i].pThis = this;
        threadArgument[i].id = i;

        if ( pthread_create(&threadHandles[i], NULL, setupThread, (void *)&threadArgument[i]) )
            std::cout << "Could NOT create thread " << i << " !"<<std::endl;
    }
}



// ****************************************************************************
//  Method: avtMassVoxelExtractor::closeThreads
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
void * avtMassVoxelExtractor::setupThread(void *arg){
    threadArguments *temp = (threadArguments *)arg;

    (temp->thisPtr)->doWork(temp->id);
}



// ****************************************************************************
//  Method: avtMassVoxelExtractor::closeThreads
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
void avtMassVoxelExtractor::doWork(int id){

    if (allPatchesProcessed == true)
        return;

    for (int i = x_Min ; i < x_Max ; i++)
        for (int j = y_Min ; j < y_Max ; j++)
        {
            double origin[4];                               // starting point where we start sampling
            double terminus[4];                             // ending point where we stop sampling
            GetSegment(i, j, origin, terminus);             // find the starting point & ending point of the ray
            SampleAlongSegment(origin, terminus, i, j);     // Go get the segments along this ray and store them in 
        }
}


// ****************************************************************************
//  Method: avtMassVoxelExtractor::closeThreads
//
//  Purpose:
//
//  Programmer: 
//  Creation:   
//
//  Modifications:
//
// ****************************************************************************
void avtMassVoxelExtractor::closeThreads(){
    enableThreads = false;
    numThreads = 0;

    for (int i=0; i<numThreads; i++)
        pthread_join(threadHandles[i],NULL);
    
    if (threadArgument != NULL)
        delete []threadArgument;
    threadArgument = NULL;

    if (threadHandles != NULL)
        delete []threadHandles;
    threadHandles = NULL;

    pthread_cond_destroy(&mutexPatchAvailable);
    pthread_mutex_destroy(&mutexPatchAvailable);

    pthread_exit(NULL);
}
