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
//                           avtRayExtractor.h                       //
// ************************************************************************* //

#ifndef AVT_RAY_EXTRACTOR_H
#define AVT_RAY_EXTRACTOR_H

#include <filters_exports.h>

#include <avtDatasetToSamplePointsFilter.h>
#include <avtDatasetToImgFilter.h>
#include <avtVolume.h>

#include <avtViewInfo.h>

#include <avtOpacityMap.h>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <utility>

#include <avtImgCommunicator.h>
#include <avtImageSource.h>

#include <deque>

class  vtkDataArray;
class  vtkDataSet;
class  vtkHexahedron;
class  vtkQuadraticHexahedron;
class  vtkPixel;
class  vtkPyramid;
class  vtkQuad;
class  vtkTetra;
class  vtkTriangle;
class  vtkVoxel;
class  vtkWedge;

class  avtHexahedronExtractor;
class  avtHexahedron20Extractor;
class  avtMassVoxelExtractor;
class  avtPointExtractor;
class  avtPyramidExtractor;
class  avtSamplePointArbitrator;
class  avtTetrahedronExtractor;
class  avtWedgeExtractor;

class  avtRayFunction;


// ****************************************************************************
//  Class: avtRayExtractor
//
//  Purpose:
//      This is a component that will take an avtDataset as an input and find
//      all of the sample points from that dataset.
//
//  Programmer: Hank Childs
//  Creation:   December 5, 2000
//
//  Modifications:
//
//    Hank Childs, Sat Jan 27 15:09:34 PST 2001
//    Added support for sending cells when doing parallel volume rendering.
//
//    Kathleen Bonnell, Sat Apr 21, 13:09:27 PDT 2001 
//    Added recursive Execute method to walk down input data tree. 
//
//    Hank Childs, Tue Nov 13 15:51:15 PST 2001
//    Remove boolean argument to Extract<Cell> calls since it is no longer
//    necessary when all of the variables are being extracted.
//
//    Hank Childs, Sun Dec 14 11:07:56 PST 2003
//    Added mass voxel extractor.
//
//    Hank Childs, Fri Nov 19 13:41:56 PST 2004
//    Added view conversion option.
//
//    Hank Childs, Sat Jan 29 13:32:54 PST 2005
//    Added 2D extractors.
//
//    Hank Childs, Sun Dec  4 19:12:42 PST 2005
//    Added support for kernel-based sampling.
//
//    Hank Childs, Sun Jan  1 10:56:19 PST 2006
//    Added RasterBasedSample and KernelBasedSample.
//
//    Hank Childs, Tue Feb 28 08:25:33 PST 2006
//    Added PreExecute.
//
//    Jeremy Meredith, Thu Feb 15 11:44:28 EST 2007
//    Added support for rectilinear grids with an inherent transform.
//
//    Hank Childs, Fri Jun  1 11:47:56 PDT 2007
//    Add method GetLoadingInfoForArrays.
//
//    Hank Childs, Thu Sep 13 14:02:40 PDT 2007
//    Added support for hex-20s.
//
//    Hank Childs, Tue Jan 15 14:17:15 PST 2008
//    Have this class set up custom sample point arbitrators, since it has
//    the most knowledge.
//
//    Hank Childs, Fri Jan  9 14:09:57 PST 2009
//    Add support for jittering.
//
// ****************************************************************************

class AVTFILTERS_API avtRayExtractor 
    : public avtDatasetToImgFilter//avtDatasetToSamplePointsFilter
{
  public:
                              avtRayExtractor(int, int, int);
    virtual                  ~avtRayExtractor();

    virtual const char       *GetType(void)
                                         { return "avtRayExtractor"; };
    virtual const char       *GetDescription(void)
                                         { return "Extracting sample points";};

    void                      SetRectilinearGridsAreInWorldSpace(bool, 
                                                   const avtViewInfo &,double);
    void                      RestrictToTile(int, int, int, int);
    void                      StopTiling(void) { shouldDoTiling = false; };

    void                      Set3DMode(bool m) { modeIs3D = m; };
    void                      SetKernelBasedSampling(bool);
    void                      SetJittering(bool);

    void                      SetUpArbitrator(std::string &name, bool min);

    void                      SetTrilinear(bool t) {trilinearInterpolation = t;  };
    void                      SetRayCastingSLIVR(bool s) {rayCastingSLIVR = s;  };
    void                      SetLighting(bool l) {lighting = l; };
    void                      SetLightPosition(double _lightPos[4]) { for (int i=0;i<4;i++) lightPosition[i]=_lightPos[i]; }
    void                      SetLightDirection(double _lightDir[3]) { for (int i=0;i<3;i++) lightDirection[i]=_lightDir[i]; }
    void                      SetMatProperties(double _matProp[4]) { for (int i=0;i<4;i++) materialProperties[i]=_matProp[i]; }
    void                      SetTransferFn(avtOpacityMap *_transferFn1D) {transferFn1D = _transferFn1D; };
    void                      SetModelViewMatrix(double _modelViewMatrix[16]) { for (int i=0;i<16;i++) modelViewMatrix[i]=_modelViewMatrix[i]; }
    void                      SetViewDirection(double *vd){ for (int i=0; i<3; i++) view_direction[i] = vd[i]; }
    void                      SetViewUp(double *vu){ for (int i=0; i<3; i++) view_up[i] = vu[i]; }
    void                      SetMeshDims(double _meshMin[3], double _meshMax[3]){ for (int i=0; i<3; i++) { meshMin[i] = _meshMin[i]; meshMax[i] = _meshMax[i];}}
    void                      SetLogicalBounds(int _l, int _w, int _h){ logicalBounds[0] = _l; logicalBounds[1] = _w; logicalBounds[2] = _h; }
    void                      SetPartitionExtents(double _currentPartitionExtents[6]){for (int i=0; i<6; i++) currentPartitionExtents[i]=_currentPartitionExtents[i]; }
    void                      SetVarName(std::string _varName){ varName = _varName;}
    void                      SetParallelOn(bool _parallelOn){parallelOn = _parallelOn;}
    void                      SetScreen(int _screen[2]){ for (int i=0;i<2;i++) screen[i]=_screen[i]; }
    void                      SetBackground(unsigned char _background[3]){ for (int i=0;i<3;i++) background[i]=_background[i]; }
    void                      SetParentChild(std::vector<int> v){parentChild = v;}
    void                      SetNumChildren(std::vector<int> v){numChildren = v;}
    void                      SetNumInEachLevel(std::vector<int> v){numInEachLevel = v;}
    void                      SetPatchLevel(std::vector<int> v){patchLevel = v;}

    // Getting image information
    int                       getTotalAssignedPatches() { return totalAssignedPatches; }              // gets the max number of patches it could have
    int                       getImgPatchSize(){ return patchCount;};                                 // gets the number of patches
    imgMetaData               getImgMetaPatch(int patchId){ return imageMetaPatchVector.at(patchId);} // gets the metadata
    void                      getnDelImgData(int patchId, imgData &tempImgData);                      // gets the image & erase its existence
    void                      getImgData(int patchId, imgData &tempImgData);
    void                      delImgData(int patchId);
    
    void                      delImgPatches();   
    avtImage_p                ExecuteRayTracer();
    avtImage_p                ExecuteRayTracerLB();

    int                       chopPartitionRT(partitionExtents parent, partitionExtents & childOne, partitionExtents & childTwo, int axisOrder[3]);
    void                      getPartitionExtents(int numDivisions, int logicalBounds[3], double minSpatialExtents[3], double maxSpatialExtents[3], double extents[6]);
    bool                      patchOverlap(float patchMinX, float patchMaxX, float patchMinY, float patchMaxY, float patchMinZ, float patchMaxZ,
    float partitionMinX, float partitionMaxX, float partitionMinY, float partitionMaxY, float partitionMinZ, float partitionMaxZ);

    std::vector<int>         getAllChildrenOfPatch(int patchId);
    std::vector<int>         getDirectChildrenOfPatch(int patchId);
    // Check if not outside; if it is not outside it has to be somewhere inside 

  protected:
    avtImgCommunicator        imgComm;
    bool                      parallelOn;
    int                       screen[2];
    unsigned char             background[3];

    int                       width, height, depth;
    int                       currentNode, totalNodes;

    bool                      shouldDoTiling;
    int                       width_min, width_max;
    int                       height_min, height_max;
    bool                      modeIs3D;
    bool                      kernelBasedSampling;
    double                    point_radius;

    bool                      shouldSetUpArbitrator;
    std::string               arbitratorVarName;
    bool                      arbitratorPrefersMinimum;
    avtSamplePointArbitrator *arbitrator;

    avtMassVoxelExtractor    *massVoxelExtractor;

    bool                      jitter;

    bool                      rectilinearGridsAreInWorldSpace;
    avtViewInfo               viewInfo;
    double                    aspect;

    int                       patchCount;
    int                       totalAssignedPatches;

    std::vector<imgMetaData>    imageMetaPatchVector;
    std::multimap<int, imgData> imgDataHashMap;
    typedef std::multimap<int, imgData>::iterator iter_t;


    // triliniear / raycastin SLIVR
    bool                      trilinearInterpolation;
    bool                      rayCastingSLIVR;

    // lighting & material
    double                    view_direction[3];
    double                    view_up[3];
    double                    modelViewMatrix[16];

    bool                      lighting;
    double                    lightPosition[4];
    double                    lightDirection[3];
    double                    materialProperties[4];
    avtOpacityMap             *transferFn1D;
    virtual void              Execute(void);
    void                      PreExecute(void);
    virtual void              PostExecute(void);
    virtual void              ExecuteTree(avtDataTree_p);
    void                      SetUpExtractors(void);
    imgMetaData               initMetaPatch(int id);    // initialize a patch

    double                    meshMin[3];
    double                    meshMax[3];
    int                       logicalBounds[3];
    double                    currentPartitionExtents[6];  // minX, minY, minZ,    maxX, maxY,maxZ
    int                       screenPartitionExtents[2];
    float                     screenPartitionDepth;
    bool                      partitionExtentsComputationDone;

    std::vector<int>          parentChild;      // parent child relationship
    std::vector<int>          numChildren;      // number of children for each patch
    std::vector<int>          numInEachLevel;   // number of patches for each level
    std::vector<int>          patchLevel;       // level of each patch

    std::string               varName;
    int                       amrLevels;

    typedef struct 
    {
      std::vector<int>                  cellDataIndex;
      std::vector<int>                  pointDataIndex;
      std::vector<int>                  cellDataSize;
      std::vector<int>                  pointDataSize;
      std::vector<vtkDataArray *>       cellArrays;
      std::vector<vtkDataArray *>       pointArrays;
      int                               nVars;
    } LoadingInfo;

    void                      RasterBasedSample(vtkDataSet *, int num = 0);

};


#endif


