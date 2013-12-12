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