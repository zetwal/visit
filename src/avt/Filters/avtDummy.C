/*****************************************************************************
*
* Copyright (c) 2000 - 2012, Lawrence Livermore National Security, LLC
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
//                           avtRayCompositer.C                              //
// ************************************************************************* //

#include <avtDummy.h>

#include <vtkImageData.h>

#include <avtImage.h>
#include <avtRayFunction.h>
#include <avtSamplePoints.h>
#include <avtVolume.h>

#include <ImproperUseException.h>




// ****************************************************************************
//  Method: avtRayCompositer constructor
//
//  Arguments:
//      rf       The ray function the compositer should use to find the
//               intensity and variable values.
//     
//  Programmer:  Hank Childs
//  Creation:    December 4, 2000
//
//  Modifications:
//     Brad Whitlock, Wed Dec 5 10:43:13 PDT 2001
//     I added code to initialize the background mode, gradient background.
//
// ****************************************************************************

avtDummy::avtDummy()
{

}



// ****************************************************************************
//  Method: avtRayCompositer::Execute
//
//  Purpose:
//      Creates the output image by compositing the rays into pixels (the
//      compositing is actually done by avtVolume).
//
//  Programmer: Hank Childs
//  Creation:   December 5, 2000
//
//  Modifications:
//
//    Hank Childs, Thu Jan 25 23:42:20 PST 2001
//    Removed section to create vtkImageData so code could be consolidated in
//    avtImageRepresentation.
//
//    Hank Childs, Mon Jan 29 20:43:38 PST 2001
//    Use the restricted portion of the screen, rather than the whole thing.
//
//    Hank Childs, Sat Feb  3 20:17:20 PST 2001
//    Pulled out pixelizers.
//
//    Brad Whitlock, Wed Dec 5 10:44:31 PDT 2001
//    Modified to support gradient backgrounds.
//
//    Hank Childs, Thu Jan  3 09:50:29 PST 2002
//    Account for case where our partition contains nothing.
//
//    Hank Childs, Wed Jan 25 12:23:59 PST 2006
//    Add error checking.
//
//    Tom Fogal, Wed Jul  2 16:32:30 EDT 2008
//    Don't assume our BG image is 3 components.
//
// ****************************************************************************

void
avtDummy::Execute(void)
{

}

