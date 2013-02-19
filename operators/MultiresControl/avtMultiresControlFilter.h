/*****************************************************************************
*
* Copyright (c) 2010, University of New Hampshire Computer Science Department
* All rights reserved.
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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  File:        avtMultiresControlFilter.h                                  //
//  Programmer:  Andrew Foulks <rafoulks@cs.unh.edu>                         //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef AVT_MultiresControl_FILTER_H
#define AVT_MultiresControl_FILTER_H

#include <avtPluginDataTreeIterator.h>
#include <MultiresControlAttributes.h>

class vtkDataSet;


//  The MultiresControlFilter does not actually filter or
//  modify the data.  Instead, it is used as a means
//  to control the current resolution.  It does this 
//  through a gui widget that comes up as part of
//  OperatorAttributes menu in VisIt.  The widget allows
//  the user to select the current resolution.  The
//  selection is passed to the MultiresData database
//  plugin using a VisIt 'contract'.
//
//  @author Andrew Foulks

class avtMultiresControlFilter : public avtPluginDataTreeIterator
{
  public:
                          avtMultiresControlFilter();
    virtual               ~avtMultiresControlFilter();
    static avtFilter*     Create();
    virtual const char*   GetType(void)  {return "avtMultiresControlFilter";};
    virtual const char*   GetDescription(void) {return "MultiresControl";};
    virtual void          SetAtts(const AttributeGroup*);
    virtual bool          Equivalent(const AttributeGroup*);
    virtual avtContract_p ModifyContract(avtContract_p contract);
    //virtual void          ExamineContract(avtContract_p contract);

  protected:
    MultiresControlAttributes   atts;
    virtual vtkDataSet   *ExecuteData(vtkDataSet *, int, std::string);
};


#endif