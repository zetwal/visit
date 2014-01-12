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

#include <avtDatasetToImgFilter.h>

#include <avtDatasetToSamplePointsFilter.h>

#include <avtDatasetExaminer.h>
#include <avtParallel.h>

#include <DebugStream.h>

#include <cstring>

#include <string>
#include <vector>

avtDatasetToImgFilter::avtDatasetToImgFilter()
{
}

avtDatasetToImgFilter::~avtDatasetToImgFilter()
{
}

void
avtDatasetToImgFilter::PreExecute(void)
{
    avtDatasetToDataObjectFilter::PreExecute();

    avtDataset_p ds = GetTypedInput();
    VarList vl;
    vl.nvars = -1;
    avtDatasetExaminer::GetVariableList(ds, vl);

    avtImage_p sp = GetTypedOutput();
    if (vl.nvars <= 0)
    {
        debug1 << "!!! Converting a dataset that has no variables to sample "
               << "points" << endl;
        vl.nvars = 0;
    }

    bool leaveAsIs = false;
    // if (vl.nvars == 0 && sp->GetNumberOfVariables() > 0)
    // {
    //     //
    //     // Someone came in and set the output so that it had more variables --
    //     // this is common practice if we are executing in parallel and we
    //     // have more processors than domains.
    //     //
    //     debug1 << "!!! The sample points already believed that it had "
    //            << "variables -- leaving as is." << endl;
    //     leaveAsIs = true;
    // }

    if (!leaveAsIs)
    {
        std::vector<std::string> varnames;
        std::vector<int>    varsize;
        int realNVars = 0;
        for (int i = 0 ; i < vl.nvars ; i++)
        {
            const char *vname = vl.varnames[i].c_str();
            if (strstr(vname, "vtk") != NULL)
                continue;
            if (strstr(vname, "avt") != NULL)
                continue;
            varnames.push_back(vl.varnames[i]);
            varsize.push_back(vl.varsizes[i]);
            realNVars++;
        }

        // Some contortions here to use existing calls in avtParallel.
        int nvars = UnifyMaximumValue((int)varnames.size());
        GetListToRootProc(varnames, nvars);
        BroadcastStringVector(varnames, PAR_Rank());

        while (varsize.size() < nvars)
            varsize.push_back(0);
        std::vector<int> varsize2(nvars);
        UnifyMaximumValue(varsize, varsize2);

        //sp->SetNumberOfVariables(varsize2, varnames);
    }
}