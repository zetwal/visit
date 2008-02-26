// ***************************************************************************
//
// Copyright (c) 2000 - 2008, Lawrence Livermore National Security, LLC
// Produced at the Lawrence Livermore National Laboratory
// LLNL-CODE-400142
// All rights reserved.
//
// This file is  part of VisIt. For  details, see https://visit.llnl.gov/.  The
// full copyright notice is contained in the file COPYRIGHT located at the root
// of the VisIt distribution or at http://www.llnl.gov/visit/copyright.html.
//
// Redistribution  and  use  in  source  and  binary  forms,  with  or  without
// modification, are permitted provided that the following conditions are met:
//
//  - Redistributions of  source code must  retain the above  copyright notice,
//    this list of conditions and the disclaimer below.
//  - Redistributions in binary form must reproduce the above copyright notice,
//    this  list of  conditions  and  the  disclaimer (as noted below)  in  the
//    documentation and/or other materials provided with the distribution.
//  - Neither the name of  the LLNS/LLNL nor the names of  its contributors may
//    be used to endorse or promote products derived from this software without
//    specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT  HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR  IMPLIED WARRANTIES, INCLUDING,  BUT NOT  LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND  FITNESS FOR A PARTICULAR  PURPOSE
// ARE  DISCLAIMED. IN  NO EVENT  SHALL LAWRENCE  LIVERMORE NATIONAL  SECURITY,
// LLC, THE  U.S.  DEPARTMENT OF  ENERGY  OR  CONTRIBUTORS BE  LIABLE  FOR  ANY
// DIRECT,  INDIRECT,   INCIDENTAL,   SPECIAL,   EXEMPLARY,  OR   CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT  LIMITED TO, PROCUREMENT OF  SUBSTITUTE GOODS OR
// SERVICES; LOSS OF  USE, DATA, OR PROFITS; OR  BUSINESS INTERRUPTION) HOWEVER
// CAUSED  AND  ON  ANY  THEORY  OF  LIABILITY,  WHETHER  IN  CONTRACT,  STRICT
// LIABILITY, OR TORT  (INCLUDING NEGLIGENCE OR OTHERWISE)  ARISING IN ANY  WAY
// OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ***************************************************************************

package llnl.visit;

import java.util.Vector;

// ****************************************************************************
// Class: HostProfileList
//
// Purpose:
//    This class contains a list of host profiles.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Mon Feb 25 15:14:55 PST 2008
//
// Modifications:
//   
// ****************************************************************************

public class HostProfileList extends AttributeSubject
{
    public HostProfileList()
    {
        super(2);

        profiles = new Vector();
        activeProfile = -1;
    }

    public HostProfileList(HostProfileList obj)
    {
        super(2);

        int i;

        // *** Copy the profiles field ***
        profiles = new Vector(obj.profiles.size());
        for(i = 0; i < obj.profiles.size(); ++i)
        {
            HostProfile newObj = (HostProfile)profiles.elementAt(i);
            profiles.addElement(new HostProfile(newObj));
        }

        activeProfile = obj.activeProfile;

        SelectAll();
    }

    public boolean equals(HostProfileList obj)
    {
        int i;

        boolean profiles_equal = (obj.profiles.size() == profiles.size());
        for(i = 0; (i < profiles.size()) && profiles_equal; ++i)
        {
            // Make references to HostProfile from Object.
            HostProfile profiles1 = (HostProfile)profiles.elementAt(i);
            HostProfile profiles2 = (HostProfile)obj.profiles.elementAt(i);
            profiles_equal = profiles1.equals(profiles2);
        }

        // Create the return value
        return (profiles_equal &&
                (activeProfile == obj.activeProfile));
    }

    // Property setting methods
    public void SetActiveProfile(int activeProfile_)
    {
        activeProfile = activeProfile_;
        Select(1);
    }

    // Property getting methods
    public Vector GetProfiles() { return profiles; }
    public int    GetActiveProfile() { return activeProfile; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
        {
            buf.WriteInt(profiles.size());
            for(int i = 0; i < profiles.size(); ++i)
            {
                HostProfile tmp = (HostProfile)profiles.elementAt(i);
                tmp.Write(buf);
            }
        }
        if(WriteSelect(1, buf))
            buf.WriteInt(activeProfile);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                {
                    int len = buf.ReadInt();
                    profiles.clear();
                    for(int j = 0; j < len; ++j)
                    {
                        HostProfile tmp = new HostProfile();
                        tmp.Read(buf);
                        profiles.addElement(tmp);
                    }
                }
                Select(0);
                break;
            case 1:
                SetActiveProfile(buf.ReadInt());
                break;
            }
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + indent + "profiles = {\n";
        for(int i = 0; i < profiles.size(); ++i)
        {
            AttributeSubject s = (AttributeSubject)profiles.elementAt(i);
            str = str + s.toString(indent + "    ");
            if(i < profiles.size()-1)
                str = str + ", ";
            str = str + "\n";
        }
        str = str + "}\n";
        str = str + intToString("activeProfile", activeProfile, indent) + "\n";
        return str;
    }

    // Attributegroup convenience methods
    public void AddProfiles(HostProfile obj)
    {
        profiles.addElement(new HostProfile(obj));
        Select(0);
    }

    public void ClearProfiles()
    {
        profiles.clear();
        Select(0);
    }

    public void RemoveProfiles(int index)
    {
        if(index >= 0 && index < profiles.size())
        {
            profiles.remove(index);
            Select(0);
        }
    }

    public int GetNumProfiles()
    {
        return profiles.size();
    }

    public HostProfile GetProfiles(int i)
    {
        HostProfile tmp = (HostProfile)profiles.elementAt(i);
        return tmp;
    }


    // Attributes
    private Vector profiles; // vector of HostProfile objects
    private int    activeProfile;
}

