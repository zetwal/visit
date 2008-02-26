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

package llnl.visit.plots;

import llnl.visit.AttributeSubject;
import llnl.visit.CommunicationBuffer;
import llnl.visit.Plugin;
import java.util.Vector;
import java.lang.Double;
import llnl.visit.ColorAttribute;

// ****************************************************************************
// Class: ParallelCoordinatesAttributes
//
// Purpose:
//    This class contains the plot attributes for the ParallelCoordinates plot.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Mon Feb 25 15:34:40 PST 2008
//
// Modifications:
//   
// ****************************************************************************

public class ParallelCoordinatesAttributes extends AttributeSubject implements Plugin
{
    public ParallelCoordinatesAttributes()
    {
        super(12);

        scalarAxisNames = new Vector();
        visualAxisNames = new Vector();
        extentMinima = new Vector();
        extentMaxima = new Vector();
        drawLines = true;
        linesColor = new ColorAttribute(128, 0, 0);
        drawContext = true;
        contextGamma = 2f;
        contextNumPartitions = 128;
        contextColor = new ColorAttribute(0, 220, 0);
        drawLinesOnlyIfExtentsOn = true;
        unifyAxisExtents = false;
    }

    public ParallelCoordinatesAttributes(ParallelCoordinatesAttributes obj)
    {
        super(12);

        int i;

        scalarAxisNames = new Vector(obj.scalarAxisNames.size());
        for(i = 0; i < obj.scalarAxisNames.size(); ++i)
            scalarAxisNames.addElement(new String((String)obj.scalarAxisNames.elementAt(i)));

        visualAxisNames = new Vector(obj.visualAxisNames.size());
        for(i = 0; i < obj.visualAxisNames.size(); ++i)
            visualAxisNames.addElement(new String((String)obj.visualAxisNames.elementAt(i)));

        extentMinima = new Vector(obj.extentMinima.size());
        for(i = 0; i < obj.extentMinima.size(); ++i)
        {
            Double dv = (Double)obj.extentMinima.elementAt(i);
            extentMinima.addElement(new Double(dv.doubleValue()));
        }

        extentMaxima = new Vector(obj.extentMaxima.size());
        for(i = 0; i < obj.extentMaxima.size(); ++i)
        {
            Double dv = (Double)obj.extentMaxima.elementAt(i);
            extentMaxima.addElement(new Double(dv.doubleValue()));
        }

        drawLines = obj.drawLines;
        linesColor = new ColorAttribute(obj.linesColor);
        drawContext = obj.drawContext;
        contextGamma = obj.contextGamma;
        contextNumPartitions = obj.contextNumPartitions;
        contextColor = new ColorAttribute(obj.contextColor);
        drawLinesOnlyIfExtentsOn = obj.drawLinesOnlyIfExtentsOn;
        unifyAxisExtents = obj.unifyAxisExtents;

        SelectAll();
    }

    public boolean equals(ParallelCoordinatesAttributes obj)
    {
        int i;

        // Create the return value
        return ((scalarAxisNames == obj.scalarAxisNames) &&
                (visualAxisNames == obj.visualAxisNames) &&
                (extentMinima == obj.extentMinima) &&
                (extentMaxima == obj.extentMaxima) &&
                (drawLines == obj.drawLines) &&
                (linesColor == obj.linesColor) &&
                (drawContext == obj.drawContext) &&
                (contextGamma == obj.contextGamma) &&
                (contextNumPartitions == obj.contextNumPartitions) &&
                (contextColor == obj.contextColor) &&
                (drawLinesOnlyIfExtentsOn == obj.drawLinesOnlyIfExtentsOn) &&
                (unifyAxisExtents == obj.unifyAxisExtents));
    }

    public String GetName() { return "ParallelCoordinates"; }
    public String GetVersion() { return "1.0"; }

    // Property setting methods
    public void SetScalarAxisNames(Vector scalarAxisNames_)
    {
        scalarAxisNames = scalarAxisNames_;
        Select(0);
    }

    public void SetVisualAxisNames(Vector visualAxisNames_)
    {
        visualAxisNames = visualAxisNames_;
        Select(1);
    }

    public void SetExtentMinima(Vector extentMinima_)
    {
        extentMinima = extentMinima_;
        Select(2);
    }

    public void SetExtentMaxima(Vector extentMaxima_)
    {
        extentMaxima = extentMaxima_;
        Select(3);
    }

    public void SetDrawLines(boolean drawLines_)
    {
        drawLines = drawLines_;
        Select(4);
    }

    public void SetLinesColor(ColorAttribute linesColor_)
    {
        linesColor = linesColor_;
        Select(5);
    }

    public void SetDrawContext(boolean drawContext_)
    {
        drawContext = drawContext_;
        Select(6);
    }

    public void SetContextGamma(float contextGamma_)
    {
        contextGamma = contextGamma_;
        Select(7);
    }

    public void SetContextNumPartitions(int contextNumPartitions_)
    {
        contextNumPartitions = contextNumPartitions_;
        Select(8);
    }

    public void SetContextColor(ColorAttribute contextColor_)
    {
        contextColor = contextColor_;
        Select(9);
    }

    public void SetDrawLinesOnlyIfExtentsOn(boolean drawLinesOnlyIfExtentsOn_)
    {
        drawLinesOnlyIfExtentsOn = drawLinesOnlyIfExtentsOn_;
        Select(10);
    }

    public void SetUnifyAxisExtents(boolean unifyAxisExtents_)
    {
        unifyAxisExtents = unifyAxisExtents_;
        Select(11);
    }

    // Property getting methods
    public Vector         GetScalarAxisNames() { return scalarAxisNames; }
    public Vector         GetVisualAxisNames() { return visualAxisNames; }
    public Vector         GetExtentMinima() { return extentMinima; }
    public Vector         GetExtentMaxima() { return extentMaxima; }
    public boolean        GetDrawLines() { return drawLines; }
    public ColorAttribute GetLinesColor() { return linesColor; }
    public boolean        GetDrawContext() { return drawContext; }
    public float          GetContextGamma() { return contextGamma; }
    public int            GetContextNumPartitions() { return contextNumPartitions; }
    public ColorAttribute GetContextColor() { return contextColor; }
    public boolean        GetDrawLinesOnlyIfExtentsOn() { return drawLinesOnlyIfExtentsOn; }
    public boolean        GetUnifyAxisExtents() { return unifyAxisExtents; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteStringVector(scalarAxisNames);
        if(WriteSelect(1, buf))
            buf.WriteStringVector(visualAxisNames);
        if(WriteSelect(2, buf))
            buf.WriteDoubleVector(extentMinima);
        if(WriteSelect(3, buf))
            buf.WriteDoubleVector(extentMaxima);
        if(WriteSelect(4, buf))
            buf.WriteBool(drawLines);
        if(WriteSelect(5, buf))
            linesColor.Write(buf);
        if(WriteSelect(6, buf))
            buf.WriteBool(drawContext);
        if(WriteSelect(7, buf))
            buf.WriteFloat(contextGamma);
        if(WriteSelect(8, buf))
            buf.WriteInt(contextNumPartitions);
        if(WriteSelect(9, buf))
            contextColor.Write(buf);
        if(WriteSelect(10, buf))
            buf.WriteBool(drawLinesOnlyIfExtentsOn);
        if(WriteSelect(11, buf))
            buf.WriteBool(unifyAxisExtents);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetScalarAxisNames(buf.ReadStringVector());
                break;
            case 1:
                SetVisualAxisNames(buf.ReadStringVector());
                break;
            case 2:
                SetExtentMinima(buf.ReadDoubleVector());
                break;
            case 3:
                SetExtentMaxima(buf.ReadDoubleVector());
                break;
            case 4:
                SetDrawLines(buf.ReadBool());
                break;
            case 5:
                linesColor.Read(buf);
                Select(5);
                break;
            case 6:
                SetDrawContext(buf.ReadBool());
                break;
            case 7:
                SetContextGamma(buf.ReadFloat());
                break;
            case 8:
                SetContextNumPartitions(buf.ReadInt());
                break;
            case 9:
                contextColor.Read(buf);
                Select(9);
                break;
            case 10:
                SetDrawLinesOnlyIfExtentsOn(buf.ReadBool());
                break;
            case 11:
                SetUnifyAxisExtents(buf.ReadBool());
                break;
            }
        }
    }

    public String toString(String indent)
    {
        String str = new String();
        str = str + stringVectorToString("scalarAxisNames", scalarAxisNames, indent) + "\n";
        str = str + stringVectorToString("visualAxisNames", visualAxisNames, indent) + "\n";
        str = str + doubleVectorToString("extentMinima", extentMinima, indent) + "\n";
        str = str + doubleVectorToString("extentMaxima", extentMaxima, indent) + "\n";
        str = str + boolToString("drawLines", drawLines, indent) + "\n";
        str = str + indent + "linesColor = {" + linesColor.Red() + ", " + linesColor.Green() + ", " + linesColor.Blue() + ", " + linesColor.Alpha() + "}\n";
        str = str + boolToString("drawContext", drawContext, indent) + "\n";
        str = str + floatToString("contextGamma", contextGamma, indent) + "\n";
        str = str + intToString("contextNumPartitions", contextNumPartitions, indent) + "\n";
        str = str + indent + "contextColor = {" + contextColor.Red() + ", " + contextColor.Green() + ", " + contextColor.Blue() + ", " + contextColor.Alpha() + "}\n";
        str = str + boolToString("drawLinesOnlyIfExtentsOn", drawLinesOnlyIfExtentsOn, indent) + "\n";
        str = str + boolToString("unifyAxisExtents", unifyAxisExtents, indent) + "\n";
        return str;
    }


    // Attributes
    private Vector         scalarAxisNames; // vector of String objects
    private Vector         visualAxisNames; // vector of String objects
    private Vector         extentMinima; // vector of Double objects
    private Vector         extentMaxima; // vector of Double objects
    private boolean        drawLines;
    private ColorAttribute linesColor;
    private boolean        drawContext;
    private float          contextGamma;
    private int            contextNumPartitions;
    private ColorAttribute contextColor;
    private boolean        drawLinesOnlyIfExtentsOn;
    private boolean        unifyAxisExtents;
}

