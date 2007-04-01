package llnl.visit;

import java.util.Vector;
import java.lang.Integer;

// ****************************************************************************
// Class: ViewerRPC
//
// Purpose:
//    This class contains the attributes for controlling the viewer.
//
// Notes:      Autogenerated by xml2java.
//
// Programmer: xml2java
// Creation:   Wed Mar 31 10:50:16 PDT 2004
//
// Modifications:
//   
// ****************************************************************************

public class ViewerRPC extends AttributeSubject
{
    // Constants
    public final static int VIEWERRPCTYPE_CLOSERPC = 0;
    public final static int VIEWERRPCTYPE_ADDWINDOWRPC = 1;
    public final static int VIEWERRPCTYPE_DELETEWINDOWRPC = 2;
    public final static int VIEWERRPCTYPE_SETWINDOWLAYOUTRPC = 3;
    public final static int VIEWERRPCTYPE_SETACTIVEWINDOWRPC = 4;
    public final static int VIEWERRPCTYPE_CLEARWINDOWRPC = 5;
    public final static int VIEWERRPCTYPE_CLEARALLWINDOWSRPC = 6;
    public final static int VIEWERRPCTYPE_OPENDATABASERPC = 7;
    public final static int VIEWERRPCTYPE_CLOSEDATABASERPC = 8;
    public final static int VIEWERRPCTYPE_ACTIVATEDATABASERPC = 9;
    public final static int VIEWERRPCTYPE_CHECKFORNEWSTATESRPC = 10;
    public final static int VIEWERRPCTYPE_CREATEDATABASECORRELATIONRPC = 11;
    public final static int VIEWERRPCTYPE_ALTERDATABASECORRELATIONRPC = 12;
    public final static int VIEWERRPCTYPE_DELETEDATABASECORRELATIONRPC = 13;
    public final static int VIEWERRPCTYPE_REOPENDATABASERPC = 14;
    public final static int VIEWERRPCTYPE_REPLACEDATABASERPC = 15;
    public final static int VIEWERRPCTYPE_OVERLAYDATABASERPC = 16;
    public final static int VIEWERRPCTYPE_OPENCOMPUTEENGINERPC = 17;
    public final static int VIEWERRPCTYPE_CLOSECOMPUTEENGINERPC = 18;
    public final static int VIEWERRPCTYPE_ANIMATIONSETNFRAMESRPC = 19;
    public final static int VIEWERRPCTYPE_ANIMATIONPLAYRPC = 20;
    public final static int VIEWERRPCTYPE_ANIMATIONREVERSEPLAYRPC = 21;
    public final static int VIEWERRPCTYPE_ANIMATIONSTOPRPC = 22;
    public final static int VIEWERRPCTYPE_TIMESLIDERNEXTSTATERPC = 23;
    public final static int VIEWERRPCTYPE_TIMESLIDERPREVIOUSSTATERPC = 24;
    public final static int VIEWERRPCTYPE_SETTIMESLIDERSTATERPC = 25;
    public final static int VIEWERRPCTYPE_SETACTIVETIMESLIDERRPC = 26;
    public final static int VIEWERRPCTYPE_ADDPLOTRPC = 27;
    public final static int VIEWERRPCTYPE_SETPLOTFRAMERANGERPC = 28;
    public final static int VIEWERRPCTYPE_DELETEPLOTKEYFRAMERPC = 29;
    public final static int VIEWERRPCTYPE_MOVEPLOTKEYFRAMERPC = 30;
    public final static int VIEWERRPCTYPE_DELETEACTIVEPLOTSRPC = 31;
    public final static int VIEWERRPCTYPE_HIDEACTIVEPLOTSRPC = 32;
    public final static int VIEWERRPCTYPE_DRAWPLOTSRPC = 33;
    public final static int VIEWERRPCTYPE_DISABLEREDRAWRPC = 34;
    public final static int VIEWERRPCTYPE_REDRAWRPC = 35;
    public final static int VIEWERRPCTYPE_SETACTIVEPLOTSRPC = 36;
    public final static int VIEWERRPCTYPE_CHANGEACTIVEPLOTSVARRPC = 37;
    public final static int VIEWERRPCTYPE_ADDOPERATORRPC = 38;
    public final static int VIEWERRPCTYPE_PROMOTEOPERATORRPC = 39;
    public final static int VIEWERRPCTYPE_DEMOTEOPERATORRPC = 40;
    public final static int VIEWERRPCTYPE_REMOVEOPERATORRPC = 41;
    public final static int VIEWERRPCTYPE_REMOVELASTOPERATORRPC = 42;
    public final static int VIEWERRPCTYPE_REMOVEALLOPERATORSRPC = 43;
    public final static int VIEWERRPCTYPE_SAVEWINDOWRPC = 44;
    public final static int VIEWERRPCTYPE_SETDEFAULTPLOTOPTIONSRPC = 45;
    public final static int VIEWERRPCTYPE_SETPLOTOPTIONSRPC = 46;
    public final static int VIEWERRPCTYPE_SETDEFAULTOPERATOROPTIONSRPC = 47;
    public final static int VIEWERRPCTYPE_SETOPERATOROPTIONSRPC = 48;
    public final static int VIEWERRPCTYPE_WRITECONFIGFILERPC = 49;
    public final static int VIEWERRPCTYPE_CONNECTTOMETADATASERVERRPC = 50;
    public final static int VIEWERRPCTYPE_ICONIFYALLWINDOWSRPC = 51;
    public final static int VIEWERRPCTYPE_DEICONIFYALLWINDOWSRPC = 52;
    public final static int VIEWERRPCTYPE_SHOWALLWINDOWSRPC = 53;
    public final static int VIEWERRPCTYPE_HIDEALLWINDOWSRPC = 54;
    public final static int VIEWERRPCTYPE_UPDATECOLORTABLERPC = 55;
    public final static int VIEWERRPCTYPE_SETANNOTATIONATTRIBUTESRPC = 56;
    public final static int VIEWERRPCTYPE_SETDEFAULTANNOTATIONATTRIBUTESRPC = 57;
    public final static int VIEWERRPCTYPE_RESETANNOTATIONATTRIBUTESRPC = 58;
    public final static int VIEWERRPCTYPE_SETKEYFRAMEATTRIBUTESRPC = 59;
    public final static int VIEWERRPCTYPE_SETPLOTSILRESTRICTIONRPC = 60;
    public final static int VIEWERRPCTYPE_SETVIEWCURVERPC = 61;
    public final static int VIEWERRPCTYPE_SETVIEW2DRPC = 62;
    public final static int VIEWERRPCTYPE_SETVIEW3DRPC = 63;
    public final static int VIEWERRPCTYPE_RESETPLOTOPTIONSRPC = 64;
    public final static int VIEWERRPCTYPE_RESETOPERATOROPTIONSRPC = 65;
    public final static int VIEWERRPCTYPE_SETAPPEARANCERPC = 66;
    public final static int VIEWERRPCTYPE_PROCESSEXPRESSIONSRPC = 67;
    public final static int VIEWERRPCTYPE_SETLIGHTLISTRPC = 68;
    public final static int VIEWERRPCTYPE_SETDEFAULTLIGHTLISTRPC = 69;
    public final static int VIEWERRPCTYPE_RESETLIGHTLISTRPC = 70;
    public final static int VIEWERRPCTYPE_SETANIMATIONATTRIBUTESRPC = 71;
    public final static int VIEWERRPCTYPE_SETWINDOWAREARPC = 72;
    public final static int VIEWERRPCTYPE_PRINTWINDOWRPC = 73;
    public final static int VIEWERRPCTYPE_RESETVIEWRPC = 74;
    public final static int VIEWERRPCTYPE_RECENTERVIEWRPC = 75;
    public final static int VIEWERRPCTYPE_TOGGLEMAINTAINVIEWMODERPC = 76;
    public final static int VIEWERRPCTYPE_TOGGLEBOUNDINGBOXMODERPC = 77;
    public final static int VIEWERRPCTYPE_TOGGLECAMERAVIEWMODERPC = 78;
    public final static int VIEWERRPCTYPE_TOGGLEPERSPECTIVEVIEWRPC = 79;
    public final static int VIEWERRPCTYPE_TOGGLESPINMODERPC = 80;
    public final static int VIEWERRPCTYPE_TOGGLELOCKTIMERPC = 81;
    public final static int VIEWERRPCTYPE_TOGGLELOCKTOOLSRPC = 82;
    public final static int VIEWERRPCTYPE_TOGGLELOCKVIEWMODERPC = 83;
    public final static int VIEWERRPCTYPE_TOGGLEFULLFRAMERPC = 84;
    public final static int VIEWERRPCTYPE_UNDOVIEWRPC = 85;
    public final static int VIEWERRPCTYPE_INVERTBACKGROUNDRPC = 86;
    public final static int VIEWERRPCTYPE_CLEARPICKPOINTSRPC = 87;
    public final static int VIEWERRPCTYPE_SETWINDOWMODERPC = 88;
    public final static int VIEWERRPCTYPE_ENABLETOOLRPC = 89;
    public final static int VIEWERRPCTYPE_COPYVIEWTOWINDOWRPC = 90;
    public final static int VIEWERRPCTYPE_COPYLIGHTINGTOWINDOWRPC = 91;
    public final static int VIEWERRPCTYPE_COPYANNOTATIONSTOWINDOWRPC = 92;
    public final static int VIEWERRPCTYPE_COPYPLOTSTOWINDOWRPC = 93;
    public final static int VIEWERRPCTYPE_CLEARCACHERPC = 94;
    public final static int VIEWERRPCTYPE_CLEARCACHEFORALLENGINESRPC = 95;
    public final static int VIEWERRPCTYPE_SETVIEWEXTENTSTYPERPC = 96;
    public final static int VIEWERRPCTYPE_CLEARREFLINESRPC = 97;
    public final static int VIEWERRPCTYPE_SETRENDERINGATTRIBUTESRPC = 98;
    public final static int VIEWERRPCTYPE_DATABASEQUERYRPC = 99;
    public final static int VIEWERRPCTYPE_POINTQUERYRPC = 100;
    public final static int VIEWERRPCTYPE_LINEQUERYRPC = 101;
    public final static int VIEWERRPCTYPE_CLONEWINDOWRPC = 102;
    public final static int VIEWERRPCTYPE_SETMATERIALATTRIBUTESRPC = 103;
    public final static int VIEWERRPCTYPE_SETDEFAULTMATERIALATTRIBUTESRPC = 104;
    public final static int VIEWERRPCTYPE_RESETMATERIALATTRIBUTESRPC = 105;
    public final static int VIEWERRPCTYPE_SETPLOTDATABASESTATERPC = 106;
    public final static int VIEWERRPCTYPE_DELETEPLOTDATABASEKEYFRAMERPC = 107;
    public final static int VIEWERRPCTYPE_MOVEPLOTDATABASEKEYFRAMERPC = 108;
    public final static int VIEWERRPCTYPE_CLEARVIEWKEYFRAMESRPC = 109;
    public final static int VIEWERRPCTYPE_DELETEVIEWKEYFRAMERPC = 110;
    public final static int VIEWERRPCTYPE_MOVEVIEWKEYFRAMERPC = 111;
    public final static int VIEWERRPCTYPE_SETVIEWKEYFRAMERPC = 112;
    public final static int VIEWERRPCTYPE_OPENMDSERVERRPC = 113;
    public final static int VIEWERRPCTYPE_ENABLETOOLBARRPC = 114;
    public final static int VIEWERRPCTYPE_HIDETOOLBARSRPC = 115;
    public final static int VIEWERRPCTYPE_HIDETOOLBARSFORALLWINDOWSRPC = 116;
    public final static int VIEWERRPCTYPE_SHOWTOOLBARSRPC = 117;
    public final static int VIEWERRPCTYPE_SHOWTOOLBARSFORALLWINDOWSRPC = 118;
    public final static int VIEWERRPCTYPE_SETTOOLBARICONSIZERPC = 119;
    public final static int VIEWERRPCTYPE_SAVEVIEWRPC = 120;
    public final static int VIEWERRPCTYPE_SETGLOBALLINEOUTATTRIBUTESRPC = 121;
    public final static int VIEWERRPCTYPE_SETPICKATTRIBUTESRPC = 122;
    public final static int VIEWERRPCTYPE_EXPORTCOLORTABLERPC = 123;
    public final static int VIEWERRPCTYPE_EXPORTENTIRESTATERPC = 124;
    public final static int VIEWERRPCTYPE_IMPORTENTIRESTATERPC = 125;
    public final static int VIEWERRPCTYPE_RESETPICKATTRIBUTESRPC = 126;
    public final static int VIEWERRPCTYPE_ADDANNOTATIONOBJECTRPC = 127;
    public final static int VIEWERRPCTYPE_HIDEACTIVEANNOTATIONOBJECTSRPC = 128;
    public final static int VIEWERRPCTYPE_DELETEACTIVEANNOTATIONOBJECTSRPC = 129;
    public final static int VIEWERRPCTYPE_RAISEACTIVEANNOTATIONOBJECTSRPC = 130;
    public final static int VIEWERRPCTYPE_LOWERACTIVEANNOTATIONOBJECTSRPC = 131;
    public final static int VIEWERRPCTYPE_SETANNOTATIONOBJECTOPTIONSRPC = 132;
    public final static int VIEWERRPCTYPE_SETDEFAULTANNOTATIONOBJECTLISTRPC = 133;
    public final static int VIEWERRPCTYPE_RESETANNOTATIONOBJECTLISTRPC = 134;
    public final static int VIEWERRPCTYPE_RESETPICKLETTERRPC = 135;
    public final static int VIEWERRPCTYPE_SETDEFAULTPICKATTRIBUTESRPC = 136;
    public final static int VIEWERRPCTYPE_CHOOSECENTEROFROTATIONRPC = 137;
    public final static int VIEWERRPCTYPE_SETCENTEROFROTATIONRPC = 138;
    public final static int VIEWERRPCTYPE_SETQUERYOVERTIMEATTRIBUTESRPC = 139;
    public final static int VIEWERRPCTYPE_SETDEFAULTQUERYOVERTIMEATTRIBUTESRPC = 140;
    public final static int VIEWERRPCTYPE_RESETQUERYOVERTIMEATTRIBUTESRPC = 141;
    public final static int VIEWERRPCTYPE_MAXRPC = 142;


    public ViewerRPC()
    {
        super(29);

        RPCType = VIEWERRPCTYPE_CLOSERPC;
        windowLayout = 1;
        windowId = 0;
        windowMode = 0;
        windowArea = new String("(null)");
        database = new String("(null)");
        programHost = new String("(null)");
        programSim = new String("(null)");
        programOptions = new Vector();
        nFrames = 0;
        stateNumber = 0;
        frameRange = new int[2];
        frameRange[0] = 0;
        frameRange[1] = 0;
        frame = 0;
        plotType = 0;
        operatorType = 0;
        variable = new String("(null)");
        activePlotIds = new Vector();
        activeOperatorIds = new Vector();
        expandedPlotIds = new Vector();
        colorTableName = new String("(null)");
        queryName = new String("(null)");
        queryPoint1 = new double[3];
        queryPoint1[0] = 0;
        queryPoint1[1] = 0;
        queryPoint1[2] = 0;
        queryPoint2 = new double[3];
        queryPoint2[0] = 0;
        queryPoint2[1] = 0;
        queryPoint2[2] = 0;
        queryVariables = new Vector();
        toolId = 0;
        boolFlag = false;
        intArg1 = 0;
        intArg2 = 0;
        intArg3 = 0;
    }

    public ViewerRPC(ViewerRPC obj)
    {
        super(29);

        int i;

        RPCType = obj.RPCType;
        windowLayout = obj.windowLayout;
        windowId = obj.windowId;
        windowMode = obj.windowMode;
        windowArea = new String(obj.windowArea);
        database = new String(obj.database);
        programHost = new String(obj.programHost);
        programSim = new String(obj.programSim);
        programOptions = new Vector(obj.programOptions.size());
        for(i = 0; i < obj.programOptions.size(); ++i)
            programOptions.addElement(new String((String)obj.programOptions.elementAt(i)));

        nFrames = obj.nFrames;
        stateNumber = obj.stateNumber;
        frameRange = new int[2];
        frameRange[0] = obj.frameRange[0];
        frameRange[1] = obj.frameRange[1];

        frame = obj.frame;
        plotType = obj.plotType;
        operatorType = obj.operatorType;
        variable = new String(obj.variable);
        activePlotIds = new Vector();
        for(i = 0; i < obj.activePlotIds.size(); ++i)
        {
            Integer iv = (Integer)obj.activePlotIds.elementAt(i);
            activePlotIds.addElement(new Integer(iv.intValue()));
        }
        activeOperatorIds = new Vector();
        for(i = 0; i < obj.activeOperatorIds.size(); ++i)
        {
            Integer iv = (Integer)obj.activeOperatorIds.elementAt(i);
            activeOperatorIds.addElement(new Integer(iv.intValue()));
        }
        expandedPlotIds = new Vector();
        for(i = 0; i < obj.expandedPlotIds.size(); ++i)
        {
            Integer iv = (Integer)obj.expandedPlotIds.elementAt(i);
            expandedPlotIds.addElement(new Integer(iv.intValue()));
        }
        colorTableName = new String(obj.colorTableName);
        queryName = new String(obj.queryName);
        queryPoint1 = new double[3];
        queryPoint1[0] = obj.queryPoint1[0];
        queryPoint1[1] = obj.queryPoint1[1];
        queryPoint1[2] = obj.queryPoint1[2];

        queryPoint2 = new double[3];
        queryPoint2[0] = obj.queryPoint2[0];
        queryPoint2[1] = obj.queryPoint2[1];
        queryPoint2[2] = obj.queryPoint2[2];

        queryVariables = new Vector(obj.queryVariables.size());
        for(i = 0; i < obj.queryVariables.size(); ++i)
            queryVariables.addElement(new String((String)obj.queryVariables.elementAt(i)));

        toolId = obj.toolId;
        boolFlag = obj.boolFlag;
        intArg1 = obj.intArg1;
        intArg2 = obj.intArg2;
        intArg3 = obj.intArg3;

        SelectAll();
    }

    public boolean equals(ViewerRPC obj)
    {
        int i;

        // Compare the frameRange arrays.
        boolean frameRange_equal = true;
        for(i = 0; i < 2 && frameRange_equal; ++i)
            frameRange_equal = (frameRange[i] == obj.frameRange[i]);

        // Compare the queryPoint1 arrays.
        boolean queryPoint1_equal = true;
        for(i = 0; i < 3 && queryPoint1_equal; ++i)
            queryPoint1_equal = (queryPoint1[i] == obj.queryPoint1[i]);

        // Compare the queryPoint2 arrays.
        boolean queryPoint2_equal = true;
        for(i = 0; i < 3 && queryPoint2_equal; ++i)
            queryPoint2_equal = (queryPoint2[i] == obj.queryPoint2[i]);

        // Create the return value
        return ((RPCType == obj.RPCType) &&
                (windowLayout == obj.windowLayout) &&
                (windowId == obj.windowId) &&
                (windowMode == obj.windowMode) &&
                (windowArea == obj.windowArea) &&
                (database == obj.database) &&
                (programHost == obj.programHost) &&
                (programSim == obj.programSim) &&
                (programOptions == obj.programOptions) &&
                (nFrames == obj.nFrames) &&
                (stateNumber == obj.stateNumber) &&
                frameRange_equal &&
                (frame == obj.frame) &&
                (plotType == obj.plotType) &&
                (operatorType == obj.operatorType) &&
                (variable == obj.variable) &&
                (activePlotIds == obj.activePlotIds) &&
                (activeOperatorIds == obj.activeOperatorIds) &&
                (expandedPlotIds == obj.expandedPlotIds) &&
                (colorTableName == obj.colorTableName) &&
                (queryName == obj.queryName) &&
                queryPoint1_equal &&
                queryPoint2_equal &&
                (queryVariables == obj.queryVariables) &&
                (toolId == obj.toolId) &&
                (boolFlag == obj.boolFlag) &&
                (intArg1 == obj.intArg1) &&
                (intArg2 == obj.intArg2) &&
                (intArg3 == obj.intArg3));
    }

    // Property setting methods
    public void SetRPCType(int RPCType_)
    {
        RPCType = RPCType_;
        Select(0);
    }

    public void SetWindowLayout(int windowLayout_)
    {
        windowLayout = windowLayout_;
        Select(1);
    }

    public void SetWindowId(int windowId_)
    {
        windowId = windowId_;
        Select(2);
    }

    public void SetWindowMode(int windowMode_)
    {
        windowMode = windowMode_;
        Select(3);
    }

    public void SetWindowArea(String windowArea_)
    {
        windowArea = windowArea_;
        Select(4);
    }

    public void SetDatabase(String database_)
    {
        database = database_;
        Select(5);
    }

    public void SetProgramHost(String programHost_)
    {
        programHost = programHost_;
        Select(6);
    }

    public void SetProgramSim(String programSim_)
    {
        programSim = programSim_;
        Select(7);
    }

    public void SetProgramOptions(Vector programOptions_)
    {
        programOptions = programOptions_;
        Select(8);
    }

    public void SetNFrames(int nFrames_)
    {
        nFrames = nFrames_;
        Select(9);
    }

    public void SetStateNumber(int stateNumber_)
    {
        stateNumber = stateNumber_;
        Select(10);
    }

    public void SetFrameRange(int[] frameRange_)
    {
        frameRange[0] = frameRange_[0];
        frameRange[1] = frameRange_[1];
        Select(11);
    }

    public void SetFrameRange(int e0, int e1)
    {
        frameRange[0] = e0;
        frameRange[1] = e1;
        Select(11);
    }

    public void SetFrame(int frame_)
    {
        frame = frame_;
        Select(12);
    }

    public void SetPlotType(int plotType_)
    {
        plotType = plotType_;
        Select(13);
    }

    public void SetOperatorType(int operatorType_)
    {
        operatorType = operatorType_;
        Select(14);
    }

    public void SetVariable(String variable_)
    {
        variable = variable_;
        Select(15);
    }

    public void SetActivePlotIds(Vector activePlotIds_)
    {
        activePlotIds = activePlotIds_;
        Select(16);
    }

    public void SetActiveOperatorIds(Vector activeOperatorIds_)
    {
        activeOperatorIds = activeOperatorIds_;
        Select(17);
    }

    public void SetExpandedPlotIds(Vector expandedPlotIds_)
    {
        expandedPlotIds = expandedPlotIds_;
        Select(18);
    }

    public void SetColorTableName(String colorTableName_)
    {
        colorTableName = colorTableName_;
        Select(19);
    }

    public void SetQueryName(String queryName_)
    {
        queryName = queryName_;
        Select(20);
    }

    public void SetQueryPoint1(double[] queryPoint1_)
    {
        queryPoint1[0] = queryPoint1_[0];
        queryPoint1[1] = queryPoint1_[1];
        queryPoint1[2] = queryPoint1_[2];
        Select(21);
    }

    public void SetQueryPoint1(double e0, double e1, double e2)
    {
        queryPoint1[0] = e0;
        queryPoint1[1] = e1;
        queryPoint1[2] = e2;
        Select(21);
    }

    public void SetQueryPoint2(double[] queryPoint2_)
    {
        queryPoint2[0] = queryPoint2_[0];
        queryPoint2[1] = queryPoint2_[1];
        queryPoint2[2] = queryPoint2_[2];
        Select(22);
    }

    public void SetQueryPoint2(double e0, double e1, double e2)
    {
        queryPoint2[0] = e0;
        queryPoint2[1] = e1;
        queryPoint2[2] = e2;
        Select(22);
    }

    public void SetQueryVariables(Vector queryVariables_)
    {
        queryVariables = queryVariables_;
        Select(23);
    }

    public void SetToolId(int toolId_)
    {
        toolId = toolId_;
        Select(24);
    }

    public void SetBoolFlag(boolean boolFlag_)
    {
        boolFlag = boolFlag_;
        Select(25);
    }

    public void SetIntArg1(int intArg1_)
    {
        intArg1 = intArg1_;
        Select(26);
    }

    public void SetIntArg2(int intArg2_)
    {
        intArg2 = intArg2_;
        Select(27);
    }

    public void SetIntArg3(int intArg3_)
    {
        intArg3 = intArg3_;
        Select(28);
    }

    // Property getting methods
    public int      GetRPCType() { return RPCType; }
    public int      GetWindowLayout() { return windowLayout; }
    public int      GetWindowId() { return windowId; }
    public int      GetWindowMode() { return windowMode; }
    public String   GetWindowArea() { return windowArea; }
    public String   GetDatabase() { return database; }
    public String   GetProgramHost() { return programHost; }
    public String   GetProgramSim() { return programSim; }
    public Vector   GetProgramOptions() { return programOptions; }
    public int      GetNFrames() { return nFrames; }
    public int      GetStateNumber() { return stateNumber; }
    public int[]    GetFrameRange() { return frameRange; }
    public int      GetFrame() { return frame; }
    public int      GetPlotType() { return plotType; }
    public int      GetOperatorType() { return operatorType; }
    public String   GetVariable() { return variable; }
    public Vector   GetActivePlotIds() { return activePlotIds; }
    public Vector   GetActiveOperatorIds() { return activeOperatorIds; }
    public Vector   GetExpandedPlotIds() { return expandedPlotIds; }
    public String   GetColorTableName() { return colorTableName; }
    public String   GetQueryName() { return queryName; }
    public double[] GetQueryPoint1() { return queryPoint1; }
    public double[] GetQueryPoint2() { return queryPoint2; }
    public Vector   GetQueryVariables() { return queryVariables; }
    public int      GetToolId() { return toolId; }
    public boolean  GetBoolFlag() { return boolFlag; }
    public int      GetIntArg1() { return intArg1; }
    public int      GetIntArg2() { return intArg2; }
    public int      GetIntArg3() { return intArg3; }

    // Write and read methods.
    public void WriteAtts(CommunicationBuffer buf)
    {
        if(WriteSelect(0, buf))
            buf.WriteInt(RPCType);
        if(WriteSelect(1, buf))
            buf.WriteInt(windowLayout);
        if(WriteSelect(2, buf))
            buf.WriteInt(windowId);
        if(WriteSelect(3, buf))
            buf.WriteInt(windowMode);
        if(WriteSelect(4, buf))
            buf.WriteString(windowArea);
        if(WriteSelect(5, buf))
            buf.WriteString(database);
        if(WriteSelect(6, buf))
            buf.WriteString(programHost);
        if(WriteSelect(7, buf))
            buf.WriteString(programSim);
        if(WriteSelect(8, buf))
            buf.WriteStringVector(programOptions);
        if(WriteSelect(9, buf))
            buf.WriteInt(nFrames);
        if(WriteSelect(10, buf))
            buf.WriteInt(stateNumber);
        if(WriteSelect(11, buf))
            buf.WriteIntArray(frameRange);
        if(WriteSelect(12, buf))
            buf.WriteInt(frame);
        if(WriteSelect(13, buf))
            buf.WriteInt(plotType);
        if(WriteSelect(14, buf))
            buf.WriteInt(operatorType);
        if(WriteSelect(15, buf))
            buf.WriteString(variable);
        if(WriteSelect(16, buf))
            buf.WriteIntVector(activePlotIds);
        if(WriteSelect(17, buf))
            buf.WriteIntVector(activeOperatorIds);
        if(WriteSelect(18, buf))
            buf.WriteIntVector(expandedPlotIds);
        if(WriteSelect(19, buf))
            buf.WriteString(colorTableName);
        if(WriteSelect(20, buf))
            buf.WriteString(queryName);
        if(WriteSelect(21, buf))
            buf.WriteDoubleArray(queryPoint1);
        if(WriteSelect(22, buf))
            buf.WriteDoubleArray(queryPoint2);
        if(WriteSelect(23, buf))
            buf.WriteStringVector(queryVariables);
        if(WriteSelect(24, buf))
            buf.WriteInt(toolId);
        if(WriteSelect(25, buf))
            buf.WriteBool(boolFlag);
        if(WriteSelect(26, buf))
            buf.WriteInt(intArg1);
        if(WriteSelect(27, buf))
            buf.WriteInt(intArg2);
        if(WriteSelect(28, buf))
            buf.WriteInt(intArg3);
    }

    public void ReadAtts(int n, CommunicationBuffer buf)
    {
        for(int i = 0; i < n; ++i)
        {
            int index = (int)buf.ReadByte();
            switch(index)
            {
            case 0:
                SetRPCType(buf.ReadInt());
                break;
            case 1:
                SetWindowLayout(buf.ReadInt());
                break;
            case 2:
                SetWindowId(buf.ReadInt());
                break;
            case 3:
                SetWindowMode(buf.ReadInt());
                break;
            case 4:
                SetWindowArea(buf.ReadString());
                break;
            case 5:
                SetDatabase(buf.ReadString());
                break;
            case 6:
                SetProgramHost(buf.ReadString());
                break;
            case 7:
                SetProgramSim(buf.ReadString());
                break;
            case 8:
                SetProgramOptions(buf.ReadStringVector());
                break;
            case 9:
                SetNFrames(buf.ReadInt());
                break;
            case 10:
                SetStateNumber(buf.ReadInt());
                break;
            case 11:
                SetFrameRange(buf.ReadIntArray());
                break;
            case 12:
                SetFrame(buf.ReadInt());
                break;
            case 13:
                SetPlotType(buf.ReadInt());
                break;
            case 14:
                SetOperatorType(buf.ReadInt());
                break;
            case 15:
                SetVariable(buf.ReadString());
                break;
            case 16:
                SetActivePlotIds(buf.ReadIntVector());
                break;
            case 17:
                SetActiveOperatorIds(buf.ReadIntVector());
                break;
            case 18:
                SetExpandedPlotIds(buf.ReadIntVector());
                break;
            case 19:
                SetColorTableName(buf.ReadString());
                break;
            case 20:
                SetQueryName(buf.ReadString());
                break;
            case 21:
                SetQueryPoint1(buf.ReadDoubleArray());
                break;
            case 22:
                SetQueryPoint2(buf.ReadDoubleArray());
                break;
            case 23:
                SetQueryVariables(buf.ReadStringVector());
                break;
            case 24:
                SetToolId(buf.ReadInt());
                break;
            case 25:
                SetBoolFlag(buf.ReadBool());
                break;
            case 26:
                SetIntArg1(buf.ReadInt());
                break;
            case 27:
                SetIntArg2(buf.ReadInt());
                break;
            case 28:
                SetIntArg3(buf.ReadInt());
                break;
            }
        }
    }


    // Attributes
    private int      RPCType;
    private int      windowLayout;
    private int      windowId;
    private int      windowMode;
    private String   windowArea;
    private String   database;
    private String   programHost;
    private String   programSim;
    private Vector   programOptions; // vector of String objects
    private int      nFrames;
    private int      stateNumber;
    private int[]    frameRange;
    private int      frame;
    private int      plotType;
    private int      operatorType;
    private String   variable;
    private Vector   activePlotIds; // vector of Integer objects
    private Vector   activeOperatorIds; // vector of Integer objects
    private Vector   expandedPlotIds; // vector of Integer objects
    private String   colorTableName;
    private String   queryName;
    private double[] queryPoint1;
    private double[] queryPoint2;
    private Vector   queryVariables; // vector of String objects
    private int      toolId;
    private boolean  boolFlag;
    private int      intArg1;
    private int      intArg2;
    private int      intArg3;
}

