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

#ifndef VECTORATTRIBUTES_H
#define VECTORATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>

#include <ColorAttribute.h>

// ****************************************************************************
// Class: VectorAttributes
//
// Purpose:
//    Attributes for the vector plot
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class VectorAttributes : public AttributeSubject
{
public:
    enum Quality
    {
        Fast,
        High
    };
    enum OriginType
    {
        Head,
        Middle,
        Tail
    };
    enum LimitsMode
    {
        OriginalData,
        CurrentPlot
    };
    enum GlyphType
    {
        Arrow,
        Ellipsoid
    };
    enum LineStem
    {
        Cylinder,
        Line
    };
    enum GlyphLocation
    {
        AdaptsToMeshResolution,
        UniformInSpace
    };

    // These constructors are for objects of this class
    VectorAttributes();
    VectorAttributes(const VectorAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    VectorAttributes(private_tmfs_t tmfs);
    VectorAttributes(const VectorAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~VectorAttributes();

    virtual VectorAttributes& operator = (const VectorAttributes &obj);
    virtual bool operator == (const VectorAttributes &obj) const;
    virtual bool operator != (const VectorAttributes &obj) const;
private:
    void Init();
    void Copy(const VectorAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectVectorColor();
    void SelectColorTableName();

    // Property setting methods
    void SetGlyphLocation(GlyphLocation glyphLocation_);
    void SetUseStride(bool useStride_);
    void SetStride(int stride_);
    void SetNVectors(int nVectors_);
    void SetLineStyle(int lineStyle_);
    void SetLineWidth(int lineWidth_);
    void SetScale(double scale_);
    void SetScaleByMagnitude(bool scaleByMagnitude_);
    void SetAutoScale(bool autoScale_);
    void SetHeadSize(double headSize_);
    void SetHeadOn(bool headOn_);
    void SetColorByMag(bool colorByMag_);
    void SetUseLegend(bool useLegend_);
    void SetVectorColor(const ColorAttribute &vectorColor_);
    void SetColorTableName(const std::string &colorTableName_);
    void SetInvertColorTable(bool invertColorTable_);
    void SetVectorOrigin(OriginType vectorOrigin_);
    void SetMinFlag(bool minFlag_);
    void SetMaxFlag(bool maxFlag_);
    void SetLimitsMode(LimitsMode limitsMode_);
    void SetMin(double min_);
    void SetMax(double max_);
    void SetLineStem(LineStem lineStem_);
    void SetGeometryQuality(Quality geometryQuality_);
    void SetStemWidth(double stemWidth_);
    void SetOrigOnly(bool origOnly_);
    void SetGlyphType(GlyphType glyphType_);

    // Property getting methods
    GlyphLocation        GetGlyphLocation() const;
    bool                 GetUseStride() const;
    int                  GetStride() const;
    int                  GetNVectors() const;
    int                  GetLineStyle() const;
    int                  GetLineWidth() const;
    double               GetScale() const;
    bool                 GetScaleByMagnitude() const;
    bool                 GetAutoScale() const;
    double               GetHeadSize() const;
    bool                 GetHeadOn() const;
    bool                 GetColorByMag() const;
    bool                 GetUseLegend() const;
    const ColorAttribute &GetVectorColor() const;
          ColorAttribute &GetVectorColor();
    const std::string    &GetColorTableName() const;
          std::string    &GetColorTableName();
    bool                 GetInvertColorTable() const;
    OriginType           GetVectorOrigin() const;
    bool                 GetMinFlag() const;
    bool                 GetMaxFlag() const;
    LimitsMode           GetLimitsMode() const;
    double               GetMin() const;
    double               GetMax() const;
    LineStem             GetLineStem() const;
    Quality              GetGeometryQuality() const;
    double               GetStemWidth() const;
    bool                 GetOrigOnly() const;
    GlyphType            GetGlyphType() const;

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string Quality_ToString(Quality);
    static bool Quality_FromString(const std::string &, Quality &);
protected:
    static std::string Quality_ToString(int);
public:
    static std::string OriginType_ToString(OriginType);
    static bool OriginType_FromString(const std::string &, OriginType &);
protected:
    static std::string OriginType_ToString(int);
public:
    static std::string LimitsMode_ToString(LimitsMode);
    static bool LimitsMode_FromString(const std::string &, LimitsMode &);
protected:
    static std::string LimitsMode_ToString(int);
public:
    static std::string GlyphType_ToString(GlyphType);
    static bool GlyphType_FromString(const std::string &, GlyphType &);
protected:
    static std::string GlyphType_ToString(int);
public:
    static std::string LineStem_ToString(LineStem);
    static bool LineStem_FromString(const std::string &, LineStem &);
protected:
    static std::string LineStem_ToString(int);
public:
    static std::string GlyphLocation_ToString(GlyphLocation);
    static bool GlyphLocation_FromString(const std::string &, GlyphLocation &);
protected:
    static std::string GlyphLocation_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;

    // User-defined methods
    bool ChangesRequireRecalculation(const VectorAttributes &obj);

    // IDs that can be used to identify fields in case statements
    enum {
        ID_glyphLocation = 0,
        ID_useStride,
        ID_stride,
        ID_nVectors,
        ID_lineStyle,
        ID_lineWidth,
        ID_scale,
        ID_scaleByMagnitude,
        ID_autoScale,
        ID_headSize,
        ID_headOn,
        ID_colorByMag,
        ID_useLegend,
        ID_vectorColor,
        ID_colorTableName,
        ID_invertColorTable,
        ID_vectorOrigin,
        ID_minFlag,
        ID_maxFlag,
        ID_limitsMode,
        ID_min,
        ID_max,
        ID_lineStem,
        ID_geometryQuality,
        ID_stemWidth,
        ID_origOnly,
        ID_glyphType,
        ID__LAST
    };

private:
    int            glyphLocation;
    bool           useStride;
    int            stride;
    int            nVectors;
    int            lineStyle;
    int            lineWidth;
    double         scale;
    bool           scaleByMagnitude;
    bool           autoScale;
    double         headSize;
    bool           headOn;
    bool           colorByMag;
    bool           useLegend;
    ColorAttribute vectorColor;
    std::string    colorTableName;
    bool           invertColorTable;
    int            vectorOrigin;
    bool           minFlag;
    bool           maxFlag;
    int            limitsMode;
    double         min;
    double         max;
    int            lineStem;
    int            geometryQuality;
    double         stemWidth;
    bool           origOnly;
    int            glyphType;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define VECTORATTRIBUTES_TMFS "ibiiiidbbdbbbasbibbiddiidbi"

#endif
