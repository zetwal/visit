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

#ifndef ISOSURFACEATTRIBUTES_H
#define ISOSURFACEATTRIBUTES_H
#include <string>
#include <AttributeSubject.h>


// ****************************************************************************
// Class: IsosurfaceAttributes
//
// Purpose:
//    Attributes for the isosurface operator
//
// Notes:      Autogenerated by xml2atts.
//
// Programmer: xml2atts
// Creation:   omitted
//
// Modifications:
//   
// ****************************************************************************

class IsosurfaceAttributes : public AttributeSubject
{
public:
    enum Select_by
    {
        Level,
        Value,
        Percent
    };
    enum Scaling
    {
        Linear,
        Log
    };

    // These constructors are for objects of this class
    IsosurfaceAttributes();
    IsosurfaceAttributes(const IsosurfaceAttributes &obj);
protected:
    // These constructors are for objects derived from this class
    IsosurfaceAttributes(private_tmfs_t tmfs);
    IsosurfaceAttributes(const IsosurfaceAttributes &obj, private_tmfs_t tmfs);
public:
    virtual ~IsosurfaceAttributes();

    virtual IsosurfaceAttributes& operator = (const IsosurfaceAttributes &obj);
    virtual bool operator == (const IsosurfaceAttributes &obj) const;
    virtual bool operator != (const IsosurfaceAttributes &obj) const;
private:
    void Init();
    void Copy(const IsosurfaceAttributes &obj);
public:

    virtual const std::string TypeName() const;
    virtual bool CopyAttributes(const AttributeGroup *);
    virtual AttributeSubject *CreateCompatible(const std::string &) const;
    virtual AttributeSubject *NewInstance(bool) const;

    // Property selection methods
    virtual void SelectAll();
    void SelectContourValue();
    void SelectContourPercent();
    void SelectVariable();

    // Property setting methods
    void SetContourNLevels(int contourNLevels_);
    void SetContourValue(const doubleVector &contourValue_);
    void SetContourPercent(const doubleVector &contourPercent_);
    void SetContourMethod(Select_by contourMethod_);
    void SetMinFlag(bool minFlag_);
    void SetMin(double min_);
    void SetMaxFlag(bool maxFlag_);
    void SetMax(double max_);
    void SetScaling(Scaling scaling_);
    void SetVariable(const std::string &variable_);

    // Property getting methods
    int                GetContourNLevels() const;
    const doubleVector &GetContourValue() const;
          doubleVector &GetContourValue();
    const doubleVector &GetContourPercent() const;
          doubleVector &GetContourPercent();
    Select_by          GetContourMethod() const;
    bool               GetMinFlag() const;
    double             GetMin() const;
    bool               GetMaxFlag() const;
    double             GetMax() const;
    Scaling            GetScaling() const;
    const std::string  &GetVariable() const;
          std::string  &GetVariable();

    // Persistence methods
    virtual bool CreateNode(DataNode *node, bool completeSave, bool forceAdd);
    virtual void SetFromNode(DataNode *node);

    // Enum conversion functions
    static std::string Select_by_ToString(Select_by);
    static bool Select_by_FromString(const std::string &, Select_by &);
protected:
    static std::string Select_by_ToString(int);
public:
    static std::string Scaling_ToString(Scaling);
    static bool Scaling_FromString(const std::string &, Scaling &);
protected:
    static std::string Scaling_ToString(int);
public:

    // Keyframing methods
    virtual std::string               GetFieldName(int index) const;
    virtual AttributeGroup::FieldType GetFieldType(int index) const;
    virtual std::string               GetFieldTypeName(int index) const;
    virtual bool                      FieldsEqual(int index, const AttributeGroup *rhs) const;


    // IDs that can be used to identify fields in case statements
    enum {
        ID_contourNLevels = 0,
        ID_contourValue,
        ID_contourPercent,
        ID_contourMethod,
        ID_minFlag,
        ID_min,
        ID_maxFlag,
        ID_max,
        ID_scaling,
        ID_variable,
        ID__LAST
    };

private:
    int          contourNLevels;
    doubleVector contourValue;
    doubleVector contourPercent;
    int          contourMethod;
    bool         minFlag;
    double       min;
    bool         maxFlag;
    double       max;
    int          scaling;
    std::string  variable;

    // Static class format string for type map.
    static const char *TypeMapFormatString;
    static const private_tmfs_t TmfsStruct;
};
#define ISOSURFACEATTRIBUTES_TMFS "id*d*ibdbdis"

#endif