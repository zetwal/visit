// ************************************************************************* //
//                            avtDuneFileFormat.h                           //
// ************************************************************************* //

#ifndef AVT_Dune_FILE_FORMAT_H
#define AVT_Dune_FILE_FORMAT_H

#include <database_exports.h>

#include <avtMTSDFileFormat.h>

#include <vector>
#include <map>
#include <visitstream.h>

// ****************************************************************************
//  Class: avtDuneFileFormat
//
//  Purpose:
//      Reads in Dune files as a plugin to VisIt.
//
//  Programmer: dslone -- generated by xml2info
//  Creation:   Thu Mar 11 09:14:32 PDT 2004
//
//    Mark C. Miller, Tue May 17 18:48:38 PDT 2005
//    Added timeState arg to PopulateDatabaseMetaData to satisfy new interface
// ****************************************************************************

using std::string;
using std::vector;
using std::map;

class avtDuneFileFormat : public avtMTSDFileFormat
{
  public:
                       avtDuneFileFormat(const char *);
    virtual           ~avtDuneFileFormat() {;};

    //
    // If you know the times and cycle numbers, overload this function.
    // Otherwise, VisIt will make up some reasonable ones for you.
    //
    virtual void        GetCycles(vector<int> &);
    virtual void        GetTimes(vector<double> &);
    //

    virtual int            GetNTimesteps(void);

    virtual const char    *GetType(void)   { return "Dune"; };
    virtual void           FreeUpResources(void); 

    virtual vtkDataSet    *GetMesh(int, const char *);
    virtual vtkDataArray  *GetVar(int, const char *);
    virtual vtkDataArray  *GetVectorVar(int, const char *);

  protected:
    // DATA MEMBERS

    string                  fname;
    int                     ntimes;
    int                     nparticles;
    vector<int>             species;
    vector<double>          radius;
    vector<double>          impulseTime;
    vector<double>          coordinates;
    vector<double>          velocities;
    vector<double>          impulseVelocities;
    vector<double>          totalVelocities;
    vector<double>          angularVelocities;

    virtual void           PopulateDatabaseMetaData(avtDatabaseMetaData *, int);

  private:
    enum fileTypes {
        UNKNOWN,
        LOADER,
        RESTART,
        TECPLOT
    };

    enum fileTypes          ftype;

    struct matInfoStruct {
        double massDensity;
        string name;
    };


    map<string, struct matInfoStruct> matInfo;

    int                      lastTimestate;
    ifstream                 ifile;
    vector<streampos>        fpos;
    vector<double>           times;
    vector<int>              cycles;
    vector<int>              num_particles;
    vector<double>           mass;
    vector<string>           matnames;

    void ReadDuneData(const int);

    inline string string_substr(const string, const string, 
                                const string, const bool);
    inline double fortranDoubleToCDouble(const string);
};


#endif
