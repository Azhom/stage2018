/**
 * @file  sedfilter.h
 * @brief Contains classes that perform calculations with SED and filter funcs
 *
 * @todo <CODE>SED::readSED</CODE> <BR> 
 *       Check if the interpolated flux value should be zero for wavelengths
 *       outside of those read in 
 *
 * @author Alex Abate and Reza Ansari
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2008
 * @date 2008
 * 
 *  Extracted from the original file, removed classes modifying SED due to absorption
 *  (dust) by Farhang Habibi (and Reza Ansari) - LAL , 2017
 */
 
#ifndef SEDFILTER_H_SEEN 
#define SEDFILTER_H_SEEN 

#include "machdefs.h"
#include "sopnamsp.h"
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>

#include "classfunc.h"
#include "tvector.h"

// FH-DEL #include "reddening.h"
// #include "slininterp.h" // FH-CHANGE "sinterp.h"
#include "sinterp.h"

#include "reddening.h"

// FH-DEL #include "igm.h"


namespace SOPHYA {

/*******************************************************************************
*                                                                              *
*                                SED CLASSES                                   *
*                                                                              *
*******************************************************************************/

class SFPathManager {
public:
    SFPathManager(string const & path) : path_(path)
    {
        size_t l = path_.length();
        if ((l>0) && (path_[l-1]!='/'))  path_ += '/';
    }
    inline string BuildFileName(string const & fname)
    {
        return path_+fname;
    }
    string path_;
};

/** @class
  * Spectral Energy Distribution (SED) class
  * 
  * SEDs are read in from a text file with 2 columns wavelength in meters and
  * flux.  The flux value at the given wavelength is returned by the class.
  * The SED can be optionally interpolated or reddened
  * 
  */
class SED : public ClassFunc1D
{
public:

    /** Default constructor */
  SED() : sed2_(NULL), sed2init_(NULL)    //  CHECK Pointers initialized to zero - Is this correct ?  Reza  April 2016 
    { isRead_=false; isRedden_=false; isInterp_=false; _test_=1.; };
  
    /** Constructor
        To read in wavelength,flux values from a file, and return flux values
        based on the linear interpolation of numbers in that file */
    SED(string& filename, double xmin, double xmax, int npt=1024);
    
    /** Copy constructor */
    SED(const SED& sed);
    virtual SED& Set(const SED& a);
      
    /** Read in SED from file 
        @param fname    name of file to read from */
    void readSED(string& fname, double xmin, double xmax, int npt=1024);
      
    /** @warning code doesn't seem to like this! */
    virtual double operator()(double x) const
            { return returnFlux(x); };
            
    /** Set up interpolation of the SED */
    void doInterp(SED* sed2,double a=1.,double b=0.);
      
    /** Set up reddening of the SED */
    void doRedden(Reddening& red, double z=0.);
    void stopRedden();
    void print_redden();
     
    /** To return the relevant flux value
        What is returned depends on isRedden and isInterp values
        Either: interpSED(), addReddening() or interpAddReddening() is called 
        @param lambda   rest-frame wavelength                                 */
    double returnFlux(double lambda) const;
    
    /** To return the reddened value 
        @param lambda   rest-frame wavelength                                 */
    double addReddening(double lambda) const;
    
    /** To return the interpolated value
        @param lambda   rest-frame wavelength                                 */
    double interpSED(double lambda) const;
    
    /** To return the reddened and interpolated value
        @param lambda   rest-frame wavelength                                 */
    double interpAddReddening(double lambda) const;
    
    /** return bool holding interpolation information*/
    bool returnInterpBool(){ return isInterp_; };
    
    /** return bool holding reddening information */
    bool returnReddenBool(){ return isRedden_; };
    
    /** return bool holding file-read information */
    bool returnReadBool(){ return isRead_; };
    
    
    double returntest(){ return _test_; };
    
protected:
    static string  sed_path_;   /**< pathname for SED definition text files   */
    SInterp1D   sed_;           /**< The SED read in and to be used           */
    bool        isRead_;        /**< true if file has been read               */
    bool        isRedden_;      /**< true if SED has been reddened            */
    bool        isInterp_;      /**< true if SED has been interpolated        */
    SED*        sed2_;          /**< 2nd SED class to interpolate between with*/
    SED*        sed2init_;      /**< a default setting of sed2                */
    double      a_,b_;          /**< how to interpolate between sed and sed2  */
    double      RvCard_,z_;
    double      _test_;
    Reddening	red_;
};


/** @class
  * SED interpolation class
  * 
  * Can linearly interpolate between 2 SED objects as specifed with the 
  * arguments. For example if you want to interpolate and make a new SED where
  * 10% of it is "sed1" and 90% of it is "sed2" set "mult1"=0.1 and "mult2"=0.9
  * 
  * @note This class is possibly redundant now this functionality is 
  * with SED class itself
  * 
  */
class SEDInterp : public ClassFunc1D
{
public:
    
    SEDInterp(SED& sed1, SED& sed2, double mult1, double mult2)
    : sed1_(sed1), sed2_(sed2), mult1_(mult1), mult2_(mult2) {
        double eps=1e-6;
        if ( std::abs(mult1+mult2-1.)>eps)
            throw ParmError("ERROR! SED factors must sum to 1");
    };
     
    /** returns interpolated SED 
        @param lambda   rest-frame wavelength                                 */
    virtual double operator()(double lambda) const
            {  return (mult1_*sed1_(lambda)+mult2_*sed2_(lambda));  };

protected:
    SED& sed1_;      /**< SED 1                                               */
    SED& sed2_;      /**< SED 2                                               */
    double mult1_;   /**< fraction of SED 1 in new SED                        */
    double mult2_;   /**< fraction of SED 2 in new SED                        */
};


/** @class SEDzFilterProd
  * Multiply SED and Filter together to create an integrand for the k-correction
  * 
  * returns: redshifted SED multiplied by the filter transmission and wavelength
  *          SED(lambda/(1+z))*Filter(lambda)*lambda
  *
  * BECAUSE for AB magnitudes, and SED in units of wavelength:
  *
  * mAB(z) = -2.5log10( int SED(lambda/(1+z))*Filter(lambda)*lambda dlambda /  
  *                                 int Filter(lambda)*lambda^-1 dlambda ) -51.6
  *
  * and when computing Kcorrection the int Filter(lambda)*lambda^-1 dlambda )-51.6
  * cancel because:
  * K12 = mAB(z) - mAB(z=0)
  */
class SEDzFilterProd : public ClassFunc1D
{
public:
    /** Constructor 
        @param f    SED
        @param g    filter
        @param z    redshift of galaxy SED                                           */
    SEDzFilterProd(ClassFunc1D& f, ClassFunc1D& g, double z=0.)
    : sed_(f) , filt_(g) , z_(z) {  };
    
    /** Returns redshifted SED multiplied by the filter transmission and wavelength,
        \f$ SED(\lambda/(1+z))*Filter(\lambda)*\lambda \f$ where \lambda is the 
        observed wavelength    
        @param lambda   observed wavelength                                   */
    virtual double operator()(double lambda) const {
        double lambdaE = lambda/(1+z_);
        double fac = 0.;
        //if (lambdaE<2.e-6)
            fac = sed_(lambdaE); //sed_(lambda*(1+z_));
        double xx = fac*filt_(lambda)*lambda/(1+z_); // we do (1+z_ division outside)

        return (xx);
        }
// *lambda if SED in units of wavelength and magnitudes have dnu/nu definition
// no lambda factor if SED in units of wavelength and magnitudes have dnu definition
// /lambda if SED in units of frequency and magnitudes have dnu/nu definition
// /lambda^2 if SED in units of frequency and magnitudes have dnu definition

protected:
  ClassFunc1D& sed_;          /**< SED                                        */
  ClassFunc1D& filt_;         /**< filter function                            */
  double z_;                  /**< redshift of galaxy SED                     */
};

    
// return the product of SED and filter myltiplied by lambda (photon count)  at rest frame
class SEDFilterProd : public ClassFunc1D
{
public:
        /** Constructor
         @param f    SED
         @param g    filter
         */

    SEDFilterProd(ClassFunc1D& f, ClassFunc1D& g)
        : sed_(f) , filt_(g) {  };
        
        /** Returns redshifted SED multiplied by the filter transmission and wavelength,
         \f$ SED(\lambda/(1+z))*Filter(\lambda)*\lambda \f$ where \lambda is the
         observed wavelength
         @param lambda   observed wavelength                                   */
        virtual double operator()(double lambda) const {
            return (sed_(lambda)*filt_(lambda)*lambda);
        }
    
    protected:
        ClassFunc1D& sed_;          /**< SED                                        */
        ClassFunc1D& filt_;         /**< filter function                            */
    };
    
    /** ReadSedList class
     *
     * Class to read in SEDs from a list of filenames
     *
     */
    class ReadSedList {
    public:
        /** Constructor, finds the file sedFile and counts number of SEDs inside
         @param sedFile  filename containing list of SEDs
         @param prt      print level, if prt>0 extra statements are printed */
        ReadSedList(string sedFile, int prt=0);
        
        /** Read environment variable $SEDLOC */
        string getSedDirEnviromentVar();
        
        /** Counts SEDs in file */
        void countSeds(string sedFileFullPath);
        
        /** Main program, reads SEDs from file and arranges them in a vector of
         pointers pointing to each SED object
         @param lmin minimum wavelength of SED in meters
         @param lmax maximum wavelength of SED in meters */
        void readSeds(double lmin=5e-8, double lmax=2.5e-6);
        
        /** If interpolating between the SEDs call this method straight after
         readSeds()
         @param nInterp number of SEDs to interpolate between each SED */
        void interpSeds(int nInterp);
        
        /** If reddening all the SEDs call this method straight after #readSeds(),
         and after #interpSeds() if interpolating.
         Currently only applies Cardelli law, and can't limit max reddening for
         elliptical galaxies
         @param nStepRed number of times to redden
         @param redMax   maximum limit to redden to (even steps between 0 and redMax */
        void reddenSeds(int nStepRed, double redMax);
        
        /** Write contents of sedArray to a file */
        void writeSpectra(string outFile, double lmin=5e-8, double lmax=2.5e-6,
                          int nl=1500);
        
        /** Return sedArray */
        vector<SED*> getSedArray()
        { return sedArray_; };
        
        /** Return total number of SEDs (could be >= to nsed depending if interpolation
         was done or reddening was applied) */
        int getNTot()
        { return ntot_; };
        
        /** Return number of SEDs read in from the initial file */
        int getNSed()
        { return nsed_; };
        
        /** Reorder SEDs so that interpolated SEDs lie between the template they
         were interpolated from */
        void reorderSEDs();
        
        /** Return a vector of the SED filenames */
        vector<string> returnSedFilenames(){ return sedFiles_; };
        
    private:
        int prt_;                   /**<  print level                             */
        string sedDir_;             /**<  path of SED list file                   */
        string sedFileFullPath_;    /**<  full path and filename to SED list      */
        int nsed_;                  /**<  number of SEDs in the file sedFile      */
        int ntot_;                  /**<  total number of SEDs after interpolation /
                                     and reddening                           */
        vector<SED*> sedArray_;     /**<  pointer to ntot_ SED objects            */
        vector<string> sedFiles_;   /**<  vector of SED file names                */
        bool isInterp_;             /**<  if SEDs are interpolated                */
        bool isRedden_;             /**<  if SEDs are reddened                    */
    };
    
    

/*******************************************************************************
*                                                                              *
*                              FILTER CLASSES                                  *
*                                                                              *
*******************************************************************************/

/** @class Filter class 
  * To load in Filter transmission functions:
  * two columns: observed wavelength (m) ; transmission probability (between 0 and 1) 
  */
class Filter : public SInterp1D
{
public:
    /** Default constructor */
    Filter(){};
    
    /** Constructor 
        @param filename     file to read transmission function from
        @param lmin         min observed wavelength value
        @param lmax         max observed wavelength value
        @param nComments    number of comment lines the filter file has at start
        @param zero_outside value of interpolation outside lmin,lmax is zero if true
        @param npt          number of points to do interpolations from        */
    Filter(string& filename, double lmin, double lmax, int nComments=0, 
                                    bool zero_outside = true, int npt = 1024);
                                      
   /** Read in filter from file 
       @param filename     file to read transmission function from
       @param lmin         min observed wavelength value
       @param lmax         max observed wavelength value
       @param nComments    number of comment lines the filter file has at start
       @param zero_outside value of interpolation outside lmin,lmax is zero if true
       @param npt          number of points to do interpolations from         */
    void readFilter(string& filename, double lmin, double lmax, int nComments=0, 
                                bool zero_outside = true, int npt=1024);
};


/** @class BlueShiftFilter
  *
  */
class BlueShiftFilter : public ClassFunc1D
{
public:
    
    /** Constructor */
    BlueShiftFilter(Filter& g, double z=0.)
    : filt_(g) , z_(z) { }
    
    /** returns the value of the filter transmission at the rest-frame 
        wavelength of the object */
    virtual double operator()(double lambda) const {
        double lambdaRF = lambda*(1+z_);
        return (filt_(lambdaRF));
        };

protected:
    Filter& filt_;    /**< class holding the filter function */
    double z_;        /**< redshift of the object */
};


/** @class FilterProd class
  * To multiply filter transmission by 1/lambda
  *
  * returns: Filter(lambda)*lambda^-1
  * BECAUSE for AB magnitudes, and SED in units of wavelength:
  * mAB(z) = -2.5log10( int SED(lambda/(1+z))*Filter(lambda)*lambda dlambda /  
  *                                    int Filter(lambda)*lambda^-1 dlambda )-51.6
  *
  * and when computing Kcorrection the int Filter(lambda)*lambda^-2 dlambda )-51.6
  * cancels because:
  *
  * K12 = mAB(z) - mAB(z=0)
  *
  * spectral flux density: W/m^2/Hz int across bandwidth but still want flux to 
  * be in units of W/m^2/Hz because magnitude in band X is defined flux density 
  * at a particular wavelength
  */
class FilterProd : public ClassFunc1D
{
public:
    
    /** Constructor */
    FilterProd(Filter& g)
    : filt_(g){ };
    
    /** returns the value of the filter transmission at #lambda, divided by 
        #lambda */
    virtual double operator()(double lambda) const {
        return (filt_(lambda)/lambda);//********* DOUBLE CHECK THIS *********//
        }  

protected:
  Filter& filt_;    /**< class holding the filter function */
};


/** @class FilterXLambda class
  * To multiply filter transmission value by the wavelength
  *
  */
class FilterXLambda : public ClassFunc1D
{
public:
    
    /** Constructor */
    FilterXLambda(ClassFunc1D& g)
    : filt_(g){ }
    
    /** returns the value of the filter transmission at #lambda, multiplied by 
        #lambda */
    virtual double operator()(double lambda) const {
        return (filt_(lambda)*lambda);//********* DOUBLE CHECK THIS *********//
        }  

protected:
  ClassFunc1D& filt_;          /**< class holding the filter function         */
};


/** @class FilterProdProd class
  * To multiply filter with 1/lambda^2
  *
  */
class FilterProdProd : public ClassFunc1D
{
public:
    FilterProdProd(Filter& g)
    : filt_(g) {};
    
    virtual double operator()(double lambda) {
        return (filt_(lambda)/(lambda*lambda));
        }  

protected:
    Filter& filt_; 
};


//--- Does int F(lambda)/lambda dlambda or SED(lambda)*F(lambda) or whatever
// Simple summing integration
class FilterIntegrator // : public ClassFunc1D
{
public:
  FilterIntegrator(ClassFunc1D& f, double xmin, double xmax, int Nstep=500) 
    : f_(f) , xmin_(xmin) , xmax_(xmax) , Nstep_(Nstep) { }
    
   /** This is defined to override the pure virtual function defined in ClassFunc1D
        otherwise FilterIntegrator is sometimes treated as an abstract class        
    virtual double operator() (double) const { };
   */
  
   double Value() const {
        if (xmin_>=xmax_){
            string emsg="FilterIntegrator::Value() integral limits don't make sense ";
            throw out_of_range(emsg);
            }
                    
            if (Nstep_ <= 0)
                Nstep_ = 500;

            double dx = (xmax_ - xmin_) / (Nstep_-1);
            double x = xmin_;
            double sum=0;
            for (int i=1; i<Nstep_; i++, x += dx)
                sum += f_(x);
            
            return sum * dx;
        }

protected:
  ClassFunc1D& f_; 
  double xmin_, xmax_;
  mutable int Nstep_;
};


/** ReadFilterList class
  * 
  * Class to read in filters from a list of filenames
  *
  */
class ReadFilterList {
public:
    /** Constructor, finds the file sedFile and counts number of filters inside 
        @param sedFile  filename containing list of filters 
        @param prt      print level, if prt>0 extra statements are printed */
    ReadFilterList(string filterFile, int prt=0);
    
    /** Read environment variable $FILTLOC */
    string getFilterDirEnviromentVar();
    
    /** Counts Filters in file */
    void countFilters(string filterFileFullPath);
    
    /** Main program, reads filters from file and arranges them in a vector of
        pointers pointing to each Filter object 
        @param lmin minimum wavelength of filter in meters
        @param lmax maximum wavelength of filter in meters */
    void readFilters(double lmin=5e-8, double lmax=2.5e-6);
    
    /** Write contents of filterArray to a file */
    void writeFilters(string outFile, double lmin=5e-8, double lmax=2.5e-6,
                                                                 int nl=1500);
                                                                 
    /** Return filterArray */
    vector<Filter*> getFilterArray()
                { return filterArray_; };
    
    /** Return total number of filters */
    int getNTot() { return ntot_; };
                
private:
    int prt_;                       /**< print level */
    string filterDir_;              /**< path of filter list file */
    string filterFileFullPath_;     /**< full path and filename to filter list */
    int ntot_;                      /**< total number of filters */
    vector<Filter*> filterArray_;   /**< pointer to ntot_ Filter objects */
};



      
}// end namespace sophya

#endif
