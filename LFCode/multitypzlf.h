/* ----
   Project   LSST/BAO/PhotoZ
   Classes to generate galaxie Absolute Magnitude and types starting  
   from Luminosity Functions (LF's)  
   F. Habibi & R.Ansari ,  LAL (CNRS/IN2P3)  & Univ. Paris Sud 
   April - July 2017                                          -------  */
    
#ifndef MULTITYPZLF_H_SEEN
#define MULTITYPZLF_H_SEEN


#include <string>
#include <vector>

#include "slininterp.h"
#include "randinterf.h"
#include "myschechter.h"

using namespace std;
using namespace SOPHYA;

//--------------------------------------------------------------------------------
// Class to handle a set of LF's corresponding to multiple types and redshifts bins
//--------------------------------------------------------------------------------
class MultiType_Z_LF
{
public:
  /* ---- the following static methods define the default values for various parameters, 
     valid for subsequently created MultiType_Z_LF objects */
  // define the Hubble constant value (km/s/Mpc), which is used to correct MStar value of the Schechter parameter
  static void set_def_H0(double H0 = 70.);
  // define the normalised to 100 km/s/Mpc Hubble constant value (H0/(100 km/s/Mpc) 
  static void set_def_h100(double h100 = 0.7);
  // define the normalised to 70 km/s/Mpc Hubble constant value (H0/(70 km/s/Mpc) 
  static void set_def_h70(double h70 = 1.);
  // return default value of the Hubble constant in km/s/Mpc
  static double get_def_H0();
  // define the redshift at which number of galaxies reaches zero
  static void set_def_zero_gal_redshift(double z = 20.);
  // define the redshift range and number of bins for magnitude distribution for all galaxies
  static void set_def_redshift_range(double zmin=0., double zmax=4., double nz=20);
  // define the Gauss-Legendre quadrature order Schechter LF integration from magMin to magMax and for incremental integration
  static void set_def_schechter_integ_glorder(int glorder=100, int incglorder=5);
  // get the Gauss-Legendre quadrature order for Schechter LF integration 
  static int get_def_schechter_integ_glorder();
  // Activate / deactivate CPU timing in class methods / constructors 
  static void activateCPUTiming(bool fgtim=true);
  
  //-------  Class constructor and methods 
  /*! \brief constructor with specification of the name of the file containing LF parameters
    \param lfparamsfilename : parameters of the LF Schechter function for elliptical, spirals, starburst and all 
      for different redshifts ranges 
    \param magMin , magMax : magnitude range 
    \param nbmag : number of magnitude points for computation of interpolated magnitude as a function of 
      integrated LF function 
    \param  useall : if true, use LF Schechter parameters for All galaxies for computing galaxy count and 
      magnitude distribution for different redshifts - if false, use Sum of LF functions for Elliptical, 
      Spiral and StarBurst 
   */
  MultiType_Z_LF (string const & lfparamsfilename, double magMin, double magMax, size_t nbmag=25, bool useall=false);

  ~MultiType_Z_LF();

  //! return minimum magnitude limit
  double getMinMagLimit() const { return magMin_; }
  //! return maximum magnitude limit
  double getMaxMagLimit() const { return magMax_; } 

  //--- return the Schechter function for each type for redshift z
  MySchechter getLF_Elliptical(double z) const ;
  MySchechter getLF_Spiral(double z)  const ;
  MySchechter getLF_StarBurst(double z)  const ;
  MySchechter getLF_All(double z)  const ;
  
  //--- return the interpolated integral of the LF for each galaxy type  at redshift z , with magMin_<=mag<=magMax_
  inline double getIntegral_Elliptical(double z) const
  { return ILF_of_z_Ell(z); }
  inline double getIntegral_Spiral(double z) const
  { return ILF_of_z_Sp(z); }
  inline double getIntegral_StarBurst(double z) const
  { return ILF_of_z_SB(z); }
  inline double getIntegral_All(double z)  const 
  { return ILF_of_z_All(z); }

  /* --- return the average total galaxy number density at redshift z , with magMin_<=mag<=magMax_
     corresponding to the interpolated Sum[Ell+Spiral+StarBurst] or from All galaxies Schechter function parameters,
     depending on useSchechall_ flag defined in the constructor call */
  double getGalaxyNumberDensity(double z)  const ;

  /* --- return the average total galaxy number density at redshift z , with magMin_<=mag<=magMax_ (same as getGalaxyNumberDensity(double z) above)
     and return also the fraction of galaxies for each type ( Elliptical, Spiral & Starburst)
     corresponding to the interpolated Sum[Ell+Spiral+StarBurst] or from All galaxies Schechter function parameters,
     depending on useSchechall_ flag defined in the constructor call */
  double getGalaxyNumberDensity(double z, double& fEll, double& fSp, double& fSB)  const ;

  /* --- return a randomly drawn magnitude and broad type for galaxies at redshift z  
       magMin <= mag <= magMax  
       type = 1 -> Elliptical ,  2 -> Spiral  ,  3 -> Star Burst   */
  void  getTypeMagnitude(double z, double& mag, int& type)  const
  { return P_getTypeMagnitude(z, mag, true, type, true);  }

  // return the absolute magnitude for redshift z 
  double getMag(double z)  const
  {
    double mag;  int type;
    P_getTypeMagnitude(z, mag, true, type, false);
    return mag;
  }
  // return a broad type from a given redshift and magnitude  (uses type fraction as a function of mag)
  int getType(double z, double mag)   const
  {
    int type;
    P_getTypeMagnitude(z, mag, false, type, true);
    return type;
  }
  
  // define the random generator to be used
  void  setRandomGenerator(RandomGeneratorInterface* rg = NULL)
  {
    if (rg) rgp_ = rg;
    else rgp_ = def_rgp_;
  }
  // return the random generator being used
  inline RandomGeneratorInterface& getRandomGenerator()
  { return *rgp_; }
  
  // prints the Schechter parameters , if fgl=true
  ostream& print(ostream& os)  const;

protected:
  double magMin_, magMax_;    // limit in magnitude
  size_t nb_mag_pts_;     // number of magnitude points for computation of mag_f_ILF_All (magnitude as a function of Integrated LF) 
  bool  useSchechall_;  // if true, use Schechter parameters for All, false: Sum[Ell+Spiral+StarBusrt]
  
  vector< pair<double,double> > zrange;       // redshift range for each LF
  bool first_z_is_zero;  // true -> first zrange starts at zero
  
  // LF per type and for given redshift bins
  vector<MySchechter> vshEll;     // Schechter parameters for Elliptical galaxies for each redshift bin 
  vector<MySchechter> vshSp;      // Schechter parameters for Spiral galaxies for each redshift bin 
  vector<MySchechter> vshSB;      // Schechter parameters for Star Burst galaxies for each redshift bin
  vector<MySchechter> vshAll;     // Schechter parameters for All galaxies for each redshift bin
    
  // Integral of LF function as a function of redshift , for each type, and for 'All' galaxies 
  SLinInterp1D  ILF_of_z_Ell;
  SLinInterp1D  ILF_of_z_Sp;
  SLinInterp1D  ILF_of_z_SB;
  SLinInterp1D  ILF_of_z_All;

  // Integral of LF function as a function of redshift, as a sum of Elliptica,+Spiral+Starburst 
  SLinInterp1D  ILF_of_z_Sum_EllSpSB;

  // Inverse of the normalised integral of LF I(m) = Integral[LF_all]_magMin^m/ Integral[LF_all]_magMin^magMax
  // computed using Schechter function for all galaxies or Sum[Ell+Spiral+StarBurst] , depending on useSchechall_ flag 
  vector<SLinInterp1D *> v_mag_f_ILF;  // m(I) = inverse of I(m) function, with 0<=I<=1

  // redshift at which number of galaxies reaches zero 
  double zero_gal_redshift_;
  // redshift range and number of number of redshifts for which we compute the Integrale[LF](mag) for random drawing of magnitudes
  double zmin_, zmax_, dz_;
  int nz_;
  
  //  The random generator used for the generation of galaxy magnitudes and types. 
  RandomGeneratorInterface* rgp_;
  //  The default random generator 
  RandomGeneratorInterface* def_rgp_;
  
  // return the bin number where a given redshift z is locating in
  size_t find_redshift_bin(double z) const ;
  // private function to draw magnitude and type , used by getMag() , getType() and getTypeMagnitude() 
  void  P_getTypeMagnitude(double z, double& mag, bool fgmag, int& type, bool fgtype) const ; 
};

/*! operator << overloading -  */
inline ostream& operator << (ostream& s, MultiType_Z_LF const & mlf)
  {  mlf.print(s);  return(s);  }

#endif  /*  Fin de MULTITYPZLF_H_SEEN */
