/* ----
   Project   LSST/BAO/PhotoZ
   Classes to compute k-corrections 
   F. Habibi & R.Ansari ,  LAL (CNRS/IN2P3)  & Univ. Paris Sud 
   January 2018                                       -------  */
    
#ifndef GALCNT_H_SEEN
#define GALCNT_H_SEEN

#include <string>
#include <vector>

#include "kcorr.h"
#include "luc.h"
#include "multitypzlf.h"
#include "reddening.h"

using namespace std;
using namespace SOPHYA;

//--------------------------------------------------------------------------------
// Class to compute galaxy number density starting from lumonosity function,
// SED's and Filters , and a Cosmology
// This class maps the three LF-type galaxies (Elliptical , Spiral StarBurst)
// into 6 SED's :  El_cww (0)  Scd_cww (1)  Sbc_cww (2)  Im_cww (3)  SB2_kin (4)  SB3_kin (5)
//  Default mapping:
//    Elliptiacl -> 100% El_cww
//    Spiral     -> 40%  Scd_cww + 40% Sbc_cww +  20% Im_cww
//    StarBurst  -> 50%  SB2_kin + 50% SB3_kin 
//--------------------------------------------------------------------------------
 
class GalaxyCountComputer {
public:
  //! Constructor
  //  LF's manager object : mlf , LF restframe filter object filt_LF and observation filter filt_Obs
  //  Cosmology (distance(redshift) su , Path for SED's ,
  //  wavelength range lambdamin,lambdamax in nanometer 
  GalaxyCountComputer(MultiType_Z_LF & mlf, Filter& filt_LF, Filter& filt_Obs,
		      SimpleUniverse & su,  SFPathManager & sedpathmgr, double lambdaSEDmin=100., double lambdaSEDmax=2500.);

  //! set/change print level
  inline void setPrintLevel(int lev=0) { prtlev_ = lev; }
  //! Define fraction of Elliptical galaxies as defined in the corresponding LF to be mapped into each of the 6 SED types
  //  Default 100% El_cww (0) 
  void setEllipticalSEDFraction(std::vector<double> & sedfrac);
  //! Define fraction of Spiral galaxies as defined in the corresponding LF to be mapped into each of the 6 SED types
  //  Default 40%  Scd_cww (1) + 40% Sbc_cww (2) +  20% Im_cww (3) 
  void setSpiralSEDFraction(std::vector<double> & sedfrac);
  //! Define fraction of Spiral galaxies as defined in the corresponding LF to be mapped into each of the 6 SED types
  //  Default 50%  SB2_kin (4) + 50% SB3_kin (5) 
  void setStarBurstSEDFraction(std::vector<double> & sedfrac);
  
  void doCompute(Reddening& red, double zmin, double zmax, double dz, double mag_obs_lim, double mag_obs_err=0,
		 double lambdamin=300., double lambdamax=1000.);
  
  inline size_t getNbRedshiftBins() const { return redshifts_.size(); }
  inline double getMinRedshift() const { return zmin_; }
  inline double getMaxRedshift() const { return zmax_; }
  inline double getDeltaRedshift() const { return dz_; }
  inline double getMagObsLimit() const { return mag_obs_limit_; }
  double getRedshift(size_t) const; 

  //! return the total and per type (Ell,Sp,SB) galaxy volume densities (ngal/Mpc^3)  (return value=total)
  double getGalDensity_Mpc3(size_t i, double& Ellcnt, double & Spcnt, double & SBcnt, double & redshift) const;
  //! return the total and per type (Ell,Sp,SB) galaxy surface densities (ngal/arcmin^2  (return value=total)
  double getGalDensity_Arcmin2(size_t i, double& Ellcnt, double & Spcnt, double & SBcnt, double & redshift) const;
  //! return the integrated (zmin<z<zmax) total and per type (Ell,Sp,SB) galaxy volume density (ngal/Mpc^3) , (return value=total)
  double getAverageGalDensity_Mpc3(double& Ellcnt, double & Spcnt, double & SBcnt) const;
  //! return the integrated (zmin<z<zmax) total and per type (Ell,Sp,SB) galaxy surface density (ngal/arcmin^2) , (return value=total)
  double getIntegratedGalDensity_Arcmin2(double& Ellcnt, double & Spcnt, double & SBcnt) const;

protected:
  double IntegrateLF(MySchechter  & sch, double magBlim, double magerr);

  int prtlev_;                // Print level 
  MultiType_Z_LF & multLF_;   // Manager class for Luminosity functions (LF's ) 
  Filter& filt_LF_;           // Filter object associated with the luminosity functions
  Filter& filt_Obs_;          // Filter object associated with the luminosity functions 
  SimpleUniverse & su_;       // Cosmology   (distance vs redshit)
  SED sedABzero_;             // Spectral energy distribution corresponding to AB zero magnitude
  std::vector<SED> sedGals_;  // SED for different galaxy types
  double lambdaSEDmin_, lambdaSEDmax_;

  std::vector<double> Ell_sed_mapping_;   // Fractional mapping for Ellipticals from LF's to the six SED's 
  std::vector<double> Sp_sed_mapping_;   // Fractional mapping for Spirals from LF's to the six SED's 
  std::vector<double> SB_sed_mapping_;   // Fractional mapping for Starburst from LF's to the six SED's 

  double zmin_, zmax_, dz_;          // redshift range for which the galaxy number density is computed
  double mag_obs_limit_, mag_obs_error_;
  double lambdamin_, lambdamax_;

  std::vector<double> redshifts_;            // redshifts for which galaxy number densities are computed 
  //----- galaxy number density , ngal / Mpc3
  std::vector<double> Ell_cnt_per_Mpc3_;       // Elliptical galaxies : ngal / Mpc3
  std::vector<double> Sp_cnt_per_Mpc3_;        // Spiral galaxies : ngal / Mpc3
  std::vector<double> SB_cnt_per_Mpc3_;        // StarBurst galaxies : ngal / Mpc3
  std::vector<double> volume_arcmin2_;         // cosmological volume / arcmin^2 / dz
};

#endif
