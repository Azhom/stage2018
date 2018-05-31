/* ----
   Project   LSST/BAO/PhotoZ
   Classes to compute k-corrections 
   F. Habibi & R.Ansari ,  LAL (CNRS/IN2P3)  & Univ. Paris Sud 
   April - Sep 2017                                       -------  */
    
#ifndef KCORR_H_SEEN
#define KCORR_H_SEEN

#include "sedfilter.h"
#include "myschechter.h"
//--------------------------------------------------------------------------------
// Class to compute k-corrections 
//--------------------------------------------------------------------------------

class KCorrector {
public:
  // constructor with specification of the name of the file containing AB magnitude reference SED = Cte/lambda^2
  //  sedABzeroFilename : filename, lambda min et max en nm
  KCorrector (string & sedABzeroFilename, double lambdaMin=100, double lambdaMax=2500)
    : sedABzero_(sedABzeroFilename,lambdaMin*1e-9, lambdaMax*1e-9), lambdaMin_(lambdaMin), lambdaMax_(lambdaMax)
  { }
  KCorrector (SED & sedABzero, double lambdaMin=100, double lambdaMax=2500)
    : sedABzero_(sedABzero), lambdaMin_(lambdaMin), lambdaMax_(lambdaMax)
  { }
  
  // it computes magnitude difference My-Mx in two filters defined by filtx and filty in rest frame
  double computeK(SED & sed, Filter & filtx, Filter &filty);
  // it computes magnitude difference Kxy=my-(Mx+DM) in two filters defined by filtx and filty redshifted to observer frame
  // DM = 5log10(dL/10pc)
  // x: filter used in rest frame, y: filter used in observer frame
  // following the relation 13 Hogs et al. 2002
  double computeKz(SED & sed, Filter & filtx, Filter &filty, double z);

  //  Computes Rest-frame magnitude limit in Filter filt_RF , starting from an observer frame magnitude
  //  limit maglimObs in Filter filt_Obs, for an object with Rest-frame SED given by sed
  //  Explain sedABzero -
  //   z is the object redshift
  static double ComputeRestFrameMagnitudeLimit(Filter & filt_RF,Filter & filt_Obs, SED &sed, SED &sedABzero,
                                               double lambdamin,double lambdamax, double dL, double maglimObs, double z, Reddening& red);
  
protected:
  SED sedABzero_;
  double lambdaMin_, lambdaMax_; // lambda values in nm
  
};



//-------------------------------------------------------------------------
// Class to represent an efficiency function as a function of magnitude
class MagEfficiencyFunc : public ClassFunc1D  {
public:
    MagEfficiencyFunc(double magcenter, double deltamag)
    : magcenter_(magcenter), magwidth_(deltamag*sqrt(2.))
    { }
    virtual double operator()(double mag) const
    {
        return 0.5*erfc((mag-magcenter_)/magwidth_);
    }

protected:
    double magcenter_, magwidth_;
};

class Eff_Schechter : public ClassFunc1D  {
public:
    Eff_Schechter(MySchechter const & sch, MagEfficiencyFunc const& eff)
    : sch_(sch), eff_(eff)
    { }
    virtual double operator()(double mag) const
    {
        return sch_(mag)*eff_(mag);
    }
protected:
    MySchechter const & sch_;
    MagEfficiencyFunc const& eff_;
};

//  original version for ComputeRestFrameMagnitudeLimit() with debug prints
double ComputeBlim(Filter & filt_B,Filter & filt_r, SED &sed, SED &sedABzero,
                   double lambdamin,double lambdamax, double dL, double magr, double z);
// if magerr > 0, use efficiency function
double ComputeNg(MySchechter & esch, double magBlim, double Xi, double deltaXi, double magerr=-1.);

#endif
