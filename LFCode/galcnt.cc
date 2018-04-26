/* ----
   Project   LSST/BAO/PhotoZ
   Classes to compute k-corrections 
   F. Habibi & R.Ansari ,  LAL (CNRS/IN2P3)  & Univ. Paris Sud 
   January 2018                                       -------  */

#include "galcnt.h"
#include "integ.h"

//--- Number of galaxy SED types 
#define  GALCNTNBSED  6
//--------------------------------------------------------------------------------
// Class to compute galaxy number density starting from lumonosity function,
// SED's and Filters , and a Cosmology
//--------------------------------------------------------------------------------
GalaxyCountComputer::GalaxyCountComputer(MultiType_Z_LF & mlf, Filter& filt_LF, Filter& filt_Obs,
					 SimpleUniverse & su,  SFPathManager & sedpathmgr,
					 double lambdaSEDmin, double lambdaSEDmax)
  : prtlev_(0), multLF_(mlf), filt_LF_(filt_LF), filt_Obs_(filt_Obs), su_(su),
    lambdaSEDmin_(lambdaSEDmin), lambdaSEDmax_(lambdaSEDmax),
    Ell_sed_mapping_(GALCNTNBSED), Sp_sed_mapping_(GALCNTNBSED), SB_sed_mapping_(GALCNTNBSED)
{
  cout << "GalaxyCountComputer()/Info Reading 6 galaxy SEDs from PATH= "<<sedpathmgr.path_<<" ..."<<endl;
  //string sedfileEl = sedpathmgr.BuildFileName("flat.txt");
  //string sedfileSp1 = sedpathmgr.BuildFileName("flat.txt");
  //string sedfileSp2 = sedpathmgr.BuildFileName("flat.txt");
  //string sedfileSp3 = sedpathmgr.BuildFileName("flat.txt");
  //string sedfileSB1 = sedpathmgr.BuildFileName("flat.txt");
  //string sedfileSB2 = sedpathmgr.BuildFileName("flat.txt");
  
  string sedfileEl = sedpathmgr.BuildFileName("function_El_cww_fix2_extend_smooth.txt");
  string sedfileSp1 = sedpathmgr.BuildFileName("function_Scd_cww_fix2_extend_smooth.txt");
  string sedfileSp2 = sedpathmgr.BuildFileName("function_Sbc_cww_fix2_extend_smooth.txt");
  string sedfileSp3 = sedpathmgr.BuildFileName("function_Im_cww_fix2_extend_smooth.txt");
  string sedfileSB1 = sedpathmgr.BuildFileName("function_SB2_kin_fix2_extend_smooth.txt");
  string sedfileSB2 = sedpathmgr.BuildFileName("function_SB3_kin_fix2_extend_smooth.txt");
  //cout<<"* Reading SED file for Ellipticals: " << sedfileEl <<endl;
  SED sedEl(sedfileEl,lambdaSEDmin*1e-9, lambdaSEDmax*1e-9);
  sedGals_.push_back(sedEl);
  //cout<<"* Reading SED file for Spirals: " << sedfileSp <<endl;
  //SED sedSp(sedfileSp,lambdaSEDmin*1e-9, lambdaSEDmax*1e-9);
  SED sedSp1(sedfileSp1,lambdaSEDmin*1e-9, lambdaSEDmax*1e-9);
  SED sedSp2(sedfileSp2,lambdaSEDmin*1e-9, lambdaSEDmax*1e-9);
  SED sedSp3(sedfileSp3,lambdaSEDmin*1e-9, lambdaSEDmax*1e-9);
  sedGals_.push_back(sedSp1);
  sedGals_.push_back(sedSp2);
  sedGals_.push_back(sedSp3);
  //cout<<"* Reading SED file for StarBursts: " << sedfileSB <<endl;
  //SED sedSB(sedfileSB,lambdaSEDmin*1e-9, lambdaSEDmax*1e-9);
  SED sedSB1(sedfileSB1,lambdaSEDmin*1e-9, lambdaSEDmax*1e-9);
  SED sedSB2(sedfileSB2,lambdaSEDmin*1e-9, lambdaSEDmax*1e-9);
  sedGals_.push_back(sedSB1);
  sedGals_.push_back(sedSB2);
  string sedfileABzero = sedpathmgr.BuildFileName("ABzero.txt");
  cout << "GalaxyCountComputer()/Info reading SED file for AB magnitude zero point: " << sedfileABzero <<endl;
  sedABzero_.readSED(sedfileABzero,lambdaSEDmin*1e-9, lambdaSEDmax*1e-9);
  //----- defining default Galaxy to SED mapping fractions
  double fracEll[GALCNTNBSED]={1.,0.,0.,0.,0.,0.};
  double fracSp[GALCNTNBSED]={0.,0.4,0.4,0.2,0.,0.};
  double fracSB[GALCNTNBSED]={0.,0.,0.,0.,0.5,0.5};
  for(size_t i=0; i<GALCNTNBSED; i++) {
    Ell_sed_mapping_[i]=fracEll[i];  Sp_sed_mapping_[i]=fracSp[i];   SB_sed_mapping_[i]=fracSB[i];
  }
}

void  GalaxyCountComputer::setEllipticalSEDFraction(std::vector<double> & sedfrac)
{
  if (sedfrac.size() != Ell_sed_mapping_.size())
    throw SzMismatchError("GalaxyCountComputer::setEllipticalSEDFraction() BAD sedfrac.size() ");
  double sum=0.;
  for(size_t i=0; i<sedfrac.size(); i++) {
    if (sedfrac[i]<0.)  throw  ParmError("GalaxyCountComputer::setEllipticalSEDFraction()  Negative fraction !");
    sum += sedfrac[i];
  }
  if (fabs(sum-1.)>1.e-19)
    throw  ParmError("GalaxyCountComputer::setEllipticalSEDFraction()  Sum[Fracs] != 1 !");
  for(size_t i=0; i<sedfrac.size(); i++)  Ell_sed_mapping_[i]=sedfrac[i];
  return;
}

void  GalaxyCountComputer::setSpiralSEDFraction(std::vector<double> & sedfrac)
{
  if (sedfrac.size() != Ell_sed_mapping_.size())
    throw SzMismatchError("GalaxyCountComputer::setSpiralSEDFraction() BAD sedfrac.size() ");
  double sum=0.;
  for(size_t i=0; i<sedfrac.size(); i++) {
    if (sedfrac[i]<0.)  throw  ParmError("GalaxyCountComputer::setSpiralSEDFraction()  Negative fraction !");
    sum += sedfrac[i];
  }
  if (fabs(sum-1.)>1.e-19)
    throw  ParmError("GalaxyCountComputer::setSpiralSEDFraction()  Sum[Fracs] != 1 !");
  for(size_t i=0; i<sedfrac.size(); i++)  Sp_sed_mapping_[i]=sedfrac[i];
  return;
}

void  GalaxyCountComputer::setStarBurstSEDFraction(std::vector<double> & sedfrac)
{
  if (sedfrac.size() != Ell_sed_mapping_.size())
    throw SzMismatchError("GalaxyCountComputer::setStarBurstSEDFraction() BAD sedfrac.size() ");
  double sum=0.;
  for(size_t i=0; i<sedfrac.size(); i++) {
    if (sedfrac[i]<0.)  throw  ParmError("GalaxyCountComputer::setStarBurstSEDFraction()  Negative fraction !");
    sum += sedfrac[i];
  }
  if (fabs(sum-1.)>1.e-19)
    throw  ParmError("GalaxyCountComputer::setStarBurstSEDFraction()  Sum[Fracs] != 1 !");
  for(size_t i=0; i<sedfrac.size(); i++)  SB_sed_mapping_[i]=sedfrac[i];
  return;
}


void GalaxyCountComputer::doCompute(double zmin, double zmax, double dz, double mag_obs_lim, double mag_obs_err, 
				    double lambdamin, double lambdamax, double ebmv)
{
  zmin_=zmin;   zmax_=zmax;   dz_=dz;
  mag_obs_limit_ = mag_obs_lim;  mag_obs_error_=mag_obs_err;
  
  redshifts_.clear();
  Ell_cnt_per_Mpc3_.clear();
  Sp_cnt_per_Mpc3_.clear();
  SB_cnt_per_Mpc3_.clear();
  volume_arcmin2_.clear();
  lambdamin_=lambdamin;   lambdamax_=lambdamax;


  size_t zcnt = 0;
  for(double z=zmin_; z<zmax_; z+=dz_)  zcnt++;
  redshifts_.resize(zcnt);
  Ell_cnt_per_Mpc3_.resize(zcnt);
  Sp_cnt_per_Mpc3_.resize(zcnt);
  SB_cnt_per_Mpc3_.resize(zcnt);
  volume_arcmin2_.resize(zcnt);

  KCorrector KC(sedABzero_, lambdaSEDmin_, lambdaSEDmax_);

  su_.SetEmissionRedShift(zmin_);
  double Xi1 = su_.LineOfSightComovDistanceMpc();
  double dL1 = su_.LuminosityDistanceMpc();

  size_t k = 0;
  for(double z=zmin_; z<zmax_; z+=dz_) {
    su_.SetEmissionRedShift(z+dz_);
    double Xi2 = su_.LineOfSightComovDistanceMpc();
    double dL2 = su_.LuminosityDistanceMpc();
    double Xi = (Xi2+Xi1)*0.5;
    double dL = (dL2+dL1)*0.5;
    double deltaXi = (Xi2-Xi1);
    Xi1 = Xi2;   dL1 = dL2;

    if (prtlev_>0)
      cout << " GalaxyCountComputer::doCompute()/Info  z="<<z<<" Xi="<<Xi<<"  dL="<<dL<<"\n ...MagLimit[0...5]= ";
    // Computing cosmological volume element for a solid angle of 1 arcmin^2
    double amin2rad = M_PI/(180.*60.);
    double deltaVol = Xi*amin2rad * Xi*amin2rad * deltaXi; // Mpc^3
    redshifts_[k] = z+0.5*dz_;
    volume_arcmin2_[k]=deltaVol;
    
    std::vector<double> vmaglims(sedGals_.size());
    for(size_t i=0; i<sedGals_.size(); i++) {
      vmaglims[i]=KC.ComputeRestFrameMagnitudeLimit(filt_LF_,filt_Obs_,sedGals_[i],sedABzero_,
						    lambdamin_,lambdamax_,dL,mag_obs_limit_,z, ebmv);
      if (prtlev_>0) cout <<i<<" : "<<vmaglims[i]<<" ;";
    }
    if (prtlev_>0) cout << endl;
    
    std::vector<double> galcnt(sedGals_.size());
    
    MySchechter schEl = multLF_.getLF_Elliptical(z);
    Ell_cnt_per_Mpc3_[k] = 0.;
    for(size_t i=0; i<sedGals_.size(); i++) {
      galcnt[i]=0.;
      if (Ell_sed_mapping_[i] > 1.e-19)  {
	galcnt[i] = IntegrateLF(schEl, vmaglims[i], mag_obs_error_);
	Ell_cnt_per_Mpc3_[k] += galcnt[i]*Ell_sed_mapping_[i];
      }
    }

    MySchechter schSp = multLF_.getLF_Spiral(z);
    Sp_cnt_per_Mpc3_[k] = 0.;
    for(size_t i=0; i<sedGals_.size(); i++) {
      galcnt[i]=0.;
      if (Sp_sed_mapping_[i] > 1.e-19)  {
	galcnt[i] = IntegrateLF(schSp, vmaglims[i], mag_obs_error_);
	Sp_cnt_per_Mpc3_[k] += galcnt[i]*Sp_sed_mapping_[i];
      }
    }

    MySchechter schSB = multLF_.getLF_StarBurst(z);
    SB_cnt_per_Mpc3_[k] = 0.;
    for(size_t i=0; i<sedGals_.size(); i++) {
      galcnt[i]=0.;
      if (SB_sed_mapping_[i] > 1.e-19)   {
	galcnt[i] = IntegrateLF(schSB, vmaglims[i], mag_obs_error_);
	SB_cnt_per_Mpc3_[k] += galcnt[i]*SB_sed_mapping_[i];
      }
    }
    k++;
  }
}

double GalaxyCountComputer::getRedshift(size_t i) const
{
  if (i >= getNbRedshiftBins())
    throw RangeCheckError("GalaxyCountComputer::getRedshift() Out of range redshift index i");
  return redshifts_[i];
}

double GalaxyCountComputer::getGalDensity_Mpc3(size_t i, double& Ellcnt, double & Spcnt, double & SBcnt, double & redshift) const
{
  if (i >= getNbRedshiftBins())
    throw RangeCheckError("GalaxyCountComputer::getGalDensity_Mpc3() Out of range redshift index i");
  Ellcnt = Ell_cnt_per_Mpc3_[i];
  Spcnt = Sp_cnt_per_Mpc3_[i];
  SBcnt = SB_cnt_per_Mpc3_[i];
  redshift=redshifts_[i];
  return (Ellcnt+Spcnt+SBcnt);
}

double GalaxyCountComputer::getGalDensity_Arcmin2(size_t i, double& Ellcnt, double & Spcnt, double & SBcnt, double & redshift) const
{
  double totcnt = getGalDensity_Mpc3(i, Ellcnt, Spcnt, SBcnt, redshift);
  double dvol = volume_arcmin2_[i];
  Ellcnt *= dvol;
  Spcnt *= dvol;
  SBcnt *= dvol;
  return (totcnt*dvol);
}

double GalaxyCountComputer::getAverageGalDensity_Mpc3(double& Ellcnt, double & Spcnt, double & SBcnt) const
{
  Ellcnt = Spcnt = SBcnt = 0.;
  for(size_t i=0; i<redshifts_.size(); i++)  {
    Ellcnt += Ell_cnt_per_Mpc3_[i];
    Spcnt += Sp_cnt_per_Mpc3_[i];
    SBcnt += SB_cnt_per_Mpc3_[i];
  }
  double fac = 1./(double)redshifts_.size();
  Ellcnt *= fac;
  Spcnt *= fac;
  SBcnt *= fac;
  return (Ellcnt+Spcnt+SBcnt);
}

double GalaxyCountComputer::getIntegratedGalDensity_Arcmin2(double& Ellcnt, double & Spcnt, double & SBcnt) const
{
  double redshift;
  Ellcnt = Spcnt = SBcnt = 0.;
  double totcnt = 0.;
  double ell, sp, sb;
  for(size_t i=0; i<redshifts_.size(); i++)  {
    totcnt += getGalDensity_Arcmin2(i, ell, sp, sb, redshift);
    Ellcnt += ell;  Spcnt += sp;  SBcnt += sb;
  }
  return totcnt;
}

double GalaxyCountComputer::IntegrateLF(MySchechter  & sch, double maglim, double magerr)
{
  double magmin=multLF_.getMinMagLimit();
  double magmax=maglim;
  bool fguseeff = false;
  double deltamag = 0.1;
  if (magerr>1.e-19)  {
      deltamag  = magerr;
      fguseeff = true;
      //DEL-DBG      cout << " -----DBG--- Using MagEfficiencyFunc , deltamag="<<deltamag<<endl;
  }
  MagEfficiencyFunc eff(maglim, deltamag);
  Eff_Schechter effsch(sch, eff);
    
  ClassFunc1D & sch4integ = sch;
  if (fguseeff)  sch4integ = effsch;
  /*
  {
    GLInteg integA(sch, magmin, magmax);  // Pourquoi 2*magBlim ?
    integA.SetOrder(MultiType_Z_LF::get_def_schechter_integ_glorder());
    GLInteg integB(effsch, magmin, magmax);  // Pourquoi 2*magBlim ?
    integB.SetOrder(MultiType_Z_LF::get_def_schechter_integ_glorder());
    cout << "--DBG--IntegrateLF() magmin="<<magmin<<" max="<<magmax<<" deltamag="<<deltamag<<" WithEff="<<integB.Value()<< " NoEff="<<integA.Value()<<endl;
  } 
  */
  if (magmax<magmin)  return 0.;
  GLInteg integ(sch4integ, magmin, magmax);  // Pourquoi 2*magBlim ?
  integ.SetOrder(MultiType_Z_LF::get_def_schechter_integ_glorder());
  return integ.Value();
}
