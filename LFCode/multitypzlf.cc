/* ----
   Project   LSST/BAO/PhotoZ
   Classes to generate galaxie Absolute Magnitude and types starting  
   from Luminosity Functions (LF's)  
   F. Habibi & R.Ansari ,  LAL (CNRS/IN2P3)  & Univ. Paris Sud 
   April - July 2017                                          -------  */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <iomanip>

#include "multitypzlf.h"
#include "pexceptions.h"
#include "array.h"   // -- for reading the ASCII input file as a TMatrix 
#include "randfmt.h"
#include "ctimer.h"

using namespace std;
using namespace SOPHYA;

//------------------------------------------------------------------------
//-----------  class MultiType_Z_LF --------------------------------------
//------------------------------------------------------------------------

//--- static variables and methods
//  renormalised value o(/70 km/s/Mpc) of the Hubble constant 
static double def_h70 = 1.;
//-- redshift value at which number of galaxies reaches zero 
static double def_zero_gal_redshift = 20.;
//-- redshift range and number of number of redshifts for which we compute the Integrale[LF](mag) for random drawing of magnitudes
static double def_zmin = 0.;
static double def_zmax = 4.;
static int def_nz = 20;
//-- Schechter integration Gauss-Legendre order, and incremental integration order
static int def_glorder = 100;
static int def_inc_glorder = 5;
//--  do CPU timing or NOT
static bool fgdo_CPUTiming = false;

// define the Hubble constant value (km/s/Mpc), which is used to correct MStar value of the Schechter parameter
void MultiType_Z_LF::set_def_H0(double H0)
{
  return set_def_h70(H0/70.);
}
// define the normalised to 100 km/s/Mpc Hubble constant value (H0/(100 km/s/Mpc) 
void MultiType_Z_LF::set_def_h100(double h100)
{
  return set_def_h70(h100/0.7);
}
// define the normalised to 70 km/s/Mpc Hubble constant value (H0/(70 km/s/Mpc) 
void MultiType_Z_LF::set_def_h70(double h70)
{
  def_h70=h70;
  cout << " MultiType_Z_LF::set_def_h70() Changing default h70 value for correcting Phistar and Mstar (magnitude) def_h70="<<def_h70<<endl;
  return;
}
// return default value of the Hubble constant in km/s/Mpc
double MultiType_Z_LF::get_def_H0()
{
  return def_h70*100.;
}


// define the redshift at which number of galaxies reaches zero
void MultiType_Z_LF::set_def_zero_gal_redshift(double z)
{
  if ((z>0.) && (z<1100.))  def_zero_gal_redshift=z;
  return;
}
// define the redshift range and number of bins for magnitude distribution for all galaxies
void MultiType_Z_LF::set_def_redshift_range(double zmin, double zmax, double nz)
{
  if ((zmin>=0.) && (zmax>zmin) && (nz>0)) {
    def_zmin=zmin;  def_zmax=zmax;  def_nz=nz;
  }
  return;  
}
// define the Gauss-Legendre quadrature order Schechter LF integration from magMin to magMax and for incremental integration
void  MultiType_Z_LF::set_def_schechter_integ_glorder(int glorder, int incglorder)
{
  if (glorder>0)  def_glorder=glorder;
  if (incglorder>0)  def_inc_glorder=incglorder;
  return;
}

int MultiType_Z_LF::get_def_schechter_integ_glorder()
{
  return def_glorder;
}


void  MultiType_Z_LF::activateCPUTiming(bool fgtim)
{
  fgdo_CPUTiming=fgtim;
  return;
}

/* --Methode-- */
MultiType_Z_LF::MultiType_Z_LF(string const & lfparamsfilename, double magMin, double magMax, size_t nbmag, bool useall)
  : magMin_(magMin), magMax_(magMax), nb_mag_pts_(nbmag), useSchechall_(useall)
{
  Timer* ptm;
  if (fgdo_CPUTiming)  ptm = new Timer("MultiType_Z_LF");
  
  TMatrix<double> lfp;
  ifstream is(lfparamsfilename.c_str());
    if (!is) {cerr << " Err: multitypzlf.cc: " << lfparamsfilename
        << " does not exist. Exiting ..." << endl; exit(0); }
  sa_size_t nr,nc;
  lfp.ReadASCII (is, nr, nc); //  char clm='#', const char *sep=" \t")
  cout << " MultiType_Z_LF[1]: read from file"<<lfparamsfilename<<" -> nb of lines="<<nr<<" number of columns=" << nc << endl;
  //DBG  lfp.Print(cout);

  // Correction MStar et PhiStar si H0 different de 70 km/s/Mpc 
  double cor_mstar=0.;   // terme additif -> +5*log10(h70) , magnitude 
  double cor_phistar=1.;  // terme multiplicatif -> phistar en 1/(Mpc/h70)^3   -> (h70)^3 
  if (fabs(def_h70-1.)>1.e-9)  {
    cor_mstar = 5.*log10(def_h70);
    cor_phistar = def_h70*def_h70*def_h70;
    cout << " MultiType_Z_LF[1.b]: Correcting Phistar and MStar magnitude from h=0.7 to h="<< def_h70*0.7 << " \n"
	 << "   ...cor_phistar=(h70)^3="<<cor_phistar<<"  cor_mstar=+5*log10(h70)="<<cor_mstar<<endl;
  }
  if (useSchechall_) {
    cout << " MultiType_Z_LF[1.c]: Phistar for Ellipticals, Spirals & StarBurst LF's will be rescaled to ensure  \n"
	 << "   that Integral[LF_all] = Integral[LF_Ell]+Integral[LF_Sp]+Integral[LF_SB] , for magMin<mag<magMax " << endl;
  }
  double z2last=0.;
  for(sa_size_t j=0; j<lfp.SizeY(); j++)  {
    double z1 = lfp(j,0);
    double z2 = lfp(j,1);
    if ((j>0)&&(fabs(z1-z2last)>1.e-9))
      throw ParmError("MultiType_Z_LF::MultiType_Z_LF: redshift interval error in file");
    pair<double,double> z12;   z12.first=z1;  z12.second=z2;
    zrange.push_back(z12);
    z2last=z2;

    sa_size_t off=2;
    MySchechter sfell(lfp(j,off)*cor_phistar, lfp(j,off+1)+cor_mstar, lfp(j,off+2));
    off+=3;
    MySchechter sfsp(lfp(j,off)*cor_phistar, lfp(j,off+1)+cor_mstar, lfp(j,off+2));
    off+=3;
    MySchechter sfsb(lfp(j,off)*cor_phistar, lfp(j,off+1)+cor_mstar, lfp(j,off+2));
    off+=3;
    MySchechter sfall(lfp(j,off)*cor_phistar, lfp(j,off+1)+cor_mstar, lfp(j,off+2));
    if (useSchechall_) {
      double phiscale=sfall.getIntegral(magMin_,magMax_,def_glorder) / 
	( sfell.getIntegral(magMin_,magMax_,def_glorder)+
	  sfsp.getIntegral(magMin_,magMax_,def_glorder) +
	  sfsb.getIntegral(magMin_,magMax_,def_glorder) ) ;
      sfell.scalePhistar(phiscale);
      sfsp.scalePhistar(phiscale);
      sfsb.scalePhistar(phiscale);
    }
    vshEll.push_back(sfell);
    vshSp.push_back(sfsp);
    vshSB.push_back(sfsb);
    vshAll.push_back(sfall);
  }

  // checking if the first range starts at zero 
  first_z_is_zero=false;
  if (fabs(zrange[0].first)<1.e-12)  {
    zrange[0].first=0.;
    first_z_is_zero=true;
  }

  zero_gal_redshift_=def_zero_gal_redshift;
  cout << " MultiType_Z_LF[2]: Computing Integral of LF's from magMin="<<magMin_<<" to magMax="<<magMax_<<" \n"
       << "  with linear interpolation as a function of redshift , with ngal=0 for redshift z="<<zero_gal_redshift_<<endl; 

  size_t vsz = (first_z_is_zero ? zrange.size()+3 : zrange.size()+4);
  size_t voff = (first_z_is_zero ? 1 : 2);
  
  vector<double> leszm(vsz); // z at the middle of the interval
  vector<double> lesilfEll(vsz); // integral of LF at a given z interval
  vector<double> lesilfSp(vsz);
  vector<double> lesilfSB(vsz);
  vector<double> lesilfAll(vsz);

  if ( first_z_is_zero ) 
    leszm[0]=zrange[0].first;
  else {
    leszm[0]=0.;
    leszm[1]=zrange[0].first;
  }
  for(size_t k=0; k<zrange.size(); k++)  {
    leszm[k+voff]=(zrange[k].first+zrange[k].second)*0.5;
    //DBG    cout << "**********DBG** k="<<k<<" z="<<leszm[k+voff]<<endl;
    MySchechter sf=vshEll[k];
    lesilfEll[k+voff]=sf.getIntegral(magMin_,magMax_,def_glorder);
    //DBG    cout<<"*DBG*Elliptical:"<<sf<<" Integral="<<lesilfEll[k+voff]<<endl;
    sf=vshSp[k];
    lesilfSp[k+voff]=sf.getIntegral(magMin_,magMax_,def_glorder);
    //DBG    cout<<"*DBG*Spiral:"<<sf<<" Integral="<<lesilfSp[k+voff]<<endl;
    sf=vshSB[k];
    lesilfSB[k+voff]=sf.getIntegral(magMin_,magMax_,def_glorder);
    //DBG    cout<<"*DBG*Starburst:"<<sf<<" Integral="<<lesilfSB[k+voff]<<endl;
    sf=vshAll[k];
    lesilfAll[k+voff]=sf.getIntegral(magMin_,magMax_,def_glorder);
    //DBG    cout<<"*DBG*All:"<<sf<<" Integral="<<lesilfAll[k+voff]<<"  Sum="<<lesilfEll[k+voff]+lesilfSp[k+voff]+lesilfSB[k+voff]
    //	<<" FracSB="<<lesilfSB[k+voff]/(lesilfEll[k+voff]+lesilfSp[k+voff]+lesilfSB[k+voff])<<endl;
  }
  
  leszm[zrange.size()+voff] = zrange[zrange.size()-1].second;
  leszm[zrange.size()+voff+1] = zero_gal_redshift_;   // maaximum z at which number of galaxies goes to zero


  if ( first_z_is_zero ) {
    lesilfEll[0]=lesilfEll[voff];
    lesilfSp[0]=lesilfSp[voff];
    lesilfSB[0]=lesilfSB[voff];
    lesilfAll[0]=lesilfAll[voff];
  }
  else {
    lesilfEll[0]=lesilfEll[1]=lesilfEll[voff];
    lesilfSp[0]=lesilfSp[1]=lesilfSp[voff];
    lesilfSB[0]=lesilfSB[1]=lesilfSB[voff];
    lesilfAll[0]=lesilfAll[1]=lesilfAll[voff];
  }
  
  lesilfEll[zrange.size()+voff]=lesilfEll[zrange.size()+voff-1];
  lesilfSp[zrange.size()+voff]=lesilfSp[zrange.size()+voff-1];
  lesilfSB[zrange.size()+voff]=lesilfSB[zrange.size()+voff-1];
  lesilfAll[zrange.size()+voff]=lesilfAll[zrange.size()+voff-1];

  lesilfEll[zrange.size()+voff+1]=0.;
  lesilfSp[zrange.size()+voff+1]=0.;
  lesilfSB[zrange.size()+voff+1]=0.;
  lesilfAll[zrange.size()+voff+1]=0.;


  ILF_of_z_Ell.DefinePoints(leszm, lesilfEll);
  ILF_of_z_Sp.DefinePoints(leszm, lesilfSp);
  ILF_of_z_SB.DefinePoints(leszm, lesilfSB);
  ILF_of_z_All.DefinePoints(leszm, lesilfAll);

  zmin_=def_zmin;
  zmax_=def_zmax;
  nz_=def_nz;
  dz_=(zmax_-zmin_)/(double)nz_;

  //  --- CPU timing 
  if (fgdo_CPUTiming)  ptm->Split("Done [1]+[2]");
			   
  cout << " MultiType_Z_LF[3]: Computing Sum_[Ell,Spiral,StarBurst]  of Integral_LF[magMin="
       <<magMin_<<","<<"magMax="<<magMax_<<"]  \n for "<<nz_+1<<" redshifts "<<zmin_<<" <=z<= "<<zmax_<<endl;
  vector<double> lesz((size_t)(nz_+1));
  vector<double> lessum((size_t)(nz_+1));
  for(size_t k=0; k<lesz.size(); k++) {
    double redshift=zmin_+(double)k*dz_;
    double sumILF=ILF_of_z_Ell(redshift)+ILF_of_z_Sp(redshift)+ILF_of_z_SB(redshift);
    lesz[k]=redshift;
    lessum[k]=sumILF;
  }
  ILF_of_z_Sum_EllSpSB.DefinePoints(lesz, lessum);
  
  //  --- CPU timing 
  if (fgdo_CPUTiming)  ptm->Split("Done [3]");
			   
  vector<double> lesmag(nb_mag_pts_+1);
  vector<double> lesILF(nb_mag_pts_+1);
  double dmag=(magMax_-magMin_)/(double)nb_mag_pts_;
  lesmag[0]=magMin_;
  double curmag=magMin_;
  for(size_t k=0; k<lesmag.size()-1; k++) {
    lesmag[k]=curmag;  curmag+=dmag;
  }
  lesmag[lesmag.size()-1]=magMax_;

  if (useSchechall_) {  // We use Schechter parameters for All galaxies
    cout << " MultiType_Z_LF[4]: Computing Integral_LF as a function of magnitude from All galaxies Schechter \n"
	 << " parameters for the corresponding " << zrange.size() << " redshift intervals ..."<<endl;
    v_mag_f_ILF.resize(zrange.size(),NULL);
    for(size_t k=0; k<zrange.size(); k++)  {
      SLinInterp1D * mag_f_ILF = new SLinInterp1D;
      double curILF = 0.;
      lesILF[0]=curILF;
      MySchechter sf=vshAll[k];
      for(size_t i=1; i<lesmag.size(); i++) {
	curILF += sf.getIntegral(lesmag[i-1],lesmag[i],def_inc_glorder);
	lesILF[i]=curILF;
      }
      for(size_t i=1; i<lesILF.size(); i++)  lesILF[i] /= curILF;
      lesILF[lesILF.size()-1]=1.;
      mag_f_ILF->DefinePoints(lesILF, lesmag);
      v_mag_f_ILF[k]=mag_f_ILF;
    }
  }
  else {  // We use Sum[Ell+Spiral+StarBurst] for computing magnitude distribution for  All galaxies
    cout << " MultiType_Z_LF[4]: Computing Integral_LF as a function of magnitude from Sum[Ell+Spiral+StarBurst] Schechters \n"
	 << " for "<<nz_+1<<" redshifts values " <<zmin_<<" <=z<= "<<zmax_<<endl;
    v_mag_f_ILF.resize((size_t)nz_,NULL);

    lesz.resize((size_t)nz_);
    for(size_t k=0; k<lesz.size(); k++) {
      double redshift=zmin_+((double)k+0.5)*dz_;
      SLinInterp1D * mag_f_ILF = new SLinInterp1D;
      double curILF = 0.;
      lesILF[0]=curILF;
      MySchechter sfE=getLF_Elliptical(redshift);
      MySchechter sfSp=getLF_Spiral(redshift);
      MySchechter sfSB=getLF_StarBurst(redshift);
      for(size_t i=1; i<lesmag.size(); i++) {
	curILF += sfE.getIntegral(lesmag[i-1],lesmag[i],def_inc_glorder);
	curILF += sfSp.getIntegral(lesmag[i-1],lesmag[i],def_inc_glorder);
	curILF += sfSB.getIntegral(lesmag[i-1],lesmag[i],def_inc_glorder);
	lesILF[i]=curILF;
      }
      for(size_t i=1; i<lesILF.size(); i++)  lesILF[i] /= curILF;
      lesILF[lesILF.size()-1]=1.;
      mag_f_ILF->DefinePoints(lesILF, lesmag);
      v_mag_f_ILF[k]=mag_f_ILF;
    }
  }
  //  --- CPU timing 
  if (fgdo_CPUTiming)  ptm->Split("Done [4]");

  //  Creating / defining default random generator 
  def_rgp_ = new FMTRandGen();
  rgp_ = def_rgp_;
  
  //  --- CPU timing 
  if (fgdo_CPUTiming)  delete ptm;
}

/* --Methode-- */
MultiType_Z_LF::~MultiType_Z_LF()
{
  delete def_rgp_ ;
  for (size_t i=0; i<v_mag_f_ILF.size(); i++)  delete v_mag_f_ILF[i];
}

/* --Methode-- */
size_t MultiType_Z_LF::find_redshift_bin(double z)  const 
{
  if (z < 0.) 
    throw ParmError("MultiType_Z_LF::find_redshift_bin()/ERROR negative redshift !");

  size_t rzb=0;
  for(size_t k=0; k<zrange.size(); k++) {
    if (z < zrange[k].first)  break;
    rzb++;
  }
  if (rzb>0)  rzb--;
  if (rzb >= zrange.size())  rzb=zrange.size()-1;
  return rzb;
}


/* --Methode-- */
MySchechter MultiType_Z_LF::getLF_Elliptical(double z)  const 
{
  size_t kbin = find_redshift_bin(z);
  MySchechter rsf = vshEll[kbin];
  // scale phistar
  rsf.scalePhistar( ILF_of_z_Ell(z)/ILF_of_z_Ell(0.5*(zrange[kbin].first+zrange[kbin].second)) );
  return rsf;
}

/* --Methode-- */
MySchechter MultiType_Z_LF::getLF_Spiral(double z)  const 
{
  size_t kbin = find_redshift_bin(z);
  MySchechter rsf = vshSp[kbin];
  // scale phistar
  rsf.scalePhistar( ILF_of_z_Sp(z)/ILF_of_z_Sp(0.5*(zrange[kbin].first+zrange[kbin].second)) );
  return rsf;
}

/* --Methode-- */
MySchechter MultiType_Z_LF::getLF_StarBurst(double z)   const 
{
  size_t kbin = find_redshift_bin(z);
  MySchechter rsf = vshSB[kbin];
  // scale phistar
  rsf.scalePhistar( ILF_of_z_SB(z)/ILF_of_z_SB(0.5*(zrange[kbin].first+zrange[kbin].second)) );
  return rsf;
}

/* --Methode-- */
MySchechter MultiType_Z_LF::getLF_All(double z)  const 
{
  size_t kbin = find_redshift_bin(z);
  MySchechter rsf = vshAll[kbin];
  // scale phistar
  rsf.scalePhistar( ILF_of_z_All(z)/ILF_of_z_All(0.5*(zrange[kbin].first+zrange[kbin].second)) );
  return rsf;
}

/* --Methode-- */
double MultiType_Z_LF::getGalaxyNumberDensity(double z)  const 
{
  if (z < 0.) 
    throw ParmError("MultiType_Z_LF::getGalaxyNumberDensity(z)/ERROR negative redshift !");

  if (useSchechall_)  return ILF_of_z_All.YInterp(z);
  return ILF_of_z_Sum_EllSpSB.YInterp(z);
}

/* --Methode-- */
double MultiType_Z_LF::getGalaxyNumberDensity(double z, double& fEll, double& fSp, double& fSB)  const 
{
  if (z < 0.) 
    throw ParmError("MultiType_Z_LF::getGalaxyNumberDensity(z, fEll, fSp, fSB)/ERROR negative redshift !");
  double galdens=0.;
  if (useSchechall_)  galdens=ILF_of_z_All.YInterp(z);
  else galdens=ILF_of_z_Sum_EllSpSB.YInterp(z);
  fEll=ILF_of_z_Ell(z);
  fSp=ILF_of_z_Sp(z);
  fSB=ILF_of_z_SB(z);
  double s=(fEll+fSp+fSB);
  if (s>1.e-39)  {
    fEll/=s;  fSp/=s;  fSB/=s;  
  }
  return galdens;
}

/* --Methode-- */
void MultiType_Z_LF::P_getTypeMagnitude(double z, double& mag, bool fgmag, int& type, bool fgtype)  const 
{
  size_t kbin = 0;
  size_t kbinm = 0;
  if (fgtype)  kbin = find_redshift_bin(z);
  if (fgmag) {  // draw magnitude at redshift z 
    if (useSchechall_)  {
      kbinm = (fgtype ? kbin : find_redshift_bin(z));
    }
    else {
      if (z < 0.) 
	throw ParmError("MultiType_Z_LF::getMag(z)/ERROR negative redshift !");
      kbinm = (z-zmin_)/dz_;
      if (kbinm >= v_mag_f_ILF.size())  kbinm=v_mag_f_ILF.size()-1;
    }
    double rnd=rgp_->Flat01();
    mag=magMin_;
    mag = v_mag_f_ILF[kbinm]->YInterp(rnd);
  }
  if (fgtype) {  // draw type for magnitude mag at redshift z 
    type=0;
    // Get rescaled Schechter functions for the three types 
    MySchechter rsfE = vshEll[kbin];
    rsfE.scalePhistar( ILF_of_z_Ell(z)/ILF_of_z_Ell(0.5*(zrange[kbin].first+zrange[kbin].second)) );
    MySchechter rsfSp = vshSp[kbin];
    rsfSp.scalePhistar( ILF_of_z_Sp(z)/ILF_of_z_Sp(0.5*(zrange[kbin].first+zrange[kbin].second)) );
    MySchechter rsfSB = vshSB[kbin];
    rsfSB.scalePhistar( ILF_of_z_SB(z)/ILF_of_z_SB(0.5*(zrange[kbin].first+zrange[kbin].second)) );
    double rE = rsfE(mag);
    double rSp = rsfSp(mag);
    double rSB = rsfSB(mag);
    double rSum = rE+rSp+rSB;
    rE /= rSum;  rSp /= rSum;  // rSB /= rSum;
    double rnd=rgp_->Flat01();
    if (rnd < rE)  type = 1;  // Elliptical
    else if (rnd < (rE+rSp))  type = 2;  // Spiral
    else type = 3;  // StarBurst  
  }
  return ;
}


ostream& MultiType_Z_LF::print(ostream& os)  const
{
  os << "# LF_Params: phistar  MStar  alpha ]   magMin="<<magMin_<<" magMax="<<magMax_<<endl;
  os << "# z-range  LF_Params_Ell  LF_Params_Spiral  LF_Params_StarBurst  LF_Params_All "<<endl;
  for(size_t k=0; k<zrange.size(); k++) {
    os << setw(8) << zrange[k].first << "  " << zrange[k].second << " \t " << vshEll[k] << " \t " << vshSp[k]
       << " \t " << vshSB[k] << " \t " << vshAll[k] << endl;
  }
  return os; 
}



