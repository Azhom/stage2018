/* ----
   Project   LSST/BAO/PhotoZ
   Classes to compute k-corrections 
   F. Habibi & R.Ansari ,  LAL (CNRS/IN2P3)  & Univ. Paris Sud 
   April - Sep 2017                                       -------  */

#include "kcorr.h"
#include "ntuple.h"
#include "integ.h"

//----------------------------------------------------------------------
//--------------------- KCorrector class code --------------------------
//----------------------------------------------------------------------
double KCorrector::computeK(SED & sed, Filter & filtx, Filter &filty)
{
  SEDFilterProd integrand_x(sed,filtx);
  SEDFilterProd integrand_y(sed,filtx);
  double Fx = FilterIntegrator(integrand_x,lambdaMin_*1e-9, lambdaMax_*1e-9).Value();
  double Fy = FilterIntegrator(integrand_y,lambdaMin_*1e-9, lambdaMax_*1e-9).Value();
  
  double Kxy = -2.5*log10(Fy/Fx);
  return(Kxy);
}

double KCorrector::computeKz(SED & sed, Filter & filtx, Filter &filty, double z)
{
  // Multiplying SED by filters and integrating them
  
  SEDFilterProd integrand_x(sed,filtx);
  SEDFilterProd integrand_x0(sedABzero_,filtx);
  // same code as above except for the next line where SEDzFilterProd replaces SEDFilterProd
  SEDzFilterProd integrand_y(sed,filty,z);
  SEDFilterProd integrand_y0(sedABzero_,filty);
  
  double Fx = FilterIntegrator(integrand_x,lambdaMin_*1e-9, lambdaMax_*1e-9).Value();
  double Fy = FilterIntegrator(integrand_y,lambdaMin_*1e-9, lambdaMax_*1e-9).Value();
  double F0x = FilterIntegrator(integrand_x0,lambdaMin_*1e-9, lambdaMax_*1e-9).Value();
  double F0y = FilterIntegrator(integrand_y0,lambdaMin_*1e-9, lambdaMax_*1e-9).Value();
  
  double Kxy = -2.5*log10(Fy/Fx*F0x/F0y);
  
  return(Kxy);
}


double KCorrector::ComputeRestFrameMagnitudeLimit(Filter & filt_RF,Filter & filt_Obs, SED &sed, SED &sedABzero, double lambdamin,double lambdamax, double dL, double maglimObs, double z, double ebmv)
{
 // B stand here for rest-frame filter , and r for observer-frame filter
    // Multiplying SED by filters and integrating them
//DBG    cout << " ComputeRestFrameMagnitudeLimit: Integrating the restframe SED with the B filter to give Flux_B: " << endl;
    
    SEDzFilterProd integrand_B(sed,filt_RF,z);
    SEDzFilterProd integrand_Babs(sed,filt_RF);
    SEDzFilterProd integrand_B0(sedABzero,filt_RF);
    // Flux emitted by the galaxy in B-band  (rest-frame)
    double FB = FilterIntegrator(integrand_B,lambdamin*1e-9, lambdamax*1e-9).Value();
    double FBabs = FilterIntegrator(integrand_Babs,lambdamin*1e-9, lambdamax*1e-9).Value();
    double F0B = FilterIntegrator(integrand_B0,lambdamin*1e-9, lambdamax*1e-9).Value();
//DBG    cout << "    FB=" << FB << endl;
//DBG    cout << "    FBabs=" << FBabs << endl;
//DBG    cout << "   KCorrector::ComputeRestFrameMagnitudeLimit: Integrating the redshifted SED with the r filter to give Flux_r at z=" << z << endl;
    SEDzFilterProd integrand_r(sed,filt_Obs,z);
    SEDzFilterProd integrand_r0(sedABzero,filt_Obs); // z=0
    double Fr = FilterIntegrator(integrand_r,lambdamin*1e-9, lambdamax*1e-9).Value();
    double F0r = FilterIntegrator(integrand_r0,lambdamin*1e-9, lambdamax*1e-9).Value();
//DBG    cout << "    Fr=" << Fr << endl;
//DBG    cout << "    F0r=" << F0r << ", F0B=" << F0B << endl;
    
    double mB = -2.5*log10(FB/F0B/4/M_PI/dL/dL);
    double mr = -2.5*log10(Fr/F0r/4/M_PI/dL/dL);
    
    // Computing the K-correction
    double KrB = mr-mB - 2.5*log10(FB/FBabs);
    
    KCorrector KC(sedABzero, lambdamin, lambdamax);
    double KClass = KC.computeKz(sed, filt_RF, filt_Obs, z);
//DBG    cout << "   KCorrector::ComputeRestFrameMagnitudeLimit:  K-correction= " << KrB << ", KClasse=" << KClass << endl;
    if(fabs(KClass-KrB)>1.e-9)
        throw ParmError("KCorrector::ComputeRestFrameMagnitudeLimitBug(): KClass-KrB ");
    
    // Computing the absolute magnitude limit for B-band according to apparen r-band mag limit
    double magb = maglimObs - KrB ;
    double magBlim = magb - 5.*log10(dL) - 25.;// + 2.5*log10(FB/FBabs)
   
    //Adding extinction sperately
    if (ebmv != 0.){
		//SEDzFilterProd integrand(sed, filt_Obs, z);
		//double F = FilterIntegrator(integrand, lambdamin*1e-9, lambdamax*1e-9).Value();
		
		sed.doRedden(ebmv, z);
		SEDzFilterProd integrand_red(sed, filt_Obs, z);
		double FRed = FilterIntegrator(integrand_red, lambdamin*1e-9, lambdamax*1e-9).Value();
		sed.stopRedden();
		double delta_m = -2.5*log10(FRed/Fr); 
		magBlim -= delta_m;
    }
    //DBG cout << magBlim << endl;
    return magBlim; 
}


double ComputeBlim(Filter & filt_B,Filter & filt_r, SED &sed, SED &sedABzero, double lambdamin,double lambdamax, double dL, double magr, double z)
{
//DEPRECATED cout << "   ComputeBlim: Reading SED from file "<< sedfile <<endl;
// DEPRECATED    SED sed(sedfile,lambdamin*1e-9, lambdamax*1e-9);
    
    
  // Multiplying SED by filters and integrating them
  cout << " ComputeBlim: Integrating the restframe SED with the B filter to give Flux_B: " << endl;
  
  SEDzFilterProd integrand_B(sed,filt_B,z);
  SEDzFilterProd integrand_Babs(sed,filt_B);
  SEDzFilterProd integrand_B0(sedABzero,filt_B);
  // Flux emitted by the galaxy in B-band
  double FB = FilterIntegrator(integrand_B,lambdamin*1e-9, lambdamax*1e-9).Value();
  double FBabs = FilterIntegrator(integrand_Babs,lambdamin*1e-9, lambdamax*1e-9).Value();
  double F0B = FilterIntegrator(integrand_B0,lambdamin*1e-9, lambdamax*1e-9).Value();
  cout << "    FB=" << FB << endl;
  cout << "    FBabs=" << FBabs << endl;
  
  cout << "   ComputeBlim: Integrating the redshifted SED with the r filter to give Flux_r at z=" << z << endl;
  SEDzFilterProd integrand_r(sed,filt_r,z);
  SEDzFilterProd integrand_r0(sedABzero,filt_r); // z=0
  double Fr = FilterIntegrator(integrand_r,lambdamin*1e-9, lambdamax*1e-9).Value();
  double F0r = FilterIntegrator(integrand_r0,lambdamin*1e-9, lambdamax*1e-9).Value();
  cout << "    Fr=" << Fr << endl;
  cout << "    F0r=" << F0r << ", F0B=" << F0B << endl;
  
  
  double mB = -2.5*log10(FB/F0B/4/M_PI/dL/dL);
  double mr = -2.5*log10(Fr/F0r/4/M_PI/dL/dL);
  
  // Computing the K-correction
  double KrB = mr-mB - 2.5*log10(FB/FBabs);
  
  KCorrector KC(sedABzero, lambdamin, lambdamax);
  double KClass = KC.computeKz(sed, filt_B, filt_r, z);
  cout << "   ComputeBlim:  K-correction= " << KrB << ", KClasse=" << KClass << endl;
  if(fabs(KClass-KrB)>1.e-9) throw ParmError("Bug: KClass-KrB ");
  
  // Computing the absolute magnitude limit for B-band according to apparen r-band mag limit
  double magb = magr - KrB ;
  double magBlim = magb - 5.*log10(dL) - 25.;// + 2.5*log10(FB/FBabs);
  
  cout << "DEB: ComputeBlim: magr=" << magr << ", magb=" << magb << ", KrB=" << KrB << ", dL=" << dL << ", magBlim=" << magBlim << endl;
  
  double magBlim_2 = magr - mr -2.5*log10(FBabs/F0B/4/M_PI) - 25.;
  if ( fabs(magBlim_2 - magBlim) > 1.e-9 )  cout << " **** BUG magBlim2="<<  magBlim_2 << " diff="<<magBlim_2 - magBlim<<endl;
  
  double magBlim_3 = magr - 5.*log10(dL) - 25 + 2.5*log10(Fr * F0B / ( FBabs * F0r) );
  if ( fabs(magBlim_3 - magBlim) > 1.e-9 )  cout << " **** BUG magBlim3="<<  magBlim_3 << " diff="<<magBlim_3 - magBlim<<endl;
  
  cout << "DEB: ComputeBlim: magBlim=" << magBlim << endl;
  
  return magBlim; 
}

double ComputeNg(MySchechter  & sch, double magBlim, double Xi, double deltaXi, double magerr)
{
  //if(magBlim>0)  throw ParmError("ComputeNg: positive magBLim given -> ERROR ");
  //if(fabs(magBlim-1111.)<1.e-6)
    //return 0.;
  
  double magmin=-24;
  //double magmin= 2.*magBlim;
  double magmax=magBlim;
  bool fguseeff = false;
  double deltamag = 0.1;
  if (magerr>0)  {
      deltamag  = magerr;
      fguseeff = true;
  }
  MagEfficiencyFunc eff(magBlim, deltamag);
  Eff_Schechter effsch(sch, eff);
    
  ClassFunc1D & sch4integ = sch;
  if (fguseeff)  sch4integ = effsch;
   
  //    TrpzInteg integ(sfh, 2.*magBlim, magBlim);  // Pourquoi 2*magBlim ?
  GLInteg integ(sch4integ, magmin, magmax);  // Pourquoi 2*magBlim ?
  integ.SetOrder(100);
  double ng = integ.Value();
  double amin2rad = M_PI/(180.*60.);
  //double V = 1.;
  //double V = 80. * Xi*amin2rad * Xi*amin2rad; //Mpc^3
  double V = Xi*amin2rad * Xi*amin2rad * deltaXi; // Mpc^3
  //cout << "   ComputeNg: Comoving volume with depth of 80 Mpc in 1 amin squared solid angle, V="
  //  << V << " Mpc^3" << endl;
  // Computing the number of galaxies in a given volume
  double Ng = ng*V;
  return Ng;
  
}
