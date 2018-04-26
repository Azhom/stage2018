/* ----
   Project   LSST/BAO/PhotoZ
   Tests de calcul de P(k)-2D et fct d'auto-correlation 
   R.Ansari - Feb 2015 
                                                     -------  */

//----- c/c++ includes 
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include <typeinfo>

//----- sophya includes 
#include "machdefs.h"
#include "sopnamsp.h"

#include "histats.h"
#include "array.h"
#include "fitsioserver.h"
#include "fiosinit.h"

// Pour les cartes spheriques
#include "skymap.h"
// pour les transformees Ylm 
#include "samba.h"

//-------------------------------------------------------------------------
//      ------------------ MAIN PROGRAM ------------------------------
//-------------------------------------------------------------------------

Vector ClMean(std::vector<Vector> cls) {
  int nsim = cls.size();
  Vector cl_mean(cls[0],false);
  cl_mean = 0;
  for(int i=0; i<cl_mean.Size(); i++){
    for(int k=0; k<nsim; k++) 
      cl_mean[i] += cls[k][i];
  }
  cl_mean /= nsim;
  return(cl_mean);
}

Vector ClStd(std::vector<Vector> cls) {
  Vector cl_mean = ClMean(cls);
  int nsim = cls.size();
  Vector cl_std(cls[0],false);
  cl_std = 0;
  for(int i=0; i<cl_std.Size(); i++){
    for(int k=0; k<nsim; k++) 
      cl_std[i] += pow(cls[k][i]-cl_mean[i],2);
  }
  cl_std /= (nsim-1);
  for(int i=0; i<cl_std.Size(); i++)
    cl_std[i] = sqrt(cl_std[i]);
  return(cl_std);
}


int main(int narg, const char* arg[])
{
  if ((narg<3)||((narg>1)&&(strcmp(arg[1],"-h")==0))) {
    cout << " Usage: clcut [options] InputClFile Out_PPFFile \n"
	 << "   InputClFile : text file  with input C(l) file \n"
	 << " options: [-N nsim] [-d dx,dy,dz] [-k kmin,kmax,nbin] [-k2d nbin2d,kmax2d] \n"
	 << "          [-s sigmaz] [-t lowval] [-r rfac] [-p lev] \n"
	 << endl;       //<< "   -N Nx,Ny,Nz : define input 3D grid size default=400x400x400 \n"
    return 1;
  }
  Timer tm("clcut");
  int rc = 0;
  try { 
    string inclname="clin.txt";
    string outppfname="toto.ppf";
	
    int prtlev=0;
    int nsim = 300;
    //----------------------------------------------------------
    // decodage arguments optionnel 
    bool fgoptarg=true;
    while (fgoptarg&&(narg>1)) {
      string fbo = arg[1];
      if (fbo=="-N")  {   // specification nombre de générations de Cl sur lesquelles moyenner
	if (narg<3) { cout << " clcut/missing/bad argument, clcut -h for help " << endl;  return 2; }
	sscanf(arg[2],"%d",&nsim);
	arg+=2; narg-=2; 
      }
      else if (fbo=="-d")  {   // specification taille de cellules (Mpc) drho/rho ou ngal 
	if (narg<3) { cout << " clcut/missing/bad argument, clcut -h for help " << endl;  return 2; }
	double dx, dy, dz;
	sscanf(arg[2],"%lf,%lf,%lf",&dx,&dy,&dz); 
	arg+=2; narg-=2; 
      }
      else fgoptarg=false;
    }
    //----------------------------------------------------------
    if (narg<3) { cout << " clcut/missing/bad argument, clcut -h for help " << endl;  return 2; }

    inclname=arg[1];
    outppfname=arg[2];

    cout << "clcut[1] : Opening output PPF file " << outppfname <<endl;
    POutPersist pof(outppfname);
    
    cout << "clcut[2] : Reading input C(l) file " << inclname <<endl;
    Vector clin;
    sa_size_t nr, nc;
    ifstream is(inclname.c_str());
    clin.ReadASCII (is, nr, nc);
    int lmax=nr-1;
    cout << " lmax from file = " << lmax <<endl; 
    cout<<clin.Transpose();
    
    // Create a set of Alm coefficients distributed according to C(l)=clin
    // Alm<double> alm(clin, 10., true);  // single-sided
    Alm<double> alm(lmax);  // il s agit d un alm symetrique; seul alm , m>=0 garde
    std::vector<Vector> clcks(nsim);
    for(int k=0; k<nsim; k++) {
      alm.GenFromCl(clin);
      //cout<<k<<" "<<alm<<endl;
      clcks[k]=alm.powerSpectrum();
    }
    Vector clck = ClMean(clcks);
    
    // Instanciate an SHT server object 
    SphericalTransformServer<r_8> shts;
    // Synthetize a HEALPix and SphereThetaPhi map from alm
    int nside = 64;  // HealPix pixelisation parameter
    SphereHEALPix<r_8> mapH(nside);
    shts.GenerateFromAlm(mapH, nside, alm); 
    cout << "clcut[3] : Computed HEALPix mapH from alm " << outppfname <<endl;
    cout<<mapH;
    // Decompose map into Ylm 
    Alm<double> almc(lmax, true);   // single-sided 
    shts.DecomposeToAlm(mapH, almc, lmax, 0.);
    cout<<"clcut[4] : Computed almc computed from HEALPix mapH"<<endl;
    // Get the power spectrum
    Vector clc=almc.powerSpectrum();
    /*
    // Synthetize also a SphereThetaPhi map from the same set of alm
    int nring = 257;  // number of rings (=theta slices) 
    SphereThetaPhi<r_8> mapT(nring);
    shts.GenerateFromAlm(mapT, nring, alm); 
    cout<<" Computed SphereThetaPhi mapT from alm "<<endl;
    cout << mapT;
    */
    pof<<PPFNameTag("clin")<<clin;
    pof<<PPFNameTag("alm")<<alm;
    pof<<PPFNameTag("clck")<<clck;
    pof<<PPFNameTag("clc")<<clc;
    pof<<PPFNameTag("mapH")<<mapH;

    cout << "clcut[5] : kept clin alm clck clc mapH ..."<<endl;
    // On applique une coupure simple angle en theta
    /* SphereHEALPix<r_8> mapHC(nside);
       mapHC=mapH;
    */
    std::vector<Vector> clcuts(nsim);
    for(int ii=0; ii<nsim; ii++) {
      alm.GenFromCl(clin);
      shts.GenerateFromAlm(mapH, nside, alm); 
      SphereHEALPix<r_8> mapHC(mapH, false);  // false -> on ne partage pas les pixels 
      for(int k=0; k<mapHC.NbPixels(); k++) {
	double tet,phi;
	mapHC.PixThetaPhi(k,tet,phi);
	if (Angle(tet).ToDegree()>30)  mapHC[k]=0.;
      }
      Alm<double> almcut(lmax, true);   // single-sided 
      shts.DecomposeToAlm(mapHC, almcut, lmax, 0.);
      //cout<<" Computed almcut computed from HEALPix mapHC"<<endl;
      // Get the power spectrum
      clcuts[ii]=almcut.powerSpectrum();
      if(ii==nsim-1) {
	pof<<PPFNameTag("mapHC")<<mapHC;
      }
    }
    Vector clcut = ClMean(clcuts);
    /*
    for(int k=0; k<nsim; k++) {
      alm.GenFromCl(clin);
      shts.GenerateFromAlm(mapH, nside, alm); 
      SphereHEALPix<r_8> mapHC(mapH, false);  // false -> on ne partage pas les pixels 
      for(int k=0; k<mapHC.NbPixels(); k++) {
	double tet,phi;
	mapHC.PixThetaPhi(k,tet,phi);
	if (Angle(tet).ToDegree()>30)  mapHC[k]=0.;
      }
      Alm<double> almcut(lmax, true);   // single-sided 
      shts.DecomposeToAlm(mapHC, almcut, lmax, 0.);
      //cout<<" Computed almcut computed from HEALPix mapHC"<<endl;
      // Get the power spectrum
      clcuts[k]=almcut.powerSpectrum();
    }
    Vector clcut2 = ClMean(clcuts);
    for(int k=0; k<nsim; k++) {
      alm.GenFromCl(clin);
      shts.GenerateFromAlm(mapH, nside, alm); 
      SphereHEALPix<r_8> mapHC(mapH, false);  // false -> on ne partage pas les pixels 
      for(int k=0; k<mapHC.NbPixels(); k++) {
	double tet,phi;
	mapHC.PixThetaPhi(k,tet,phi);
	if (Angle(tet).ToDegree()>30)  mapHC[k]=0.;
      }
      Alm<double> almcut(lmax, true);   // single-sided 
      shts.DecomposeToAlm(mapHC, almcut, lmax, 0.);
      //cout<<" Computed almcut computed from HEALPix mapHC"<<endl;
      // Get the power spectrum
      clcuts[k]=almcut.powerSpectrum();
    }
    Vector clcut3 = ClMean(clcuts);
    */

    pof<<PPFNameTag("clcut")<<clcut;
    //pof<<PPFNameTag("clcut2")<<clcut2;
    //pof<<PPFNameTag("clcut3")<<clcut3;
    //pof<<PPFNameTag("mapHC")<<mapHC;
    cout << "clcut[6] : kept clcut mapHC ..."<<endl;
    Vector clratio;// = clcut.Div(clin);
    clcut.DivElt(clin,clratio);
    cout<< clratio<<" "<<clcut<<endl;
    clratio[lmax]=0;
    pof<<PPFNameTag("clratio")<<clratio;

  }  // End of try bloc 
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " clcut.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " clcut.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " clcut.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
  }
  cout << " ==== End of clcut.cc program  Rc= " << rc << endl;
  return rc;    
}
