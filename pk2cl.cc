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

#include "interpcosmoc.h"
#include "slininterp.h"
#include "integ.h"  // include pour les integrateurs de SOPHYA

//------ boost includes
#include <boost/math/special_functions/bessel.hpp>



//-------------------------------------------------------------------------
//      ------------------ MAIN PROGRAM ------------------------------
//-------------------------------------------------------------------------


class CosmoCoord {
public:
  CosmoCoord(SimpleUniverse &su);
  inline double getLOS(double z) { return isu_.LineOfSightComovDistanceMpc(z); } 
  SimpleUniverse su_;
  InterpCosmoCalc isu_;
  SLinInterp1D  iisu_;
};

CosmoCoord::CosmoCoord(SimpleUniverse &su)
  : su_(su), isu_(su_)
{
  size_t npts=1000;
  vector<double> leslos(npts);
  vector<double> lesz(npts);
  for(size_t i=0; i<npts; i++) {
    lesz[i]=i*10./(double)npts;
    leslos[i]=isu_.LineOfSightComovDistanceMpc(lesz[i]);
  }
  iisu_.DefinePoints(leslos, lesz);
}

class PkApprox : public ClassFunc1D {
public:
  PkApprox(double mu=log(0.1), double sigma=1.4, double A=24000.)
    :mu_(mu), sigma_(sigma), A_(A) { }
  virtual double operator()(double k) const {
    double arg = (-mu_ + log(k))/sigma_;
    arg *= arg;
    arg *= 0.5;
    if(arg>200.) return 0.;
    return A_*exp(-arg)/(sqrt(2.*M_PI)*k*sigma_);
  }
  double mu_, sigma_, A_;
};


class WGauss : public ClassFunc1D {
public:
  WGauss(double z0, double sigma, double A=1.)
    :z0_(z0), sigma_(sigma), A_(A) { }
  virtual double operator()(double z)  const {
    // on code une gaussienne
    z=(z-z0_)/sigma_;
    return A_*exp(-0.5*z*z); 
  }
  double z0_, sigma_, A_;
};

class WTopHat : public ClassFunc1D {
public:
  WTopHat(double zmin, double zmax, double A=1.)
    : zmin_(zmin), zmax_(zmax), A_(A) { }
  virtual double operator()(double z)  const {
    // on code une fenetre
    if(z > zmin_ and z < zmax_) {
      return A_;
    } else {
      return 0.0;
    }
  }
  double zmin_, zmax_, A_;
};

class ClIntegrandNoW : public ClassFunc1D {
public:
  ClIntegrandNoW(int l, double r, SLinInterp1D Pk)
    : l_(l), r_(r), Pk_(Pk) { }
  virtual double operator()(double k)  const {
    // d'après 1307.1459 equation 2.16
    return 2*k*k*Pk_(k)*pow(boost::math::sph_bessel(l_,r_*k),2)/Pi;
  }
  double l_, r_;
  SLinInterp1D Pk_;
};

class ClIntegrandNoWPkApprox : public ClassFunc1D {
public:
  ClIntegrandNoWPkApprox(int l, double r, PkApprox PkApprox)
    : l_(l), r_(r), PkApprox_(PkApprox) { }
  virtual double operator()(double k)  const {
    // d'après 1307.1459 equation 2.16
    return 2*k*k*PkApprox_(k)*pow(boost::math::sph_bessel(l_,r_*k),2)/Pi;
  }
  double l_, r_;
  PkApprox PkApprox_;
};


int main(int narg, const char* arg[])
{
  if ((narg<3)||((narg>1)&&(strcmp(arg[1],"-h")==0))) {
    cout << " Usage: pk2cl [options] InputClFile Out_PPFFile \n"
	 << "   InputClFile : text file  with input C(l) file \n"
	 << " options: [-N nsim] [-d dx,dy,dz] [-k kmin,kmax,nbin] [-k2d nbin2d,kmax2d] \n"
	 << "          [-s sigmaz] [-t lowval] [-r rfac] [-p lev] \n"
      	 << "   -N Nx,Ny,Nz : define input 3D grid size default=400x400x400 \n"
	 << endl;
    return 1;
  }
  Timer tm("pk2cl");
  int rc = 0;
  try { 
    string inpkname="Pk_z_1_0.txt";
    string outppfname="clout.ppf";
	
    int prtlev=0;
    int nsim = 300;
    //----------------------------------------------------------
    // decodage arguments optionnel 
    bool fgoptarg=true;
    while (fgoptarg&&(narg>1)) {
      string fbo = arg[1];
      if (fbo=="-N")  {   // specification nombre de générations de Cl sur lesquelles moyenner
	if (narg<3) { cout << " pk2cl/missing/bad argument, pk2cl -h for help " << endl;  return 2; }
	sscanf(arg[2],"%d",&nsim);
	arg+=2; narg-=2; 
      }
      else if (fbo=="-d")  {   // specification taille de cellules (Mpc) drho/rho ou ngal 
	if (narg<3) { cout << " pk2cl/missing/bad argument, pk2cl -h for help " << endl;  return 2; }
	double dx, dy, dz;
	sscanf(arg[2],"%lf,%lf,%lf",&dx,&dy,&dz); 
	arg+=2; narg-=2; 
      }
      else fgoptarg=false;
    }
    //----------------------------------------------------------
    if (narg<3) { cout << " pk2cl/missing/bad argument, pk2cl -h for help " << endl;  return 2; }

    inpkname=arg[1];
    outppfname=arg[2];

    cout << "pk2cl[1] : Opening output C(l) file " << outppfname <<endl;
    POutPersist pof(outppfname);
    
    cout << "pk2cl[2] : Reading input P(k) file " << inpkname <<endl;
    TMatrix<r_8> pkin;
    sa_size_t nr, nc;
    ifstream is(inpkname.c_str());
    pkin.ReadASCII (is, nr, nc);
    double kmin=pkin(0,0);
    double kmax=pkin(nr-1,0);
    kmin = 0.0001;
    kmax = 10;
    //kmax=5;
    cout << " kmin from file = " << kmin<<endl; 
    cout << " kmax from file = " << kmax <<endl; 
    pkin.Print(cout);
    
    pof<<PPFNameTag("pkin")<<pkin;
    cout << "pk2cl[3] : kept pkin ..."<<endl;

    // Define cosmology and P(k)
    SimpleUniverse su;
    InterpCosmoCalc isu(su);
    CosmoCoord csu(su);
    Vector ks = pkin.Column(0);
    Vector Pks = pkin.Column(1);
    vector<double> vx = ks.ConvertTostdvec();
    vector<double> vy = Pks.ConvertTostdvec();
    //SLinInterp1D Pk(vx,vy);
    double redshift = 1.0;
    cout << "pk2cl[4] : Cosmology and P(k) function defined ..."<<endl;
    PkApprox Pk(log(0.1),1.4,24000);
    cout << "pk2cl[4.2] : Approximate analytic P(k) function defined ..."<<endl;

    // On calcule l'integrale par la methode des Trapezes
    double LOS = csu.getLOS(redshift);
    LOS = 3600;
    int npt=(int)((kmax-kmin)*LOS*2/Pi); // 2 points per period (enough)
    TrpzInteg trpinteg;
    trpinteg.NStep(npt);
    cout << "Trapezoidal integral with npts="<<npt<<endl;
    int lmax = 512;
    Vector cl(lmax);
    for(int l=0; l<lmax; l++) {
      trpinteg.NStep(npt);
      //ClIntegrandNoW f(l,LOS,Pk);
      ClIntegrandNoWPkApprox f(l,LOS,Pk);
      trpinteg.SetFunc(f);
      cl(l) = trpinteg.ValueBetween(kmin, kmax);
      /*
	double kmin_l = max(0.95*l/LOS,kmin);
      double kmax_l = min(kmax,200*l/LOS);
      if(l<40) {
	kmin_l = kmin; kmax_l = kmax;
      }
      int npt_l=(int)((kmax_l-kmin_l)*LOS*2/Pi); // 2 points per period (enough)
      trpinteg.NStep(npt_l);
      double cl_test = 0;//trpinteg.ValueBetween(kmin_l, kmax_l);
      */
      cout << " ---- C(l="<<l<<",z="<<redshift<<") = "<<cl(l)<<endl;//" "<<cl_test<<" "<<(cl(l)-cl_test)/cl(l)<<" "<<kmin_l<<" "<<kmax_l<<" "<<npt_l<<endl;
    }
    pof<<PPFNameTag("cl")<<cl;
    GLInteg glinteg;
    int order = 4000;
    glinteg.SetOrder(order);
    cout << "Gauss-Legendre integral with order="<<order<<endl;
    Vector cl2(lmax);
    for(int l=0; l<lmax; l++) {
      //ClIntegrandNoW f(l,LOS,Pk);
      ClIntegrandNoWPkApprox f(l,LOS,Pk);
      glinteg.SetFunc(f);
      cl2(l) = glinteg.ValueBetween(kmin, kmax);
      cout << " ---- C(l="<<l<<",z="<<redshift<<") = "<<cl2(l)<<endl;
    }
    pof<<PPFNameTag("cl2")<<cl2;
    cout << "pk2cl[4] : kept cl ..."<<endl;
    
  }  // End of try bloc 
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " pk2cl.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " pk2cl.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " pk2cl.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
  }
  cout << " ==== End of pk2cl.cc program  Rc= " << rc << endl;
  return rc;    
}
