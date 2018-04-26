/* ----
   Project   LSST/BAO/PhotoZ
   R.Ansari , Oct 2016 
                                                     -------  */

#ifndef GPKUTIL_SEEN
#define GPKUTIL_SEEN

#include <iostream>
#include <fstream>
#include <string>

#include <vector>

//----- sophya includes 
#include "machdefs.h"
#include "sopnamsp.h"
#include <slininterp.h>

using namespace SOPHYA;
using namespace std;

// Quelques classes utilitaires - pour transmettre les arguments, etc ...
class GridCenter {
public:
  GridCenter(double rc=3000., double t0=0., double p0=0.)
  { r_center=rc; theta0=t0;  phi0=p0;  }
  GridCenter(GridCenter const &a)
  { r_center=a.r_center; theta0=a.theta0;  phi0=a.phi0; }
  double r_center;
  double theta0, phi0;
};

class GridDef {
public:
  GridDef()
  { Nx=Ny=Nz=256; dx=dy=dz=6.; r_center=3000.; theta0=phi0=0.; }
  GridDef(long int nx, long int ny, long int nz, double ddx, double ddy, double ddz, double rc, double t0=0., double p0=0.)
  { Set(nx,ny,nz,ddx,ddy,ddz,rc,t0,p0); }
  void Set(long int nx, long int ny, long int nz, double ddx, double ddy, double ddz, double rc, double t0=0., double p0=0.)
  { Nx=nx; Ny=ny; Nz=nz; dx=ddx; dy=ddy; dz=ddz; r_center=rc; theta0=t0;  phi0=p0;  return; }
  GridDef(GridDef const& a)
  { Nx=a.Nx; Ny=a.Ny; Nz=a.Nz; dx=a.dx; dy=a.dy; dz=a.dz;  r_center=a.r_center; theta0=a.theta0;  phi0=a.phi0;  }
  GridDef& operator = (GridDef const& a)
  { Nx=a.Nx; Ny=a.Ny; Nz=a.Nz; dx=a.dx; dy=a.dy; dz=a.dz;  r_center=a.r_center; theta0=a.theta0;  phi0=a.phi0; return *this; }
  void SetCenter(GridCenter const &a)
  { r_center=a.r_center; theta0=a.theta0;  phi0=a.phi0;  }
  
  double GetCellVolume() const
  { return (dx*dy*dz); }
  double GetVolume() const
  { return (double)(Nx*Ny*Nz)*(dx*dy*dz); }
  
  long int Nx,Ny,Nz;
  double dx, dy, dz;
  double r_center;
  double theta0, phi0;
};
class HPkDef {
public:
  HPkDef()
  { nkbin=100; kmin=0.005; kmax=1.005; delta_k=0.01; }
  HPkDef(int nb, double min, double max)
  { Set(nb, min, max); }
  HPkDef(HPkDef const & a) :
    nkbin(a.nkbin), kmin(a.kmin), kmax(a.kmax)  { }
  HPkDef& operator= (HPkDef const & a)
  { nkbin=a.nkbin; kmin=a.kmin;  kmax=a.kmax;  return *this; }
  void Set(int nb, double min, double max) 
  { nkbin=nb; kmin=min; kmax=max;  delta_k=(kmax-kmin)/(double)nkbin; return; }
  double  Getk(int i) { return (kmin+((double)i+0.5)*delta_k); }
  int nkbin;
  double kmin, kmax;
  double delta_k;
};
class GFkParam {
public:
  GFkParam()
  {
    sigmaz_=0.; fgtrunclowval_=false; lowval_=-1.; fgrfac_=false; rfac_=1.; fgpoisson_=false; numberdensity_=1.;
    SetB();
    SetPrintLevel();
  }
  GFkParam(double sigmaz, bool fgtrunc, double lowval, bool fgrfac, double rfac, bool fgpoiss, double numberdens)
  {
    Set(sigmaz, fgtrunc, lowval, fgrfac, rfac, fgpoiss, numberdens);
    SetB();
    SetPrintLevel();
  }
  GFkParam(GFkParam const& a)
  {
    Set(a.sigmaz_, a.fgtrunclowval_, a.lowval_, a.fgrfac_, a.rfac_, a.fgpoisson_, a.numberdensity_); 
    SetB(a.fgsmneg_, a.fglognormal_, a.fgshotnoiseonly_);
    SetPrintLevel(a.prtlevel_);
  }
  GFkParam& operator=(GFkParam const& a)
  {
    Set(a.sigmaz_, a.fgtrunclowval_, a.lowval_, a.fgrfac_, a.rfac_, a.fgpoisson_, a.numberdensity_); 
    SetB(a.fgsmneg_, a.fglognormal_, a.fgshotnoiseonly_);
    SetPrintLevel(a.prtlevel_);
    return *this;
  }    
  void Set(double sigmaz, bool fgtrunc, double lowval, bool fgrfac, double rfac, bool fgpoiss, double numberdens)
  { sigmaz_=sigmaz; fgtrunclowval_=fgtrunc; lowval_=lowval; fgrfac_=fgrfac; rfac_=rfac; fgpoisson_=fgpoiss; numberdensity_=numberdens; }
  void SetB(bool fgsn=false, bool fgln=false, bool fgsno=false)
  { fgsmneg_=fgsn; fglognormal_=fgln; fgshotnoiseonly_=fgsno; }
  void SetPrintLevel(int lev=0) { prtlevel_=lev; }
  
  bool fgtrunclowval_, fgrfac_, fgpoisson_;
  double sigmaz_,lowval_,rfac_,numberdensity_;
  bool fgsmneg_, fglognormal_;
  bool fgshotnoiseonly_;
  int prtlevel_;
};

class GPkArgDecoder {
public:
  GPkArgDecoder();
  int DecodeArgs(int narg, const char* arg[], int na, const char* msg);
  int UsageOptions();
  int ReadSimLSSPkFile(string const & filename);
  int ReadSelectionFunctionFile(string const & filename);

  long int Nx,Ny,Nz;
  double dx,dy,dz;
  double center;
  double t0, p0;
  long int inNx,inNy,inNz;
  double indx,indy,indz;
  double incenter;
  double in_t0, in_p0;
  long int outNx,outNy,outNz;
  double outdx,outdy,outdz;
  double outcenter;
  double out_t0, out_p0;

  int nkbin;
  double kmin, kmax;
  int nkbin2d;
  double kmax2d;
  bool fgtrunclowval;
  double lowval;
  bool fgsmneg;
  bool fglognormal;
  bool fgpoisson;
  double numberdensity;
  double sigmaz;
  bool fgrfac;
  double rfac;
  int prtlev;
  bool fgdebug;
  int debuglevel;
  int NLoop; 
  bool fg_rand_autoinit_;
  vector<string> lastargs;

  GridDef grid_;
  GridDef ingrid_;
  GridDef outgrid_;
  vector<GridDef> voutgrids_; 
  HPkDef  hpkdef_; 
  GFkParam gfkparm_;
  
  SLinInterp1D fpk_;         // Fonction interpolee P(k) 
  SLinInterp1D fpknosc_;     // Fonction interpolee P(k)-No-Oscillation
  bool fgnosc_;  // true -> use fpknosc_, fpk_ if false
  
  SLinInterp1D selfunc_;     // Fonction de selection interpolee  eta(radial_distance_en_Mpc)
  bool fgselfunc;            // If true -> has selection function 
  bool fgshotnoiseonly;     // If true, compute for shot-noise only (input d rho/rho grid = 0)
  //  bool fgpkdoboth_;  // true -> perform computation for fpk_ and fpknosc_
  bool fgpkfits;         // true -> write individual computed P(k) as a 3D array to FITS file
  string pkfitsfilename;  // FITS filename for 3D array with individual computed P(k) 
};

#endif
