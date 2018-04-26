#ifndef JREGRID_H_SEEN
#define JREGRID_H_SEEN
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
#include "fitsioserver.h"
#include "fiosinit.h"

#include "ctimer.h"
#include "interpcosmoc.h"
#include "spherethetaphi.h"

//---- simlss includes
#include "pkspectrum.h"

//---- this modules include files 
#include "gpkspec.h"

//--- Change this to r_8 if one needs to work with double precision arrays 
#define TF  r_4 


class CosmoCoord {
public:
  CosmoCoord(SimpleUniverse &su);
  void getAngleRedshift(double dX, double dY, double dZ, double & theta, double & phi, double & los, double & redshift);
  inline double getRedshift(double los) { return iisu_.YInterp(los); } // iisu_(los); 
  inline double getLOS(double z) { return isu_.LineOfSightComovDistanceMpc(z); } 
  SimpleUniverse su_;
  InterpCosmoCalc isu_;
  SLinInterp1D  iisu_;
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

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
/*class PowerSpec : public ClassFunc1D {
 public:
 PowerSpec(double zref=0, double kmin=0.001, double kmax=100) : zref_(zref), kmin_(kmin), kmax_(kmax) {}
 PowerSpec(string inpkname, double zref, double kmin, double kmax);
 PowerSpec(SimpleUniverse su, double zref, double kmin, double kmax, int npt);
  virtual ~PowerSpec(void) {};
  double Max() const {return pmax_;}
  double kPeak() const {return kpeak_;}
  double Inverse(double pkval);
  virtual void   SetZ(double z) {zref_ = z;}
  virtual double GetZ(void) {return zref_;}
  virtual double operator()(double k) const { return Pk_(k); }
  virtual SLinInterp1D& GetPknosc() { return Pknosc_; }
  virtual SLinInterp1D& GetPk() { return Pk_; }
  virtual NTuple& GetNt()  { return nt_; }
 protected:
  SLinInterp1D Pk_;
  double kmin_, kmax_, kpeak_;
  double zref_;
  double pmin_, pmax_;
  double h100, Om0, Ob0, Ol0, w0;
  double T_CMB;
  double ns, as;
  SLinInterp1D Pknosc_;
  SimpleUniverse su_;
  double lkmin_, lkmax_, dlk_;
  int npt_;
  double scale_;
  NTuple nt_;
};
*/


class PowerSpecBase : public SLinInterp1D {
 public:
  PowerSpecBase(double zref=0, double kmin=0.001, double kmax=100);
  virtual ~PowerSpecBase(); 
  double Max() const {return pmax_;}
  double kPeak() const {return kpeak_;}
  double Inverse(double pkval);
  virtual void   SetZ(double z) {zref_ = z;}
  virtual double GetZ(void) {return zref_;}
  virtual double operator()(double k) const { return Pk_(k); }
  virtual SLinInterp1D& GetPk() { return Pk_; }
  virtual NTuple& GetNt()  { return nt_; }
 protected:
  double kmin_, kmax_, kpeak_;
  SLinInterp1D Pk_;
  double zref_;
  double pmin_, pmax_;
  NTuple nt_;
};

/*! 
  Ici fonction de test pour P(k)
  mu = Log[0.1]; sigma = 1.4;
  PkApprox[x_] :=
  N[24000*Exp[-(-mu + Log[x])^2/(2*sigma^2)]/(Sqrt[2*Pi]*x*sigma),
  70]
*/
class PowerSpecToy : public PowerSpecBase {
 public:
  PowerSpecToy(double mu=0.1, double sigma=1.4, double norm=24000, double zref=0., double kmin=0.001, double kmax=100);
  virtual double operator()(double k) const {
    double arg = (-mu_ + log(k))/sigma_;
    arg *= arg;
    arg *= 0.5;
    if(arg>300.) return 0.;
    return norm_*exp(-arg)/(sqrt(2.*M_PI)*k*sigma_);
  }
 private:
  double mu_, sigma_, norm_;
};

/*
  Load a power spectrum from a two column k P(k) ascii file
*/
class PowerSpecFile : public PowerSpecBase {
 public:
  PowerSpecFile(string inpkname, double zref, double kmin, double kmax);
 protected:
  string inpkname_;
  double lkmin_, lkmax_, dlk_;
  int npt_;
};

/*
  Compute P(k) given a cosmological model.
  Set osc to false to get power spectra without oscillations. Default=true.
*/
class PowerSpecCalc : public PowerSpecBase {
 public:
  PowerSpecCalc(SimpleUniverse su, bool osc, double zref, double kmin, double kmax, int npt);
  virtual bool GetOsc() const { return osc_; }
 protected:
  bool osc_;
  double h100, Om0, Ob0, Ol0, w0;
  double T_CMB;
  double ns, as;
  SimpleUniverse su_;
  double lkmin_, lkmax_, dlk_;
  int npt_;
  double scale_;
};



//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
class ProjectGridOnShells {
public:
  // Constructor
  ProjectGridOnShells(TArray< TF > & mgrid, CosmoCoord csu, double zcenter);
  virtual ~ProjectGridOnShells() { } 

  // Give the central index of the axis from its length
  double getCentralIndex(size_t npts);

  // Set the grid cell size (in Mpc) 
  void SetGridCellSize(double dx=1., double dy=1., double dz=1., bool fgprt=false) { dx_=dx; dy_=dy;  dz_=dz; };

  // Define the print level 0,1,2 ...
  inline int SetPrtLevel(int lev=0, int prtmod=10) 
  { int olev=prtlev_; prtlev_=lev; prtmodulo_=prtmod; return olev; }

  // return the grid cell size 
  inline double getdX() const {  return dx_; } 
  inline double getdY() const {  return dy_; } 
  inline double getdZ() const {  return dz_; }
  
  // Set the number of thin shells of depth dz
  void SetNShells() {nshells_=mgal_grid_.Size(2); };
  void SetShells();
  void FillShells();
  // Return the spherical map shell index given a redshift or a distance
  size_t GetShellbyLOS(double los, size_t ishell_start);
  //  inline size_t GetShellbyLOS(double los) { return (los-ShellsLOS0_)/ShellsdLOS_; } 
  size_t GetShellbyRedshift(double redshift, size_t ishell_start) ; 

protected:
  // member attribute
  TArray< TF > mgal_grid_;          // drho/rho or galaxy count array (grid)
  double dx_, dy_, dz_;             // mgal_grid_ cell size
  double LOScenter_;                // LOS distance of the cube center
  double iX0_, iY0_, iZ0_;          // central index of the cube
  sa_size_t nshells_;                  // number of shells of size dz
  CosmoCoord csu_;                  // class to convert cube coordinates in cosmo coordinates
  std::vector< SphereThetaPhi<r_8> > Shells_;
  std::vector< SphereThetaPhi<r_8> > ShellsN_; // auxiliary maps
  std::vector<double> ShellsLOS_;
  std::vector<double> ShellsRedshifts_;
  // double ShellsLOS0_, ShellsdLOS_;
  int prtlev_;
  int prtmodulo_;
  double ShellBadVal;
};


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
class ReProjectGrid {
public:
  // Constructor
  ReProjectGrid(TArray< TF > & mgrid, CosmoCoord csu, double zcenter);
  virtual ~ReProjectGrid() { }
  void UpdateInputGrid(TArray< TF > & mgrid);

  // Give the central index of the axis from its length
  inline double getCentralIndex(size_t npts)
  {
    return 0.5*(double)npts; 
    /* equivalent a 
    if (npts%2 == 0) {
      return(npts/2);
    } else {
      return(npts/2+0.5);
    }
    */
  }

  // Set the grid cell size (in Mpc) 
  void SetInputGridCellSize(double dx=1., double dy=1., double dz=1.)
  { odx_=dx_=dx;  ody_=dy_=dy;  odz_=dz_=dz; };

  // Define output grid
  void DefineOutputGrid(sa_size_t nx, sa_size_t ny, sa_size_t nz, double dx=-1., double dy=-1., double dz=-1.);
  void DefineOutputGridRedshift();

  // Define the print level 0,1,2 ...
  inline int SetPrtLevel(int lev=0, int prtmod=10) 
  { int olev=prtlev_; prtlev_=lev; prtmodulo_=prtmod; return olev; }

  // return the grid cell size 
  inline double getdX() const {  return dx_; } 
  inline double getdY() const {  return dy_; } 
  inline double getdZ() const {  return dz_; }
  // return the Output grid cell size 
  inline double getOdX() const {  return odx_; } 
  inline double getOdY() const {  return ody_; } 
  inline double getOdZ() const {  return odz_; }

  inline TArray< TF >& getOutputGrid() { return outgrid_; }
  inline TArray< TF >& getOutputGridRedshift() { return outgridz_; }
  void FillOutputGrid();
  void FillOutputGridRedshift();
  TArray< TF > AddPhotoZ(TArray< TF > ingrid, double phoz_err=0.02);
  // void FillOutputGridLOS();
 /*
  // Return the spherical map shell index given a redshift or a distance
  size_t GetShellbyLOS(double los, size_t ishell_start);
  //  inline size_t GetShellbyLOS(double los) { return (los-ShellsLOS0_)/ShellsdLOS_; } 
  size_t GetShellbyRedshift(double redshift, size_t ishell_start) ; 
  */
protected:
  // member attribute
  TArray< TF > mgal_grid_;          // drho/rho or galaxy count array (grid)
  double dx_, dy_, dz_;             // mgal_grid_ cell size
  //double dtheta_, dphi_, dzredshift_;    // mgal_grid cell size
  double LOScenter_;                // LOS distance of the cube center
  double zcenter_, zmin_, zmax_;       // Redshift cube parameters
  double iX0_, iY0_, iZ0_;          // central index of the cube
  TArray< TF > outgrid_;          //   Reprojected grid
  TArray< TF > outgridN_;          //  cell count in reprojected grid
  TArray< TF > outgridz_;          //   Reprojected redshift grid
  TArray< TF > outgridzN_;          //  cell count in redshift grid
  double odx_, ody_, odz_;             // reprojected grid cell size
  double odzredshift_;               // reprojected grid cell size in redshift
  double ioX0_, ioY0_, ioZ0_;          // central index of the reprojected cube
  double iozX0_, iozY0_, iozZ0_;          // central index of the redshift cube
  CosmoCoord csu_;                  // class to convert cube coordinates in cosmo coordinates
  int prtlev_;
  int prtmodulo_;
};


#endif /*  fin de JREGRID_H_SEEN */
