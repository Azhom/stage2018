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
void CosmoCoord::getAngleRedshift(double dX, double dY, double dZ, double & theta, double & phi, double & los, double & redshift)
{
  los = sqrt(dX*dX+dY*dY+dZ*dZ);
  theta=atan2(sqrt(dX*dX+dY*dY),dZ);
  phi=atan2(dX,dY);
  redshift=getRedshift(los);
}


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

  // Give the central index of the axis from its length
  double getCentralIndex(size_t npts);

  // Set the grid cell size (in Mpc) 
  void SetInputGridCellSize(double dx=1., double dy=1., double dz=1.)
  { odx_=dx_=dx;  ody_=dy_=dy;  odz_=dz_=dz; };

  // Define output grid
  void DefineOutputGrid(sa_size_t nx, sa_size_t ny, sa_size_t nz, double dx=-1., double dy=-1., double dz=-1.)
  {
    sa_size_t siz[3];  siz[0]=nx; siz[1]=ny; siz[2]=nz;
    outgrid_.SetSize(3,siz);
    if (dx>0.)  odx_=dx;
    if (dy>0.)  ody_=dy;
    if (dz>0.)  odz_=dz; 
  }

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
  
  void FillOutputGrid();
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
  TArray< TF > outgrid_;          //   Reprojected grid
  TArray< TF > outgridN_;          //  cell count in reprojected grid
  double odx_, ody_, odz_;             // reprojected grid cell size
  double LOScenter_;                // LOS distance of the cube center
  double iX0_, iY0_, iZ0_;          // central index of the cube
  CosmoCoord csu_;                  // class to convert cube coordinates in cosmo coordinates
  int prtlev_;
  int prtmodulo_;
};

