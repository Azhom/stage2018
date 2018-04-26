/* ----
   Project   LSST/BAO/PhotoZ
   Class to compute power spectrum from from mass density of 
   galaxy count number grid (array) 
   R.Ansari for A. Choyer , Feb 2015 
                                                     -------  */

#ifndef GPKSPEC_SEEN
#define GPKSPEC_SEEN

#include "machdefs.h"      
#include "sopnamsp.h"       
#include <math.h>
#include <iostream>
#include <vector>
#include <string> 

#include "classfunc.h"  
#include "array.h"         
#include "histats.h"       
#include "fftwserver.h"    
#include "randinterf.h"      
#include <rotation3d.h>

#include "gpkutil.h"
#include "corfunc.h"


//--- Change this to r_8 if one needs to work with double precision arrays 
#define TF  r_4 

// -- GFour3DPk class :  3D fourier amplitudes and power spectrum 
class GFour3DPk {
public:
// Constructor
  GFour3DPk(TArray< TF > & mgrid);
  virtual ~GFour3DPk(); 

  // Defining a specific RandomGenerator to be used 
  inline void SetRandomGenerator(RandomGeneratorInterface* rgp)
  { rgp_=rgp; }
  
  // Set the grid cell size (in Mpc) 
  void SetGridCellSize(double dx=1., double dy=1., double dz=1., bool fgprt=false);

  // Define the print level 0,1,2 ...
  inline int SetPrtLevel(int lev=0, int prtmod=10) 
  { int olev=prtlev_; prtlev_=lev; prtmodulo_=prtmod; return olev; }

// Return the total number of modes , and modes along each axis  
  inline sa_size_t TotN_k() const { return fourAmp.Size(); }
  inline sa_size_t N_kX()   const { return fourAmp.SizeX(); }
  inline sa_size_t N_kY()   const { return fourAmp.SizeY(); }
  inline sa_size_t N_kZ()   const { return fourAmp.SizeZ(); }

  // return the cell/step of Fourier modes 
  inline double getdkX() const {  return dkx_; } 
  inline double getdkY() const {  return dky_; } 
  inline double getdkZ() const {  return dkz_; } 

// Return the number of cells alons each direction  
  inline sa_size_t Nx()   const { return mgal_grid_.SizeX(); }
  inline sa_size_t Ny()   const { return mgal_grid_.SizeY(); }
  inline sa_size_t Nz()   const { return mgal_grid_.SizeZ(); }

// return the cell size along each direction
  inline double getdX() const {  return dx_; } 
  inline double getdY() const {  return dy_; } 
  inline double getdZ() const {  return dz_; } 

// Return the fourier amplitude matrix  
  TArray< complex<TF> > GetFourierAmp()  { return fourAmp; }

// Return the reconstructed power spectrum as a profile histogram   
  HProf ComputePk(int nbin=100, double kmin=0., double kmax=-1., bool fgmodcnt=false);

  // Return the reconstructed power spectrum as a profile histogram   
  Histo2D ComputePk2D(int nbin=100, double kmax=-1.);

  // Fills a data table from the computed P(k) profile histogram and mode count 
  Histo FillPkDataTable(DataTable& dt, double rfac=1.);
  inline HProf& getPk() { return *hp_pk_p_; }
  // un peu dangereux , l'histo peut ne pas etre alloue 
  inline Histo& getModCount() { return *hmcnt_p_; }

  // perform the FFT of the input array, it is called by ComputePk() if not called already
  void doFFT();
  // perform the inverse FFT, computing the mass grid array from Fourier coefficients
  //  if fgNOk0 true, force (0,0) fourier coefficient to zero
  void doInverseFFT(bool fgNOk0=true);
  // Generate Fourier coefficients according to an isotropic power spectrum
  void generateFourierAmp(ClassFunc1D & pk, double sigma_z=0., double pkfac=1.);
  
  // Changes mass/ngal grid cells < thr to thr 
  size_t CleanNegatives(TF thr=-1.);
  // Changes mass grid cells d rho/rho < -0.9 to exp(-3.5 *x*x*x*x) 
  size_t SmoothCleanNegatives();

  // Convert delta rho/rho to galaxy number count 
  void ConvertToGalaxyDensity(double galdens, bool fgpoiss, bool fgrenorm=true);

  void  HisPkCumul();

protected:

  // member attribute
  RandomGeneratorInterface* rgp_;
  TArray< TF > mgal_grid_;          // drho/rho or galaxy count array (grid)
  TArray< complex<TF> > fourAmp;    // complex array of fourier coefficients
  double dx_, dy_, dz_;       // mgal_grid_ cell size
  double dkx_, dky_, dkz_;    // fourier mode steps 
  int prtlev_;
  int prtmodulo_;
  // Profile histograms for power spectrum and number of modes 
  HProf* hp_pk_p_;
  Histo* hmcnt_p_;
  Histo* hmcntok_p_;
  // 2D histograms for 2D power spectrum and number of modes 
  Histo2D* hpk2_p_;
  Histo2D* h2mcnt_p_;
  // double s2cut_;    for later
};

// -- ReprojGrid class: Reprojecting an input grid to an output grid  
class ReprojGrid {
public:
// Constructor
  ReprojGrid(GridDef& ing, GridDef& outg, TArray< TF > & ingrid); 
  //	     double incenterdist, sa_size_t onx, sa_size_t ony, sa_size_t onz, double outcenterdist=-1.);
  virtual ~ReprojGrid();
  
  // Project the input grid to the output grid
  void Project();
  void ProjectSameCenter();

  // Project the input grid to the output grid, through galaxies distributed in the input grid cells 
  void ProjectWithGalaxies(double galdens, double sigmaR=-1., bool fgpoiss=false, SLinInterp1D* selfuncp=NULL);
  
  // Project the input grid to the output grid, assuming spherical geometry for the output grid
  void ProjectSpherical(bool fgisoangle=true);
  // Project the input grid to the output grid, through galaxies distributed in the input grid cells 
  void ProjectSphericalWithGalaxies(double galdens, bool fgisoangle=true, bool fgpoiss=false);

  // For debugging - creates and saves an NTuple with R,theta,phi coordinates of all input array cells after rotation
  void Grid2RThetaPhiNTuple(string& ppfname);
  
  // Define the print level 0,1,2 ...
  inline int SetPrtLevel(int lev=0, int prtmod=10) 
  { int olev=prtlev_; prtlev_=lev; prtmodulo_=prtmod; return olev; }

  // Access to input / output grids
  inline TArray< TF > & getInGrid()  { return in_grid; }
  inline TArray< TF > & getOutGrid()  { return out_grid; }

  // Access to input / output grids
  inline GridDef & getInGridDef()  { return ing_; }
  inline GridDef & getOutGridDeg()  { return outg_; }

  // Return the number of cells alons each direction  
  inline sa_size_t InNx()   const { return in_grid.SizeX(); }
  inline sa_size_t InNy()   const { return in_grid.SizeY(); }
  inline sa_size_t InNz()   const { return in_grid.SizeZ(); }

  // return the cell size along each direction
  inline double getIndX() const {  return ing_.dx; } 
  inline double getIndY() const {  return ing_.dy; } 
  inline double getIndZ() const {  return ing_.dz; } 

  // Return the number of cells alons each direction  
  inline sa_size_t OutNx()   const { return out_grid.SizeX(); }
  inline sa_size_t OutNy()   const { return out_grid.SizeY(); }
  inline sa_size_t OutNz()   const { return out_grid.SizeZ(); }

  // return the cell size along each direction
  inline double getOutdX() const {  return outg_.dx; } 
  inline double getOutdY() const {  return outg_.dy; } 
  inline double getOutdZ() const {  return outg_.dz; } 

  inline Vector3d In_XYZ_RThetaPhi(double x, double y, double z, Vector3d& v)
  {
    z+=ing_.r_center;
    v=Vector3d(x,y,z);
    return in_rot.RotateBack(v);
  }
  inline bool In_XYZ_2OutIJK_NoRot(double x, double y, double z, sa_size_t& ix, sa_size_t& jy, sa_size_t& kz)
  {
    ix=(sa_size_t)floor(x/outg_.dx+ocx_);
    jy=(sa_size_t)floor(y/outg_.dy+ocy_);
    kz=(sa_size_t)floor((z+ing_.r_center-outg_.r_center)/outg_.dz+ocz_);
    if ((ix<0) || (ix>=outg_.Nx) || (jy<0) || (jy>=outg_.Ny) || (kz<0) || (kz>=outg_.Nz))  return false;
    return true;
  }

  inline bool In_XYZ_2OutIJK(double x, double y, double z, sa_size_t& ix, sa_size_t& jy, sa_size_t& kz)
  {
    z+=ing_.r_center;
    Vector3d vout=out_rot.Rotate( in_rot.RotateBack(Vector3d(x,y,z)) );
    ix=(sa_size_t)floor(vout.X()/outg_.dx+ocx_);
    jy=(sa_size_t)floor(vout.Y()/outg_.dy+ocy_);
    kz=(sa_size_t)floor((vout.Z()-outg_.r_center)/outg_.dz+ocz_);
    if ((ix<0) || (ix>=outg_.Nx) || (jy<0) || (jy>=outg_.Ny) || (kz<0) || (kz>=outg_.Nz))  return false;
    return true;
  }
  inline bool In_XYZ_2OutIJK_SR(RandomGeneratorInterface& rg, double sigma_R, double x, double y, double z, 
				sa_size_t& ix, sa_size_t& jy, sa_size_t& kz)
  {
    z+=ing_.r_center;
    double r=sqrt(x*x+y*y+z*z);
    double sc=(r+rg.Gaussian(sigma_R))/r;
    x *= sc;   y *= sc;   z *= sc;
    Vector3d vout=out_rot.Rotate( in_rot.RotateBack(Vector3d(x,y,z)) );
    ix=(sa_size_t)floor(vout.X()/outg_.dx+ocx_);
    jy=(sa_size_t)floor(vout.Y()/outg_.dy+ocy_);
    kz=(sa_size_t)floor((vout.Z()-outg_.r_center)/outg_.dz+ocz_);
    if ((ix<0) || (ix>=outg_.Nx) || (jy<0) || (jy>=outg_.Ny) || (kz<0) || (kz>=outg_.Nz))  return false;
    return true;
  }
  
  protected:
  GridDef& ing_;              // input grid definition 
  GridDef& outg_;              // input grid definition 
  TArray< TF > in_grid;       // input grid 
  TArray< TF > out_grid;       // output grid

  EulerRotation3D in_rot;
  EulerRotation3D out_rot;
  double icx_, icy_, icz_;  // normalized coordinates (cell size=1) of the central cell for the input grid 
  double ocx_, ocy_, ocz_;  // normalized coordinates (cell size=1) of the central cell for the output grid 
  
  int prtlev_;
  int prtmodulo_;
};

#endif 
