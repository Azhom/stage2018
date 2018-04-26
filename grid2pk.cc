/*  ------------------------ Projet BAO/PhotoZ/LSST -------------------- 
  Programme de calcul du spectre de puissance (3D) a partir d'un 
  cube de donnees (delta rho/rho ou NGal )
    R. Ansari (pour Adline Choyer) - Feb 2015 
  Usage: grid2pkgrid2pk [options] In3DMap_FitsName OutPk_TextFile [OutPk_PPFFile]
         Options: [-d dx,dy,dz] [-k kmin,kmax,nbin] [-t lowval] [-r rfac] [-p lev]
---------------------------------------------------------------  */

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
#include "randfmt.h"

#include "ctimer.h"

//---- this modules include files 
#include "gpkspec.h"

//-------------------------------------------------------------------------
//      ------------------ MAIN PROGRAM ------------------------------
//-------------------------------------------------------------------------
int main(int narg, const char* arg[])
{
  if ((narg<3)||((narg>1)&&(strcmp(arg[1],"-h")==0))) {
    cout << " Usage: grid2pk [options] In3DMap_FitsName OutPk_TextFile [OutPk_PPFFile] \n"
	 << " options: [-in3dmap] [-zgrid] [-d dx,dy,dz] [-N Nx,Ny,Nz] [-k kmin,kmax,nbin]  \n"
	 << "          [-t lowval] [-sneg ] [-r rfac] [-p gal_density] [-rgi] [-prt lev] \n"
	 << "   -N Nx,Ny,Nz : define input 3D grid size default=400x400x400 \n"
	 << "   -d dx,dy,dz : define input 3D grid cell size (in Mpc) default=5x5x5 Mpc\n"
	 << "   -zgrid: Start with a grid (Nx*Ny*Nz) put to zero (evaluate white noise) \n"
	 << "   -k nbin,kmin,kmax : define number of bins and k-range for P(k) computation \n"
	 << "   -t lowval : set input map cells with val<lowval to lowval before computing P(k) \n"
      	 << "   -sneg: set input map cells with val<-0.9 smoothly to zero (exp(-3.5 x^4)) before computing P(k) \n"
	 << "   -r rfac : P(k) renormalisation factor \n"
      	 << "   -p gal_density : Generate Poisson randoms with the specified cell density  \n"
      	 << "   -rgi : Auto_initialize Random number generator \n"
	 << "   -prt lev : define print level (0,1,2..) \n" << endl;
    return 1;
  }
  Timer tm("grid2pk");
  int rc = 0;
  try { 
    string inpfitsfile;
    string outtextname;
    string outppfname;

    // Taille du cube 
    long int Nx,Ny,Nz;
    Nx=Ny=Nz=400;
    bool fgsetcellsize=false;
    double dx=5.,dy=5.,dz=5.;
    int nkbin=100;
    double kmin=0.001, kmax=1.001;
    bool fgtrunclowval=false;
    double lowval=-9.e19;
    bool fgsmneg=false;
    bool fgpoisson=false;
    double galdensity=1.;
    bool fgrfac=false;
    double rfac=1.;
    bool fgzerogrid=false;
    bool fgrgautoinit=false;
    
    int prtlev=0;
    //----------------------------------------------------------
    // decodage arguments optionnel 
    bool fgoptarg=true;
    while (fgoptarg&&(narg>1)) {
      string fbo = arg[1];
      if (fbo=="-zgrid")  { 
	if (narg<2) { cout << " grid2pk/missing/bad argument, grid2pk -h for help " << endl;  return 2; }
	fgzerogrid=true; arg++; narg--; 
      }
      else if (fbo=="-d")  {   // specification taille de cellules (Mpc) drho/rho ou ngal 
	if (narg<3) { cout << " grid2pk/missing/bad argument, grid2pk -h for help " << endl;  return 2; }
	sscanf(arg[2],"%lf,%lf,%lf",&dx,&dy,&dz); 
	fgsetcellsize=true; arg+=2; narg-=2; 
      }
      else if (fbo=="-N")  {   // specification taille de cellules (Mpc) drho/rho ou ngal 
	if (narg<3) { cout << " grid2pk/missing/bad argument, tpk2d -h for help " << endl;  return 2; }
	sscanf(arg[2],"%ld,%ld,%ld",&Nx,&Ny,&Nz);   arg+=2; narg-=2; 
      }
      else if (fbo=="-k")  {   // specification taille de cellules (Mpc) drho/rho ou ngal 
	if (narg<3) { cout << " grid2pk/missing/bad argument, grid2pk -h for help " << endl;  return 2; }
	sscanf(arg[2],"%d,%lf,%lf",&nkbin,&kmin,&kmax);  arg+=2; narg-=2;
      }
      else if (fbo=="-t")  { 
	if (narg<3) { cout << " grid2pk/missing/bad argument, grid2pk -h for help " << endl;  return 2; }
	lowval=atof(arg[2]);  fgtrunclowval=true; arg+=2; narg-=2; 
      }
      else if (fbo=="-sneg")  { 
	if (narg<2) { cout << " grid2pk/missing/bad argument, grid2pk -h for help " << endl;  return 2; }
	fgsmneg=true; arg++; narg--; 
      }
      else if (fbo=="-r")  { 
	if (narg<3) { cout << " grid2pk/missing/bad argument, grid2pk -h for help " << endl;  return 2; }
	rfac=atof(arg[2]);  fgrfac=true; arg+=2; narg-=2; 
      }
      else if (fbo=="-p")  { 
	if (narg<3) { cout << " grid2pk/missing/bad argument, grid2pk -h for help " << endl;  return 2; }
	galdensity=atof(arg[2]);  fgpoisson=true; arg+=2; narg-=2; 
      }
      else if (fbo=="-rgi")  { 
	if (narg<2) { cout << " grid2pk/missing/bad argument, grid2pk -h for help " << endl;  return 2; }
	fgrgautoinit=true; arg++; narg--; 
      }
      else if (fbo=="-prt")  { 
	if (narg<3) { cout << " grid2pk/missing/bad argument, grid2pk -h for help " << endl;  return 2; }
	prtlev=atoi(arg[2]);  fgrfac=true; arg+=2; narg-=2; 
      }
      else fgoptarg=false;
    }

    //----------------------------------------------------------
    if (narg<3) { cout << " grid2pk/missing/bad argument, grid2pk -h for help " << endl;  return 2; }

    if (fgsmneg) fgtrunclowval=false;
    inpfitsfile=arg[1];
    outtextname=arg[2];
    if (narg>3) outppfname=arg[3];

    FMTRandGen *RandGen = new FMTRandGen;
    RandGen->SelectGaussianAlgo(C_Gaussian_RandLibSNorm);
    if (fgrgautoinit)  RandGen->AutoInit(2);
    RandomGeneratorInterface::SetGlobalRandGenP(RandGen);

    TArray<r_4> ingrid;
    InterpPk* interpk_p=NULL;
    double mean, sigma;

    if (fgzerogrid) {
      cout << "tpk2d[1]  Setting input grid with size= " <<Nx<<" x "<<Ny<<" x "<<Nz<<" to zero"<<endl;
      sa_size_t sz[3]; sz[0]=Nx; sz[1]=Ny;  sz[2]=Nz;
      ingrid.ReSize(3,sz);
      ingrid = (TF)0.;
    }
    else {
      cout << "grid2pk[1] : reading 3D map from fits file " << inpfitsfile << endl;
      FitsInOutFile fin(inpfitsfile, FitsInOutFile::Fits_RO);
      fin >> ingrid;
      MeanSigma(ingrid, mean, sigma);
      cout << "grid2pk[1.b] Input grid sizes " << ingrid.InfoString() << endl;
      ingrid.Show(); 
      cout << "... Input grid Mean=" << mean << " Sigma=" << sigma << endl;
    }
    
    tm.Split(" After read/init input grid ");
    
    GFour3DPk  gpkc(ingrid);
    gpkc.SetPrtLevel(prtlev);
    gpkc.SetGridCellSize(dx,dy,dz);
      
    if (fgsmneg) {  // smoothly setting values below -0.9 to zero 
      cout << "grid2pk[2] : calling SmoothCleanNegatives() ... " << endl;    
      gpkc.SmoothCleanNegatives();
      MeanSigma(ingrid, mean, sigma);
      cout << "... After SmoothCleanNegatives grid Mean=" << mean << " Sigma=" << sigma << endl;
      tm.Split(" After SmoothCleanNegatives ");
    }
    else if (fgtrunclowval) {  // truncating grid value below threshold 
      cout << "grid2pk[2] : calling CleanNegatives() ... " << endl;    
      gpkc.CleanNegatives(lowval);
      MeanSigma(ingrid, mean, sigma);
      cout << "... After CleanNegatives grid Mean=" << mean << " Sigma=" << sigma << endl;
      tm.Split(" After CleanNegatives ");
    }
    if (fgpoisson) {  // replacing array with galaxy numbers generarted according to Poisson
      cout << "grid2pk[2.b] : calling ConvertToGalaxyDensity("<<galdensity<<") ... " << endl;    
      gpkc.ConvertToGalaxyDensity(galdensity,true,true);
      tm.Split(" After ConvertToGalaxyDensity ");
    }
    cout << "grid2pk[3] : Computing Fourier coefficients ... " << endl;    
    gpkc.doFFT();
    tm.Split(" After doFFT ");

    cout << "calcpk[4] : computing power spectrum ... " << endl;
    gpkc.ComputePk(nkbin,kmin,kmax,true);
    for(int ii=0; ii<50; ii++) {
      cout << " ----- Loop ii="<<ii<<endl;
      ingrid=(r_4)0.;
      gpkc.ConvertToGalaxyDensity(galdensity,true,true);
      gpkc.doFFT();
      gpkc.HisPkCumul();
    }
    HProf hpk = gpkc.getPk();
    DataTable dtpk; 
    Histo hrap = gpkc.FillPkDataTable(dtpk, rfac);
    tm.Split(" After ComputePk ");
    if (prtlev>0)  { 
      dtpk.SetShowMinMaxFlag(true);
      cout << dtpk; 
    }
    {   // Saving computed P(k) to text (ascii) file 
      cout << "calcpk[5] : wrting P(k) to text file " << outtextname << endl;
      ofstream tof(outtextname.c_str());
      dtpk.WriteASCII(tof);
    }
    if (outppfname.length()>0)  {   // Saving computed P(k) to PPF file 
      cout << "calcpk[6] : writing profile histo P(k) (hpk) and DataTable P(k) dtpk to PPF file  " << outppfname << endl;
      POutPersist po(outppfname);
      po<<PPFNameTag("hpk")<<hpk;
      po<<PPFNameTag("dtpk")<<dtpk;
    }
  }  // End of try bloc 
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " grid2pk.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " grid2pk.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " grid2pk.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
  }
  cout << " ==== End of grid2pk.cc program  Rc= " << rc << endl;
  return rc;    
}
