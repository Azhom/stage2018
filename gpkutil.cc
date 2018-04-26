
//----- c/c++ includes 
#include <iostream>
#include <fstream>
#include <string>

//----- sophya includes 
#include "machdefs.h"
#include "sopnamsp.h"

#include "histats.h"
#include "array.h"

#include "gpkutil.h"


GPkArgDecoder::GPkArgDecoder()
  : grid_(), ingrid_(), outgrid_(), hpkdef_(), gfkparm_(), fg_rand_autoinit_(false), prtlev(0), pkfitsfilename("")
{
  Nx=Ny=Nz=400;  dx=dy=dz=6.;  center=3000.;  t0=0.; p0=0.;
  inNx=inNy=inNz=401;  indx=indy=indz=6.; incenter=3000.;  in_t0=0.; in_p0=0.;
  outNx=outNy=outNz=251;  outdx=outdy=outdz=6.;  outcenter=3000.;  out_t0=0.; out_p0=0.;

  nkbin=200;  kmin=0.;   kmax=1.;
  nkbin2d=100;  kmax2d=1.;
  fgtrunclowval=false;  lowval=-1.;
  fgsmneg=false;  fglognormal=false;
  fgpoisson=false;  numberdensity=1.;
  sigmaz=0.;
  fgrfac=false;  rfac=1.;
  prtlev=0;   fg_rand_autoinit_=false;
  fgdebug=false;  debuglevel=0;
  NLoop=1; 
  fgnosc_=false;
  fgselfunc=false;
  fgshotnoiseonly=false;
  fgpkfits=false;   
}

int GPkArgDecoder::DecodeArgs(int narg, const char* arg[], int na, const char* msg)
{
  // Taille du cube
  Nx=Ny=Nz=400;
  dx=dy=dz=6.;
  center=3000.;
  t0=0.; p0=0.;
  inNx=inNy=inNz=401;
  indx=indy=indz=6.;
  incenter=3000.;
  in_t0=0.; in_p0=0.;
  outNx=outNy=outNz=251;
  outdx=outdy=outdz=6.;
  outcenter=3000.;
  out_t0=0.; out_p0=0.;

  nkbin=200;
  kmin=0.;   kmax=1.;
  nkbin2d=100;
  kmax2d=1.;
  fgtrunclowval=false;
  lowval=-1.;
  fgsmneg=false;  fglognormal=false;
  fgpoisson=false;
  numberdensity=1.;
  sigmaz=0.;
  fgrfac=false;
  rfac=1.;
  prtlev=0;
  fgdebug=false;  debuglevel=0;
  NLoop=1;
  fg_rand_autoinit_=false;
  fgnosc_=false;
  fgselfunc=false;
  fgshotnoiseonly=false;
  vector<GridCenter> vgc;
  
  // decodage arguments optionnel 
  bool fgoptarg=true;
  while (fgoptarg&&(narg>1)) {
    string fbo = arg[1];
    //DBG    cout << " **DBG**A / Decoding , narg="<<narg<<" fbo="<<fbo<<endl;
    if (fbo=="-grid")  {   // specification taille grille/cellules drho/rho ou ngal 
      if (narg<na+2) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1; }
      sscanf(arg[2],"%ld,%ld,%ld,%lg,%lg,%lg",&Nx,&Ny,&Nz, &dx,&dy,&dz);   arg+=2; narg-=2; 
    }
    else if (fbo=="-in")  {   // specification taille grille/cellules drho/rho ou ngal  grille input 
      if (narg<na+2) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1; }
      sscanf(arg[2],"%ld,%ld,%ld,%lg,%lg,%lg",&inNx,&inNy,&inNz, &indx,&indy,&indz);   arg+=2; narg-=2; 
    }
    else if (fbo=="-out")  {   // specification taille grille/cellules drho/rho ou ngal  grille output 
      if (narg<na+2) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1; }
      sscanf(arg[2],"%ld,%ld,%ld,%lg,%lg,%lg",&outNx,&outNy,&outNz, &outdx,&outdy,&outdz);   arg+=2; narg-=2; 
    }
    else if (fbo=="-center")  { 
      if (narg<na+2) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1; }
      sscanf(arg[2],"%lg,%lg,%lg",&center,&t0,&p0); arg+=2; narg-=2; 
    }
    else if (fbo=="-ic")  { 
      if (narg<na+2) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1; }
      sscanf(arg[2],"%lg,%lg,%lg",&incenter,&in_t0,&in_p0);  arg+=2; narg-=2;
       in_t0 *= (M_PI/180.);    in_p0 *= (M_PI/180.); 
    }
    else if (fbo=="-oc")  { 
      if (narg<na+2) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1; }
      sscanf(arg[2],"%lg,%lg,%lg",&outcenter,&out_t0,&out_p0); arg+=2; narg-=2;
      out_t0 *= (M_PI/180.);    out_p0 *= (M_PI/180.); 
    }
    else if (fbo=="-moc")  { 
      if (narg<na+2) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1; }
      GridCenter gc;
      sscanf(arg[2],"%lg,%lg,%lg",&(gc.r_center),&(gc.theta0),&(gc.phi0));   arg+=2; narg-=2;
      gc.theta0 *= (M_PI/180.);    gc.phi0 *= (M_PI/180.);
      vgc.push_back(gc);
    }
    else if (fbo=="-k")  {   // specification nb-bin kmin,kmax pour calcul P(k)  
      if (narg<na+2) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1; }
      sscanf(arg[2],"%d,%lf,%lf",&nkbin,&kmin,&kmax);   arg+=2; narg-=2;
    }
    else if (fbo=="-k2d")  {   // specification nb-bin kmax pour calcul 2D-P(k)  
      if (narg<na+2) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1; }
      sscanf(arg[2],"%d,%lf",&nkbin2d,&kmax2d);  arg+=2; narg-=2; 
    }
    else if (fbo=="-N")  { 
      if (narg<na+2) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1; }
      NLoop=atoi(arg[2]);  arg+=2; narg-=2; 
    }
    else if (fbo=="-thr")  { 
      if (narg<na+2) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1; }
      lowval=atof(arg[2]);  arg+=2; narg-=2; 
    }
    else if (fbo=="-neg")  { 
      if (narg<na+1) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1;  }
      fgtrunclowval=true; arg++; narg--;    fgsmneg=fglognormal=false;
    }
    else if (fbo=="-sneg")  { 
      if (narg<na+1) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1;  }
      fgsmneg=true; arg++; narg--;   fgtrunclowval=fglognormal=false;
    }
    else if (fbo=="-lognormal")  { 
      if (narg<na+1) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1;}
      fgsmneg=true; arg++; narg--;   fgtrunclowval=fgsmneg=false;
    }
    else if (fbo=="-s")  { 
      if (narg<na+2) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1; }
      sigmaz=atof(arg[2]);  arg+=2; narg-=2; 
    }
    else if (fbo=="-ngal")  { 
      if (narg<na+2) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1; }
      numberdensity=atof(arg[2]);  arg+=2;  narg-=2; 
    }
    else if (fbo=="-selfunc")  { 
      if (narg<na+2) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1; }
      ReadSelectionFunctionFile(string(arg[2]));  arg+=2;  narg-=2; 
    }      
    else if (fbo=="-poiss")  { 
      if (narg<na+1) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1; }
      fgpoisson=true; arg++;  narg--; 
    }      
    else if (fbo=="-snonly")  { 
      if (narg<na+1) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1; }
      fgshotnoiseonly=true; arg++;  narg--; 
    }      
    else if (fbo=="-r")  { 
      if (narg<na+2) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1; }
      rfac=atof(arg[2]);  fgrfac=true; arg+=2; narg-=2; 
    }
    else if (fbo=="-rgi")  { 
      if (narg<na+1) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1;}
      fg_rand_autoinit_=true; arg++; narg--; 
    }
    else if (fbo=="-prt")  { 
      if (narg<na+2) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1; }
      prtlev=atoi(arg[2]);  fgrfac=true; arg+=2; narg-=2; 
    }
    else if (fbo=="-dbg")  { 
      if (narg<na+2) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1; }
      debuglevel=atoi(arg[2]);  fgdebug=true; arg+=2; narg-=2; 
    }
    else if (fbo=="-nosc")  { 
      if (narg<na+1) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1; }
      fgnosc_=true; arg++; narg--; 
    }
    else if (fbo=="-spfits")  { 
      if (narg<na+2) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1; }
      fgpkfits=true;  pkfitsfilename=arg[2]; arg+=2; narg-=2; 
    }
    
    /*   DEL 
    else if (fbo=="-wnosc")  { 
      if (narg<na+1) { cout << msg << " missing/bad argument, -h for help " << endl;  return -1; }
      fgpkdoboth_=true; arg++; narg--; 
    }
    */
    else fgoptarg=false;
  }

  for(int k=1; k<narg; k++) lastargs.push_back(arg[k]);
  
  grid_.Set(Nx,Ny,Nz,dx,dy,dz,center,t0,p0);
  ingrid_.Set(inNx,inNy,inNz,indx,indy,indz,incenter,in_t0,in_p0);
  outgrid_.Set(outNx,outNy,outNz,outdx,outdy,outdz,outcenter,out_t0,out_p0);
  for(size_t jj=0; jj<vgc.size(); jj++) {
    GridDef cgrid(outgrid_);
    cgrid.SetCenter(vgc[jj]);
    voutgrids_.push_back(cgrid);
  }
  hpkdef_.Set(nkbin,kmin,kmax);
  gfkparm_.Set(sigmaz,fgtrunclowval,lowval,fgrfac,rfac,fgpoisson,numberdensity);
  gfkparm_.SetB(fgsmneg,fglognormal,fgshotnoiseonly);
  gfkparm_.SetPrintLevel(prtlev);
  return 0;
}


int GPkArgDecoder::UsageOptions()
{
  cout << " options: [-grid Nx,Ny,Nz,dx,dy,dz] [-in Nx,Ny,Nz,dx,dy,dz] [-out Nx,Ny,Nz,dx,dy,dz] \n"
       << "          [-center center_dist,tet0,phi0] [-ic center_dist,tet0,phi0]  \n"
       << "          [-oc center_dist,tet0,phi0] [-moc center_dist,tet0,phi0] \n"
       << "          [-k kmin,kmax,nbin] [-k2d nbin,kmax] \n"
       << "          [-neg] [-sneg] [-lognormal] [-ltn lowval] [-r rfac] \n"
       << "          [-s sigmaz] [-ngal ngal_per_cell] [-selfunc sfname] \n"
       << "          [-poiss] [-snonly] [-nosc] [-spfits fitsfilename] \n"
       << "          [-N nloop] [-thr NThreads] [-rgi] [-prt lev] [-dbg lev] \n"
       << "   -grid Nx,Ny,Nz,dx,dy,dz : define 3D grid  default=400x400x400,6,6,6  (6 MPc)\n"
       << "   -in Nx,Ny,Nz,dx,dy,dz : define input 3D grid  default=401x401x401,6,6,6  (6 MPc)\n"
       << "   -out Nx,Ny,Nz,dx,dy,dz : define out 3D grid  default=251x251x251,8,8,8 (8 Mpc)\n"
       << "   -center center_dist,theta0,phi0: define grid center distance(in Mpc) and direction default=3000,0.,0.\n"
       << "   -ic center_dist,theta0,phi0: define input grid center distance(in Mpc) and direction default=3000,0.,0.\n"
       << "   -oc center_dist,theta0,phi0: define output grid center distance(in Mpc) and direction default=3000,0.,0.\n"
       << "   -moc center_dist,theta0,phi0: definition of multiple output grids \n"
       << "   -k nbin,kmin,kmax : define number of bins and k-range for P(k) computation (def=200,0.,1.) \n"
       << "   -k2d nbin,kmax : define number of bins and k-max for 2D computation (def=100,1.) \n"
       << "   -ltn lowval : set input map cells with val<lowval to lowval before computing P(k) (def=-1)\n"
       << "   -neg: set input map cells with val<lowval to lowval (def=-1) before computing P(k) - default behaviour \n"
       << "   -sneg: set input map cells with val<-0.9 smoothly to -1 (exp(-3.5 x^4)) before computing P(k) \n"
       << "   -lognormal: replace input map cells with exp(cell_val) \n"
       << "   -r rfac : P(k) renormalisation factor \n"
       << "   -s sigmaz : gaussian smearing along z (in Mpc units) default=0.\n"
       << "   -selfunc sfname : read selection function from file sfname \n"
       << "   -ngal density : Specify galaxy density (number of galaxies / cell) \n"
       << "   -poiss : Generate Poisson randoms using the specified per cell galaxy density \n"
       << "   -snonly : shot-noise only (input d rho/rho grid = 0.) \n"
       << "   -nosc : use Pk-no-oscillation from the SimLSS output  \n"
       << "   -spfits filename : write individual computed spectra as a 3D array in FITS format to filename \n" 
       << "   -N nloop : number of different sky cubes generated (default=1) \n"
       << "   -thr NThreads: number of parallel threads used for computation \n"
       << "   -rgi : Auto_initialize Random number generator \n"
       << "   -dbg lev : activate debug flag and defines debug level \n"
       << "   -prt lev : define print level (0,1,2..) \n" 

       << endl;
//DEL       << "   -wnosc : use both Pk , with and without oscillation from the SimLSS output  \n"
  
  return 0;
}

int GPkArgDecoder::ReadSimLSSPkFile(string const & filename)
{
  // Read a file defining P(k) , P-noosc(k) and P-current(k) , produced by SimLSS
  // an text file with four space/tab separated values on each line: k  P(k)  P-noosc(k)
  //  lines with # as the first character are considered as comments and are ignored 
  cout<<"---GPkArgDecoder::ReadSimLSSPkFile("<<filename<<")"<<endl;
  const char * names[4] = {"k", "pk", "pknos","cepk"};
   // NTuple (Table) creation with 4 columns (double precision)
   NTuple  nt(4, names);
   nt.FillFromASCIIFile (filename);
   vector< r_8 > lesk, lespk, lespknos;
   nt.GetColumn (0, lesk);
   nt.GetColumn (1, lespk);
   nt.GetColumn (2, lespknos);
   fpk_.DefinePoints (lesk, lespk);
   fpknosc_.DefinePoints (lesk, lespknos);
   return 0;
}

int GPkArgDecoder::ReadSelectionFunctionFile(string const & filename)
{
  //  Read a file containing containing a selection function eta(r) where r is the radial distance in Mpc , 0<=eta<=1
  //  a text file with two space/tab separated values on each line:  r  eta(r)
  //  lines with # as the first character are considered as comments and are ignored 
  cout<<"---GPkArgDecoder::ReadSelectionFunctionFile("<<filename<<")"<<endl;
  const char * names[4] = {"r", "eta"};
   // NTuple (Table) creation with 2 columns (double precision)
   NTuple  nt(2, names);
   nt.FillFromASCIIFile (filename);
   vector< r_8 > lesr, leseta;
   nt.GetColumn (0, lesr);
   nt.GetColumn (1, leseta);
   selfunc_.DefinePoints (lesr, leseta);
   fgselfunc=true;
   return 0;
}
