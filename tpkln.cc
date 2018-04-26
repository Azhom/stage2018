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
#include "randfmt.h"

//---- this modules include files
// #include "corfunc.h"
// #include "myinteg2d.h"
#include "gpkutil.h"
#include "hsplfit.h"
#include "gpkspec.h"

//---- Fonctions de ce fichier
vector<HProf> DoCompPk(ClassFunc1D& interpk, GridDef& grid, GFkParam& gfkparm, HPkDef& hpkparm, NTuple& ntms, int prtlev=0);
void VHP_Update(int k, HPkDef& hpkparm, vector<HProf>& vhpk, vector<HProf>& vh, POutPersist& pof);
DataTable HPk2Dt(vector<HProf>& vhpk, GridDef& grid, ClassFunc1D& interpk);

//-------------------------------------------------------------------------
//      ------------------ MAIN PROGRAM ------------------------------
//-------------------------------------------------------------------------
int main(int narg, const char* arg[])
{
  int rc = 0;
  try {
    GPkArgDecoder decoder;
    if ((narg<3)||((narg>1)&&(strcmp(arg[1],"-h")==0))) {
      cout << " Usage: Usage: tpkln [options] InputPkFile Out_PPFFile  \n"
	   << "   InputPkFile : text file with pairs of values k  Pk on each line \n"
	   << " options: [-N Nx,Ny,Nz] [-d dx,dy,dz] [-k kmin,kmax,nbin] [-k2d nbin2d,kmax2d] \n"
	   << "            [-neg] [-sneg] [-lognormal] [-thr lowval] [-r rfac] [-rgi] [-prt lev]"<<endl;
      return decoder.UsageOptions();
    }
    Timer tm("tpkln");
    string inpkname=".";
    string outppfname="toto.ppf";

    //----------------------------------------------------------
    decoder.DecodeArgs(narg, arg, 2, "tpkln");
    if (decoder.lastargs.size()<2)
      { cout << " tpkln/missing/bad argument, tpkln -h for help " << endl;  return 2; }

    inpkname=decoder.lastargs[0];
    outppfname=decoder.lastargs[1];

    FMTRandGen *RandGen = new FMTRandGen;
    RandGen->SelectGaussianAlgo(C_Gaussian_RandLibSNorm);
    if (decoder.fg_rand_autoinit_)  RandGen->AutoInit(2);
    RandomGeneratorInterface::SetGlobalRandGenP(RandGen);
 
    cout << "tpksph[1.a]: reading input power spectrum from file "<<inpkname<<endl;
    decoder.ReadSimLSSPkFile(inpkname);
    SLinInterp1D& interpk = decoder.fpk_;
    if (decoder.fgnosc_)  interpk = decoder.fpknosc_;

    cout << "tpkln[1.b]: Opening output PPF file " << outppfname <<endl;
    POutPersist pof(outppfname);
    string outppfvhpk=outppfname;
    size_t pdfn=outppfname.find_last_of('.');
    if (pdfn < outppfname.length()) outppfvhpk=outppfname.substr(0,pdfn);
    outppfvhpk+="_vhpk.ppf";
    cout << "tpkprj[1.c]: Opening output PPF file " << outppfvhpk <<" for individual P(k) HProf's" <<endl;
    POutPersist pofvhpk(outppfvhpk);
    /*
    if ((decoder.gfkparm_.fgtrunclowval_==false) && (decoder.gfkparm_.fgsmneg_==false) && (decoder.gfkparm_.fglognormal_==false))  {
      cout << " tpkln / Forcing TruncLowVal, d rho/rho < -1 -> -1 ..."<<endl;
      decoder.gfkparm_.fgtrunclowval_=true;
      decoder.gfkparm_.lowval_=-1.;
    }
    */
    //    cout << " *DBG***B** nkbin="<<decoder.hpkdef_.nkbin<<" kmin="<<decoder.hpkdef_.kmin<<" kmax="<<decoder.hpkdef_.kmax<<endl;

    const char * ntmsnames[4]={"gin_mean","gin_sigma","gsneg_mean","gsneg_sigma"};
    NTuple ntms(4,ntmsnames);
    
    vector<HProf>  vhpk;
    
    for(int k=0; k<decoder.NLoop; k++) {
      cout << "=============tpkln[2."<<k+1<<"] Calling DoCompPk() in loop, k="<<k<<" / NLoop="<<decoder.NLoop<<endl;
      vector<HProf> vh=DoCompPk(interpk, decoder.grid_, decoder.gfkparm_, decoder.hpkdef_, ntms, decoder.prtlev);      
      VHP_Update(k, decoder.hpkdef_, vhpk, vh, pofvhpk);
    }
    cout << " tpkln[3.a] saving PkIn, recPk to PPF file "<<outppfname<<endl;
    pof<<PPFNameTag("PkIn")<<vhpk[0];
    pof<<PPFNameTag("recPk")<<vhpk[1];
    DataTable dt = HPk2Dt(vhpk, decoder.grid_, interpk);
    cout << " tpksph[3.b] recPk DataTable (dtpk) to "<<outppfname<<endl;
    pof<<PPFNameTag("dtpk")<<dt;
    cout << " tpksph[3.c] grid mean-sigma NTuple (ntms) saved to "<<outppfname<<endl;
    cout << ntms;
    pof<<PPFNameTag("ntms")<<ntms;
  }  // End of try bloc 
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " tpkln.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " tpkln.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " tpkln.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
  }
  cout << " ==== End of tpkln.cc program  Rc= " << rc << endl;
  return rc;    
}


vector<HProf> DoCompPk(ClassFunc1D& interpk, GridDef& grid, GFkParam& gfkparm, HPkDef& hpkparm, NTuple& ntms, int prtlev)
{
  Timer tm("DoCompPk");
  int nkbin=hpkparm.nkbin;
  double kmin=hpkparm.kmin;
  double kmax=hpkparm.kmax;

  vector<HProf> rvhpk;

  cout << "DoCompPk[1]: Creating 3D map Nx Ny Nz= " << grid.Nx<<" x "<<grid.Ny<<" x "<<grid.Nz<<" and GFour3DPk"<<endl;
  TArray<r_4> ingrid(grid.Nx,grid.Ny,grid.Nz);
  
  GFour3DPk  gpkc(ingrid);
  gpkc.SetPrtLevel(prtlev);
  cout << "tpkln[1.c]: setting grid cell size: "<< grid.dx<<" x "<<grid.dy<<" x "<<grid.dz<<endl;    
  gpkc.SetGridCellSize(grid.dx,grid.dy,grid.dz);
  
  /*
    if (fgtrunclowval) {  // truncating grid value below threshold 
      cout << "tpkln[2] : calling CleanNegatives() ... " << endl;    
      gpkc.CleanNegatives(lowval);
      MeanSigma(ingrid, mean, sigma);
      cout << "... After CleanNegatives grid Mean=" << mean << " Sigma=" << sigma << endl;
      tm.Split(" After CleanNegatives ");
    }
    */
  
  cout << "DoCompPk[2.a]: Generating Fourier coefficients from power spectrum ... " << endl;
  if (gfkparm.sigmaz_>0.) cout << " Smearing along the z-axis with sigma="<<gfkparm.sigmaz_<<endl;
  gpkc.generateFourierAmp(interpk, gfkparm.sigmaz_, gfkparm.rfac_);
  tm.Split(" After generateFourierAmp ");
  cout << "DoCompPk[2.b]: Computing 1D power spectrum ... nbin="<<nkbin<<" kmin="<<kmin<<" kmax="<<kmax<<endl;
  HProf recPkIn=gpkc.ComputePk(nkbin, kmin, kmax, true);
  rvhpk.push_back(recPkIn);
  tm.Split(" After ComputePk");


  double xnt[4];
  
  cout << "DoCompPk[2.c]: computing d rho/rho-bar grid gpkc.doInverseFFT() ..."<<endl;
  gpkc.doInverseFFT();
  double mean, sigma;
  MeanSigma(ingrid, mean, sigma);
  xnt[0]=mean;  xnt[1]=sigma;
  cout << "DoCompPk[2.d]: d rho/rho-bar, mean="<<mean<<" sigma="<<sigma<<endl;
  tm.Split(" After doInverseFFT + mean-sigma");
  
  if (gfkparm.fgsmneg_) {
    cout << "DoCompPk[3]: calling SmoothCleanNegatives()  ..."<<endl;
    gpkc.SmoothCleanNegatives();
    MeanSigma(ingrid, mean, sigma);
    xnt[2]=mean;  xnt[3]=sigma;
    cout << "... After SmoothCleanNegatives grid Mean=" << mean << " Sigma=" << sigma << endl;
    tm.Split(" After SmoothCleanNegatives ");
  }
  else if (gfkparm.fgtrunclowval_) {
    cout << "DoCompPk[3]:: calling CleanNegatives() ... " << endl;    
    gpkc.CleanNegatives(gfkparm.lowval_);
    MeanSigma(ingrid, mean, sigma);
    xnt[2]=mean;  xnt[3]=sigma;
    cout << "... After CleanNegatives grid Mean=" << mean << " Sigma=" << sigma << endl;
    tm.Split(" After CleanNegatives ");
  }
  else if (gfkparm.fglognormal_) { 
    cout << "DoCompPk[3]: Computing the exp(d rho/rho) ..."<<endl;
    MathArray<r_4> ma;
    ma.ApplyFunctionInPlace(ingrid, exp);
    MeanSigma(ingrid, mean, sigma);
    xnt[2]=mean;  xnt[3]=sigma;
    cout << "... Exp[d rho/rho-bar], mean="<<mean<<" sigma="<<sigma<<endl;
    ingrid -= (r_4)mean;   ingrid /= (r_4)mean;
    MeanSigma(ingrid, mean, sigma);
    cout << "... (Exp[d rho/rho-bar]-mean)/mean -> mean="<<mean<<" sigma="<<sigma<<endl;
  }
  ntms.Fill(xnt);
  
  if (gfkparm.fgpoisson_) {
    cout << "DoCompPk[3]: Converting to Galaxy density with Poisson, GalDens="<< gfkparm.numberdensity_<< "/Mpc^3..."<<endl;
    gpkc.ConvertToGalaxyDensity(gfkparm.numberdensity_, true, true);
  }
  else if (gfkparm.numberdensity_ > 0.1) {
    cout << "DoCompPk[3]: Converting to Galaxy density TruncInteger GalDens="<< gfkparm.numberdensity_<< "/Mpc^3..."<<endl;
    gpkc.ConvertToGalaxyDensity(gfkparm.numberdensity_, false, true);
  }
  cout << "DoCompPk[5.a]: Computing Fourier coefficients ..."<<endl;
  gpkc.doFFT();
  cout << "DoCompPk[5.b]: Computing 1D power spectrum ... nbin="<<nkbin<<" kmin="<<kmin<<" kmax="<<kmax<<endl;
  HProf recPk=gpkc.ComputePk(nkbin, kmin, kmax, true);
  rvhpk.push_back(recPk);

  return rvhpk;
}

/*-- Fonction --*/
void VHP_Update(int k, HPkDef& hpkparm, vector<HProf>& vhpk, vector<HProf>& vh, POutPersist& pof)
{
  if (k==0)
    for(int i=0; i<vh.size(); i++)
      vhpk.push_back(HProf(hpkparm.kmin, hpkparm.kmax, hpkparm.nkbin));
  
  for(int i=0; i<vhpk.size(); i++) {
    HProf& hpa=vhpk[i];
    HProf& hpb=vh[i];
    for(int_4 ib=0; ib<hpa.NBins(); ib++) {
      hpa.Add(hpb.BinCenter(ib), hpb(ib));
    }
  }
  char buff[64];
  sprintf(buff,"genpk_%d",k);
  pof << PPFNameTag(buff) << vh[0];
  sprintf(buff,"recpk_%d",k);
  pof << PPFNameTag(buff) << vh[1];

  return;
}

/*-- Fonction --*/
DataTable HPk2Dt(vector<HProf>& vhpk, GridDef& grid, ClassFunc1D& interpk)
{
  const char* nomcol[6] = {"k","inpk","genpk","recpk","sig_genpk","sig_recpk"};
  DataTable dt;
  for(int i=0; i<6; i++)   dt.AddDoubleColumn(nomcol[i]);    

  double vol = grid.Nx*grid.Ny*grid.Nz*grid.dx*grid.dy*grid.dz;
  HProf& hp=vhpk[0];
  HProf& hpr=vhpk[1];

  DataTableRow 	dtr = dt.EmptyRow();
  double k;
  for(int_4 ib=0; ib<hp.NBins(); ib++) {
    dtr[0]=k=hp.BinCenter(ib);
    dtr[1]=interpk(k);
    dtr[2]=hp(ib)*vol;
    dtr[3]=hpr(ib)*vol;
    // les erreurs/sigmas 
    dtr[4]=hp.Error(ib)*vol;
    dtr[5]=hpr.Error(ib)*vol;

    dt.AddRow(dtr);
  }
  dt.SetShowMinMaxFlag(true);
  cout << dt;
  return dt; 
}
