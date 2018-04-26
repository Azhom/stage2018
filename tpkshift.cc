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
vector<HProf> DoCompPk(ClassFunc1D& interpk, GridDef& ingrid, GridDef& outgrid, sa_size_t shiftnz, double deltaXYZ, r_4 shiftweight,
		       RandomGeneratorInterface& rg, GFkParam& gfkparm, HPkDef& hpkparm, NTuple& ntms, int prtlev=0);
void VHP_Update(int k, HPkDef& hpkparm, vector<HProf>& vhpk, vector<HProf>& vh, POutPersist& pof);
DataTable HPk2Dt(vector<HProf>& vhpk, GridDef& ingrid, GridDef& outgrid, ClassFunc1D& interpk);

//-------------------------------------------------------------------------
//      ------------------ MAIN PROGRAM ------------------------------
//-------------------------------------------------------------------------
int main(int narg, const char* arg[])
{
  int rc = 0;
  try {
    GPkArgDecoder decoder;
    if ((narg<3)||((narg>1)&&(strcmp(arg[1],"-h")==0))) {
      cout << " Usage: Usage: tpkshift [options] InputPkFile Out_PPFFile ShiftNZ DeltaXYZ ShiftWeight  \n"
	   << "   InputPkFile : text file with pairs of values k  Pk on each line \n"<<endl;
      return decoder.UsageOptions();
    }
    Timer tm("tpkshift");
    string inpkname=".";
    string outppfname="toto.ppf";

    //----------------------------------------------------------
    decoder.DecodeArgs(narg, arg, 2, "tpkshift");
    if (decoder.lastargs.size()<5)
      { cout << " tpkshift/missing/bad argument, tpkshift -h for help " << endl;  return 2; }

    inpkname=decoder.lastargs[0];
    outppfname=decoder.lastargs[1];
    sa_size_t shiftnz = atoi(decoder.lastargs[2].c_str());
    double deltaxyz = atof(decoder.lastargs[3].c_str());
    r_4 shiftweight =  atof(decoder.lastargs[4].c_str());
    FMTRandGen randg;
    randg.SelectGaussianAlgo(C_Gaussian_RandLibSNorm);
    if (decoder.fg_rand_autoinit_)  randg.AutoInit(2);
    RandomGeneratorInterface::SetGlobalRandGenP(&randg);
 
    cout << "tpksph[1.a]: reading input power spectrum from file "<<inpkname<<endl;
    decoder.ReadSimLSSPkFile(inpkname);
    SLinInterp1D& interpk = decoder.fpk_;
    if (decoder.fgnosc_)  interpk = decoder.fpknosc_;

    cout << "tpkshift[1.b]: Opening output PPF file " << outppfname <<endl;
    POutPersist pof(outppfname);
    string outppfvhpk=outppfname;
    size_t pdfn=outppfname.find_last_of('.');
    if (pdfn < outppfname.length()) outppfvhpk=outppfname.substr(0,pdfn);
    outppfvhpk+="_vhpk.ppf";
    cout << "tpkprj[1.c]: Opening output PPF file " << outppfvhpk <<" for individual P(k) HProf's" <<endl;
    POutPersist pofvhpk(outppfvhpk);
    /*
    if ((decoder.gfkparm_.fgtrunclowval_==false) && (decoder.gfkparm_.fgsmneg_==false) && (decoder.gfkparm_.fglognormal_==false))  {
      cout << " tpkshift / Forcing TruncLowVal, d rho/rho < -1 -> -1 ..."<<endl;
      decoder.gfkparm_.fgtrunclowval_=true;
      decoder.gfkparm_.lowval_=-1.;
    }
    */
    //    cout << " *DBG***B** nkbin="<<decoder.hpkdef_.nkbin<<" kmin="<<decoder.hpkdef_.kmin<<" kmax="<<decoder.hpkdef_.kmax<<endl;

    const char * ntmsnames[4]={"gin_mean","gin_sigma","gsneg_mean","gsneg_sigma"};
    NTuple ntms(4,ntmsnames);
    
    vector<HProf>  vhpk;
    
    for(int k=0; k<decoder.NLoop; k++) {
      cout << "=============tpkshift[2."<<k+1<<"] Calling DoCompPk() in loop, k="<<k<<" / NLoop="<<decoder.NLoop<<endl;
      vector<HProf> vh=DoCompPk(interpk, decoder.ingrid_, decoder.outgrid_, shiftnz, deltaxyz, shiftweight, 
				randg, decoder.gfkparm_, decoder.hpkdef_, ntms, decoder.prtlev);      
      VHP_Update(k, decoder.hpkdef_, vhpk, vh, pofvhpk);
    }
    cout << " tpkshift[3.a] saving PkIn, recPk to PPF file "<<outppfname<<endl;
    pof<<PPFNameTag("PkIn")<<vhpk[0];
    pof<<PPFNameTag("recPk")<<vhpk[1];
    DataTable dt = HPk2Dt(vhpk, decoder.ingrid_, decoder.outgrid_, interpk);
    cout << " tpksph[3.b] recPk DataTable (dtpk) to "<<outppfname<<endl;
    pof<<PPFNameTag("dtpk")<<dt;
    cout << " tpksph[3.c] grid mean-sigma NTuple (ntms) saved to "<<outppfname<<endl;
    cout << ntms;
    pof<<PPFNameTag("ntms")<<ntms;
  }  // End of try bloc 
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " tpkshift.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " tpkshift.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " tpkshift.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
  }
  cout << " ==== End of tpkshift.cc program  Rc= " << rc << endl;
  return rc;    
}

vector<HProf> DoCompPk(ClassFunc1D& interpk, GridDef& igrid, GridDef& ogrid, sa_size_t shiftnz, double deltaXYZ, r_4 shiftweight,
		       RandomGeneratorInterface& rg, GFkParam& gfkparm, HPkDef& hpkparm, NTuple& ntms, int prtlev)
{
  Timer tm("DoCompPk");
  int nkbin=hpkparm.nkbin;
  double kmin=hpkparm.kmin;
  double kmax=hpkparm.kmax;

  vector<HProf> rvhpk;

  cout << "DoCompPk[1]: Creating 3D map Nx Ny Nz= " << igrid.Nx<<" x "<<igrid.Ny<<" x "<<igrid.Nz<<" and GFour3DPk"<<endl;
  TArray<r_4> ingrid(igrid.Nx,igrid.Ny,igrid.Nz);

  GFour3DPk  gpkc(ingrid);
  gpkc.SetPrtLevel(prtlev);
  cout << "tpkshift[1.c]: setting grid cell size: "<< igrid.dx<<" x "<<igrid.dy<<" x "<<igrid.dz<<endl;    
  gpkc.SetGridCellSize(igrid.dx,igrid.dy,igrid.dz);
  
  /*
    if (fgtrunclowval) {  // truncating grid value below threshold 
      cout << "tpkshift[2] : calling CleanNegatives() ... " << endl;    
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

  cout << "DoCompPk[5.a]: Filling output grid["<<ogrid.Nx<<"x"<<ogrid.Ny
       <<"x"<<ogrid.Nz<<"]"<<endl;
  TArray<r_4> gridS(ogrid.Nx,ogrid.Ny,ogrid.Nz);
  //  Range xrng((igrid.Nx-ogrid.Nx)/2, -1, ogrid.Nx, 1);
  //  Range yrng((igrid.Ny-ogrid.Ny)/2, -1, ogrid.Ny, 1);
  Range xrng((igrid.Nx-ogrid.Nx)/2, ogrid.Nx-1);
  Range yrng((igrid.Ny-ogrid.Ny)/2, ogrid.Ny-1);

  if (deltaXYZ < 1e-6) {
    cout << "DoCompPk[5.b] NO XYZ smearing ShiftOnly: ShiftNZ="<<shiftnz<<" ShiftWeight="<<shiftweight<<endl;
    gridS = (r_4)(1.-shiftweight)*ingrid(xrng, yrng, Range(shiftnz,-1,gridS.SizeZ(),1));
    gridS += shiftweight*ingrid(xrng, yrng, Range(0,-1,gridS.SizeZ(),1));
  }
  else {
    cout << "DoCompPk[5.b] XYZ smearing with flat deltaXYZ="<<deltaXYZ<<" ShiftNZ="<<shiftnz<<" ShiftWeight="<<shiftweight<<endl;
    gridS = (r_4)(1.-shiftweight)*ingrid(xrng, yrng, Range(shiftnz,-1,gridS.SizeZ(),1));
    sa_size_t deltax = (igrid.Nx-ogrid.Nx)/2;
    sa_size_t deltay = (igrid.Ny-ogrid.Ny)/2;
    for(sa_size_t kz=0; kz<ingrid.SizeZ(); kz++) {
      for(sa_size_t ky=0; ky<ingrid.SizeY(); ky++) {
	for(sa_size_t kx=0; kx<ingrid.SizeX(); kx++) {
	  sa_size_t okz = kz+shiftnz+deltaXYZ*rg.Flatpm1();
	  if ((okz<0)||(okz>=gridS.SizeZ()))  continue;
	  sa_size_t oky = ky+deltaXYZ*rg.Flatpm1();
	  if ((oky<0)||(oky>=gridS.SizeY()))  continue;
	  sa_size_t okx = kx+deltaXYZ*rg.Flatpm1();
	  if ((okx<0)||(okx>=gridS.SizeX()))  continue;
	  gridS(okx,oky,okz) += shiftweight*ingrid(kx,ky,kz);
	}
      }
    }
  }
  MeanSigma(gridS, mean, sigma);
  cout << "... Shifted grid -> mean="<<mean<<" sigma="<<sigma<<" -> gridS-=mean"<<endl;
  gridS -= (r_4)mean;

  cout << "DoCompPk[6.a]: Computing Fourier coefficients on input grid["<<ingrid.SizeX()<<"x"<<ingrid.SizeY()
       <<"x"<<ingrid.SizeZ()<<"]"<<endl;
  gpkc.doFFT();
  cout << "DoCompPk[6.b]: Computing 1D power spectrum on input grid ... nbin="<<nkbin<<" kmin="<<kmin<<" kmax="<<kmax<<endl;
  HProf recPk=gpkc.ComputePk(nkbin, kmin, kmax, true);
  rvhpk.push_back(recPk);

  GFour3DPk  gpkcS(gridS);
  gpkcS.SetPrtLevel(prtlev);
  gpkcS.SetGridCellSize(ogrid.dx,ogrid.dy,ogrid.dz);
  cout << "DoCompPk[7.b]: Computing Fourier coefficients on output + shift grid"<<endl;
  gpkcS.doFFT();
  cout << "DoCompPk[7.b]: Computing 1D power spectrum on output + shift grid"<<endl;
  HProf recPkS=gpkcS.ComputePk(nkbin, kmin, kmax, true);
  rvhpk.push_back(recPkS);

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
  if (vh.size()>1) 
    pof << PPFNameTag(buff) << vh[1];
  sprintf(buff,"recpkshift_%d",k);
  if (vh.size()>2) 
    pof << PPFNameTag(buff) << vh[2];

  return;
}

/*-- Fonction --*/
DataTable HPk2Dt(vector<HProf>& vhpk, GridDef& ingrid, GridDef& outgrid, ClassFunc1D& interpk)
{
  const char* nomcol[8] = {"k","inpk","genpk","recpk","recpkS","sig_genpk","sig_recpk","sig_recpkS"};
  DataTable dt;
  for(int i=0; i<8; i++)   dt.AddDoubleColumn(nomcol[i]);    

  double volin = ingrid.Nx*ingrid.Ny*ingrid.Nz*ingrid.dx*ingrid.dy*ingrid.dz;
  double volout = outgrid.Nx*outgrid.Ny*outgrid.Nz*outgrid.dx*outgrid.dy*outgrid.dz;

  HProf& hp=vhpk[0];
  HProf& hpr=vhpk[1];
  HProf& hprS=vhpk[2];

  DataTableRow 	dtr = dt.EmptyRow();
  double k;
  for(int_4 ib=0; ib<hp.NBins(); ib++) {
    dtr[0]=k=hp.BinCenter(ib);
    dtr[1]=interpk(k);
    dtr[2]=hp(ib)*volin;
    dtr[3]=hpr(ib)*volin;
    dtr[4]=hprS(ib)*volout;
    // les erreurs/sigmas 
    dtr[5]=hp.Error(ib)*volin;
    dtr[6]=hpr.Error(ib)*volin;
    dtr[7]=hprS.Error(ib)*volout;

    dt.AddRow(dtr);
  }
  dt.SetShowMinMaxFlag(true);
  cout << dt;
  return dt; 
}
