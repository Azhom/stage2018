/* ----
   Project   LSST/BAO/PhotoZ
   Tests de calcul de P(k) avec reprojection de grille 
   R.Ansari - Octobre 2016 
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
#include "corfunc.h"
#include "myinteg2d.h"
#include "hsplfit.h"
#include "gpkspec.h"
#include "gpkutil.h"

//---- Fonctions de ce fichier
vector<HProf> DoReprojPk(ClassFunc1D& interpk, GridDef& ing, GridDef& outg, GFkParam& gfkparm, HPkDef& hpkparm, int prtlev=0);
void VHP_Update(int k, HPkDef& hpkparm, vector<HProf>& vhpk, vector<HProf>& vh);
DataTable HPk2Dt(vector<HProf>& vhpk, GridDef& ing, GridDef& outg);
		
//-------------------------------------------------------------------------
//      ------------------ MAIN PROGRAM ------------------------------
//-------------------------------------------------------------------------
int main(int narg, const char* arg[])
{
  int rc = 0;
  try {
    GPkArgDecoder decoder;
    if ((narg<3)||((narg>1)&&(strcmp(arg[1],"-h")==0))) {
      cout << " Usage: tpksph [options] InputPkFile Out_PPFFile \n"
	   << "   InputPkFile : text file with pairs of values k  Pk on each line" << endl;
      return decoder.UsageOptions();
    }
    Timer tm("tpksph");
    string inpkname=".";
    string outppfname="toto.ppf";
	
    //----------------------------------------------------------
    decoder.DecodeArgs(narg, arg, 2, "tpksph");
    if (decoder.lastargs.size()<2)
      { cout << " tpksph/missing/bad argument, tpksph -h for help " << endl;  return 2; }

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
    cout << "tpksph[1.b]: Opening output PPF file " << outppfname <<endl;
    POutPersist pof(outppfname);
    
    vector<HProf>  vhpk;
    
    for(int k=0; k<decoder.NLoop; k++) {
      cout << "=============tpksph[2."<<k+1<<"] Calling DoReprojPk() in loop, k="<<k<<" / NLoop="<<decoder.NLoop<<endl;
      vector<HProf> vh=DoReprojPk(interpk, decoder.ingrid_, decoder.outgrid_, decoder.gfkparm_, decoder.hpkdef_, decoder.prtlev);
      VHP_Update(k, decoder.hpkdef_, vhpk, vh);
    }
    cout << " tpksph[3.a] saving recPkIn, reprojPk, reprojsphPk to PPF file "<<outppfname<<endl;
    pof<<PPFNameTag("recPkIn")<<vhpk[0];
    pof<<PPFNameTag("reprojPk")<<vhpk[1];
    pof<<PPFNameTag("reprojsphAPk")<<vhpk[2];
    pof<<PPFNameTag("reprojsphBPk")<<vhpk[3];
    DataTable dt = HPk2Dt(vhpk, decoder.ingrid_, decoder.outgrid_);
    cout << " tpksph[3.b] recPk DataTable (dtpk) to "<<outppfname<<endl;
    pof<<PPFNameTag("dtpk")<<dt;

  }  // End of try bloc 
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " tpksph.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " tpksph.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " tpksph.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
  }
  cout << " ==== End of tpksph.cc program  Rc= " << rc << endl;
  return rc;    
}



/*-- Fonction --*/
vector<HProf> DoReprojPk(ClassFunc1D& interpk, GridDef& ing, GridDef& outg, GFkParam& gfkparm, HPkDef& hpkparm, int prtlev)
{
  Timer tm("DoReprojPk");

  // Decodage des parametres 
  long int inNx,inNy,inNz;
  inNx=ing.Nx; inNy=ing.Ny; inNz=ing.Nz;
  double indx,indy,indz;
  indx=ing.dx; indy=ing.dy; indz=ing.dz;
  double incenter=ing.r_center;
  long int outNx,outNy,outNz;
  outNx=outg.Nx; outNy=outg.Ny; outNz=outg.Nz;
  double outdx,outdy,outdz;
  outdx=outg.dx; outdy=outg.dy; outdz=outg.dz;
  double outcenter=outg.r_center;
  int nkbin=hpkparm.nkbin;
  double kmin=hpkparm.kmin;
  double kmax=hpkparm.kmax;

  vector<HProf> rvhpk;
  
  cout << "DoReprojPk[1]: Creating 3D map Nx Ny Nz= " << inNx<<" x "<<inNy<<" x "<<inNz<<endl;
  TArray<r_4> ingrid_g(inNx,inNy,inNz);
  {
    TArray<r_4>& ingrid=(ingrid_g);
    GFour3DPk  gpkc(ingrid);
    gpkc.SetPrtLevel(prtlev);
    cout << "DoReprojPk[2]: GFour3DPk created, setting setting grid cell size: "<< indx<<" x "<<indy<<" x "<<indz<<endl;    
    gpkc.SetGridCellSize(indx,indy,indz);
    
    cout << "DoReprojPk[3.a]: Generating Fourier coefficients from power spectrum ... " << endl;
    if (gfkparm.sigmaz_>0.) cout << " Smearing along the z-axis with sigma="<<gfkparm.sigmaz_<<endl;
    gpkc.generateFourierAmp(interpk, gfkparm.sigmaz_, gfkparm.rfac_);
    tm.Split(" After generateFourierAmp ");
    gpkc.doInverseFFT();
    cout << "DoReprojPk[3.b]: Computing 1D power spectrum ... nbin="<<nkbin<<" kmin="<<kmin<<" kmax="<<kmax<<endl;
    HProf recPkIn=gpkc.ComputePk(nkbin, kmin, kmax, true);
    rvhpk.push_back(recPkIn);
    tm.Split(" end-of[3]");
  }
  {
    cout << "DoReprojPk[4]: Creating ReprojGrid... with Output Grid Nx Ny Nz= " << outNx<<" x "<<outNy<<" x "<<outNz<<endl;
    cout << " ... Input grid center="<<incenter<<" Output grid center="<<outcenter<<" Mpc"<<endl;
    cout << " ... Input grid cell size: "<<indx<<" x "<<indy<<" x "<<indz<<" Output cell size: "
	 <<outdx<<" x "<<outdy<<" x "<<outdz<<" Mpc^3"<<endl;   
    TArray<r_4>& ingrid=(ingrid_g);
    ReprojGrid reproj(ing, outg, ingrid);
    reproj.SetPrtLevel(prtlev);
    cout<<"tpksph[4.a]: calling reproj.Project() ..."<<endl;
    reproj.Project();
    double mean, sigma;
    MeanSigma(reproj.getOutGrid(), mean, sigma);
    cout<<" reprojected grid mean="<<mean<<" sigma="<<sigma<<endl;
    tm.Split(" end-of[4.a]-Project");
    {
      cout << "tpksph[4.b]: Computing reprojected grid power spectrum using GFour3DPk..."<<endl; 
      TArray<r_4> prjgrid=reproj.getOutGrid();
      GFour3DPk  gpkc(prjgrid);
      gpkc.SetPrtLevel(prtlev);
      gpkc.SetGridCellSize(outdx,outdy,outdz);
      HProf reprojPk=gpkc.ComputePk(nkbin, kmin, kmax, true);
      rvhpk.push_back(reprojPk);
      tm.Split(" end-of[4]");
    }
    cout<<"tpksph[5.a]: calling reproj.ProjectSpherical(true) - reprojection with spherical geometry in cells with SAME solid angle"<<endl;
    reproj.ProjectSpherical(true);
    MeanSigma(reproj.getOutGrid(), mean, sigma);
    cout<<" reprojected grid mean="<<mean<<" sigma="<<sigma<<endl;
    tm.Split(" end-of[5.a]-ProjectSpherical");
    {
      cout << "tpksph[5.b]: Computing reprojected spherical(SameOmega) grid power spectrum using GFour3DPk..."<<endl; 
      TArray<r_4> prjgrid=reproj.getOutGrid();
      GFour3DPk  gpkc(prjgrid);
      gpkc.SetPrtLevel(prtlev);
      gpkc.SetGridCellSize(outdx,outdy,outdz);
      HProf reprojsphAPk=gpkc.ComputePk(nkbin, kmin, kmax, true);
      rvhpk.push_back(reprojsphAPk);
      tm.Split(" end-of[5]");
    }
    cout<<"tpksph[6.a]: calling reproj.ProjectSpherical(true) - reprojection with spherical geometry in cells with SAME SIZE (Mpc)"<<endl;
    reproj.ProjectSpherical(false);
    MeanSigma(reproj.getOutGrid(), mean, sigma);
    cout<<" reprojected grid mean="<<mean<<" sigma="<<sigma<<endl;
    tm.Split(" end-of[6.a]-ProjectSpherical");
    {
      cout << "tpksph[6.b]: Computing reprojected spherical(SameL) grid power spectrum using GFour3DPk..."<<endl; 
      TArray<r_4> prjgrid=reproj.getOutGrid();
      GFour3DPk  gpkc(prjgrid);
      gpkc.SetPrtLevel(prtlev);
      gpkc.SetGridCellSize(outdx,outdy,outdz);
      HProf reprojsphBPk=gpkc.ComputePk(nkbin, kmin, kmax, true);
      rvhpk.push_back(reprojsphBPk);
      tm.Split(" end-of[6]");
    }

  }
  return rvhpk;
}

/*-- Fonction --*/
void VHP_Update(int k, HPkDef& hpkparm, vector<HProf>& vhpk, vector<HProf>& vh)
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
  return;
}

/*-- Fonction --*/
DataTable HPk2Dt(vector<HProf>& vhpk, GridDef& ing, GridDef& outg)
{
  const char* nomcol[9] = {"k","inpk","prjpk","sphApk","sphBpk","sig_inpk","sig_prjpk","sig_sphApk","sig_sphBpk"};
  DataTable dt;
  for(int i=0; i<9; i++)   dt.AddDoubleColumn(nomcol[i]);    

  double volin = ing.Nx*ing.Ny*ing.Nz*ing.dx*ing.dy*ing.dz;
  double volout = outg.Nx*outg.Ny*outg.Nz*outg.dx*outg.dy*outg.dz;
  HProf& hp=vhpk[0];
  HProf& hpr=vhpk[1];
  HProf& hpA=vhpk[2];
  HProf& hpB=vhpk[3];

  DataTableRow 	dtr = dt.EmptyRow();
  for(int_4 ib=0; ib<hp.NBins(); ib++) {
    dtr[0]=hp.BinCenter(ib);
    dtr[1]=hp(ib)*volin;
    dtr[2]=hpr(ib)*volout;
    dtr[3]=hpA(ib)*volout;
    dtr[4]=hpB(ib)*volout;
    // les erreurs/sigmas 
    dtr[5]=hp.Error(ib)*volin;
    dtr[6]=hpr.Error(ib)*volout;
    dtr[7]=hpA.Error(ib)*volout;
    dtr[8]=hpB.Error(ib)*volout;

    dt.AddRow(dtr);
  }
  dt.SetShowMinMaxFlag(true);
  cout << dt;
  return dt; 
}
