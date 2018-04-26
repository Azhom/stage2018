/* ----
   Project   LSST/BAO/PhotoZ
   Tests de calcul de P(k) avec projection ds grilles 
   R.Ansari - Decembre 2016 
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

#include "zthread.h"

//---- this modules include files
#include "corfunc.h"
#include "myinteg2d.h"
#include "hsplfit.h"
#include "gpkspec.h"
#include "gpkutil.h"

//---- Fonctions de ce fichier
vector<HProf> DoReprojPk(ClassFunc1D& interpk, GridDef& ing, GridDef& outg, GFkParam& gfkparm, HPkDef& hpkparm, SLinInterp1D* selfuncp=NULL);
vector<HProf> DoReprojPkMG(ClassFunc1D& interpk, GridDef& ing, vector<GridDef>& voutg, GFkParam& gfkparm, HPkDef& hpkparm, SLinInterp1D* selfuncp=NULL);
void VHP_Update(int k, HPkDef& hpkparm, vector<HProf>& vhpk, vector<HProf>& vh, POutPersist& pof, TArray< r_8 >* pkspectrasetp);
DataTable HPk2Dt(vector<HProf>& vhpk, GridDef& ing, GridDef& outg, ClassFunc1D& interpk, int noutgrid=1);
void SetDT_Pk_GridInfo(DataTable& dt, GPkArgDecoder& decoder);
		
//-------------------------------------------------------------------------
//      ------------------ MAIN PROGRAM ------------------------------
//-------------------------------------------------------------------------
int main(int narg, const char* arg[])
{
  int rc = 0;
  try {
    GPkArgDecoder decoder;
    if ((narg<3)||((narg>1)&&(strcmp(arg[1],"-h")==0))) {
      cout << " Usage: tpkprj [options] InputPkFile Out_PPFFile \n"
	   << "   InputPkFile : text file with pairs of values k  Pk on each line" << endl;
      return decoder.UsageOptions();
    }
    Timer tm("tpksph");
    string inpkname=".";
    string outppfname="toto.ppf";
	
    //----------------------------------------------------------
    decoder.DecodeArgs(narg, arg, 2, "tpkprj");
    if (decoder.lastargs.size()<2)
      { cout << " tpksph/missing/bad argument, tpkprj -h for help " << endl;  return 2; }

    inpkname=decoder.lastargs[0];
    outppfname=decoder.lastargs[1];
    
    FMTRandGen *RandGen = new FMTRandGen;
    RandGen->SelectGaussianAlgo(C_Gaussian_RandLibSNorm);
    if (decoder.fg_rand_autoinit_)  RandGen->AutoInit(2);
    RandomGeneratorInterface::SetGlobalRandGenP(RandGen);

    if (decoder.fgdebug &&(decoder.debuglevel > 1))   {   /*   For debugging */
      TArray<r_4> ingrid_g(decoder.ingrid_.Nx,decoder.ingrid_.Ny,decoder.ingrid_.Nz);
      ReprojGrid reproj(decoder.ingrid_, decoder.outgrid_, ingrid_g);
      reproj.Grid2RThetaPhiNTuple(outppfname);
      return(0);
    }

    cout << "tpkprj[1.a]: reading input power spectrum from file "<<inpkname<<endl;
    decoder.ReadSimLSSPkFile(inpkname);
    SLinInterp1D& interpk = decoder.fpk_;
    if (decoder.fgnosc_)  interpk = decoder.fpknosc_;
    cout << "tpkprj[1.b]: Opening output PPF file " << outppfname <<endl;
    POutPersist pof(outppfname);
    string outppfvhpk=outppfname;
    size_t pdfn=outppfname.find_last_of('.');
    if (pdfn < outppfname.length()) outppfvhpk=outppfname.substr(0,pdfn);
    outppfvhpk+="_vhpk.ppf";
    cout << "tpkprj[1.c]: Opening output PPF file " << outppfvhpk <<" for individual P(k) HProf's" <<endl;
    POutPersist pofvhpk(outppfvhpk);

    if (decoder.fgshotnoiseonly) cout << "tpkprj[1.d]: Shot Noise Only run ..."<<endl;

    // Pour garder les spectres individuels dans un tableau 
    TArray< r_8 >* pkspectrasetp=NULL;
    if (decoder.fgpkfits) {
      cout << "tpkprj[1.e]: individual spectra P(k)-gen , P(k)-prj-rec will be saved to FITS file "<< decoder.pkfitsfilename << endl;
      pkspectrasetp = new TArray< r_8 >(decoder.hpkdef_.nkbin, decoder.NLoop, 2);
      pkspectrasetp->Info()["CONTENT"]="kz=1,2:  P(k) before/after projection"; 
      pkspectrasetp->Info()["Nbkbin"]=decoder.hpkdef_.nkbin;
      pkspectrasetp->Info()["kmin"]=decoder.hpkdef_.kmin;
      pkspectrasetp->Info()["kmax"]=decoder.hpkdef_.kmax;
      pkspectrasetp->Info()["dk"]=(decoder.hpkdef_.kmax-decoder.hpkdef_.kmin)/(double)decoder.hpkdef_.nkbin;
    }
    

    SLinInterp1D* selfuncp=NULL;
    if (decoder.fgselfunc) {
      cout << "tpkprj[1.d]: run with selection function ..."<<endl;
      selfuncp=&(decoder.selfunc_);
    }
    
    vector<HProf>  vhpk;
    string funcnameinloop=((decoder.voutgrids_.size()>0)?"DoReprojPkMG()":"DoReprojPk()");
    for(int k=0; k<decoder.NLoop; k++) {
      vector<HProf> vh;
      cout << "=============tpkprj[2."<<k+1<<"] Calling "<<funcnameinloop<<" in loop, k="<<k<<" / NLoop="<<decoder.NLoop<<endl;
      if (decoder.voutgrids_.size()>0)  
	vh=DoReprojPkMG(interpk, decoder.ingrid_, decoder.voutgrids_, decoder.gfkparm_, decoder.hpkdef_, selfuncp);
      else 
	vh=DoReprojPk(interpk, decoder.ingrid_, decoder.outgrid_, decoder.gfkparm_, decoder.hpkdef_, selfuncp);
      VHP_Update(k, decoder.hpkdef_, vhpk, vh, pofvhpk, pkspectrasetp);
    }
    cout << " tpkprj[3.a] saving recPkIn, reprojPk, reprojsphPk to PPF file "<<outppfname<<endl;
    pof<<PPFNameTag("recPkIn")<<vhpk[0];
    pof<<PPFNameTag("reprojPk")<<vhpk[1];
    DataTable dt = HPk2Dt(vhpk, decoder.ingrid_, decoder.outgrid_, interpk);
    SetDT_Pk_GridInfo(dt, decoder);
    cout << " tpksph[3.b] recPk DataTable (dtpk) to "<<outppfname<<endl;
    pof<<PPFNameTag("dtpk")<<dt;
    if (decoder.fgpkfits) {
      cout << " tpksph[3.b] saving individual spectra array to FITS file "<<decoder.pkfitsfilename<<endl;
      GridDef& ing=decoder.ingrid_;
      double volin = ing.Nx*ing.Ny*ing.Nz*ing.dx*ing.dy*ing.dz;
      GridDef& outg=decoder.outgrid_;
      double volout = outg.Nx*outg.Ny*outg.Nz*outg.dx*outg.dy*outg.dz;
      TArray< r_8 >& pkspecsa=*pkspectrasetp;
      pkspecsa(Range::all(), Range::all(), Range(0,0)) *= (r_8)volin;
      pkspecsa(Range::all(), Range::all(), Range(1,1)) *= (r_8)volout;
      pkspecsa.Info().Merge(dt.Info());
      pkspecsa.Show(cout);
      FitsInOutFile fos(decoder.pkfitsfilename, FitsInOutFile::Fits_Create);
      fos<<pkspecsa;
      cout << " tpksph[3.c] saving recPk DataTable to 2nd HDU FITS file "<<decoder.pkfitsfilename<<endl;
      fos.SetNextExtensionName("dtpk");
      fos<<dt;
    }	  

  }  // End of try bloc 
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " tpkprj.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " tpkprj.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " tpkprj.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
  }
  cout << " ==== End of tpkprj.cc program  Rc= " << rc << endl;
  return rc;    
}



/*-- Fonction --*/
vector<HProf> DoReprojPk(ClassFunc1D& interpk, GridDef& ing, GridDef& outg, GFkParam& gfkparm, HPkDef& hpkparm, SLinInterp1D* selfuncp)
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

  int prtlev=gfkparm.prtlevel_;
  
  vector<HProf> rvhpk;
  double mean, sigma;

  cout << "DoReprojPk[1]: Creating 3D map Nx Ny Nz= " << inNx<<" x "<<inNy<<" x "<<inNz<<endl;
  TArray<r_4> ingrid_g(inNx,inNy,inNz);
  if (gfkparm.fgshotnoiseonly_)  {
    cout << "DoReprojPk[2]: Shot Noise only ->  ingrid_g=0."<<endl;
    ingrid_g=(TF)0.;
    cout << "DoReprojPk[3.e]: Null power spectrum (HProf) for recPkIn: "<<nkbin<<" kmin="<<kmin<<" kmax="<<kmax<<endl;
    rvhpk.push_back(HProf(kmin, kmax, nkbin));
  }
  else  {
    TArray<r_4>& ingrid=(ingrid_g);
    GFour3DPk  gpkc(ingrid);
    gpkc.SetPrtLevel(prtlev);
    cout << "DoReprojPk[2]: GFour3DPk created, setting setting grid cell size: "<< indx<<" x "<<indy<<" x "<<indz<<endl;    
    gpkc.SetGridCellSize(indx,indy,indz);
    
    cout << "DoReprojPk[3.a]: Generating Fourier coefficients from power spectrum ... " << endl;
    gpkc.generateFourierAmp(interpk, -1., gfkparm.rfac_);
    tm.Split(" After generateFourierAmp ");
    cout << "DoReprojPk[3.b]: computing d rho/rho-bar grid gpkc.doInverseFFT() ... " << endl;
    gpkc.doInverseFFT();
    MeanSigma(ingrid, mean, sigma);
    cout << "... d rho/rho-bar, mean="<<mean<<" sigma="<<sigma<<endl;
    if (gfkparm.fgtrunclowval_) {
      cout << "DoReprojPk[3.c]:: calling CleanNegatives() ... " << endl;    
      gpkc.CleanNegatives(gfkparm.lowval_);
      MeanSigma(ingrid, mean, sigma);
      cout << "... After CleanNegatives grid Mean=" << mean << " Sigma=" << sigma << endl;
      cout << "DoReprojPk[3.d]: Recomputing Fourier coefficients ... " << endl;
      gpkc.doFFT();
      tm.Split(" After CleanNegatives ");
    }

    cout << "DoReprojPk[3.e]: Computing 1D power spectrum ... nbin="<<nkbin<<" kmin="<<kmin<<" kmax="<<kmax<<endl;
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
    cout<<"tpksph[6.a]: calling reproj.Project() ..."<<endl;
    if ((!gfkparm.fgpoisson_)&&(gfkparm.numberdensity_ < 1.))  reproj.Project();
    else reproj.ProjectWithGalaxies(gfkparm.numberdensity_, gfkparm.sigmaz_, gfkparm.fgpoisson_, selfuncp);

    double mean, sigma;
    MeanSigma(reproj.getOutGrid(), mean, sigma);
    cout<<" reprojected grid mean="<<mean<<" sigma="<<sigma<<endl;
    tm.Split(" end-of[6.a]-Project");
    {
      cout << "tpksph[6.b]: Computing reprojected grid power spectrum using GFour3DPk..."<<endl; 
      TArray<r_4> prjgrid=reproj.getOutGrid();
      GFour3DPk  gpkc(prjgrid);
      gpkc.SetPrtLevel(prtlev);
      gpkc.SetGridCellSize(outdx,outdy,outdz);
      HProf reprojPk=gpkc.ComputePk(nkbin, kmin, kmax, true);
      rvhpk.push_back(reprojPk);
      tm.Split(" end-of[6]");
    }
  }
  return rvhpk;
}

//----------------------------------------------------------
//------ Classe de thread pour projection en parallele
//----------------------------------------------------------

class PkPrjThread : public ZThread {
public:
  PkPrjThread(GridDef& ing, GridDef& outg, TArray< TF > & ingrid, GFkParam& gfkparm, HPkDef& hpkparm, SLinInterp1D* selfuncp=NULL)
    : ing_(ing), outg_(outg), ingrid_(ingrid), gfkparm_(gfkparm), hpkparm_(hpkparm), selfuncp_(selfuncp), prtlev_(gfkparm.prtlevel_)
  {
  }
  PkPrjThread(PkPrjThread const& a)
    : ing_(a.ing_), outg_(a.outg_), ingrid_(a.ingrid_), gfkparm_(a.gfkparm_), hpkparm_(a.hpkparm_),
      selfuncp_(a.selfuncp_), prtlev_(a.prtlev_), pkprj_(a.pkprj_) 
  {
  }
  PkPrjThread operator=(PkPrjThread const& a)
  {
    ing_=a.ing_; outg_=a.outg_; ingrid_=a.ingrid_;
    gfkparm_=a.gfkparm_;  hpkparm_=a.hpkparm_;
    prtlev_=a.prtlev_;   selfuncp_=a.selfuncp_;
    pkprj_=a.pkprj_;
    return *this;
  }
  virtual void run()
  {
    ReprojGrid reproj(ing_, outg_, ingrid_);
    reproj.SetPrtLevel(prtlev_);
    cout<<"tpksph[6.a]: calling reproj.Project() ..."<<endl;
    if ((!gfkparm_.fgpoisson_)&&(gfkparm_.numberdensity_ < 1.))  reproj.Project();
    else reproj.ProjectWithGalaxies(gfkparm_.numberdensity_, gfkparm_.sigmaz_, gfkparm_.fgpoisson_, selfuncp_);

    TArray<r_4> prjgrid=reproj.getOutGrid();
    GFour3DPk  gpkc(prjgrid);
    gpkc.SetPrtLevel(prtlev_);
    gpkc.SetGridCellSize(outg_.dx,outg_.dy,outg_.dz);
    pkprj_=gpkc.ComputePk(hpkparm_.nkbin, hpkparm_.kmin, hpkparm_.kmax, true);

    return;
  }

  inline HProf& getHPk() { return pkprj_; }
  
  GridDef& ing_;
  GridDef& outg_;
  TArray< TF > & ingrid_;
  GFkParam& gfkparm_;
  HPkDef& hpkparm_;
  SLinInterp1D* selfuncp_;
  int prtlev_;
  HProf pkprj_;
};


vector<HProf> DoReprojPkMG(ClassFunc1D& interpk, GridDef& ing, vector<GridDef>& voutg, GFkParam& gfkparm, HPkDef& hpkparm, SLinInterp1D* selfuncp)
{
  Timer tm("DoReprojPkMG");

  vector<HProf> rvhpk;
  double mean, sigma;

  int prtlev=gfkparm.prtlevel_;
  if (prtlev>0) cout << "DoReprojPkMG[1]: Creating 3D map Nx Ny Nz= " << ing.Nx<<" x "<<ing.Ny<<" x "<<ing.Nz<<endl;
  TArray<r_4> ingrid_g(ing.Nx,ing.Ny,ing.Nz);
  if (gfkparm.fgshotnoiseonly_)  {
    if (prtlev>0)
      cout << "DoReprojPkMG[2]: Shot Noise only run ->  ingrid_g=0."<<endl;
    ingrid_g=(TF)0.;
    if (prtlev>0)
      cout << "DoReprojPkMG[3]: Null power spectrum (HProf) for recPkIn "<<endl;
    rvhpk.push_back(HProf(hpkparm.kmin, hpkparm.kmax, hpkparm.nkbin));
  }
  else  {
    TArray<r_4>& ingrid=(ingrid_g);
    GFour3DPk  gpkc(ingrid);
    gpkc.SetPrtLevel(prtlev);
    if (prtlev>0)
      cout << "DoReprojPkMG[2]: GFour3DPk created, setting setting grid cell size: "<< ing.dx<<" x "<<ing.dy<<" x "<<ing.dz<<endl;    
    gpkc.SetGridCellSize(ing.dx,ing.dy,ing.dz);
    
    if (prtlev>0) cout << "DoReprojPkMG[3.a]: Generating Fourier coefficients from power spectrum ... " << endl;
    gpkc.generateFourierAmp(interpk, -1., gfkparm.rfac_);
    if (prtlev>1) tm.Split(" After generateFourierAmp ");
    if (prtlev>0) cout << "DoReprojPkMG[3.b]: computing d rho/rho-bar grid gpkc.doInverseFFT() ... " << endl;
    gpkc.doInverseFFT();
    if (prtlev>0) {
      MeanSigma(ingrid, mean, sigma);
      cout << "... d rho/rho-bar, mean="<<mean<<" sigma="<<sigma<<endl;
    }
    if (gfkparm.fgtrunclowval_) {
      if (prtlev>0)  cout << "DoReprojPkMG[3.c]:: calling CleanNegatives() ... " << endl;    
      gpkc.CleanNegatives(gfkparm.lowval_);
      if (prtlev>0) {
	MeanSigma(ingrid, mean, sigma);
	cout << "... After CleanNegatives grid Mean=" << mean << " Sigma=" << sigma << endl;
	cout << "DoReprojPkMG[3.d]: Recomputing Fourier coefficients ... " << endl;
      }
      gpkc.doFFT();
      if (prtlev>1) tm.Split(" After CleanNegatives ");
    }

    if (prtlev>0)
      cout << "DoReprojPkMG[3.e]: Computing 1D power spectrum ... nbin="<<hpkparm.nkbin
	   <<" kmin="<<hpkparm.kmin<<" kmax="<<hpkparm.kmax<<endl;
    HProf recPkIn=gpkc.ComputePk(hpkparm.nkbin, hpkparm.kmin, hpkparm.kmax, true);
    rvhpk.push_back(recPkIn);
    if (prtlev>1) tm.Split(" end-of[3]");
  }
  vector<PkPrjThread> vthr;
  for(size_t jj=0; jj<voutg.size(); jj++)
    vthr.push_back(PkPrjThread(ing, voutg[jj], ingrid_g, gfkparm, hpkparm, selfuncp) );
  for(size_t jj=0; jj<vthr.size(); jj++)  {
    cout << "DoReprojPkMG[4.a] Starting thread for outgrid No "<<jj<<" / "<<vthr.size()<<endl;
    vthr[jj].start();
  }
  cout << "DoReprojPkMG[4.b] Waiting for "<< vthr.size() << " threads to finish ..."<<endl;
  for(size_t jj=0; jj<vthr.size(); jj++)   vthr[jj].join();

  HProf pkprj(hpkparm.kmin, hpkparm.kmax, hpkparm.nkbin);
  
  for(size_t jj=0; jj<vthr.size(); jj++)  pkprj += vthr[jj].getHPk();
  rvhpk.push_back(pkprj);

  return rvhpk;
}


/*-- Fonction --*/
void VHP_Update(int k, HPkDef& hpkparm, vector<HProf>& vhpk, vector<HProf>& vh, POutPersist& pof, TArray< r_8 >* pkspectrasetp)
{
  if (k==0)  {
    for(int i=0; i<vh.size(); i++)
      vhpk.push_back(HProf(hpkparm.kmin, hpkparm.kmax, hpkparm.nkbin));
  }
  
  for(int i=0; i<vhpk.size(); i++) {
    HProf& hpa=vhpk[i];
    HProf& hpb=vh[i];
    for(int_4 ib=0; ib<hpa.NBins(); ib++) {
      hpa.Add(hpb.BinCenter(ib), hpb(ib));
    }
  }
  char buff[64];
  sprintf(buff,"recpkin_%d",k);
  pof << PPFNameTag(buff) << vh[0];
  sprintf(buff,"prjpk_%d",k);
  pof << PPFNameTag(buff) << vh[1];

  if ( pkspectrasetp != NULL ) { // on veut garder les spectres P(k) individuels
    TArray< r_8 >& pkspecsa = *pkspectrasetp;
    HProf& hpb=vh[0];
    HProf& hpc=vh[1];
    for(int_4 ib=0; ib<hpb.NBins(); ib++)  {
      pkspecsa(ib, k, 0)=hpb(ib);
      pkspecsa(ib, k, 1)=hpc(ib);
    }
  }
  return;
}

/*-- Fonction --*/
DataTable HPk2Dt(vector<HProf>& vhpk, GridDef& ing, GridDef& outg, ClassFunc1D& interpk, int noutgrid)
{
  const char* nomcol[9] = {"k","inpk","recpkin","prjpk","sig_recpkin","sig_prjpk","thsig_inpk","thsig_recpkin","thsig_prjpk"};
  DataTable dt;
  for(int i=0; i<9; i++)   dt.AddDoubleColumn(nomcol[i]);    

  double volin = ing.Nx*ing.Ny*ing.Nz*ing.dx*ing.dy*ing.dz;
  double volout = outg.Nx*outg.Ny*outg.Nz*outg.dx*outg.dy*outg.dz;
  HProf& hpa=vhpk[0];
  HProf& hpb=vhpk[1];

  double deltak=hpa.BinWidth(); 
  double CstSigmaPkin=2.*M_PI/sqrt(deltak*volin);
  double CstSigmaPkout=2.*M_PI/sqrt(deltak*volout*(double)noutgrid);

  DataTableRow 	dtr = dt.EmptyRow();
  double k,pkin,recpkin,prjpk;
  for(int_4 ib=0; ib<hpa.NBins(); ib++) {
    dtr[0]=k=hpa.BinCenter(ib);
    dtr[1]=pkin=interpk(k);
    dtr[2]=recpkin=hpa(ib)*volin;
    dtr[3]=prjpk=hpb(ib)*volout;
    // les erreurs/sigmas 
    dtr[4]=hpa.Error(ib)*volin;
    dtr[5]=hpb.Error(ib)*volout;
    // les erreurs theoriques 
    dtr[6]=(CstSigmaPkin/k)*pkin;
    dtr[7]=(CstSigmaPkin/k)*recpkin;
    dtr[8]=(CstSigmaPkout/k)*prjpk;
    dt.AddRow(dtr);
  }
  dt.SetShowMinMaxFlag(true);
  cout << dt;
  return dt;
  
}/*-- Fonction --*/
void SetDT_Pk_GridInfo(DataTable& dt, GPkArgDecoder& decoder)
{
  dt.Info()["InNx"]=(int_8)decoder.ingrid_.Nx;   dt.Info()["InNy"]=(int_8)decoder.ingrid_.Ny;   dt.Info()["InNz"]=(int_8)decoder.ingrid_.Nz;
  dt.Info()["Indx"]=(r_8)decoder.ingrid_.dx;   dt.Info()["Indy"]=(r_8)decoder.ingrid_.dy;   dt.Info()["Indz"]=(r_8)decoder.ingrid_.dz;
  dt.Info()["InCenter"]=(r_8)decoder.ingrid_.r_center;  dt.Info()["InTheta0"]=(r_8)decoder.ingrid_.theta0;   dt.Info()["InPhi0"]=(r_8)decoder.ingrid_.phi0;
  if (decoder.voutgrids_.size() > 0)  
    dt.Info()["NOutGrid"]=(int_8)decoder.voutgrids_.size();
  dt.Info()["NOutGrid"]=(int_8)1;
  dt.Info()["OutNx"]=(int_8)decoder.outgrid_.Nx;   dt.Info()["OutNy"]=(int_8)decoder.outgrid_.Ny;   dt.Info()["OutNz"]=(int_8)decoder.outgrid_.Nz;
  dt.Info()["Outdx"]=(r_8)decoder.outgrid_.dx;   dt.Info()["Outdy"]=(r_8)decoder.outgrid_.dy;   dt.Info()["Outdz"]=(r_8)decoder.outgrid_.dz;
  dt.Info()["OutCenter"]=(r_8)decoder.outgrid_.r_center;  dt.Info()["OutTheta0"]=(r_8)decoder.outgrid_.theta0;   dt.Info()["OutPhi0"]=(r_8)decoder.outgrid_.phi0;
  dt.Info()["InGalDens"]=(r_8)decoder.numberdensity;   dt.Info()["FgPoisson"]=std::string(decoder.fgpoisson?"True":"False");
  dt.Info()["SigmaR"]=(r_8)decoder.sigmaz;   dt.Info()["NLoop"]=(int_8)decoder.NLoop;
  return;
}
