/* ----
   Project   LSST/BAO/PhotoZ
   Tests de calcul de P(k) avec erreur associe 
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
#include "fitkbao.h"
#include "randfmt.h"

#include "fiosinit.h"


//---- this modules include files
#include "corfunc.h"
#include "myinteg2d.h"
#include "hsplfit.h"
#include "gpkspec.h"
#include "gpkutil.h"
		


//-------------------------------------------------------------------------
//      ------------------ MAIN PROGRAM ------------------------------
//-------------------------------------------------------------------------
int main(int narg, const char* arg[])
{
  int rc = 0;
  try {
    GPkArgDecoder decoder;
    if ((narg<3)||((narg>1)&&(strcmp(arg[1],"-h")==0))) {
      cout << " Usage: terrpk [options] InputPkFile Out_PPFFile \n"
	   << "   options: [-grid Nx,Ny,Nz,dx,dy,dz] [-k kmin,kmax,nbin] [-p density] \n"
	   << "   InputPkFile : text file with Pk as output by SimLSS" << endl;
      return decoder.UsageOptions();
    }
    Timer tm("terrpk");
    string inpkname="pksimlss.txt";
    string outppfname="toto.ppf";
	
    //----------------------------------------------------------
    decoder.DecodeArgs(narg, arg, 2, "tpksph");
    if (decoder.lastargs.size()<2)
      { cout << " terrpk/missing/bad argument, terrpk -h for help " << endl;  return 2; }

    inpkname=decoder.lastargs[0];
    outppfname=decoder.lastargs[1];

    cout << "terrpk[1]: reading input power spectrum from file "<<inpkname<<endl;
    decoder.ReadSimLSSPkFile(inpkname);

    double Vsurv=decoder.grid_.GetVolume();
    double deltak=decoder.hpkdef_.delta_k;
    double Pnoise=decoder.grid_.GetCellVolume()/decoder.numberdensity;
    cout << " SurveyVolume="<<Vsurv<<" Mpc^3 deltak="<<deltak<<" Mpc^-1 \n"
	 << "   GalDens="<<decoder.numberdensity/decoder.grid_.GetCellVolume()
	 << "  PNoise="<<Pnoise<<" Mpc^-3"<<endl;
    double CstSigmaPk=2.*M_PI/sqrt(deltak*Vsurv);
    
    const char * names[6] = {"k", "pk","pknosc","pkerr","rap_pk","rap_pk_err"};
    NTuple  nt(6, names);
    double xnt[10];

    int N =decoder.hpkdef_.nkbin;    
    int nkok=0;
    for(int i=0; i<decoder.hpkdef_.nkbin; i++) {
      double k=xnt[0]=decoder.hpkdef_.Getk(i);
      double Pk=xnt[1]=decoder.fpk_(k);
      double Pknosc=xnt[2]=decoder.fpknosc_(k);
      double sigmaPk=xnt[3]=(CstSigmaPk/k)*(Pk+Pnoise);
      xnt[4]=Pk/Pknosc;     xnt[5]=sigmaPk/Pknosc;  
      nt.Fill(xnt);
      if ((k>0.0199)&&(k<0.221))  nkok++;
    }
    cout << nt;
    cout << " terrpk[2] saving Pk NTuple to PPF file "<<outppfname<<endl;
    POutPersist pof(outppfname);
    pof<<nt;

    cout << " Nb k OK = " << nkok << endl;
    FMTRandGen rg;  // Random generator
    if (decoder.fg_rand_autoinit_) rg.AutoInit();

    GeneralFitData mGdata(1, N);

    for(int i=0; i<decoder.hpkdef_.nkbin; i++)  {
      double k=decoder.hpkdef_.Getk(i);
      if (!((k>0.0199)&&(k<0.221))) continue;
      double Pk=decoder.fpk_(k);
      double Pknosc=decoder.fpknosc_(k);
      double sigmaPk=(CstSigmaPk/k)*(Pk+Pnoise);
      double rapp=(Pk+rg.Gaussian(sigmaPk))/Pknosc;
      double errrap=sigmaPk/Pknosc;
      mGdata.AddData1(k,rapp,errrap); // Fill x, y and error on y
      
    }

    mGdata.PrintStatus();
    double A,errA;
    double sdamp,errsdamp;
    double sbao,errsbao;
    double xi2r;
    int rcfit;
    rcfit=doFitkBAO_A(mGdata, xi2r, A,errA, sdamp,errsdamp, sbao,errsbao, 2);
    rcfit=doFitkBAO_B(mGdata, xi2r, A,errA, sdamp,errsdamp, sbao,errsbao, 2);

  }  // End of try bloc 
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " terrpk.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
	 << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " terrpk.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " terrpk.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
  }
  cout << " ==== End of terrpk.cc program  Rc= " << rc << endl;
  return rc;    
}





