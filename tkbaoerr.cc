/* ----
   Project   LSST/BAO/PhotoZ
   Calcul des erreurs sur k_BAO avec le rapport P(k)/P_nos(k)  avec erreur associe 
   R.Ansari - Mars 2018 
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
#include "luc.h"

#include "fiosinit.h"


//---- this modules include files
#include "corfunc.h"
#include "myinteg2d.h"
#include "hsplfit.h"
#include "gpkspec.h"
#include "gpkutil.h"
		


//----- Damping function xsi(k) 
double damping_xsi(double k, double sigmaR) 
{
  double rk=k*sigmaR;
  return (1./sqrt(1+rk*rk));
}

// IQR values in percent from paper Fig 3 
#define NIQR 32
double lesz_[NIQR] = { 0.1,0.3,0.4,0.48,0.5,0.53,0.6,0.7,0.8,0.85,0.9,
		       1.0,1.1,1.15,1.2,1.28,1.38,1.5,1.6,1.65,1.7,
		       1.76,1.83,1.9,1.95,2.,2.1,2.15,2.2,2.3,2.4,2.6};
double lesIQR_[NIQR] = { 2.2,2.3,3.2,3.5,3.,2.4,2.2,2.,2.2,2.5,2.,
			 2.15,3.,3.1,3.,2.8,3.3,3.,3.7,4.,4.8,
			 4.8,4.2,5.4,6.,5.8,6.2,6.3,6.,5.5,5.,5.};

//-------------------------------------------------------------------------
//      ------------------ MAIN PROGRAM ------------------------------
//-------------------------------------------------------------------------
int main(int narg, const char* arg[])
{
  int rc = 0;
  try {
    GPkArgDecoder decoder;
    if ((narg<3)||((narg>1)&&(strcmp(arg[1],"-h")==0))) {
      cout << " Usage: tkbaoerr [options] InputPkFile Out_FITS_File [sigmaZ=0.] [z1,z2=0.7,1.2] [galDensMpc=1e-2] \n" 
	   <<"                  [kminfit,kmaxfit=0.0199,0.221] [OmegaSurvey=10000 deg^2] [f_out,P_out=0.,1000] \n"
	   << "   options: [-k nbin,kmin,kmax] [-N nloop] [-prt lev] [-rgi] \n"
	   << "   InputPkFile : text file with Pk as output by SimLSS \n" 
	   << "   sigmaZ : photoZ gaussian smearing std-dev (if sigmaZ=A compute it)\n"
	   << "   z1,z2 : redshift range  \n"
	   << "   galDensMpc : galaxy number density (n_gal / Mpc^3) \n"
	   << "   kminfit,kmaxfit : spatial wave number k-range used for fit \n"
	   << "   OmegaSurvey: survey sky surface in deg^2 \n"
	   << "   f_out,P_out : Outlier fraction and corresponding flat power spectrum level (1/Mpc^3) \n" << endl; 
      return 1;
    }
    Timer tm("tkbaoerr");
    string inpkname="pksimlss.txt";
    string outfitsname="!toto.fits";
	
    //----------------------------------------------------------
    decoder.DecodeArgs(narg, arg, 2, "tkbaoerr");
    if (decoder.lastargs.size()<2)
      { cout << " tkbaoerr/missing/bad argument, tkbaoerr -h for help " << endl;  return 2; }

    inpkname=decoder.lastargs[0];
    outfitsname=decoder.lastargs[1];

    double sigmaZ=0.0;
    bool fgautosigmaZ=false;
    if (decoder.lastargs.size()>2) {
      if (decoder.lastargs[2]=="A")  fgautosigmaZ=true;
      else fgautosigmaZ=atof(decoder.lastargs[2].c_str());
    } 
    double z1=0.7;  
    double z2=1.2;
    if (decoder.lastargs.size()>3) 
      sscanf(decoder.lastargs[3].c_str(),"%lg,%lg",&z1,&z2);
    double galdensMpc3=1.e-2;
    if (decoder.lastargs.size()>4)   galdensMpc3=atof(decoder.lastargs[4].c_str());
    double kminfit=0.0199, kmaxfit=0.221;
    if (decoder.lastargs.size()>5) 
      sscanf(decoder.lastargs[5].c_str(),"%lg,%lg",&kminfit,&kmaxfit);    
    double OmegaSurv=10000.;
    if (decoder.lastargs.size()>6)   OmegaSurv=atof(decoder.lastargs[6].c_str());
    double f_out=0.0, P_out=1000.;
    if (decoder.lastargs.size()>7) 
      sscanf(decoder.lastargs[7].c_str(),"%lg,%lg",&f_out,&P_out);    

    cout << "tkbaoerr[1]: reading input power spectrum from file "<<inpkname<<endl;
    decoder.ReadSimLSSPkFile(inpkname);

    SimpleUniverse su;
    cout << "tkbaoerr[2]: Creating interpolated sigmaZ=f(z) "<<endl;
    vector<double> lesz(NIQR), lessigmaz(NIQR);
    for(size_t i=0; i<lesz.size(); i++) {
      lesz[i]=lesz_[i];  lessigmaz[i]=lesIQR_[i]*(1.+lesz[i])*0.01/1.35;
    }
    SLinInterp1D lin_sigmaz_z(lesz, lessigmaz);
    {
      const char * namesZ[4] = {"z","H_z","sigmaZ","sigmaR"};
      NTuple  ntZ(4, namesZ);
      double xntZ[5];
      for(double z=0.3;z<2.5; z+=0.05) {
	su.SetEmissionRedShift(z);
	xntZ[0]=z;  xntZ[1]=su.H(z);
	xntZ[2]=lin_sigmaz_z(z);  xntZ[3]=(3.e5/xntZ[1])*xntZ[2];
	ntZ.Fill(xntZ);
      }
      cout << ntZ;
      POutPersist pof("sigmazr.ppf");  pof<<ntZ;
      cout << "tkbaoerr[2.b]:  ntZ ntuple saved to file sigmazr.ppf"<<endl;
    }
    
    su.SetEmissionRedShift(z1);
    double V1=su.VolumeMpc3();
    su.SetEmissionRedShift(z2);
    double V2=su.VolumeMpc3();
    double zcenter=0.5*(z1+z2);
    su.SetEmissionRedShift(zcenter);
    double Hz = su.H(zcenter);
    double sigmaR = (3.e5/Hz)*sigmaZ;
    if (fgautosigmaZ) {
      cout << " Computing sigmaR from interpolated sigmaZ = f(z) ... ==> sigmaR= ";
      int mcnt=0;  double msigmar=0.;
      for(double z=z1; z<=z2; z+=0.05) {
	su.SetEmissionRedShift(z);
	msigmar += (3.e5/su.H(z))*lin_sigmaz_z(z);  mcnt++;
      }
      sigmaR = msigmar/(double)mcnt;
      cout << sigmaR << " Mpc"<<endl;
    }

    double deg2 = (M_PI/180.)*(M_PI/180.);
    double OmegaTotdeg2 = 4.*M_PI/deg2;
    double Vsurv=(V2-V1)*OmegaSurv/OmegaTotdeg2;

    double deltak=decoder.hpkdef_.delta_k;
    double Pnoise=1./galdensMpc3;

    cout << " OmegaSurv="<<OmegaSurv<<" deg^2 -> "<<OmegaSurv/OmegaTotdeg2*100.<<" % of 4Pi SurveyVolume="<<Vsurv<<" Mpc^3 \n"
	 << " deltak="<<deltak<<" Mpc^-1 "<< "   GalDens="<<galdensMpc3<< "/Mpc^3  PNoise="<<Pnoise<<" Mpc^-3"<<endl;
    cout << " Redshift range: "<<z1<<" < z < "<<z2<<" zc="<<zcenter<<" sigmaZ="<<sigmaZ<<" ->sigmaR="<<sigmaR<<" Mpc"
	 <<" ( H(zc)= "<<Hz<<" km/s/Mpc)"<<endl;
    cout << " k range and binning kmin="<<decoder.hpkdef_.kmin<<" kmax= "<<decoder.hpkdef_.kmax<<" nkbin="<<decoder.hpkdef_.nkbin<<endl;
    cout << " fit range for k: "<<kminfit<<" < k < "<<kmaxfit<<" Outliers: f_out="<<f_out<<" P_out="<<P_out<<endl;
    double CstSigmaPk=2.*M_PI/sqrt(deltak*Vsurv);
    
    const char * names[8] = {"rcfit","xi2red","A", "errA","sdamp","errsdamp","sbao","errsbao"};
    NTuple  nt(8, names);
    double xnt[10];

    FMTRandGen rg;  // Random generator
    if (decoder.fg_rand_autoinit_) rg.AutoInit();
    int nokfit = 0;
    double mean_sbao=0., sig_sbao=0.;
    double mean_errsbao=0., sig_errsbao=0.;
    for(int ll=0; ll<decoder.NLoop; ll++) {  // Loop over  
      GeneralFitData mGdata(1, decoder.hpkdef_.nkbin);

      for(int i=0; i<decoder.hpkdef_.nkbin; i++)  {
	double k=decoder.hpkdef_.Getk(i);
	if (!((k>kminfit)&&(k<kmaxfit))) continue;
	double Pk=decoder.fpk_(k);
	Pk *= ((1.-f_out)*(1.-f_out));
	double Pknosc=decoder.fpknosc_(k);
	Pknosc *= ((1.-f_out)*(1.-f_out));
	if (sigmaR>1.)  { 
	  double fdamp=damping_xsi(k,sigmaR);
	  Pk *= fdamp;
	  Pknosc *= fdamp;
	}
	double sigmaPk=(CstSigmaPk/k)*(Pk+Pnoise+(P_out*f_out*f_out));
	double rapp=(Pk+rg.Gaussian(sigmaPk))/Pknosc;
	double errrap=sigmaPk/Pknosc;
	mGdata.AddData1(k,rapp,errrap); // Fill x, y and error on y	
      }

      double A,errA;
      double sdamp,errsdamp;
      double sbao,errsbao;
      double xi2r;
      int rcfit=doFitkBAO_B(mGdata, xi2r, A,errA, sdamp,errsdamp, sbao,errsbao, decoder.prtlev);
      if (rcfit > 0) {
	xnt[0]=rcfit;  xnt[1]=xi2r;
	xnt[2]=A;  xnt[3]=errA;
	xnt[4]=sdamp;  xnt[5]=errsdamp;
	xnt[6]=sbao;  xnt[7]=errsbao;
	nt.Fill(xnt);    nokfit ++;
	mean_sbao+=sbao;  sig_sbao+=(sbao*sbao);
	mean_errsbao+=errsbao;  sig_errsbao+=(errsbao*errsbao);

      }
      if (decoder.prtlev > 0) 
	cout << " done k_BAO/s_BAO fit in loop l="<<ll<<" -> rcfit="<<rcfit<<" sbao="<<sbao<<" +/- "<<errsbao<<endl;;
    }
    cout << nt;
    cout << "-------------------------------------------------------------"<<endl;
    cout << "--------   NOkFit="<<nokfit<<" / N="<<decoder.NLoop<<endl;
    if (nokfit>0) {
      double oon=1./(double)nokfit;
      mean_sbao*=oon;   sig_sbao*=oon;  sig_sbao=sqrt(sig_sbao-mean_sbao*mean_sbao);  
      mean_errsbao*=oon;   sig_errsbao*=oon;  sig_errsbao=sqrt(sig_errsbao-mean_errsbao*mean_errsbao);  
    }
    cout << "    sbao (mean+/-sigma) = " << mean_sbao << " +/- " << sig_sbao << " ("<<100.*sig_sbao/mean_sbao<<" %)"<<endl;
    cout << " errsbao (mean+/-sigma) = " << mean_errsbao << " +/- " << sig_errsbao << endl;
    cout << "-------------------------------------------------------------"<<endl;

    cout << " Sauvegarde NTuple ds le fichier FITS " << outfitsname << endl;
    FitsInOutFile fos(outfitsname, FitsInOutFile::Fits_Create);
    fos << nt;

  }  // End of try bloc 
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " tkbaoerr.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " tkbaoerr.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " tkbaoerr.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
  }
  cout << " ==== End of tkbaoerr.cc program  Rc= " << rc << endl;
  return rc;    
}





