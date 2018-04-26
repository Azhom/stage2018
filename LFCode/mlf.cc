/* ----
   Project   LSST/BAO/PhotoZ
   Test program for  MultiType_Z_LF classes which generate galaxie Absolute  
   Magnitude and types starting  from Luminosity Functions (LF's)  
   F. Habibi & R.Ansari ,  LAL (CNRS/IN2P3)  & Univ. Paris Sud 
   April - July 2017                                          -------  */

#include <iostream>

using namespace std;

#include "multitypzlf.h"
#include "pexceptions.h"
#include "array.h"
#include "histats.h"
#include "ctimer.h"
#include "randfmt.h"

using namespace SOPHYA;

/* To perform the tests and check the results 
  (a) Execution speed   ( random generation of 5 x 10^6 galaxies ) 
csh> time ./Objs/mlf lfparamRamosAll.txt . -25 -17 1000000
  (b) execution and checking results 
csh> ./Objs/mlf lfparamRamosAll.txt ckmlf.ppf -25 -17
  run piapp and display results
csh> spiapp -term
setaxesatt 'font=helvetica,bold,20 fixedfontsize minorticks'
delobjs *
openppf ckmlf.ppf
nt2d galcount z nGal - - - - 'line=solid,2 cpts nsta black'
nt2d galcount z nAll - - - - 'line=dashed,2 cpts nsta same navyblue'
nt2d galcount z nEll - - - - 'line=solid,1 cpts nsta same red'
nt2d galcount z nSp - - - - 'line=solid,1 cpts nsta same green'
nt2d galcount z nSB - - - - 'line=solid,1 cpts nsta same blue'
sleep 3
newh1d hAll -25 -17 80 
newh1d hEll -25 -17 80 
newh1d hSp -25 -17 80 
newh1d hSB -25 -17 80 
projh1d hAll galaxies mag 1 
projh1d hEll galaxies mag 1 fabs(type-1)<0.01 
projh1d hSp galaxies mag 1 fabs(type-2)<0.01 
projh1d hSB galaxies mag 1 fabs(type-3)<0.01 
disp hAll 'logy black'
disp hEll 'same red'
disp hSp 'same green'
disp hSB 'same blue'

*/

int main(int narg, char* arg[])
{
  int rc = 0;
  
  cout<< endl;
  cout << "\n ---- mlf.cc Test of MultiType_Z_LF class ------ " << endl;
  if ((narg<3)||((narg>1)&&(strcmp(arg[1],"-h")==0))||((narg>1)&&(narg<4)&&(strcmp(arg[1],"-useall")==0)) ) {
    cout << " Test of MultiType_Z_LF class  Usage: \n"
	 << " mlf [-useall] input_schech_param_file outppf_file [magMin=-25] [magMax=-17] [ngal=200000] [H0=70 km/s/Mpc] \n"
	 << "  o input_schech_param_file: input text file name containing Schechter LF parameters \n" 
    	 << "  o outppf_file: input text output PPF file name  \n"
	 << "    specify - (dash) as outppf_file to avoid creating output PPF file \n"
	 << "    ( Use - for outppf_file to check galaxy generation speed with larger ngal) \n" << endl;
    return 1;
  }

  try  {
    Timer tm("mlf.cc");
    int off=1;
    bool useall=false;
    if (strcmp(arg[1],"-useall")==0)  { useall=true; off=2; }
    string inlfname = arg[off];
    string outppfname = arg[off+1];
    double magMin=-25. , magMax=-17;
    double H0=100.;
    if (narg>off+2)  magMin=atof(arg[off+2]);
    if (narg>off+3)  magMax=atof(arg[off+3]);
    int ngal=200000;
    if (narg>off+4)  ngal=atoi(arg[off+4]);
    if (narg>off+5)  {
      H0=atof(arg[off+5]);
      MultiType_Z_LF::set_def_H0(H0);
    }
    
    // Checking if we have to save NTuples to an output PPF file 
    bool fgoutppf=false;
    if (isalpha(outppfname[0]))  fgoutppf=true;
    
    cout << "mlf.cc/Info: Input file containing Scheshter parameters: " << inlfname << " UseAll="<<(useall?"TRUE":"FALSE")<<endl;
    cout << "mlf.cc/Info: Absolute B-mag min/max: " << magMin << ", " << magMax << "  Ngal="<<ngal<<endl;
    if (fgoutppf)   cout << "mlf.cc/Info: NTuples will be saved to output PPF file: " << outppfname << endl;
    else  cout << "mlf.cc/Info: NO output PPF file will be created"<<endl;
      
    cout << "[1] Creating MultiType_Z_LF multLF(inlfname, MagMin, MagMax, useall) ..."<<endl;
    int nbmagpts=25;
    MultiType_Z_LF multLF(inlfname, magMin, magMax, nbmagpts, useall);
    cout << "... Creating and AutoInit of Random Generator ..."<<endl;
    FMTRandGen rg;
    rg.AutoInit(1);
    multLF.setRandomGenerator(&rg);
    cout << multLF;
    
    cout << "[2] Checking interpolated Shchechter parameters..." << endl;
    cout << " z, multLF.getLF_Elliptical(z) / getLF_Spiral(z) / getLF_StarBurst(z) / multLF.getLF_All(z) -> Phistar Mstar alpha " << endl;
    for(double z=0.1; z<0.7; z+=0.05) {
      cout << " z="<<z<<" Ell: " << multLF.getLF_Elliptical(z) << " Sp: "<< multLF.getLF_Spiral(z)
	   << " SB: "<< multLF.getLF_StarBurst(z) << " All: "<< multLF.getLF_All(z) << endl;
      MySchechter lfEll = multLF.getLF_Elliptical(z);
    }
    tm.Split("Done [2]");
    cout << "[3] Checking interpolated number of galaxies ..." << endl;
    const char * ntnamesa[13]={"z","nEll","nSp","nSB","nAll","nGal",
                               "LumDensEll","LumDensSp","LumDensSB","LumDensAll","fEll","fSp","fSB"};
// Pour tester scaling des schechter- fait 14/12/2017  "nEllSch","nSpSch","nSBSch","nAllSch"};
    NTuple nta(13,ntnamesa,128);
    double xnt[20];
    double fEll, fSp, fSB;
    for(int i=0; i<100; i++) {
      double z=(double)i*0.05;
      xnt[0]=z;
      xnt[1]=multLF.getIntegral_Elliptical(z);
      xnt[2]=multLF.getIntegral_Spiral(z);
      xnt[3]=multLF.getIntegral_StarBurst(z);
      xnt[4]=multLF.getIntegral_All(z);
      xnt[5]=multLF.getGalaxyNumberDensity(z, fEll, fSp, fSB);
      xnt[6]=multLF.getLF_Elliptical(z).getLumWeightedIntegral(magMin, magMax);
      xnt[7]=multLF.getLF_Spiral(z).getLumWeightedIntegral(magMin, magMax);
      xnt[8]=multLF.getLF_StarBurst(z).getLumWeightedIntegral(magMin, magMax);
      xnt[9]=multLF.getLF_All(z).getLumWeightedIntegral(magMin, magMax);
      xnt[10]=fEll;   xnt[11]=fSp;   xnt[12]=fSB;   
/*  Pour verifier que le scaling des Schechter se fait correctement - 
    Test fait le 14 decembre 2017 avec Farhang
      xnt[6]=multLF.getLF_Elliptical(z).getIntegral(magMin, magMax);
      xnt[7]=multLF.getLF_Spiral(z).getIntegral(magMin, magMax);
      xnt[8]=multLF.getLF_StarBurst(z).getIntegral(magMin, magMax);
      xnt[9]=multLF.getLF_All(z).getIntegral(magMin, magMax);
 */
      if (i%5==0) {
        cout << "z="<<z<<" ILF: Ell="<<xnt[1]<<" Sp="<<xnt[2]<<" SB="<<xnt[3]
             <<" All="<<xnt[4]<<" NGal=?SumOrAll="<<xnt[5]<<endl;
        cout << "  Luminosity Density: Schech.getLumWeightedIntegral(): Ell="<<xnt[6]<<" Sp="<<xnt[7]<<" SB="<<xnt[8]
             << " All="<<xnt[9]<<endl;
      }
      nta.Fill(xnt);
    }

    tm.Split("Done [3]");
    cout << nta;

    const char * ntnamesb[2]={"mag","type"};
    NTuple ntb(2,ntnamesb,1024,false);  // false -> using float as column content 

    if (fgoutppf) {
      double z = 0.5+rg.Flat01()*2.5;
      cout<<"[4.a] Generating " << ngal << " galaxies at resdhift z="<<z<<" and saving them as NTuple in PPF file"<<endl;
      double mag;
      int type;
      float ynt[5];
      for(int i=0; i<ngal; i++) {
	multLF.getTypeMagnitude(z,mag,type);
	ynt[0]=mag;  ynt[1]=type;
	ntb.Fill(ynt);
      }
      cout << ntb;
      cout << "[4.b] Opening output PPF file: " << outppfname << endl;
      POutPersist pof(outppfname);
      cout << "[4.c] Saving NTuple nta as galcount in the PPF file " << endl;
      pof << PPFNameTag("galcount")<<nta;
      cout << "[4.d] Saving NTuple ntb as galaxies in the PPF file " << endl;
      pof << PPFNameTag("galaxies")<<ntb;
    }
    else {
      for(int k=0; k<5; k++) {
	double z = 0.5+rg.Flat01()*2.5;
	cout<<"[4."<<k+1<<"] Generating " << ngal << " galaxies at resdhift z="<<z<<endl;
	double mag;
	int type;
	for(int i=0; i<ngal; i++) {
	  multLF.getTypeMagnitude(z,mag,type);
	}
      }
    }
    tm.Split("Done [4]");
  } // end try
    
    
  catch (PThrowable& exc) {
    cerr << " mlf.cc catched Exception " << exc.Msg() << endl;
    rc = 77;
  }
  catch (std::exception& sex) {
    cerr << "\n mlf.cc std::exception :"
	 << (string)typeid(sex).name() << "\n msg= "
	 << sex.what() << endl;
  }
  catch (...) {
    cerr << " mlf.cc catched unknown (...) exception  " << endl;
    rc = 78; 
  } 
  
  cout << endl;
  cout << " --------------------- End of mlf.cc (RC="<<rc<<") ----------------------- " << endl;
  cout << endl;
  return rc;
}
