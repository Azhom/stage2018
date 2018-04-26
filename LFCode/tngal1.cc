
/*
 Project   LSST/BAO/PhotoZ
 Test program for computing galaxy number counts starting from MultiType_Z_LF , 
 and set of SED's and Filters 
 F. Habibi & R.Ansari ,  LAL (CNRS/IN2P3)  & Univ. Paris Sud
 Novembre 2017                                          -------  */

/*
 The code computes the limit on the absolute B mag for ellipiticals, spirals and star bursts galaxies,
 for a given apparent mag limit on r-band, at different redshifts.
 
 Then, the number density of the elliptical galexis are computed throuh integrating the Scheshter function
 with parameters given by Dahlen et al. 2005
 
 
 R. Ansari & F. Habibi Juillet 2016, Tir 1395, LAL

*/


#include <unistd.h>

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <cmath>

#include "pexceptions.h"
#include "array.h"
#include "histats.h"
#include "ctimer.h"
#include "randfmt.h"
#include "randr48.h"

#include "multitypzlf.h"
#include "sedfilter.h"
#include "kcorr.h"

#include "luc.h"
#include "slininterp.h"
#include "integ.h"

#include "fitsioserver.h"





int main(int narg, char* arg[]) {
    
    cout << endl;
    cout << " ---------- tngal.cc --------- " << endl;
    cout << endl;
    
    string inlfname = "lfparamDahlenAll.txt";
    string outppfname = "otngal.ppf";
    string sedpath = "../SEDs/";
    string filtpath = "../Filters/";
    bool fguseall = false;
    bool fgLFfilB = true;
    double magMin=-25. , magMax=-17;
    int nbmagpts=25;
    double zmin = 0.1, zmax=1., dz=0.05;
    double lambdamin=300., lambdamax=1000.;  // in nanometre
    double lambdaSEDmin=100.;  // on nanometre, SED wavelength interval
    double lambdaSEDmax=2500.;
    double maglim = 25.3;  // Magnitude limit in observation band (LSST-i-band for example)
    
    
    if ((narg>1)&&(strcmp(arg[1],"-h")==0)) {
        cout << " Test computation of ngal using MultiType_Z_LF class  Usage: \n"
        << " tngal [-i input_schech_param_file] [-magMLF magMin,magMax,nbpts] [-useall] \n"
        << "       [-sedp SEDFilePath] [-filtp FilterPath] [-LFB] [-LFiS] \n"
        << "       [-limag magLimit] [-z zmin,zmax,dz] \n"
	     << "       [-lambda min,max] [-o outfilename] \n" << endl;
	cout << "  -i input_schech_param_file: input text file name containing Schechter LF parameters \n" << endl;
        cout << " -sedp sedpath: path for the SEDs directory. By default sedpath=" << sedpath << endl;
        cout << " -filtp filtpath: path for the filters directory. By default filtpath=" << filtpath << endl;
        cout << " -magMLF magMin,magMax,nbmagpts: minimum, maximum and step magnitude for the absolut magnitudes in LF function. By default magMin=-25., magMax=-17 and nbmagpts=25." << endl;
        cout << " -limag maglim: apparent magnitude limit in i-band. By default maglim=25.3" << endl;
        cout << " -z zmin,zmax,dz: minumum, maximum and step redshifts. By default zmin=0.1, zmax=1. and dz=0.05" << endl;
        cout << " -lambda lambdamin,lambdamax: minimum and maximum eawelwngth to integrate the SEDs. By default lambdamin=300., lambdamax=1000 nm." << endl;
        cout << " -o outppfname: output file name in ppf format. Bydefault outppfname = otngal.ppf " << endl;
        cout << " -useall: reads the Scheshter parameters for total galaxies from the input file. False by default." << endl;
        cout << " -LFB: if the Scheshter parameters are given in B-band. True by default" << endl;
        cout << " -LFiS: if the Scheshter parameters are given in iS-band. False by default" << endl;
        cout << endl;
        
        return 1;
    }
    
    
    bool fgoptarg=true;
    while (fgoptarg&&(narg>1))
    {
        string fbo = arg[1];
        
        if (fbo=="-i")  {
            if (narg<3) { cout << " tngal.cc: missing/bad argument, -h for help " << endl;  return -1; }
            inlfname = arg[2];   arg+=2; narg-=2;
        }
        else if (fbo=="-sedp")  {
            if (narg<3) { cout << " tngal.cc: missing/bad argument, -h for help " << endl;  return -1; }
            sedpath = arg[2];   arg+=2; narg-=2;
        }
        else if (fbo=="-filtp")  {
            if (narg<3) { cout << " tngal.cc: missing/bad argument, -h for help " << endl;  return -1; }
            filtpath=arg[2];   arg+=2; narg-=2;
        }
        else if (fbo=="-magMLF")  {
            if (narg<3) { cout << " tngal.cc: missing/bad argument, -h for help " << endl;  return -1; }
            sscanf(arg[2],"%lg,%lg,%d",&magMin,&magMax,&nbmagpts);   arg+=2; narg-=2;
        }
        else if (fbo=="-limag")  {
            if (narg<3) { cout << " tngal.cc: missing/bad argument, -h for help " << endl;  return -1; }
            maglim = atof(arg[2]);     arg+=2; narg-=2;
        }
        else if (fbo=="-z")  {
            if (narg<3) { cout << " tngal.cc: missing/bad argument, -h for help " << endl;  return -1; }
            sscanf(arg[2],"%lg,%lg,%lg",&zmin,&zmax,&dz);   arg+=2; narg-=2;
        }
        else if (fbo=="-lambda")  {
            if (narg<3) { cout << " tngal.cc: missing/bad argument, -h for help " << endl;  return -1; }
            sscanf(arg[2],"%lg,%lg",&lambdamin,&lambdamax);   arg+=2; narg-=2;
        }
        /*
        else if (fbo=="-lambdaSED")  {
            if (narg<3) { cout << " tngal.cc: missing/bad argument, -h for help " << endl;  return -1; }
            sscanf(arg[2],"%lg,%lg",&lambdamin,&lambdamax);   arg+=2; narg-=2;
        }
         */
        else if (fbo=="-o")  {
            if (narg<3) { cout << " tngal.cc: missing/bad argument, -h for help " << endl;  return -1; }
            outppfname=arg[2];   arg+=2; narg-=2;
        }
        else if (fbo=="-useall")  {
            fguseall=true; arg++; narg--;
        }
        else if (fbo=="-LFB")  {
            fgLFfilB=true; arg++; narg--;
        }
        else if (fbo=="-LFiS")  {
            fgLFfilB=false; arg++; narg--;
        }
        else {
            cout << " tngal.cc wrong option "<<fbo<<" use -h to see list of options"<<endl;
            return 9;
        }
    }
    
    cout << " Inpute LF file: " << inlfname << endl;
    cout << " Output file: " << outppfname << endl;
    cout << " Reading SDEs from: " << sedpath << endl;
    cout << " Reading Filters from: " << filtpath << endl;
    cout << " MultiType_Z_LF magnitude range magMin=" << magMin << ", magMax=" << magMax << endl;
    cout << " z-min=" << zmin << ", z-max=" << zmax << " z-step=" << dz << endl;
    // setting the lambda's interval for the spectra
    cout << " Wavelength interval(nm): [" << lambdamin << "," << lambdamax << "]" << endl;
    cout << " SED Wavelength interval(nm): [" << lambdaSEDmin << "," << lambdaSEDmax << "]" << endl;
    cout << " MagnitudeLimit (Observation band) = " << maglim << endl;
    cout << " Luminosity Functions filter = " << (fgLFfilB?"B-Filter":"iS-Filter") << endl;
    
    // make sure SOPHYA modules are initialized
    SophyaInit();
    FitsIOServerInit();
    
    SFPathManager sedpathmgr(sedpath);
    SFPathManager filtpathmgr(filtpath);

    ThSDR48RandGen rg;
    rg.AutoInit();
    
    int rc = 0;
    
    double m5i = 24.+2.5;
    //DEL    double x = pow(10,.4*(maglim-m5i));
    //DEL    double gammai = 0.039;
    // error bar in for mag=maglim in i-band (science book)
    double delmagi = 0.; //sqrt( (.04-gammai)*x + gammai*x*x );
    try
    {
        Timer tm("tngal.cc");
        cout << "[1] Creating MultiType_Z_LF multLF(inlfname, MagMin, MagMax, nbmagpts, useall) ..."<<endl;
        MultiType_Z_LF multLF(inlfname, magMin, magMax, nbmagpts, fguseall);
        cout << " ... Creating and AutoInit of Random Generator ..."<<endl;
        FMTRandGen rg;
        rg.AutoInit(1);
        multLF.setRandomGenerator(&rg);
        cout << multLF;
        

        double delmag = delmagi;
        
        cout << "* magnitude limit=" << maglim << ", delmag=" << delmagi << endl;
   
        cout << "[2] Creating SimpleUniverse su(h,OmegaT,OmgaM,OmegaL) ..."<<endl;
      
        double h = .697; //.7; //Hubble parameter
        double OmegaT = 1.; // Omega total
        double OmegaM = 0.3065; //.3; // Omega Matter
        double OmegaL = OmegaT - OmegaM; // Omega Lambda
        
        // Importing the cosmological parameters
        SimpleUniverse su(h,OmegaT);
        su.SetOmegaMatter(OmegaM);
        su.SetOmegaLambda(OmegaL);

       
        cout << "[3] Reading filters ..."<<endl;

        // importing the filter shapes
 // LF definition filters
        string filterfile_iS = filtpathmgr.BuildFileName("iS_cfht.txt");  // i' (iS) filter of CFHT telescope (Ramos et al. 2011, LF)
        cout << "* Reading Filter from file " << filterfile_iS << endl;
        Filter filt_iS(filterfile_iS, lambdamin*1e-9, lambdamax*1e-9);
        
        string filterfile_B = filtpathmgr.BuildFileName("B_wfi_eso.txt");  // ../filters/g_SDSS.txt");
        cout << "* Reading Filter from file " << filterfile_B << endl;
        Filter filt_B(filterfile_B, lambdamin*1e-9, lambdamax*1e-9);

        // Observation filters
        string filterfile_r = filtpathmgr.BuildFileName("r_lsst_etc.txt");  // ../filters/g_SDSS.txt");
        cout << "* Reading Filter from file " << filterfile_r <<endl;
        Filter filt_r(filterfile_r, lambdamin*1e-9, lambdamax*1e-9);

        string filterfile_i = filtpathmgr.BuildFileName("i_lsst_etc.txt");  // ../filters/g_SDSS.txt");
        cout << "* Reading Filter from file " << filterfile_i <<endl;
        Filter filt_i(filterfile_i, lambdamin*1e-9, lambdamax*1e-9);
        
        Filter & filt_LF = filt_B;  // this is for Dahlen
        if (!fgLFfilB)  filt_LF = filt_iS; // this is for Ramos

        Filter & filt_Obs = filt_i;  // we choose the observation filter in which magLimit is set
        /*
        string filterfile_i_atm = filtpathmgr.BuildFileName("LSST_i_atm.txt");
        //string filterfile_i_atm("test.txt");
        cout << "* Reading Filter from file " << filterfile_i_atm <<endl;
        Filter filt_i_atm(filterfile_i_atm, lambdamin*1.e-9, lambdamax*1.e-9);
        Filter & filter = filt_i_atm;

        */

        cout << "[4] Reading SEDs ..."<<endl;

        string sedfileABzero = sedpathmgr.BuildFileName("ABzero.txt");
        string sedfileEl = sedpathmgr.BuildFileName("function_El_cww_fix2_extend_smooth.txt");
        string sedfileSp1 = sedpathmgr.BuildFileName("function_Scd_cww_fix2_extend_smooth.txt");
        string sedfileSp2 = sedpathmgr.BuildFileName("function_Sbc_cww_fix2_extend_smooth.txt");
        string sedfileSp3 = sedpathmgr.BuildFileName("function_Im_cww_fix2_extend_smooth.txt");
        string sedfileSB1 = sedpathmgr.BuildFileName("function_SB2_kin_fix2_extend_smooth.txt");
        string sedfileSB2 = sedpathmgr.BuildFileName("function_SB3_kin_fix2_extend_smooth.txt");

        
        //cout<<"* Reading SED file for Ellipticals: " << sedfileEl <<endl;
        SED sedEl(sedfileEl,lambdaSEDmin*1e-9, lambdaSEDmax*1e-9);
        
        //cout<<"* Reading SED file for Spirals: " << sedfileSp <<endl;
        //SED sedSp(sedfileSp,lambdaSEDmin*1e-9, lambdaSEDmax*1e-9);
        SED sedSp1(sedfileSp1,lambdaSEDmin*1e-9, lambdaSEDmax*1e-9);
        SED sedSp2(sedfileSp2,lambdaSEDmin*1e-9, lambdaSEDmax*1e-9);
        SED sedSp3(sedfileSp3,lambdaSEDmin*1e-9, lambdaSEDmax*1e-9);
        vector<SED> sp;
        sp.push_back(sedSp1);
        sp.push_back(sedSp2);
        sp.push_back(sedSp3);
        
        //cout<<"* Reading SED file for StarBursts: " << sedfileSB <<endl;
        //SED sedSB(sedfileSB,lambdaSEDmin*1e-9, lambdaSEDmax*1e-9);
        SED sedSB1(sedfileSB1,lambdaSEDmin*1e-9, lambdaSEDmax*1e-9);
        SED sedSB2(sedfileSB2,lambdaSEDmin*1e-9, lambdaSEDmax*1e-9);
        SED sedSB3 = sedSp3;
        vector<SED> sb;
        sb.push_back(sedSB1);
        sb.push_back(sedSB2);
        sb.push_back(sedSB3);

        
        cout<<"* Reading SED file for AB zero point: " << sedfileABzero <<endl;
        SED sedABzero(sedfileABzero,lambdaSEDmin*1e-9, lambdaSEDmax*1e-9);
        
        // Saving the redshoft and the number of galaxies in a file
        string outname = "./Ng-z.dat";
        cout << "* Saving the z and Ng in " << outname << endl;
        //ofstream outp(outname, ios::app | ios::binary);
        ofstream outp(outname.c_str(), ios::out | ios::binary);
        outp << "# z \t NgElliptical \t NgSpiral \t NgSB"<<endl;

        
        const char * names[8] = {"z", "ngel", "ngsp","ngsb","ngall","MBel","MBsp","MBsb"};
        // NTuple (Table) creation with 4 columns (double precision)
        NTuple  nt(8, names);
        double xnt[10];

        double totNgEl=0.;
        double totNgSp=0.;
        double totNgSB=0.;
        double totNgAll=0.; // sum of all types given by the catalogue on all redshifts
        double totNgAllEllCut=0.; // totNgAll with magnitude limit for early types
        double totNgAllSpCut=0.; // same for late types
        double totNgAllSBCut=0.; // same for star bursts

        vector<MySchechter> vsfh;

        // filt_LF : rest-frame -> Lum. functiuons are defined in rest-frame magnitudes
        // filt_Obs : observer-frame
        KCorrector KC(sedABzero, lambdaSEDmin, lambdaSEDmax);

        for(double z=zmin; z<zmax; z+=dz) {
            su.SetEmissionRedShift(z);
            cout << "* z=" << z << endl;
            double dL = su.LuminosityDistanceMpc();
            cout << "* Luminosity distance dL=" << dL << " Mpc" << endl;
            double Xi = su.LineOfSightComovDistanceMpc();
            cout << "* Comoving distance  Xi=" << Xi << " Mpc" << endl;
            su.SetEmissionRedShift(z+dz);
            double deltaXi = su.LineOfSightComovDistanceMpc()-Xi;

            // Computing the absolute B mag limit according to the SED and the apparent magnitude limit in r-band
            cout << "* Computing M_LF (rest-frame Abs mag limit) for the early types in rest-frame/LF filter ... " << endl;
            double magLFlimEl = KC.ComputeRestFrameMagnitudeLimit(filt_LF,filt_Obs,sedEl,sedABzero,lambdamin,lambdamax,dL,maglim,z);

            cout << "* Computing M_LF (rest-frame Abs mag limit) for the late types ... " << endl;

            double magLFlimSp1 = KC.ComputeRestFrameMagnitudeLimit(filt_LF,filt_Obs,sp[0],sedABzero,lambdamin,lambdamax,dL,maglim,z);
            double magLFlimSp2 = KC.ComputeRestFrameMagnitudeLimit(filt_LF,filt_Obs,sp[1],sedABzero,lambdamin,lambdamax,dL,maglim,z);
            double magLFlimSp3 = KC.ComputeRestFrameMagnitudeLimit(filt_LF,filt_Obs,sp[2],sedABzero,lambdamin,lambdamax,dL,maglim,z);
           
            double magLFlimSp = magLFlimSp2; //(magLFlimSp1+magLFlimSp2)/2.;
            cout << "* Computing M_LF (rest-frame Abs mag limit) for the Star bursts ... " << endl;

            double magLFlimSB1 = KC.ComputeRestFrameMagnitudeLimit(filt_LF,filt_Obs,sb[0],sedABzero,lambdamin,lambdamax,dL,maglim,z);
            /*
            double magLFlimcheck = ComputeBlim(filt_LF,filt_Obs,sb[0],sedABzero,lambdamin,lambdamax,dL,maglim,z);
            if (magLFlimcheck != magLFlimSB1)  {
                cout << " **** BUG magLFlimcheck != magLFlimSB1: "<<magLFlimcheck<<" ,"<<magLFlimSB1<<endl;
                throw ParmError(" **** BUG magLFlimcheck != magLFlimSB1");
            }
             */
            double magLFlimSB2 = KC.ComputeRestFrameMagnitudeLimit(filt_LF,filt_Obs,sb[1],sedABzero,lambdamin,lambdamax,dL,maglim,z);
            double magLFlimSB3 = KC.ComputeRestFrameMagnitudeLimit(filt_LF,filt_Obs,sb[2],sedABzero,lambdamin,lambdamax,dL,maglim,z);

            //double magBlimSB = (magBlimSB1+magBlimSB2)/2.;
            double magLFlimSB = magLFlimSB2; //(magLFlimSB1+magLFlimSB2+magLFlimSB3)/3.;
            //double magBlimSB = -13;
            
            cout << "* z="<<z<<" Limit on absolute mag B" << endl;
            cout << "       " << magLFlimEl << " (Early types)" << endl;
            cout << "       " << magLFlimSp << " (Late types)" << endl;
            cout << "       " << magLFlimSB << " (Star Bursts)" << endl;
            
            bool fgeffic = true;
            MySchechter schEl = multLF.getLF_Elliptical(z);
            double NgEl = ComputeNg(schEl,magLFlimEl,Xi,deltaXi,(fgeffic?delmag:-1));
            totNgEl+=NgEl;

            MySchechter schSp = multLF.getLF_Spiral(z);
            double NgSp = ComputeNg(schSp,magLFlimSp,Xi,deltaXi,(fgeffic?delmag:-1));
            totNgSp+=NgSp;

            MySchechter schSB = multLF.getLF_StarBurst(z);
            double NgSB = ComputeNg(schSB, magLFlimSB,Xi,deltaXi,(fgeffic?delmag:-1));
            totNgSB+=NgSB;
            
            MySchechter schAll = multLF.getLF_All(z);
            double NgAllEllCut = ComputeNg(schAll, magLFlimEl,Xi,deltaXi,(fgeffic?delmag:-1));
            double NgAllSpCut = ComputeNg(schAll, magLFlimSp,Xi,deltaXi,(fgeffic?delmag:-1));
            double NgAllSBCut = ComputeNg(schAll, magLFlimSB,Xi,deltaXi,(fgeffic?delmag:-1));
            double somme = NgEl+NgSp+NgSB;
            double NgAll = NgAllEllCut; //(NgAllEllCut*NgEl + NgAllSpCut*NgSp + NgAllSBCut*NgSB)/somme;
            totNgAll+=NgAll;
            
            totNgAllEllCut += NgAllEllCut;
            totNgAllSpCut += NgAllSpCut;
            totNgAllSBCut += NgAllSBCut;

            xnt[0]=z;  xnt[1]=NgEl;  xnt[2]=NgSp;  xnt[3]=NgSB;
            //if(fguseall)
                xnt[4]=NgAll;
            //else xnt[4]=NgEl+NgSp+NgSB;
            xnt[5]=magLFlimEl; xnt[6]=magLFlimSp; xnt[7]=magLFlimSB;
            nt.Fill(xnt);

            cout << "* z="<<z<<" Number of galaxies  Ng:" << endl;
            cout << "       " << NgEl << " (Early types), magBlimEl=" << magLFlimEl << endl;
            cout << "       " << NgSp << " (Late types), magBlimSp=" << magLFlimSp << endl;
            cout << "       " << NgSB << " (Star Bursts), magBlimSB=" << magLFlimSB << endl;
            
        
            outp << z << " \t " << NgEl << " \t " << NgSp << " \t " << NgSB << endl;
        }  // Fin boucle sur les redshifts z
        cout<<nt;
        POutPersist po(outppfname);
        cout<<" Saving NTuple to file "<< outppfname<<endl;
        po << PPFNameTag("Ng") << nt;
        
        cout << " NgEl=" << totNgEl << ", NgSp=" << totNgSp << ", NgSB=" << totNgSB << endl;
        double sum = totNgEl+totNgSp+totNgSB;
        cout << " NgAll=" << totNgAll << " per arcmin^2" << endl;
        cout << " Sum of all types=" << sum << " per arcmin^2" << endl;
        cout << "NgAllEllCut=" << totNgAllEllCut << ", NgAllSpCut=" << totNgAllSpCut << ", NgAllSBCut=" << totNgAllSBCut << endl;
       cout << " Fel=" << totNgEl/sum << ", Fsp=" << totNgSp/sum << ", Fsb=" << totNgSB/sum << endl;
        
        cout << "* magnitude limit=" << maglim << ", delmag=" << delmagi << endl;
    
    }  // End of try bloc
    
    
    catch (PThrowable & exc) {  // catching SOPHYA exceptions
        cerr << " tngal.cc: Catched Exception (PThrowable)"
        << (string)typeid(exc).name()
        << "\n --> exc.Msg= " << exc.Msg() << endl;
        rc = 99;
    }
    catch (std::exception & e) {  // catching standard C++ exceptions
        cerr << " tngal.cc: Catched std::exception "  << " - what()= "
        << e.what() << endl;
        rc = 98;
    }
    catch (...) {  // catching other exceptions
        cerr << " tngal.cc: some other exception (...) was caught ! "
        << endl;
        rc = 97;
    }
    cout << endl;
    cout << " rc=" << rc << endl;
    cout << " ---------- End ---------- " << endl;
    return rc;	
}

