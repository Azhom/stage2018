/* ----
   Project   LSST/BAO/PhotoZ
   R.Ansari , Mars 2018 
                                                     -------  */
//----- c/c++ includes 
#include <iostream>
#include <fstream>
#include <string>

#include "fitkbao.h"

/*-- Fonction --*/
int doFitkBAO_A(GeneralFitData& gdata, double & xi2r, 
		double & A, double & errA, double & sdamp, double & errsdamp, double & sbao, double & errsbao, int prtlev)
{
  A=9999.;  errA=-99.;
  sdamp=9999.;  errsdamp=-99.;
  sbao=9999.;  errsbao=-99.;
  xi2r=-99999.;

  if (prtlev>1) {
    cout << " ------------------------------------------------------------------------" << endl;
    cout << " ------- doFitkBAO_A() 1+A*k*Exp(-pow(k*sdump,alpha))*sin(k*sbao)  ------" << endl;
  }
  MykBAOOScA mFunction;
  GeneralFit mFit(&mFunction);  // create a fitter for the choosen function
  mFit.SetData(&gdata);        // connect data to the fitter
  // Set and initialise the parameters (that's non-linear fitting!)
  //           (parameter number, parameter name, initial guess, step, [limits min and max])
  //    mFit.SetParam(0,"phistar",phistar*0.2,0.005*phistar);
  // setting the itteration for xi2 minimising
  mFit.SetMaxStep(1000);
  //mFit.SetLambda_Fac(.0001);
  mFit.SetParam(0,"A",5,0.1,0.1,50.);
  mFit.SetParam(1,"SDump",10.,0.5,0.05,100.);
  //  mFit.SetFix (1,20.);
  mFit.SetParam(2,"alpha",1.5,0.05,0.2,3.);
  mFit.SetFix (2,1.4);
  mFit.SetParam(3,"SBAO",145.,0.5,50.,250.);
  
  // Fit and print result
  int rcfit = mFit.Fit();
  if (prtlev>2)  mFit.PrintFit();
  
  if(rcfit>0) {
    if (prtlev>0) cout<< "doFitkBAO_A() Reduce_Chisquare = " << mFit.GetChi2Red()<< " nstep="<<mFit.GetNStep() << " rc="<<rcfit<<endl;
    // Get the result of fitted paraneters with fits' errors
    TVector<r_8> ParResult = mFit.GetParm();
    A=mFit.GetParm(0);   errA=mFit.GetParmErr(0);
    sdamp=mFit.GetParm(1);   errsdamp=mFit.GetParmErr(1);
    sbao=mFit.GetParm(3);  errsbao=mFit.GetParmErr(3);
    xi2r=mFit.GetChi2Red();
    if (prtlev>1) {
      cout << " A= "<< mFit.GetParm(0)<<" +/- "<<mFit.GetParmErr(0);
      cout << " SDamp= "<< mFit.GetParm(1)<<" +/- "<<mFit.GetParmErr(1);
      //      cout << " alpha= "<< mFit.GetParm(2)<<" Err:"<<mFit.GetParmErr(2) << endl; 
      cout << " SBAO= "<< mFit.GetParm(3)<<" +/- "<<mFit.GetParmErr(3) << endl; 
    }
  }
  else {
    if (prtlev>0) cout << "doFitkBAO_A() Fit_Error, rc = " << rcfit << "  nstep="<<mFit.GetNStep()<<endl;
    if (prtlev>1) mFit.PrintFitErr(rcfit);
  } 
  return rcfit;
}

/*-- Fonction --*/
int doFitkBAO_B(GeneralFitData& gdata, double & xi2r, 
		double & A, double & errA, double & sdamp, double & errsdamp, double & sbao, double & errsbao, int prtlev)
{
  A=9999.;  errA=-99.;
  sdamp=9999.;  errsdamp=-99.;
  sbao=9999.;  errsbao=-99.;
  xi2r=-99999.;

  if (prtlev>1) {
    cout << " ------------------------------------------------------------------------" << endl;
    cout << " ------- doFitkBAO_B() 1+A*k*Exp(-pow(k*sdump,alpha))*sin(k*sbao)  ------" << endl;
  }
  MykBAOOScB mFunction;
  GeneralFit mFit(&mFunction);  // create a fitter for the choosen function
  mFit.SetData(&gdata);        // connect data to the fitter
  // Set and initialise the parameters (that's non-linear fitting!)
  //           (parameter number, parameter name, initial guess, step, [limits min and max])
  //    mFit.SetParam(0,"phistar",phistar*0.2,0.005*phistar);
  // setting the itteration for xi2 minimising
  mFit.SetMaxStep(1000);
  //mFit.SetLambda_Fac(.0001);
  mFit.SetParam(0,"A",0.1,0.005,0.005,1.);
  mFit.SetParam(1,"SDump",10.,0.5,0.05,100.);
  //    mFit.SetFix (1,10.);
  mFit.SetParam(2,"SBAO",145.,0.5,50.,250.);
  
  // Fit and print result
  int rcfit = mFit.Fit();
  if (prtlev>2)  mFit.PrintFit();
  
  if(rcfit>0) {
    if (prtlev>0) cout<< "doFitkBAO_B() Reduce_Chisquare = " << mFit.GetChi2Red()<< " nstep="<<mFit.GetNStep() << " rc="<<rcfit<<endl;
    // Get the result of fitted paraneters with fits' errors
    TVector<r_8> ParResult = mFit.GetParm();
    A=mFit.GetParm(0);   errA=mFit.GetParmErr(0);
    sdamp=mFit.GetParm(1);   errsdamp=mFit.GetParmErr(1);
    sbao=mFit.GetParm(2);  errsbao=mFit.GetParmErr(2);
    xi2r=mFit.GetChi2Red();
    if (prtlev>1) {
      cout << " A= "<< mFit.GetParm(0)<<" +/- "<<mFit.GetParmErr(0);
      cout << " SDamp= "<< mFit.GetParm(1)<<" +/- "<<mFit.GetParmErr(1);
      cout << " SBAO= "<< mFit.GetParm(2)<<" +/- "<<mFit.GetParmErr(2) << endl; 
    }
  }
  else {
    if (prtlev>0) cout << "doFitkBAO_B() Fit_Error, rc = " << rcfit << "  nstep="<<mFit.GetNStep()<<endl;
    if (prtlev>1) mFit.PrintFitErr(rcfit);
  } 
  return rcfit;
}
