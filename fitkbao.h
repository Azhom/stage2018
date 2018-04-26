/* ----
   Project   LSST/BAO/PhotoZ
   R.Ansari , Mars 2018 
                                                     -------  */
#ifndef FITKBAO_SEEN
#define FITKBAO_SEEN

#include "machdefs.h"
#include "sopnamsp.h"
#include "generalfit.h"

int doFitkBAO_A(GeneralFitData& gdata, double & xi2r, 
		double & A, double & errA, double & sdamp, double & errsdamp, double & sbao, double & errsbao, int prtlev=1);

int doFitkBAO_B(GeneralFitData& gdata, double & xi2r, 
		double & A, double & errA, double & sdamp, double & errsdamp, double & sbao, double & errsbao, int prtlev=1);

////------- Pour le fit :
class MykBAOOScA : public GeneralFunction
{
public:
  // We define here a function of a single variable (M) with three parameters (phistar, Mstar and alpha)
  // a Scheshter function
  MykBAOOScA() : GeneralFunction(1,4)
  {
  }

  // Return the value of the function for xp[] and parameters parm[]
  virtual double Value (double const xp[], double const *parm)
  {
      double A = parm[0];
      double sdump = parm[1];
      double alpha = parm[2];
      double sbao = parm[3];
      double k = xp[0];
      
      double fdump=exp(-pow(k*sdump,alpha));
      double sik=sin(k*sbao);
      return 1.+A*k*fdump*sik;
  }

  // return the value of the function and its derivatives for the point xp[], the parameters parm 
  virtual double Val_Der (double const xp[], double const *parm, double *DgDpar)
  {
      double A = parm[0];
      double sdump = parm[1];
      double alpha = parm[2];
      double sbao = parm[3];
      double k = xp[0];
 
      double fdump=exp(-pow(k*sdump,alpha));
      double sik=sin(k*sbao);
      double cok=cos(k*sbao);

      // Derivatives
      DgDpar[0]= k*fdump*sik; // with respect to A
      DgDpar[1]= A*k*fdump*(-alpha*k*pow(k*sdump,alpha-1.))*sik;  // with respect to sdump 
      DgDpar[2]= A*k*fdump*(-pow(k*sdump,alpha)*log(k*sdump))*sik;  // with respect to alpha 
      DgDpar[3]= A*k*fdump*cok*k;  // with respect to sbao 

      return 1.+A*k*fdump*sik;
  }
};

////------- Pour le fit 2 :
class MykBAOOScB : public GeneralFunction
{
public:
  // We define here a function of a single variable (M) with three parameters (phistar, Mstar and alpha)
  // a Scheshter function
  MykBAOOScB() : GeneralFunction(1,3)
  {
  }

  // Return the value of the function for xp[] and parameters parm[]
  virtual double Value (double const xp[], double const *parm)
  {
      double A = parm[0];
      double sdump = parm[1];
      double sbao = parm[2];
      double k = xp[0];
      
      double fdump=exp(-sdump*k);
      double sik=sin(k*sbao);
      return 1.+A*fdump*sik;
  }

  // return the value of the function and its derivatives for the point xp[], the parameters parm 
  virtual double Val_Der (double const xp[], double const *parm, double *DgDpar)
  {
      double A = parm[0];
      double sdump = parm[1];
      double sbao = parm[2];
      double k = xp[0];
 
      double fdump=exp(-sdump*k);
      double sik=sin(k*sbao);
      double cok=cos(k*sbao);

      // Derivatives
      DgDpar[0]= fdump*sik; // with respect to A
      DgDpar[1]= -A*fdump*k*sik;  // with respect to sdump 
      DgDpar[2]= A*fdump*cok*k;  // with respect to sbao 

      return 1.+A*fdump*sik;
  }
};


#endif
