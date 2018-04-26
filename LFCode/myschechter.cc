/* ----
   Project   LSST/BAO/PhotoZ
   Classes to represent Schechter functions - 
   would be integrated in SOPHYA later 
   F. Habibi & R.Ansari ,  LAL (CNRS/IN2P3)  & Univ. Paris Sud 
   April - July 2017                                          -------  */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <iostream>
#include <iomanip>

#include "myschechter.h"
#include "integ.h"

using namespace std;
using namespace SOPHYA;

//------------------------------------------------------------------------
//-----------  class MySchechter -----------------------------------------
//------------------------------------------------------------------------

double MySchechter::getIntegral(double magmin, double magmax, int glorder)
{
  GLInteg integ((*this), magmin, magmax);  
  integ.SetOrder(glorder);
  
  return integ.Value();
}

class MyLWSchechter : public ClassFunc1D
{
public:
  MyLWSchechter(MySchechter const & mysch) : mysch_(mysch) { }
  // return the LF value weighted by luminosity for a given abs-mag Mag
  virtual double operator()(double Mag) const
  {
    return mysch_.LumWeightedValue(Mag);
  }
  MySchechter const & mysch_;
};

double MySchechter::getLumWeightedIntegral(double magmin, double magmax, int glorder)
{
  MyLWSchechter lwsc((*this));
  GLInteg integ(lwsc, magmin, magmax);
  integ.SetOrder(glorder);
  return integ.Value();
}

std::ostream& MySchechter::print(std::ostream& os, bool fgl)  const
{
  if (fgl)  os << "phistar="<<phistar_<<" Mstar="<< Mstar_<<" alpha="<<alpha_;
  else os << setw(12) << phistar_ << '\t' << setw(12) << Mstar_ << '\t' << setw(12) << alpha_ << " ";
  return os;
}



