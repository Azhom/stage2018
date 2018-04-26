/* ----
   Project   LSST/BAO/PhotoZ
   Classes to represent Schechter functions - 
   would be integrated in SOPHYA later 
   F. Habibi & R.Ansari ,  LAL (CNRS/IN2P3)  & Univ. Paris Sud 
   April - July 2017                                          -------  */
    
#ifndef MYSCHECHTER_H_SEEN
#define MYSCHECHTER_H_SEEN

#include <string>
#include <vector>

#include <iostream> 
#include "classfunc.h"

using namespace std;
using namespace SOPHYA;

//--------------------------------------------------------------------------------
// We need a class to handle Schechter functions - Maybe rename this class ....
//--------------------------------------------------------------------------------

class MySchechter : public ClassFunc1D  {
public:
  MySchechter()
    : phistar_(0.01), Mstar_(-20), alpha_(-1)
    {
        A = 0.4*log(10.);
    }

  MySchechter(double phistar, double Mstar, double alpha)
    : phistar_(phistar), Mstar_(Mstar), alpha_(alpha)
  {
    A = 0.4*log(10.);  
  }
  
  MySchechter( MySchechter const & a )
    : phistar_(a.phistar_), Mstar_(a.Mstar_), alpha_(a.alpha_), A(a.A)
  {
    
  }
    
  MySchechter& operator=(MySchechter const & a )
  {
    phistar_=a.phistar_; Mstar_=a.Mstar_;
    alpha_=a.alpha_; A=a.A;
    return *this;
  }

  // return the LF value (number density per abs-mag) for a given abs-mag Mag
  virtual double operator()(double Mag) const
  {
    return Value(Mag);
  }
  
  // computes the LF
  inline double Value(double Mag) const
  {
    double X=pow(10., 0.4*(Mstar_-Mag));
    double phiret= A * phistar_ * pow(X,alpha_+1);// * exp(-X) ; // The Shceshter function
    return phiret;
  }
  // computes the LF
  inline double LumWeightedValue(double Mag) const
  {
    // luminosity for a magnitude Mag
    double Lum=pow(10., -0.4*Mag);
    double X=pow(10., 0.4*(Mstar_-Mag));
    double valret= Lum* A * phistar_ * pow(X,alpha_+1);// * exp(-X) ;
    return valret;
  }
  
  // to change the phiStar of LF
  inline void setPhistar(double phis)
  { phistar_=phis; }
  
  // to scale the phiStar of LF
  inline void scalePhistar(double scale)
  { phistar_*=scale; }
  
  // return the phiStar of LF
  inline double getPhistar()
  {  return phistar_;  }
  
  // return the MStar of LF
  inline double getMstar()
  { return Mstar_; }

  // return the Alpha of LF
  inline double getAlpha()
  {  return alpha_;  }
  
// return the integral of the LF function for a given abs-mag interval
  double getIntegral(double magmin, double magmax, int glorder=100);
// return the integral of the luminosity weighted LF function for a given abs-mag interval
  double getLumWeightedIntegral(double magmin, double magmax, int glorder=100);
// prints the Schechter parameters , if fgl=true
  std::ostream& print(std::ostream& os, bool fgl=false)  const;

protected:
  double phistar_, Mstar_, alpha_, A; 
};

/*! operator << overloading - Prints Schechter parameters in string format on \b os*/
inline std::ostream& operator << (std::ostream& s, MySchechter const & sch)
  {  sch.print(s);  return(s);  }


#endif  /*  Fin de MYSCHECHTER_H_SEEN */
