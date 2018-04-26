/* ----   
   Exemple d'integration numerique avec SOPHYA/NTools 
   pour Jeremy -   Fichier jinteg.h 
   execution par runcxx :
   runcxx -inc jinteg.h -f jinteg.cc
*/

#include <math.h>
#include "integ.h"  // include pour les integrateurs de SOPHYA

// On definit une classe heritant de ClassFunc1D
class myFuncPoly : public ClassFunc1D {
public:
  myFuncPoly(double a, double b, double c=0.)
    : a_(a), b_(b), c_(a) { }
  // Il faut coder cet operateur ()(double x)
  virtual double operator()(double x)  const
  {
    // on code un polynome a x^2 + b x + c
    
    return a_*x*x+b_*x+c_; 
  }
  double a_, b_, c_; 
};

// On definit une 2eme classe heritant de ClassFunc1D
class myFuncGauss : public ClassFunc1D {
public:
  myFuncGauss(double x0, double sigma, double A=1.)
    : x0_(x0), sigma_(sigma), A_(A) { }
  // Il faut coder cet operateur ()(double x)
  virtual double operator()(double x)  const
  {
    // on code une gaussienne
    x=(x-x0_)/sigma_;
    return A_*exp(-0.5*x*x); 
  }
  double x0_, sigma_, A_;
};
