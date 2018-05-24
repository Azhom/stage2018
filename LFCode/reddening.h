#ifndef REDDENING_H_SEEN 
#define REDDENING_H_SEEN 

#include <math.h>
#include <iostream>
#include <vector>

// sophya
#include "machdefs.h"
#include "sopnamsp.h"
#include "pexceptions.h"
#include "cspline.h"

class Reddening
{
public:
	//Constructor
	Reddening(double rv=3.1);
	//Absoprtion in lambda wavelength using Cardelli, Clayton, Mathis 1989 extinction law
	double cardelli_A_lambda(double lambda, double EBmV);
	//Absoprtion in lambda wavelength using  Fitzpatrick 2007 extinction law
	double fitzpatrick_A_lambda(double lambda, double ebmv);
	//change value of R_V
	void update_rv(double new_rv);
protected:
	double rv_;
};

#endif
