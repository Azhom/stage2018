#ifndef REDDENING_H_SEEN 
#define REDDENING_H_SEEN 

#include <math.h>
#include <iostream>
#include <vector>
#include <string>

// sophya
#include "machdefs.h"
#include "sopnamsp.h"
#include "pexceptions.h"
#include "cspline.h"
#include "randfmt.h"
#include "randr48.h"


class Reddening
{
public:
	//Constructor
	Reddening(double rv=3.1);
	Reddening(const Reddening &red);
	double A_lambda(double lambda, bool use_random = false);
	double randomize(double sigma=.0, int nb_points=5);
	//Absoprtion in lambda wavelength using Cardelli, Clayton, Mathis 1989 extinction law
	double cardelli_A_lambda(double lambda);
	//Absoprtion in lambda wavelength using  Fitzpatrick 1999 extinction law
	double fitzpatrick_A_lambda(double lambda);
	//change value of R_V
	void update_rv(double new_rv);
	void update_ebv(double new_ebv);
protected:
	double rv_;
	double ebv_;
	string law_;
	CSpline modifier_;
};

#endif
