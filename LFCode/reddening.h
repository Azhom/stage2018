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
	Reddening(double rv=3.1);
	double compute_A_lambda(double lambda, double EBmV);
	double fitzpatrick_A_lambda(double lambda, double ebmv);
protected:
	double rv_;
};

#endif
