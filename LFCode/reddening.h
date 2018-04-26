#ifndef REDDENING_H_SEEN 
#define REDDENING_H_SEEN 

#include <math.h>
#include <iostream>

// sophya
#include "machdefs.h"
#include "sopnamsp.h"
#include "pexceptions.h"

class Reddening
{
public:
	Reddening();
	double compute_A_lambda(double lambda, double EBmV);
};

#endif
