#include "reddening.h"

Reddening::Reddening(){};

double Reddening::compute_A_lambda(double lambda, double EBmV)
{
	double A_V;
	double A_lambda;
	double x, y;
	double a, b;
	double F_a, F_b;
	
	x = 1/(lambda*1e6);
	y = x - 1.82;
	
	if (x <= 1.1){
		a= 0.574 * pow(x, 1.61);
		b=-0.527 * pow(x, 1.61);
	}
	else if (x > 1.1 && x <= 3.3){
		a = 1 + 0.17699*y - 0.50447*y*y - 0.02427*pow(y,3)+ 0.72085*pow(y,4) + 0.01979*pow(y,5) - 0.77530*pow(y,6) + 0.32999*pow(y,7);
		
		b = 1.41338*y + 2.28305*y*y + 1.07233 *pow(y,3) - 5.38434*pow(y,4) - 0.62251*pow(y,5) + 5.30260*pow(y,6) - 2.09002*pow(y,7);
	}
	else if (x > 3.3 && x <= 8){
		if (x <= 5.9){
			F_a = 0;
			F_b = 0;
		}
		else {
			F_a = -0.04473*pow(x-5.9, 2) - 0.009779*pow(x-5.9, 3);
			F_b = 0.2130*pow(x-5.9, 2) + 0.1207*pow(x-5.9, 3);
		}
		
		a = 1.752 - 0.316*x - 0.104/(pow(x-4.67, 2) + 0.341) + F_a;
		b = -3.090 + 1.825*x + 1.206/(pow(x-4.62, 2) + 0.263) + F_b;
	}
	else if (x > 8) {
		a = -1.073 - 0.628*(x-8) + 0.137*pow(x-8, 2) - 0.070*pow(x-8, 3);
		b = 13.670 + 4.257*(x-8) - 0.420*pow(x-8, 2) + 0.374*pow(x-8, 3);
	}
	else
	{
		cout << "Error: x is not a number." << endl;
	}
	A_lambda = (a * 3.1 + b) * EBmV;
	return A_lambda;
};
