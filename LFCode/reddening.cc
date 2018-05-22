#include "reddening.h"

Reddening::Reddening(double rv){
	rv_ = rv;
};

double Reddening::compute_A_lambda(double lambda, double EBmV)
{
	double A_V;
	double A_lambda;// In meters
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
	A_lambda = (a * rv_ + b) * EBmV;
	return A_lambda;
};

double Reddening::fitzpatrick_A_lambda(double lambda, double ebmv){
	
	double A_lambda;
	double x = 1/(lambda*1e6); //in um-1
	
	if (lambda < 2700e-6){
		//if lambda < 2700 e-6 m we use FM90 law
		
		double c2 = -0.824 + 4.717/rv_;
		double c1 = 2.030 - 3.007*c2;
		double x0 = 4.596e-6, gamma=0.99e-6, c3=3.23, c4=0.41;
		
		double D = x*x / (pow(x*x-x0*x0) + x*x+gamma*gamma);
		double F = 0.;
		if(x >= 5.9){
			F = 0.5392*pow(x-5.9, 2) + 0.05644*pow(x-5.9, 3);
		}
		
		double k = c1 + c2*x +c3*D+c4*F;
		A_lambda = k - rv_;
	}
	else{
		vector<int> x_points = [0, 0.377, 0.820, 1.667, 1.828, 2.141, 2.433, 3.704, 3.846];
		vector<int> y_points = [0, 0.265, 0.829, -0.426+1.0044*rv_, -0.05+1.0016*rv_, 0.701+1.0016*rv_, -1.208+1.0032*rv_-0.00033*rv_*rv_, 6.265, 6.591]
		CSpline optical_IR(x_points.size(), x_points, y_points);
		A_lambda = optical_IR.CSplineInt(x);
	}
	return A_lambda;
}
