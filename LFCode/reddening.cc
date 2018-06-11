#include "reddening.h"

Reddening::Reddening(double rv){
	rv_ = rv;
	law_ = "f99";
	ebv_ = 0;
};

Reddening::Reddening(const Reddening &red)
{
	modifier_ = red.modifier_;
	rv_ = red.rv_;
	law_ = red.law_;
	ebv_ = red.ebv_;
};

double Reddening::A_lambda(double lambda, bool use_random)
{
	double A_lambda;
	if(law_ == "ccm"){
		A_lambda = cardelli_A_lambda(lambda);
	}
	else if(law_ == "f99"){
		A_lambda = fitzpatrick_A_lambda(lambda);
	}
	if(use_random){
		A_lambda += modifier_.CSplineInt(1/lambda*1e-6)*A_lambda;
	}
	if (A_lambda<0){return 0;}
	return A_lambda;
};

double Reddening::randomize(double sigma, int nb_points){
	ThSDR48RandGen rng;
	rng.AutoInit();
	//rng.SetSeed(123456789);
	int nb_pts = 5;
	double mod_x[nb_pts+1];
	double mod_values[nb_pts+1];
	double xmin=0., xmax=10., curr_x;
	mod_x[0] = -1e-6;
	mod_values[0] = 0.0;
	for(int ii=1; ii<nb_pts+1; ii++){
		curr_x = xmin + ii*(xmax-xmin)/nb_pts;
		mod_x[ii] = curr_x;
		mod_values[ii] = rng.Gaussian(sigma);
		cout << mod_values[ii] << endl;
		//if(mod_values[ii]<0){mod_values[ii] = 0;}
	}
	modifier_.SetNewTab(nb_pts+1, mod_x, mod_values, true);
	modifier_.ComputeCSpline();
};

double Reddening::cardelli_A_lambda(double lambda)
{
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
		a = 1 + 0.17699*y - 0.50447*y*y - 0.02427*pow(y,3) + 0.72085*pow(y,4) + 0.01979*pow(y,5) - 0.77530*pow(y,6) + 0.32999*pow(y,7);
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
	A_lambda = (a * rv_ + b) * ebv_;
	return A_lambda;
};

double Reddening::fitzpatrick_A_lambda(double lambda)
{
	double A_lambda;
	double x = 1/(lambda*1e6); //in um-1
	
	if(ebv_==0.){
		return 0.;
	}
	
	if (x>=3.703){
		//if lambda < 2700 e-6 m we use FM90 law
		
		double c2 = -0.824 + 4.717/rv_;
		double c1 = 2.030 - 3.007*c2;
		double x0 = 4.596, gamma=0.99, c3=3.23, c4=0.41;
		
		double D = x*x / (pow(x*x - x0*x0, 2) + x*x*gamma*gamma);
		double F = 0.;
		if(x >= 5.9){
			F = 0.5392*pow(x-5.9, 2) + 0.05644*pow(x-5.9, 3);
		}
		
		double k = c1 + c2*x + c3*D + c4*F;
		A_lambda = (k+rv_)*ebv_;
	}
	else{
		double x_points1[9] = {0., 0.377, 0.820, 1.667, 1.828, 2.141, 2.433, 3.704, 3.846};
		double y_points1[9] = {0., 0.265, 0.829, -0.426+1.0044*rv_, -0.05+1.0016*rv_, 0.701+1.0016*rv_, 1.208+1.0032*rv_-0.00033*rv_*rv_, this->fitzpatrick_A_lambda(2700e-10)/ebv_, this->fitzpatrick_A_lambda(2600e-10)/ebv_};
		CSpline optical_IR(9, x_points1, y_points1, 0., 0., 3, false);
		optical_IR.ComputeCSpline();
		A_lambda = optical_IR.CSplineInt(x)*ebv_;
	}
	return A_lambda;
}

void Reddening::update_rv(double new_rv){
	rv_ = new_rv;
}

void Reddening::update_ebv(double new_ebv){
	ebv_ = new_ebv;
}

void Reddening::print_law(int suffix){
	ofstream file;
	stringstream ss;
	ss<<"Objs/extlaw" <<suffix<<".txt";
	file.open(ss.str().c_str());
	double x_min = 0.0;
	double x_max = 8.45;
	double dx = 0.02;
	double a;
	
	for(double m = x_min; m<x_max; m+=dx){
		a = fitzpatrick_A_lambda(1/m*1e-6);
		file << m << " " <<  a;
		for(int e=0; e<1; e++){
			file << " " << a + modifier_.CSplineInt(m)*a << " " << modifier_.CSplineInt(m);
		}
		file << endl;
	}
	file.close();
}
