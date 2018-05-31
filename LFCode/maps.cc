#include "maps.h"

SphereHEALPix<T> galcnt_map(SphereHEALPix<T> &ebmv_map, SInterp1D &ngal_interpolated){
	
	SphereHEALPix<r_4> diff_map;
	for(int ii = 0; ii<ebmv_map.NbPixels(); ii++) {
		if (ebmv_map[ii]<7.35){
			totcnt_reddened 	= ngal_interpolated(ebmv_map[ii]);
			diff_map[ii] 		= totcnt_reddened;
		}
		else{
			diff_map[ii] = 0.;
		}
	}
	return diff_map;
}

SphereHEALPix<T> galcnt_corrected(SphereHEALPix<T> &ebmv_map, SInterp1D &true_interpolation, SInterp1D &corr_interpolation){
	return diff_map;
}

void create_interpolated_points(string &filename, double ebmv_max=7.35, double ebmv_min=0., int ebmv_steps=1000){
	vector<double> ngal_vector;
	vector<double> ebmv_vector;
	double ebmv;
	
	ofstream ngal_points;
	ngal_points.open(filename);
	
	for(int ii=0; ii<ebmv_steps; ii++){
		ebmv = (ebmv_max-ebmv_min)/ebmv_steps*ii+ebmv_min;
		galcntc.doCompute(zmin, zmax, dz, maglim, magerr, lambdamin,lambdamax, ebmv);
		ngal_vector.push_back(galcntc.getIntegratedGalDensity_Arcmin2(Ellcnt, Spcnt, SBcnt));
		ebmv_vector.push_back(ebmv);
		ngal_points << ebmv_vector[ii] << " " << ngal_vector[ii] << endl;
		cout << ebmv_vector[ii] << endl;
	}
	ngal_points.close();
}
