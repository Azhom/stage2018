#include "maps.h"

template <class T>
TVector<T> errcl(SphereHEALPix<T> &diff_ebv, SphereHEALPix<T> &ebv1, SphereHEALPix<T> &ebv2, vector<int> &mask_indices, SLinInterp1D& ngal_interpolated, double angle_cut)
{
	double theta, phi;
	for(int e=0; e<ebv1.NbPixels(); e++){
		diff_ebv.PixThetaPhi(e, theta, phi);
		if(theta>(90-angle_cut)*M_PI/180.0 and theta<(90+angle_cut)*M_PI/180.0){
			diff_ebv[e] = 0;
		}
	}

	int nside = ebv1.SizeIndex();
	int lmax = 3*nside-1;
	//cout << "Computing Cls" << endl;
	SphericalTransformServer<r_8> sts;
	SphereHEALPix<r_8> err_ebv(nside, true);
	SphereHEALPix<r_8> mod_ebv(nside, true);
	TVector<double> err_cl = sts.DecomposeToCl(diff_ebv, lmax, 0.);
	sts.GenerateFromCl(err_ebv, nside, err_cl, 0.);
	
	//cout << "Computing mod_ebv" << endl;
	for(int ii=0; ii<mod_ebv.NbPixels(); ii++){
		mod_ebv[ii] = abs((ebv2[ii]+ebv1[ii])/2 - err_ebv[ii]/2);
	}
	
	SphereHEALPix<r_8> ngal_map(nside);
	//cout << "Computing galmap" << endl;
	for(int ii = 0; ii<mod_ebv.NbPixels(); ii++){
		ngal_map[ii] = ngal_interpolated(mod_ebv[ii]) / ngal_interpolated((ebv2[ii]+ebv1[ii])/2);
	}
	//cout << "Masking" << endl;
	for(int e=0; e<mask_indices.size(); e++){
		ngal_map[mask_indices[e]] = 0; //mask
	}
	
	//cout << "Decomposing" << endl;
	for(int e=0; e<ebv1.NbPixels(); e++){
		ngal_map.PixThetaPhi(e, theta, phi);
		if(theta>70*M_PI/180.0 and theta<110*M_PI/180.0){
			ngal_map[e] = 0;
		}
	}
	TVector<double> cls = sts.DecomposeToCl(ngal_map, lmax, 0.);
	return cls;
};
