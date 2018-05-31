#ifndef MAPS_H_SEEN
#define MAPS_H_SEEN

#include <iostream>
#include <fstream>
#include <string>

using namespace std;
using namespace SOPHYA;


SphereHEALPix<T> galcnt_map(SphereHEALPix<T> &ebmv_map);
SphereHEALPix<T> galcnt_corrected(SphereHEALPix<T> &ebmv_map, SInterp1D &true_interpolation, SInterp1D &corr_interpolation);
void create_interpolated_points(string &filename);
	 
#endif
