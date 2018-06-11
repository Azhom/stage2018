#ifndef MAPS_H_SEEN
#define MAPS_H_SEEN

//#include "machdefs.h"
//#include "sopnamsp.h"
#include <iostream>
#include <fstream>
#include <string>

#include "spherehealpix.h"
#include "tvector.h"
#include "slininterp.h"
#include "sphericaltransformserver.h"
//#include "samba.h"

using namespace std;
using namespace SOPHYA;

template <class T>
TVector<T> errcl(SphereHEALPix<T> &diff_ebv, SphereHEALPix<T> &ebv1, SphereHEALPix<T> &ebv2, vector<int> &mask_indices, SLinInterp1D& ngal_interpolated, double angle_cut=20.);

#endif
