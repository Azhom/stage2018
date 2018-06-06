#ifndef MAPS_H_SEEN
#define MAPS_H_SEEN

#include "machdefs.h"
#include "sopnamsp.h"
#include <iostream>
#include <fstream>
#include <string>

#include "spherehealpix.h"
#include "tvector.h"
#include "samba.h"

using namespace std;
using namespace SOPHYA;

SphereHEALPix<T> errcl(SphereHEALPix& diff_ebv, SphereHEALPix& ebv1, SphereHEALPix& ebv2, double angle_cut=20.);
	 
#endif
