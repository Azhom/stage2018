/*  
in piapp 
c++import SkyMap
c++import Samba
delobjs *
c++execfrf tsht.cc
*/
// Generate an input spectra  a + b* l + c * gaussienne(l, 50, 20)
int lmax = 192;
Vector clin(lmax+1);
for(int l=0; l<=lmax; l++) {
  double xx = (l-50.)/10.;
  clin(l) = 1.e-2*(1.-(double)l/(double)lmax) + 0.1*exp(-xx*xx);
}
// Create a set of Alm coefficients distributed according to C(l)=clin
// Alm<double> alm(clin, 10., true);  // single-sided
Alm<double> alm(lmax);  // il s agit d un alm symmetrique; seul alm , m>=0 garde
alm.GenFromCl(clin);
cout<<alm<<endl;
Vector clck=alm.powerSpectrum();
// Instanciate an SHT server object 
SphericalTransformServer<r_8> shts;
// Synthetize a HEALPix and SphereThetaPhi map from alm
int nside = 64;  // HealPix pixelisation parameter
SphereHEALPix<r_8> mapH(nside);
shts.GenerateFromAlm(mapH, nside, alm); 
cout<<" Computed HEALPix mapH from alm "<<endl;
cout<<mapH;
// Decompose map into Ylm 
Alm<double> almc(lmax, true);   // single-sided 
shts.DecomposeToAlm(mapH, almc, lmax, 0.);
cout<<" Computed almc computed from HEALPix mapH"<<endl;
// Get the power spectrum
Vector clc=almc.powerSpectrum();
/*
// Synthetize also a SphereThetaPhi map from the same set of alm
int nring = 257;  // number of rings (=theta slices) 
SphereThetaPhi<r_8> mapT(nring);
shts.GenerateFromAlm(mapT, nring, alm); 
cout<<" Computed SphereThetaPhi mapT from alm "<<endl;
cout << mapT;
*/
KeepObj(clin);
KeepObj(alm);
KeepObj(clck);
KeepObj(clc);
KeepObj(mapH);
cout << " kept clin alm clck clc mapH ..."<<endl;
// On applique une coupure simple angle en theta
/* SphereHEALPix<r_8> mapHC(nside);
   mapHC=mapH;
*/
SphereHEALPix<r_8> mapHC(mapH, false);  // false -> on ne partage pas les pixels 
for(int k=0; k<mapHC.NbPixels(); k++) {
  double tet,phi;
  mapHC.PixThetaPhi(k,tet,phi);
  if (Angle(tet).ToDegree()>30)  mapHC[k]=0.;
}
Alm<double> almcut(lmax, true);   // single-sided 
shts.DecomposeToAlm(mapHC, almcut, lmax, 0.);
cout<<" Computed almcut computed from HEALPix mapHC"<<endl;
// Get the power spectrum
Vector clcut=almcut.powerSpectrum();
KeepObj(clcut);
KeepObj(mapHC);
cout << " kept clcut mapHC ..."<<endl;


// Keepobj(mapT);

