/* ----
   Project   LSST/BAO/PhotoZ
   Class to compute power spectrum from from mass density of 
   galaxy count number grid (array) 
   R.Ansari for A. Choyer , Feb 2015 
                                                     -------  */


#include "gpkspec.h"
#include "randr48.h"      
#include "ctimer.h"      
#include "fftwserver.h"
#include "randfmt.h"

#include "corfunc.h"

#define DeuxPI 2.*M_PI 

//---------------------------------------------------------------
// -- GFour3DPk class : power spectrum computation from a 3D-grid 
//---------------------------------------------------------------
// Constructeur a partir du tableau d rho/rho ou n-gal 
GFour3DPk::GFour3DPk(TArray< TF > & mgrid)
  : mgal_grid_(mgrid), rgp_(RandomGeneratorInterface::GetGlobalRandGenP())
{
  SetPrtLevel();
  // Get the d rho/rho , n_gal array cell size
  double dx = mgal_grid_.Info().GetD("DX",1.);
  double dy = mgal_grid_.Info().GetD("DY",1.);
  double dz = mgal_grid_.Info().GetD("DZ",1.);
  cout<<"GFour3DPk() - calling SetGridCellSize() with cell size from the input grid Info() object ..."<<endl;
  SetGridCellSize(dx, dy, dz, true);
  hp_pk_p_=NULL;  hmcnt_p_=NULL;  hmcntok_p_=NULL;
  hpk2_p_=h2mcnt_p_=NULL;
}

// Destructor
GFour3DPk::~GFour3DPk()
{
  if (hp_pk_p_) delete hp_pk_p_;
  if (hmcnt_p_) delete hmcnt_p_;
  if (hmcntok_p_) delete hmcntok_p_;
}

void GFour3DPk::SetGridCellSize(double dx, double dy, double dz, bool fgprt)
{
  dx_=dx; dy_=dy;  dz_=dz;
  dkx_ = DeuxPI/(mgal_grid_.SizeX()*dx);
  dky_ = DeuxPI/(mgal_grid_.SizeY()*dy);
  dkz_ = DeuxPI/(mgal_grid_.SizeZ()*dz);
  if ((prtlev_>0)||(fgprt))  {
    cout<<"GFour3DPk::SetGridCellSize() dxyz="<<dx_<<"x"<<dy_<<"x"<<dz_<<" --> dkxyz="
	<<dkx_<<"x"<<dky_<<"x"<<dkz_<<endl;
  }
  return;
}

// Generate mass field Fourier Coefficient
void GFour3DPk::doFFT()
{
  FFTWServer ffts;
  // ATTENTION : il faut qu'on soit normalize a true pour avoir la bonne normalisation (c'est le mode par defaut)
  ffts.setNormalize(true);
  ffts.FFTForward(mgal_grid_, fourAmp);
  if (prtlev_>1)
    cout << " GFour3DPk::doFFT done ..." << endl;
  return;
}

// Generate mass field from Fourier Coefficient
void GFour3DPk::doInverseFFT(bool fgNOk0)
{
  if (fgNOk0) fourAmp(0,0,0)=complex<TF>(0.,0.);
  FFTWServer ffts;
  // ATTENTION : il faut qu'on soit normalize a true pour avoir la bonne normalisation (c'est le mode par defaut)
  ffts.setNormalize(true);
  ffts.FFTBackward(fourAmp, mgal_grid_, true);
  if (prtlev_>1)
    cout << " GFour3DPk::doInverseFFT done ..." << endl;
  return;
}

// Generate mass field Fourier Coefficient
void GFour3DPk::generateFourierAmp(ClassFunc1D & pk, double sigma_z, double pkfac)
{
  if (! fourAmp.IsAllocated())  doFFT();
  RandomGeneratorInterface& rg(*rgp_);
  bool fgsmg=false;  // True, gaussian smoothing 
  double dsig2=0.;
  double Asmooth=1.;
  if (sigma_z>1.e-19) {
    dsig2 = 0.5*sigma_z*sigma_z;
    Asmooth=1./sqrt(2.*M_PI);
    fgsmg=true;
  }
  // fourAmp represent 3-D fourier transform of a real input array. 
  // The second half of the array along Y and Z contain negative frequencies
  double kxx, kyy, kzz;
  // sa_size_t is large integer type
  cout << " Four3DPk::generateFourierAmp/Info : generating fourier coefficients - sigma_z="<<sigma_z
       <<" pkfac="<<pkfac<<endl;
  ProgressBar pgb(fourAmp.SizeZ());
  if (prtlev_ > 0)  pgb.setMode(ProgBarM_Percent);
  else pgb.setMode(ProgBarM_None);
  
  // To compute the conversion of P(k) for continuous P(k) for each spatial frequency bin spanning 1/V = 1/(dk_x dk_y dk_z) 
  double Vol = (dx_*mgal_grid_.SizeX())*(dy_*mgal_grid_.SizeY())*(dz_*mgal_grid_.SizeZ());
  double invVol=pkfac/Vol;
  // We ignore 0th term in all frequency directions ...
  for(sa_size_t kz=0; kz<fourAmp.SizeZ(); kz++) {
    kzz =  (kz > fourAmp.SizeZ()/2) ? (double)(fourAmp.SizeZ()-kz)*dkz_ : (double)kz*dkz_; 
    for(sa_size_t ky=0; ky<fourAmp.SizeY(); ky++) {
      kyy =  (ky > fourAmp.SizeY()/2) ? (double)(fourAmp.SizeY()-ky)*dky_ : (double)ky*dky_; 
      for(sa_size_t kx=0; kx<fourAmp.SizeX(); kx++) {  // ignore the 0th coefficient (constant term)
	kxx=(double)kx*dkx_;
	complex<TF> za = fourAmp(kx, ky, kz);
	//	if (za.real()>8.e19) continue;
	double wk = sqrt(kxx*kxx+kyy*kyy+kzz*kzz);
        double amp = sqrt(pk(wk)*invVol/2.);
	if (fgsmg)  amp *= (Asmooth*exp(-(kzz*kzz*dsig2)));
	fourAmp(kx,ky,kz) = complex<TF>(rg.Gaussian(amp), rg.Gaussian(amp));
      }
    }
    pgb.update(kz);
  }
  return;
}

// Compute power spectrum as a function of wave number k 
// cells with amp^2=re^2+im^2>s2cut are ignored
// Output : power spectrum (profile histogram)
HProf GFour3DPk::ComputePk(int nbin, double kmin, double kmax, bool fgmodcnt)
{
  // The second half of the array along Y (matrix rows) contain
  // negative frequencies
  //  int nbh = sqrt(fourAmp.SizeX()*fourAmp.SizeX()+fourAmp.SizeY()*fourAmp.SizeY()/4.+fourAmp.SizeZ()*fourAmp.SizeY()/4.);
  // The profile histogram will contain the mean value of FFT amplitude
  // as a function of wave-number k = sqrt((double)(kx*kx+ky*ky))
  if ((kmax<0.)||(kmax<kmin)) {
    kmin=0.;
    double maxx=fourAmp.SizeX()*dkx_;
    double maxy=fourAmp.SizeY()/2*dky_;
    double maxz=fourAmp.SizeZ()/2*dkz_;
    kmax=sqrt(maxx*maxx+maxy*maxy+maxz*maxz);
  }
  if (nbin<2) nbin=100;
  if (hp_pk_p_) delete hp_pk_p_;
  hp_pk_p_ = new HProf(kmin, kmax, nbin);
  hp_pk_p_->SetErrOpt(false);
  if (fgmodcnt) {
    if (hmcnt_p_) delete hmcnt_p_;
    hmcnt_p_ = new Histo(kmin, kmax, nbin);
    if (hmcntok_p_) delete hmcntok_p_;
    hmcntok_p_ = new Histo(kmin, kmax, nbin);
  }
  if (! fourAmp.IsAllocated())  doFFT();
  if (prtlev_>1)  
    cout << " GFour3DPk::ComputePk()  calling HisPkCumul() kmin,kmax,nbin="
	 <<kmin<<","<<kmax<<","<<nbin<<" ..."<<endl;
  HisPkCumul();

  return *hp_pk_p_;
}

Histo2D GFour3DPk::ComputePk2D(int nbin, double kmax)
{
  if (kmax<1.e-19)  {
    double maxx=fourAmp.SizeX()*dkx_;
    double maxy=fourAmp.SizeY()/2*dky_;
    double maxz=fourAmp.SizeZ()/2*dkz_;
    kmax=sqrt((maxx*maxx+maxy*maxy+maxz*maxz)*0.5);
  }
  if (hpk2_p_) delete hpk2_p_;
  hpk2_p_ = new Histo2D(-kmax,kmax,nbin,-kmax,kmax,nbin);
  if (h2mcnt_p_)  delete h2mcnt_p_;
  h2mcnt_p_ = new Histo2D(-kmax,kmax,nbin,-kmax,kmax,nbin);
  
  uint_8 nmodeok=0;
  // fourAmp represent 3-D fourier transform of a real input array. 
  // The second half of the array along Y and Z contain negative frequencies
  double kxx, kyy, kzz;
  // sa_size_t is large integer type  
  // We ignore 0th term in all frequency directions ...
  for(sa_size_t kz=0; kz<fourAmp.SizeZ(); kz++) {
    kzz =  (kz > fourAmp.SizeZ()/2) ? -(double)(fourAmp.SizeZ()-kz)*dkz_ : (double)kz*dkz_; 
    for(sa_size_t ky=0; ky<fourAmp.SizeY(); ky++) {
      kyy =  (ky > fourAmp.SizeY()/2) ? -(double)(fourAmp.SizeY()-ky)*dky_ : (double)ky*dky_; 
      for(sa_size_t kx=0; kx<fourAmp.SizeX(); kx++) {  // ignore the 0th coefficient (constant term)
	kxx=(double)kx*dkx_;
	complex<TF> za = fourAmp(kx, ky, kz);
	//	if (za.real()>8.e19) continue;
	double ktrans = sqrt(kxx*kxx+kyy*kyy);
	double klong = kzz;
	double amp2 = norm(za); 
	if (h2mcnt_p_) {
	  h2mcnt_p_->Add(ktrans, klong);
	  h2mcnt_p_->Add(-ktrans, klong);  // pour les kx negatifs qu'on n'a pas F(-kx) = conj(F(kx)) 
	}
	hpk2_p_->Add(ktrans, klong, amp2);
	hpk2_p_->Add(-ktrans, klong, amp2);
	nmodeok++;
      }
    }
  }
  (*hpk2_p_) /= (*h2mcnt_p_);
  //  if ((prtlev_>1)||((prtlev_>0)&&(s2cut_>1.e-9))) {
  if (prtlev_>0) {
    cout << " Four3DPk::ComputePk2D/Info : NModeOK=" << nmodeok << " / NMode=" << fourAmp.Size() 
	 << " -> " << 100.*(double)nmodeok/(double)fourAmp.Size() << "%" << endl;
  }
  return (*hpk2_p_);
}

// Compute power spectrum as a function of wave number k 
// Cumul dans hp - cells with amp^2=re^2+im^2>s2cut are ignored
void GFour3DPk::HisPkCumul()
{
  /*  Normalisation fait lors du remplissage du DataTable 
  // On normalise le spectre avec la taille de la grille (Volume V) - L'unite de P(k) est en (Mpc^3)
  double V=Nx()*getdX()*Ny()*getdY()*Nz()*getdZ();
  if (prtlev_>0)  cout<<cout << " GFour3DPk::HisPkCumul() - Normalising power spectrum -> x Volume V="<<V<<" Mpc^3"<<endl;
  */  
  uint_8 nmodeok=0;
  // fourAmp represent 3-D fourier transform of a real input array. 
  // The second half of the array along Y and Z contain negative frequencies
  double kxx, kyy, kzz;
  // sa_size_t is large integer type  
  // We ignore 0th term in all frequency directions ...
  for(sa_size_t kz=0; kz<fourAmp.SizeZ(); kz++) {
    kzz =  (kz > fourAmp.SizeZ()/2) ? (double)(fourAmp.SizeZ()-kz)*dkz_ : (double)kz*dkz_; 
    for(sa_size_t ky=0; ky<fourAmp.SizeY(); ky++) {
      kyy =  (ky > fourAmp.SizeY()/2) ? (double)(fourAmp.SizeY()-ky)*dky_ : (double)ky*dky_; 
      for(sa_size_t kx=0; kx<fourAmp.SizeX(); kx++) {  // ignore the 0th coefficient (constant term)
	kxx=(double)kx*dkx_;
	complex<TF> za = fourAmp(kx, ky, kz);
	//	if (za.real()>8.e19) continue;
	double wk = sqrt(kxx*kxx+kyy*kyy+kzz*kzz);
	double amp2 = norm(za); // *V;   (za.real()*za.real()+za.imag()*za.imag());
	if (hmcnt_p_) hmcnt_p_->Add(wk);
	//	if ((s2cut_>1.e-9)&&(amp2>s2cut_))  continue;
	if (hmcntok_p_) hmcntok_p_->Add(wk);
	hp_pk_p_->Add(wk, amp2);
	nmodeok++;
      }
    }
  }
  
  //  if ((prtlev_>1)||((prtlev_>0)&&(s2cut_>1.e-9))) {
  if (prtlev_>0) {
    cout << " Four3DPk::HisPkCumul/Info : NModeOK=" << nmodeok << " / NMode=" << fourAmp.Size() 
	 << " -> " << 100.*(double)nmodeok/(double)fourAmp.Size() << "%"
	 << " hpk.NEntries()="<<hp_pk_p_->NEntries() << endl;
  }
  return;
}

size_t GFour3DPk::CleanNegatives(TF seuil)
{
  size_t nneg = 0.;
  for(sa_size_t kz=0; kz<mgal_grid_.SizeZ(); kz++) 
    for(sa_size_t ky=0; ky<mgal_grid_.SizeY(); ky++) 
      for(sa_size_t kx=0; kx<mgal_grid_.SizeX(); kx++) 
        if (mgal_grid_(kx, ky, kz) < seuil)  {
          nneg++; mgal_grid_(kx, ky, kz)=seuil;
        }
  cout << " GFour3DPk::CleanNegatives " << nneg << " cells <" << seuil << " changed to" << seuil << endl;
  return nneg;
}

inline double sm_cln_neg(double x) {
  double x2=x*x;
  return exp(-3.5 * x2*x2);
}

size_t GFour3DPk::SmoothCleanNegatives()
{
  double seuil = -0.9;
  double val=0.;
  size_t nneg = 0.;
  for(sa_size_t kz=0; kz<mgal_grid_.SizeZ(); kz++) 
    for(sa_size_t ky=0; ky<mgal_grid_.SizeY(); ky++) 
      for(sa_size_t kx=0; kx<mgal_grid_.SizeX(); kx++) 
        if ((val=(double)mgal_grid_(kx, ky, kz)) < seuil)  {
          nneg++; mgal_grid_(kx, ky, kz)=(TF)sm_cln_neg(val);
        }
  cout << " GFour3DPk::SmoothCleanNegatives " << nneg << " cells <" << seuil << " changed ->"
       << 100*nneg/mgal_grid_.Size() << " % of total cell count"<<endl;
  return nneg;
}

void GFour3DPk::ConvertToGalaxyDensity(double galdens, bool fgpoiss, bool fgrenorm)
{
  RandomGeneratorInterface& rg(*rgp_);
  size_t nneg = 0.;
  ProgressBar pgb(mgal_grid_.SizeZ());
  if (prtlev_ > 0)  pgb.setMode(ProgBarM_Percent);
  else pgb.setMode(ProgBarM_None);
  for(sa_size_t kz=0; kz<mgal_grid_.SizeZ(); kz++) {
    // double fac=0.1+0.9*(1.-(double)kz/(double)mgal_grid_.SizeZ());
    //    double fac=1.;
    for(sa_size_t ky=0; ky<mgal_grid_.SizeY(); ky++) {
      for(sa_size_t kx=0; kx<mgal_grid_.SizeX(); kx++)  {
	// if (fgpoiss)  mgal_grid_(kx, ky, kz) = (TF)(rg.Poisson(fac*galdens*(mgal_grid_(kx,ky,kz)+1.))/fac);
	if (fgpoiss)  mgal_grid_(kx, ky, kz) = (TF)(rg.Poisson(galdens*(mgal_grid_(kx,ky,kz)+1.)));
	else mgal_grid_(kx, ky, kz) = (TF)((int)(galdens*(mgal_grid_(kx,ky,kz)+1.)+0.5));
      }
    }
    pgb.update(kz);
  }
  double mean, sigma;
  MeanSigma(mgal_grid_, mean, sigma);
  cout << " GFour3DPk::ConvertToGalaxyDensity() before renormalizing- mean="<<mean<< " sigma="<<sigma<<endl;
  if (fgrenorm) {
    mgal_grid_ -= (TF)galdens;  mgal_grid_ *= (1./(TF)galdens);
    MeanSigma(mgal_grid_, mean, sigma);
    cout << "  ... AFTER renormalizing- mean="<<mean<< " sigma="<<sigma<<endl;
  }  
  return;
}



// Fills a data table from the computed P(k) profile histogram and mode count 
Histo GFour3DPk::FillPkDataTable(DataTable& dt, double rfac)
{
  if (hp_pk_p_==NULL) throw ParmError("Four3DPk::FillPkDataTable P(k) NOT computed");
  const char* nomcol[5] = {"k","Pk","nmode","nmodok","fracmodok"};
  dt.Clear();
  dt.AddDoubleColumn(nomcol[0]);    
  dt.AddDoubleColumn(nomcol[1]);    
  double Vol = (dx_*mgal_grid_.SizeX())*(dy_*mgal_grid_.SizeY())*(dz_*mgal_grid_.SizeZ());
  if (prtlev_>0) 
    cout << " GFour3DPk::FillPkDataTable()/Info: normalizing P(k) by V=dx*dy*dz="<<Vol<<" Mpc^3 x fac="<<rfac<<endl;
  Vol *= rfac;

  bool fgokmodcnt=true;
  if ((hmcnt_p_==NULL)||(hmcntok_p_==NULL)) {
    cout << " GFour3DPk::FillPkDataTable()/Warning Mode count histos NOT filled, using only P(k) ..."<<endl; 
    fgokmodcnt=false;
    HProf& hp=(*hp_pk_p_);
    DataTableRow 	dtr = dt.EmptyRow();
    for(int_4 ib=0; ib<hp.NBins(); ib++) {
      dtr[0]=hp.BinCenter(ib);
      dtr[1]=hp(ib)*Vol;
      dt.AddRow(dtr);
    }
    Histo fracmodok(hp.XMin(), hp.XMax(), hp.NBins());
    return fracmodok;
  } 
  HProf& hp=(*hp_pk_p_);
  Histo& hmcnt=(*hmcnt_p_);
  Histo& hmcntok=(*hmcntok_p_);
  Histo fracmodok=hmcntok/hmcnt;
  dt.AddIntegerColumn(nomcol[2]);    
  dt.AddIntegerColumn(nomcol[3]);    
  dt.AddFloatColumn(nomcol[4]);    
  DataTableRow 	dtr = dt.EmptyRow();
  for(int_4 ib=0; ib<hp.NBins(); ib++) {
    dtr[0]=hp.BinCenter(ib);
    dtr[1]=hp(ib)*Vol;
    dtr[2]=hmcnt(ib);
    dtr[3]=hmcntok(ib);
    dtr[4]=fracmodok(ib);
    dt.AddRow(dtr);
  }
return fracmodok;
}


//---------------------------------------------------------------
// -- ReprojGrid class: Reprojecting an input grid to an output grid  
//---------------------------------------------------------------

ReprojGrid::ReprojGrid(GridDef& ing, GridDef& outg, TArray< TF > & ingrid)
  : ing_(ing), outg_(outg), in_grid(ingrid), out_grid(outg_.Nx, outg_.Ny, outg_.Nz),
    in_rot((fabs(ing_.theta0)>1.e-9)?ing_.phi0+Angle::PioTwoCst():0., ing_.theta0, 0.) ,
    out_rot((fabs(outg_.theta0)>1.e-9)?outg_.phi0+Angle::PioTwoCst():0., outg_.theta0, 0.)
{
  if ( (InNx() != ing_.Nx) || (InNy() != ing_.Ny) || (InNz() != ing_.Nz) )
    throw ParmError("ReprojGrid::ReprojGrid()/ERROR Incompatible input grid definition and array size !");
  icx_= (double)InNx()*0.5;
  icy_= (double)InNy()*0.5;
  icz_= (double)InNz()*0.5;
  ocx_= (double)OutNx()*0.5;
  ocy_= (double)OutNy()*0.5;
  ocz_= (double)OutNz()*0.5;
}

ReprojGrid::~ReprojGrid()
{
}

void ReprojGrid::Project()
{
  //  if ((fabs(ing_.theta0-outg_.theta0)<1.e-9) && (fabs(ing_.phi0-outg_.phi0)<1.e-9)) return ProjectSameCenter();
  TArray<uint_2> cnt(OutNx(), OutNy(), OutNz());

  cout << "---- ReprojGrid::Project() :"<<endl;
  cout <<"Input grid Size="<<in_grid.Size()<<" ("<<InNx()<<"x"<<InNy()<<"x"<<InNz()<<") "
       <<" r_center="<<ing_.r_center<<" Theta0="<<ing_.theta0*180./M_PI<< " Phi0="<<ing_.phi0*180./M_PI<<" (deg)"<<endl;
  cout<<"Output grid Size="<<out_grid.Size()<<" ("<<OutNx()<<"x"<<OutNy()<<"x"<<OutNz()<<") "
       <<" r_center="<<outg_.r_center<<" Theta0="<<outg_.theta0*180./M_PI<< " Phi0="<<outg_.phi0*180./M_PI<<" (deg)"<<endl;
  if (prtlev_>1) {
    cout << " ingrid_rotation matrix:" << in_rot <<endl;
    cout << " outgrid_rotation matrix:" << out_rot <<endl;
  }
  
  size_t okcount=0;

  ProgressBar pgb(in_grid.SizeZ());
  if (prtlev_ > 0)  pgb.setMode(ProgBarM_Percent);
  else pgb.setMode(ProgBarM_None);

  out_grid=(TF)0.;
  for(sa_size_t kz=0; kz<in_grid.SizeZ(); kz++) {
    double z=((double)kz+0.5-icz_)*ing_.dz;
    for(sa_size_t ky=0; ky<in_grid.SizeY(); ky++) {
      double y=((double)ky+0.5-icy_)*ing_.dy;
      for(sa_size_t kx=0; kx<in_grid.SizeX(); kx++) {
	double x=((double)kx+0.5-icx_)*ing_.dx;
	sa_size_t okx, oky, okz;
	bool fgok=In_XYZ_2OutIJK(x,y,z,okx,oky,okz);
       	if (fgok) {   // inside the output grid 
	  cnt(okx,oky,okz)++;  out_grid(okx,oky,okz)+=in_grid(kx, ky, kz);
	  okcount++;
	}
      }
    }
    pgb.update(kz);
  }
  cout << "ReprojGrid::Project(),End of Loop over input: "<<okcount<<" cells within output grid ( "
       <<(double)okcount*100./(double)in_grid.Size()<<" %)"<<endl;
  size_t ookcnt=0;
  for(sa_size_t kz=0; kz<out_grid.SizeZ(); kz++) 
    for(sa_size_t ky=0; ky<out_grid.SizeY(); ky++) 
      for(sa_size_t kx=0; kx<out_grid.SizeX(); kx++) {
	if (cnt(kx,ky,kz)>0)  { out_grid(kx,ky,kz) /= (TF)cnt(kx,ky,kz);  ookcnt++; }
      }
  cout << "ReprojGrid::Project(), reprojection count renormalisation-NonZero count="<<ookcnt
       <<  " cells / total"<<out_grid.Size()<<" Filled="<<(100*ookcnt)/out_grid.Size()<<" %"<<endl;
  return;
}

void ReprojGrid::ProjectSameCenter()
{
  TArray<uint_2> cnt(OutNx(), OutNy(), OutNz());
  double icx, icy, icz;
  icx=(double)in_grid.SizeX()*0.5-0.5;
  icy=(double)in_grid.SizeY()*0.5-0.5;
  icz=(double)in_grid.SizeZ()*0.5-0.5;
  double ocx, ocy, ocz;
  ocx=(double)out_grid.SizeX()*0.5-0.5;
  ocy=(double)out_grid.SizeY()*0.5-0.5;
  ocz=(double)out_grid.SizeZ()*0.5-0.5;
  
  cout << "ReprojGrid::ProjectSameCenter(), Input grid with "<<in_grid.Size()<<" cells distance to center"
       << " InGrid:"<<ing_.r_center<<" OutGrid:"<<outg_.r_center<<" Mpc"<<endl;
  size_t okcount=0;
  ProgressBar pgb(in_grid.SizeZ());
  if (prtlev_ > 0)  pgb.setMode(ProgBarM_Percent);
  else pgb.setMode(ProgBarM_None);

  out_grid=(TF)0.;
  for(sa_size_t kz=0; kz<in_grid.SizeZ(); kz++) {
    sa_size_t okz=(sa_size_t)(ocz+(((double)kz-icz)*ing_.dz+ing_.r_center-outg_.r_center)/outg_.dz);
    if ((okz<0)||(okz>=OutNz()))  {
      pgb.update(kz);
      continue;
    }
    /*DBG
    cout<<"*DEBUG* kz="<<kz<<" okz="<<okz<<" --- icz="<<icz<<" ocz="<<ocz<<" indz="<<in_dz<<" outdz="<<out_dz<<endl;
    cout<<"*DBG*  icx="<<icx<<" ocx="<<ocx<<" indx="<<in_dx<<" outdx="<<out_dx
	<<" icy="<<icy<<" ocy="<<ocy<<" indy="<<in_dy<<" outdy="<<out_dy<<endl;
    */
    for(sa_size_t ky=0; ky<in_grid.SizeY(); ky++) {
      sa_size_t oky=(sa_size_t)(ocy+((double)ky-icy)*ing_.dy/outg_.dy);
      if ((oky<0)||(oky>=OutNy()))  continue;   // out of the output grid 
      for(sa_size_t kx=0; kx<in_grid.SizeX(); kx++) {
	sa_size_t okx=(sa_size_t)(ocx+((double)kx-icx)*ing_.dx/outg_.dx);
	if ((okx<0)||(okx>=OutNx()))  continue;   // out of the output grid 
	cnt(okx,oky,okz)++;  out_grid(okx,oky,okz)+=in_grid(kx, ky, kz);
	okcount++;
      }
    }
    pgb.update(kz);
  }
  cout << "ReprojGrid::ProjectSameCenter(),End of Loop over input: "<<okcount<<" cells within output grid ( "
       <<(double)okcount*100./(double)in_grid.Size()<<" %)"<<endl;
  size_t ookcnt=0;
  for(sa_size_t kz=0; kz<out_grid.SizeZ(); kz++) 
    for(sa_size_t ky=0; ky<out_grid.SizeY(); ky++) 
      for(sa_size_t kx=0; kx<out_grid.SizeX(); kx++) {
	if (cnt(kx,ky,kz)>0)  { out_grid(kx,ky,kz) /= (TF)cnt(kx,ky,kz);  ookcnt++; }
      }
  cout << "ReprojGrid::ProjectSameCenter(), reprojection count renormalisation-NonZero count="<<ookcnt
       <<  " cells / total"<<out_grid.Size()<<" Filled="<<(100*ookcnt)/out_grid.Size()<<" %"<<endl;
  /*
  POutPersist po("prj.ppf");
  po<<PPFNameTag("cnt")<<cnt;
  po<<PPFNameTag("grid")<<out_grid;
  */
  return;
}

void ReprojGrid::ProjectWithGalaxies(double galdens, double sigmaR, bool fgpoiss, SLinInterp1D* selfuncp)
{
  cout << "ReprojGrid::ProjectWithGalaxies() GalaxyDensity="<<galdens
       <<" /cell ="<<galdens/(ing_.dx*ing_.dy*ing_.dz)<<" /Mpc^3"
       <<" sigmaR="<<sigmaR<<" Mpc"<<endl;       
  cout <<"Input grid Size="<<in_grid.Size()<<" ("<<InNx()<<"x"<<InNy()<<"x"<<InNz()<<") "
       <<" r_center="<<ing_.r_center<<" Theta0="<<ing_.theta0*180./M_PI<< " Phi0="<<ing_.phi0*180./M_PI<<" (deg)"<<endl;
  cout<<"Output grid Size="<<out_grid.Size()<<" ("<<OutNx()<<"x"<<OutNy()<<"x"<<OutNz()<<") "
       <<" r_center="<<outg_.r_center<<" Theta0="<<outg_.theta0*180./M_PI<< " Phi0="<<outg_.phi0*180./M_PI<<" (deg)"<<endl;
  if (prtlev_>1) {
    cout << " ingrid_rotation matrix:" << in_rot <<endl;
    cout << " outgrid_rotation matrix:" << out_rot <<endl;
  }
  if (selfuncp) {
    cout<<"---> Using Selection function, eta(r_center_in)= "<<(*selfuncp)(ing_.r_center)<<endl;
  }
  bool fgsr=false;
  if (sigmaR>1e-3)  fgsr=true;
  
  FMTRandGen rg; 
  rg.SelectGaussianAlgo(C_Gaussian_RandLibSNorm);
  rg.AutoInit(0);
  
  size_t okcount=0, totcount=0;
  ProgressBar pgb(in_grid.SizeZ());
  if (prtlev_ > 0)  pgb.setMode(ProgBarM_Percent);
  else pgb.setMode(ProgBarM_None);

  out_grid=(TF)0.;
  for(sa_size_t kz=0; kz<in_grid.SizeZ(); kz++) {
    double z0=((double)kz-icz_)*ing_.dz;
    for(sa_size_t ky=0; ky<in_grid.SizeY(); ky++) {
      double y0=((double)ky-icy_)*ing_.dy;
      for(sa_size_t kx=0; kx<in_grid.SizeX(); kx++) {
	double x0=((double)kx-icx_)*ing_.dx;
	double ng0=(in_grid(kx, ky, kz)+1.)*galdens;
	if (ng0<0.) continue;
	double etasel=1., weight=1.;
	if (selfuncp) {
	  double za=z0+ing_.r_center;
	  etasel=(*selfuncp)(sqrt(x0*x0+y0*y0+za*za));
	  weight=1./etasel; ng0*=etasel;
	}
	int ng=(fgpoiss ? rg.PoissonSimple(ng0) : (int)(ng0+0.5));
	sa_size_t okx, oky, okz;
	double x,y,z;
	bool fgok=false;
	for(int i=0; i<ng; i++) {
	  x=x0+(rg.Flat01())*ing_.dx;
	  y=y0+(rg.Flat01())*ing_.dy;
	  z=z0+(rg.Flat01())*ing_.dz;
	  /*   DEBUG 
	  sa_size_t okxb, okyb, okzb;
	  bool fgokb=In_XYZ_2OutIJK_NoRot(x,y,z,okxb,okyb,okzb);
	  fgok=In_XYZ_2OutIJK(x,y,z,okx,oky,okz);
	  if ((fgokb!=fgok)||(okxb!=okx)||(okyb!=oky)||(okyb!=oky)) {
	    cout<<"*BUG* x,y,z="<<x<<","<<y<<","<<z<<"  okx,oky,okz="<<okx<<","<<oky<<","<<okz
		<<" b="<<okxb<<","<<okyb<<","<<okzb<<endl;
	  }
	  fgok=In_XYZ_2OutIJK_NoRot(x,y,z,okx,oky,okz);
	  -- END DEBUG */
	  if (fgsr)  fgok=In_XYZ_2OutIJK_SR(rg, sigmaR, x,y,z,okx,oky,okz);
	  else fgok=In_XYZ_2OutIJK(x,y,z,okx,oky,okz);
	  if (fgok) {   // inside the output grid 
	    out_grid(okx,oky,okz)+=weight;
	    okcount++;
	  }
	  totcount++;
	}
      }
    }
    pgb.update(kz);
  }
  cout << "ReprojGrid::Project(),End of Loop over input: "<<okcount<<" gals within output grid "
       <<(double)okcount*100./((double)totcount)<<" %)"<<endl;
  double mean, sigma;
  MeanSigma(out_grid, mean, sigma);
  double galdensout=galdens*(outg_.dx*outg_.dy*outg_.dz)/(ing_.dx*ing_.dy*ing_.dz);
  if (selfuncp) {
    cout << "ReprojGrid::Project(), output grid after reprojection mean="<<mean<<" sigma="<<sigma
	 <<"  (SelFunction-> renormalized by mean)"<<endl; 
    out_grid -= (TF)mean;
    out_grid /= (TF)mean;
  }
  else {
    cout << "ReprojGrid::Project(), output grid after reprojection mean="<<mean<<" sigma="<<sigma
	 <<"   renormalized by galdensout="<<galdensout<<endl; 
    out_grid -= (TF)galdensout;
    out_grid /= (TF)galdensout;
  }
  return;
}


void ReprojGrid::ProjectSpherical(bool fgisoangle)
{
  TArray<uint_2> cnt(OutNx(), OutNy(), OutNz());
  double icx, icy, icz;
  icx=(double)in_grid.SizeX()*0.5-0.5;
  icy=(double)in_grid.SizeY()*0.5-0.5;
  icz=(double)in_grid.SizeZ()*0.5-0.5;
  double ocx, ocy, ocz;
  ocx=(double)out_grid.SizeX()*0.5-0.5;
  ocy=(double)out_grid.SizeY()*0.5-0.5;
  ocz=(double)out_grid.SizeZ()*0.5-0.5;

  if (fgisoangle) cout << "ReprojGrid::ProjectSpherical() - projecting in cells with UNIFORM SOLID ANGLE ..."<<endl;
  else cout << "ReprojGrid::ProjectSpherical() - projecting in cells with UNIFORM TRANSVERSE SIZE (in Mpc) ..."<<endl;
  cout << "ReprojGrid::ProjectSpherical(), Input grid with "<<in_grid.Size()<<" cells distance to center"
       << " InGrid:"<<ing_.r_center<<" OutGrid:"<<outg_.r_center<<" Mpc"<<endl;
  /*
  cout<<"*DEBUG* --- icz="<<icz<<" ocz="<<ocz<<" indz="<<in_dz<<" outdz="<<out_dz<<endl;
  cout<<"*DBG*  icx="<<icx<<" ocx="<<ocx<<" indx="<<in_dx<<" outdx="<<out_dx
	<<" icy="<<icy<<" ocy="<<ocy<<" indy="<<in_dy<<" outdy="<<out_dy<<endl;
  */
  size_t okcount=0;
  ProgressBar pgb(in_grid.SizeZ());
  if (prtlev_ > 0)  pgb.setMode(ProgBarM_Percent);
  else pgb.setMode(ProgBarM_None);

  out_grid=(TF)0.;
  for(sa_size_t kz=0; kz<in_grid.SizeZ(); kz++) {
    double z=((double)kz-icz)*ing_.dz+ing_.r_center;
    for(sa_size_t ky=0; ky<in_grid.SizeY(); ky++) {
      double y=((double)ky-icy)*ing_.dy;
      //DBG      sa_size_t kyc=(sa_size_t)((double)ky-icy);
      for(sa_size_t kx=0; kx<in_grid.SizeX(); kx++) {
	//DBG	sa_size_t kxc=(sa_size_t)((double)kx-icx);
	double x=((double)kx-icx)*ing_.dx;
	double r=sqrt(x*x+y*y+z*z);
	sa_size_t okz=(sa_size_t)(ocz+(r-outg_.r_center)/outg_.dz);
	if ((okz<0)||(okz>=OutNz()))   continue;   // out of the output grid
	double phi=atan2(y,x);
	double theta=acos(z/r);
	double rt=((fgisoangle)?outg_.r_center*theta:r*theta);
	sa_size_t okx=(sa_size_t)(ocx+rt*cos(phi)/outg_.dx);
	/*
	if ((kxc==0)||(kyc==0)) {
	  sa_size_t oky2=(sa_size_t)(ocy+rt*sin(phi)/out_dy);
	  cout<<"*DBG*B* kx,y,z="<<kx<<","<<ky<<","<<kz<<" -> r="<<r<<" okz="<<okz<<" theta="<<theta<<endl
	      <<" phi="<<phi<<" rt="<<rt<<" okx="<<okx<<" oky="<<oky2<<endl;
	}
	*/
	if ((okx<0)||(okx>=OutNx()))  continue;   // out of the output grid 
	sa_size_t oky=(sa_size_t)(ocy+rt*sin(phi)/outg_.dy);
	if ((oky<0)||(oky>=OutNy()))  continue;   // out of the output grid 
	cnt(okx,oky,okz)++;  out_grid(okx,oky,okz)+=in_grid(kx, ky, kz);
	okcount++;
      }
    }
    pgb.update(kz);
  }
  cout << "ReprojGrid::ProjectSpherical(),End of Loop over input: "<<okcount<<" cells within output grid ( "
       <<(double)okcount*100./(double)in_grid.Size()<<" %)"<<endl;
  size_t ookcnt=0;
  for(sa_size_t kz=0; kz<out_grid.SizeZ(); kz++) 
    for(sa_size_t ky=0; ky<out_grid.SizeY(); ky++) 
      for(sa_size_t kx=0; kx<out_grid.SizeX(); kx++) {
	if (cnt(kx,ky,kx)>0)  { out_grid(kx,ky,kz) /= (TF)cnt(kx,ky,kx);  ookcnt++; }
      }
  cout << "ReprojGrid::ProjectSpherical(), reprojection count renormalisation-NonZero count="<<ookcnt
       <<  " cells / total"<<out_grid.Size()<<" Filled="<<(100*ookcnt)/out_grid.Size()<<" %"<<endl;

  /*
  POutPersist po("prjsph.ppf");
  po<<PPFNameTag("cntsph")<<cnt;
  po<<PPFNameTag("gridsph")<<out_grid;
  */
  return;
  
}

void ReprojGrid::ProjectSphericalWithGalaxies(double galdens, bool fgisoangle, bool fgpoiss)
{
}


void ReprojGrid::Grid2RThetaPhiNTuple(string& ppfname)
{
  cout << "----- ReprojGrid::Grid2RThetaPhiNTuple()------ \n"
       <<"Input grid Size="<<in_grid.Size()<<" ("<<InNx()<<"x"<<InNy()<<"x"<<InNz()<<") "
       <<" r_center="<<ing_.r_center<<" Theta0="<<ing_.theta0*180./M_PI<< " Phi0="<<ing_.phi0*180./M_PI<<" (deg)"<<endl;
  cout<<"Output grid Size="<<out_grid.Size()<<" ("<<OutNx()<<"x"<<OutNy()<<"x"<<OutNz()<<") "
       <<" r_center="<<outg_.r_center<<" Theta0="<<outg_.theta0*180./M_PI<< " Phi0="<<outg_.phi0*180./M_PI<<" (deg)"<<endl;
  cout << " in_rot:" << in_rot <<endl;
  cout << " out_rot:" << out_rot <<endl;

  vector<Vector3d> vv;
  vv.push_back( Vector3d(1.,0.,0.) );
  vv.push_back( Vector3d(0.,1.,0.) );
  vv.push_back( Vector3d(0.,0.,1.) );
  vv.push_back( Vector3d(0.1,0.,sqrt(1.-0.1*0.1) ) );
  vv.push_back( Vector3d(0,0.1,sqrt(1.-0.1*0.1) ) );
  vv.push_back( Vector3d(0.1,0.1,sqrt(1.-0.2*0.1) ) );
  cout << "==========> v2= in_rot.Rotate(v1); "<<endl;
  for(size_t jj=0; jj<vv.size(); jj++) {
    Vector3d v1=vv[jj];
    Vector3d v2= in_rot.Rotate(v1);
    cout << "v1="<<v1<<" (T,P="<<v1.GetThetaPhi()<<") \n      -> v2="<<v2<<" (T,P="<<v2.GetThetaPhi()<<")"<<endl;
  }
  cout << "==========> v2= in_rot.RotateBack(v1); "<<endl;
  for(size_t jj=0; jj<vv.size(); jj++) {
    Vector3d v1=vv[jj];
    Vector3d v2= in_rot.RotateBack(v1);
    cout << "v1="<<v1<<" (T,P="<<v1.GetThetaPhi()<<") \n      -> v2="<<v2<<" (T,P="<<v2.GetThetaPhi()<<")"<<endl;
  }

  cout << " Checking Index computation  In_XYZ_2OutIJK_NoRot() - without rotation ... " << endl;
  // For debugging
  double icx, icy, icz;
  icx=(double)in_grid.SizeX()*0.5-0.5;
  icy=(double)in_grid.SizeY()*0.5-0.5;
  icz=(double)in_grid.SizeZ()*0.5-0.5;
  double ocx, ocy, ocz;
  ocx=(double)out_grid.SizeX()*0.5-0.5;
  ocy=(double)out_grid.SizeY()*0.5-0.5;
  ocz=(double)out_grid.SizeZ()*0.5-0.5;


  cout<<"DBG icx,icy,icz  ="<<icx<<","<<icy<<","<<icz<<" ocx,ocy,ocz  ="<<ocx<<","<<ocy<<","<<ocz<<endl;
  cout<<"DBG icx,icy,icz_  ="<<icx_<<","<<icy_<<","<<icz_<<" ocx,ocy,ocz_ ="<<ocx_<<","<<ocy_<<","<<ocz_<<endl;

  size_t okcount=0;
  size_t okcounta=0;  
  size_t errcount=0;
  sa_size_t okx,oky,okz;
    
  for(sa_size_t kz=0; kz<in_grid.SizeZ(); kz++) {
    double z=((double)kz+0.5-icz_)*ing_.dz;
    for(sa_size_t ky=0; ky<in_grid.SizeY(); ky++) {
      double y=((double)ky+0.5-icy_)*ing_.dy;
      for(sa_size_t kx=0; kx<in_grid.SizeX(); kx++) {
	double x=((double)kx+0.5-icx_)*ing_.dx;
	sa_size_t okx, oky, okz;
	bool fgok=In_XYZ_2OutIJK_NoRot(x,y,z,okx,oky,okz); 
	sa_size_t okza=(sa_size_t)(ocz+(((double)kz-icz)*ing_.dz+ing_.r_center-outg_.r_center)/outg_.dz);  
	sa_size_t okya=(sa_size_t)(ocy+((double)ky-icy)*ing_.dy/outg_.dy);
	sa_size_t okxa=(sa_size_t)(ocx+((double)kx-icx)*ing_.dx/outg_.dx);
	bool fgoka=true;
	if ((okza<0)||(okza>=OutNz())||(okxa<0)||(okxa>=OutNx())||(okya<0)||(okya>=OutNy()))  fgoka=false;
	if ((fgoka!=fgok)||(okxa!=okx)||(okya!=oky)||(okya!=oky)) {
	  cout<<"*BUG* kx,ky,kz="<<kx<<","<<ky<<","<<kz<<"  x,y,z="<<x<<","<<y<<","<<z
	      <<"  okx,oky,okz="<<okx<<","<<oky<<","<<okz
	      <<" a="<<okxa<<","<<okya<<","<<okza<<endl;
	  errcount++;
	}
      }
    }
  }
  if (errcount==0)  cout<< "Check  In_XYZ_2OutIJK_NoRot() OK (errcount=0)";
  else  cout<< "ERROR Checking  In_XYZ_2OutIJK_NoRot() -> errcount="<<errcount;
  cout << " [ okcount="<<okcount<<" okcounta="<<okcounta<<" ]"<<endl;  

  const char* names[16]={"r1","t1","p1","r2","t2","p2","ix","jy","kz","fg","x","y","z","ox","oy","oz"};
  NTuple nt(16,names,2048,false);
  float xnt[20];
  for(sa_size_t kz=0; kz<in_grid.SizeZ(); kz++) {
    double z=((double)kz+0.5-icz_)*ing_.dz;
    for(sa_size_t ky=0; ky<in_grid.SizeY(); ky++) {
      double y=((double)ky+0.5-icy_)*ing_.dy;
      for(sa_size_t kx=0; kx<in_grid.SizeX(); kx++) {
	double x=((double)kx+0.5-icx_)*ing_.dx;
	Vector3d vc;
	Vector3d vr=In_XYZ_RThetaPhi(x,y,z,vc);
	bool fg=In_XYZ_2OutIJK(x,y,z,okx,oky,okz);
	xnt[0]=vc.Norm();
	LongitudeLatitude ll=vc.GetThetaPhi();
	xnt[1]=ll.Theta();
	xnt[2]=ll.Phi();
	xnt[3]=vr.Norm();
	ll=vr.GetThetaPhi();
	xnt[4]=ll.Theta();
	xnt[5]=ll.Phi();
	xnt[6]=okx;  xnt[7]=oky;  xnt[8]=okz; xnt[9]=(fg?1.:0.);
	xnt[10]=x-icx_;  xnt[11]=y-icy_;  xnt[12]=z-icz_;
	xnt[13]=((double)okx+0.5)*outg_.dx-ocx_;
	xnt[14]=((double)oky+0.5)*outg_.dy-ocy_;
	xnt[15]=((double)okz+0.5)*outg_.dz-ocz_;
	nt.Fill(xnt);
      }
    }
  }
  cout<<"ReprojGrid::Grid2RThetaPhiNTuple() Saving NTuple to PPF file "<<ppfname<<endl;
  cout<<nt;
  POutPersist po(ppfname);
  po<<nt;
}

