/*  ------------------------ Projet BAO/PhotoZ/LSST -------------------- 
  Programme de calcul du spectre de puissance (3D) a partir d'un 
  cube de donnees (delta rho/rho ou NGal )
    R. Ansari (pour Adline Choyer) - Feb 2015 
  Usage: grid2pkgrid2pk [options] In3DMap_FitsName OutPk_TextFile [OutPk_PPFFile]
         Options: [-d dx,dy,dz] [-k kmin,kmax,nbin] [-t lowval] [-r rfac] [-p lev]
---------------------------------------------------------------  */

#include "jregrid.h"

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

CosmoCoord::CosmoCoord(SimpleUniverse &su)
  : su_(su), isu_(su_)
{
  size_t npts=1000;
  vector<double> leslos(npts);
  vector<double> lesz(npts);
  for(size_t i=0; i<npts; i++) {
    lesz[i]=i*10./(double)npts;
    leslos[i]=isu_.LineOfSightComovDistanceMpc(lesz[i]);
  }
  iisu_.DefinePoints(leslos, lesz);
}
void CosmoCoord::getAngleRedshift(double dX, double dY, double dZ, double & theta, double & phi, double & los, double & redshift)
{
  los = sqrt(dX*dX+dY*dY+dZ*dZ);
  theta=atan2(sqrt(dX*dX+dY*dY),dZ);
  phi=atan2(dY,dX);
  redshift=getRedshift(los);
}


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------
PowerSpecBase::PowerSpecBase(double zref, double kmin, double kmax)
  :  zref_(zref), kmin_(kmin), kmax_(kmax)
{
}

PowerSpecBase::~PowerSpecBase()
{
}

double PowerSpecBase::Inverse(double pkval) {
  //Return k of power spectrum maximum between kpeak and kmax.
  double kmin = kpeak_;
  double kmax = kmax_;
  double kval = (kmax_+kmin_)*0.5;
  double diff = operator()(kval) - pkval;
  size_t iter=0;
  size_t Maxiter=100;
  double eps = 1e-9;
  while(fabs(diff)>eps && iter<Maxiter){
    iter++;
    if(diff>0){
      kmin = kval;
    } else {
      kmax = kval;
    }
    kval = (kmax+kmin)*0.5;
    diff = operator()(kval) - pkval;
  }
  if(iter>=Maxiter-1) std::cout << "Warning iter large" << std::endl;
  return kval;
}

PowerSpecToy::PowerSpecToy(double mu, double sigma, double norm, double zref, double kmin, double kmax)
  : mu_(mu), sigma_(sigma), norm_(norm)
{
  zref_=zref; kmin_=kmin; kmax_=kmax;
  mu_ = log(mu);
  pmax_ = operator()(exp(mu_-sigma_*sigma_));
  kpeak_=  exp(mu_-sigma_*sigma_);
}


PowerSpecFile::PowerSpecFile(string inpkname, double zref=0., double kmin=0.001, double kmax=100) {
  zref_=zref; kmin_=kmin; kmax_=kmax;
  cout << " Reading input Pk from file " << inpkname << " at zref="<<zref_<< endl;
  TMatrix<r_8> pkin;
  sa_size_t nr, nc;
  ifstream is(inpkname.c_str());
  pkin.ReadASCII (is, nr, nc);
  Vector ks = pkin.Column(0);
  Vector Pks = pkin.Column(1);
  vector<double> vx;// = ks.ConvertTostdvec();
  vector<double> vy;// = Pks.ConvertTostdvec();
  for(size_t i=0;i<nr;i++){
    if(ks(i)<kmin_||ks(i)>kmax_) continue;
    vx.push_back(ks(i));
    vy.push_back(Pks(i));
  }
  Pk_.DefinePoints(vx,vy);
  //kmin_=pkin(0,0);
  //kmax_=pkin(nr-1,0);
  npt_ = nr; lkmin_ = log10(kmin); lkmax_ = log10(kmax); dlk_=(lkmax_-lkmin_)/npt_;
  //cout << "   kmin from file = " << kmin_<<endl; 
  //cout << "   kmax from file = " << kmax_<<endl;
  Pks.MinMax(pmin_,pmax_);
  kpeak_ = Inverse(pmax_);
  //----------- Fill Pk function --------
  const int n = 2;
  const char *vname[n] = { "k","pk" };
  NTuple nt(n,vname); nt_ = nt;
  double xnt[n];     
  for(double lk=lkmin_;lk<lkmax_+dlk_/2.;lk+=dlk_) {
   double k = pow(10.,lk);
   double pkzk = Pk_(k);
   xnt[0] = k;
   xnt[1] = pkzk;
   nt_.Fill(xnt);
  }
}

PowerSpecCalc::PowerSpecCalc(SimpleUniverse su, bool osc=true, double zref=0., double kmin=0.001, double kmax=1, int npt=200) :
  su_(su), osc_(osc) {
  //---------- Initialisation ------------
  zref_=zref; kmin_=kmin; kmax_=kmax;
  scale_ = 2.54499e+07;  // normalisation du spectre a z=0 selon SDSS
  lkmin_ = log10(kmin); lkmax_ = log10(kmax); dlk_=(lkmax_-lkmin_)/npt;
  h100 = su.h();
  Om0 = su.OmegaMatter();
  Ob0 = su.OmegaBaryon();
  Ol0 = su.OmegaLambda();
  w0 = su.wDE();
  ns = su.Ns();
  as =  1.;
  T_CMB = su.T_CMB();
  cout<<"h100="<<h100<<" Om0="<<Om0<<" Ob0="<<Ob0<<" Ocdm="<<Om0-Ob0<<" Ol0="<<Ol0<<" w0="<<w0<<endl;
  cout<<"zref="<<zref<<" lkmin="<<lkmin_<<" lkmax="<<lkmax_<<" npt="<<npt<<" dlk="<<dlk_<<endl;
  //---------- Computation ------------
  InitialPowerLaw pkini(ns,as);
  pkini.SetNorm(scale_);
  TransfertEisenstein tf(h100,Om0-Ob0,Ob0,T_CMB,false);
   double tf0 = tf(0);
   kpeak_ = tf.KPeak();
   cout<<"kpeak="<<kpeak_<<endl;
  GrowthEisenstein d1(Om0,Ol0);
  cout<<"GrowthFactor: D1(0)="<<d1(0.)<<" D1("<<zref<<") = "<<d1(zref)<<endl;
  PkSpecCalc pkz(pkini,tf,d1,zref);

  TransfertEisenstein tfnosc2(tf); tfnosc2.SetNoOscEnv(2);
    double tfnosc20 = tfnosc2(0);
    cout<<"tf(0)="<<tf0<<" tfnosc2(0)="<<tfnosc20<<endl;
  PkSpecCalc pkznosc2(pkini,tfnosc2,d1,zref);
  /*
  TransfertEisenstein tfnosc1(tf); tfnosc1.SetNoOscEnv(1);
   double tfnosc10 = tfnosc1(0);
  TransfertEisenstein tfnob(h100,Om0,0.,T_CMB_Par,true);
   double tfnob0 = tfnob(0);
  PkSpecCalc pkznosc1(pkini,tfnosc1,d1,zval);
  PkSpecCalc pkznob(pkini,tfnob,d1,zval);
   cout<<"tfnosc1(0)="<<tfnosc10<<" tfnob(0)="<<tfnob0<<endl;
  */
  //----------- Fill Pk function --------
  const int n = 2;
  const char *vname[n] = { "k","pk" };
  NTuple nt(n,vname); nt_ = nt;
  double xnt[n];     
  vector<double> ks;
  vector<double> Pks;
  //vector<double> Pknosc2;
  double max=-1;
  for(double lk=lkmin_;lk<lkmax_+dlk_/2.;lk+=dlk_) {
   double k = pow(10.,lk);
   ks.push_back(k);
   double pkzk = 0;
   if(osc) {
     pkzk = pkz(k);
   } else {
     pkzk = pkznosc2(k);
   }
   Pks.push_back(pkzk);
   //Pknosc2.push_back(pkznosc2(k));
   xnt[0] = k;
   xnt[1] = pkzk;
   //xnt[2] = pkznosc2(k);
   if(pkzk>max) max=pkzk;
   nt_.Fill(xnt);
  }
  Pk_.DefinePoints(ks,Pks);
  //Pknosc_.DefinePoints(ks,Pknosc2);
  pmax_ = max;
}


//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

ProjectGridOnShells::ProjectGridOnShells(TArray< TF > & mgrid, CosmoCoord csu, double zcenter=0.5)
  : mgal_grid_(mgrid), csu_(csu) {
  ShellBadVal = -10;
  SetPrtLevel();
  // Get the d rho/rho , n_gal array cell size
  double dx = mgal_grid_.Info().GetD("DX",1.);
  double dy = mgal_grid_.Info().GetD("DY",1.);
  double dz = mgal_grid_.Info().GetD("DZ",1.);
  cout<<"ProjectGridOnShells() - calling SetGridCellSize() with cell size from the input grid Info() object ..."<<endl;
  SetGridCellSize(dx, dy, dz, true);
  // Get the center of the grid
  iX0_ = getCentralIndex(mgal_grid_.Size(0));
  iY0_ = getCentralIndex(mgal_grid_.Size(1));
  iZ0_ = getCentralIndex(mgal_grid_.Size(2));
  // Set the number of thin shells of depth dz
  LOScenter_ = csu_.getLOS(zcenter);
}

double ProjectGridOnShells::getCentralIndex(size_t npts) {
  if(npts % 2 == 0) {
    return(npts/2);
  } else {
    return(npts/2+0.5);
  }
}

void ProjectGridOnShells::SetShells() {
  SetNShells();
  Shells_.resize(nshells_);
  ShellsN_.resize(nshells_);
  ShellsLOS_.resize(nshells_);
  ShellsRedshifts_.resize(nshells_);
  for(size_t i=0; i<nshells_; i++) {
    double ShellLOS = LOScenter_ + (i-iZ0_)*dz_;
    double reso = atan2(min(dx_,dy_),ShellLOS);
    int M = SphereThetaPhi<r_8>::ResolToSizeIndex(reso);
    SphereThetaPhi<r_8> sphtp(M);
    sphtp = ShellBadVal;
    Shells_[i] = sphtp;
    ShellsLOS_[i] = ShellLOS;
    ShellsRedshifts_[i] = csu_.getRedshift(ShellLOS);
    SphereThetaPhi<r_8> sphn(M);
    sphn = 0.;
    ShellsN_[i] = sphn;
  }
}

size_t ProjectGridOnShells::GetShellbyLOS(double los, size_t ishell_start=0) {
  size_t ishell = 0;
  for(size_t i=ishell_start;i<nshells_;i++) {
    if(ShellsLOS_[i]-dz_/2. < los and ShellsLOS_[i]+dz_/2 > los){
      ishell = i;
      break;
    }
  }
  return(ishell);
}

size_t ProjectGridOnShells::GetShellbyRedshift(double redshift, size_t ishell_start=0) {
  size_t ishell = 0;
  double zdiff = ShellsRedshifts_[ishell_start]-redshift;
  for(size_t i=ishell_start+1;i<nshells_;i++) {
    double zdiff_tmp = ShellsRedshifts_[i]-redshift;
    if(zdiff_tmp < zdiff){
      zdiff = zdiff_tmp;
    } else {
      ishell = i-1;
      break;
    }
  }
  return(ishell);
}

void ProjectGridOnShells::FillShells() {
  for(size_t i=0; i<mgal_grid_.SizeX(); i++) {
    for(size_t j=0; j<mgal_grid_.SizeY(); j++) {
      size_t ishell_start = 0;
      for(size_t k=0; k<mgal_grid_.SizeZ(); k++) {
	double dX = (i-iX0_)*dx_;
	double dY = (j-iY0_)*dy_;
	double dZ = (k-iZ0_)*dz_;
	double theta, phi, los, redshift;
	csu_.getAngleRedshift(dX, dY, dZ, theta, phi, los, redshift);
	size_t ishell = GetShellbyLOS(los,ishell_start);
	if (ishell>=Shells_.size())  continue ; // out of range LOS
	Shells_[ishell](theta,phi) += mgal_grid_(i,j,k);
	ShellsN_[ishell](theta,phi) += 1.;
	  /*
	ishell_start = ishell;
	if(Shells_[ishell](theta,phi) == ShellBadVal) {
	  Shells_[ishell](theta,phi) = mgal_grid_(i,j,k);
	} else {
	  Shells_[ishell](theta,phi) += mgal_grid_(i,j,k);
	  ShellsN_[ishell](theta,phi) += 1;
	}
	  */
      }
    }
  }
  for(size_t i=0;i<nshells_;i++){
    Shells_[i] /= ShellsN_[i];  // je pense que lq division n'est pas protegee 
      //    Shells_[i] /= const_cast<const SphereThetaPhi<r_8>& (ShellsN_[i]);
  }
}

//-----------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------

ReProjectGrid::ReProjectGrid(TArray< TF > & mgrid, CosmoCoord csu, double zcenter=0.5)
  : mgal_grid_(mgrid), csu_(csu) {
  SetPrtLevel();
  // Get the d rho/rho , n_gal array cell size
  double dx = mgal_grid_.Info().GetD("DX",1.);
  double dy = mgal_grid_.Info().GetD("DY",1.);
  double dz = mgal_grid_.Info().GetD("DZ",1.);
  cout<<"ReProjectGrid() - calling SetInputGridCellSize() with cell size from the input grid Info() object ..."<<endl;
  SetInputGridCellSize(dx, dy, dz);
  // Get the center of the grid
  iX0_ = getCentralIndex(mgal_grid_.Size(0));
  iY0_ = getCentralIndex(mgal_grid_.Size(1));
  iZ0_ = getCentralIndex(mgal_grid_.Size(2));
  // Set the spherical projection parameters
  LOScenter_ = csu_.getLOS(zcenter);
  zcenter_ = zcenter;
  cout<<"ReProjectGrid() - LOScenter="<< LOScenter_<<endl;
}

void ReProjectGrid::UpdateInputGrid(TArray< TF > & mgrid) {
  mgal_grid_ = mgrid;
  // Get the d rho/rho , n_gal array cell size
  double dx = mgrid.Info().GetD("DX",1.);
  double dy = mgrid.Info().GetD("DY",1.);
  double dz = mgrid.Info().GetD("DZ",1.);
  SetInputGridCellSize(dx, dy, dz);
  // Get the center of the grid
  iX0_ = getCentralIndex(mgrid.Size(0));
  iY0_ = getCentralIndex(mgrid.Size(1));
  iZ0_ = getCentralIndex(mgrid.Size(2));
}

void ReProjectGrid::DefineOutputGrid(sa_size_t nx, sa_size_t ny, sa_size_t nz, double dx, double dy, double dz)
{
  sa_size_t siz[3];  siz[0]=nx; siz[1]=ny; siz[2]=nz;
  outgrid_.SetSize(3,siz);
  outgridN_.SetSize(3,siz);
  if (dx>0.)  odx_=dx;
  if (dy>0.)  ody_=dy;
  if (dz>0.)  odz_=dz;
  ioX0_ = getCentralIndex(outgrid_.Size(0));
  ioY0_ = getCentralIndex(outgrid_.Size(1));
  ioZ0_ = getCentralIndex(outgrid_.Size(2));
}

void ReProjectGrid::FillOutputGrid() {
  outgrid_=0.;
  outgridN_=0;
  size_t okcnt=0;
  for(size_t i=0; i<mgal_grid_.SizeX(); i++) {
    double dX = (i-iX0_)*dx_;
    sa_size_t io=(dX/odx_)+ioX0_;
    if ((io<0)||(io>=outgrid_.SizeX()))  continue;
    for(size_t j=0; j<mgal_grid_.SizeY(); j++) {
      double dY = (j-iY0_)*dy_;
      sa_size_t jo=(dY/ody_)+ioY0_;
      if ((jo<0)||(jo>=outgrid_.SizeY()))  continue;
      for(size_t k=0; k<mgal_grid_.SizeZ(); k++) {
	double dZ = (k-iZ0_)*dz_;
	sa_size_t ko=(dZ/odz_)+ioZ0_;
	if ((ko<0)||(ko>=outgrid_.SizeZ()))  continue;
	outgrid_(io,jo,ko) += mgal_grid_(i,j,k);
	outgridN_(io,jo,ko) += 1.;
	okcnt++;
      }
    }
  }
  
  cout << "ReProjectGrid::FillOutputGrid() okcount="<<okcnt<<" / Total="<<mgal_grid_.Size()<<endl;
  outgrid_.DivElt(outgridN_, outgrid_, false, true);  // division par zero protegee
}

void ReProjectGrid::DefineOutputGridRedshift()
{
  //odtheta_ = atan2(odx_,LOScenter_-odz_*iZ0_); // dtheta_ is set at the closest edge of the reprojected cude from observer
  //odphi_ = atan2(ody_,iX0_*odx_);// dphi_ is set at the closest edge of the reprojected cude from observer
  sa_size_t siz[3];  siz[0]=outgrid_.SizeX(); siz[1]=outgrid_.SizeY();
  siz[2]=outgrid_.SizeZ();
  zmin_ = csu_.getRedshift(LOScenter_-odz_*ioZ0_);
  zmax_ = csu_.getRedshift(LOScenter_+odz_*ioZ0_);
  odzredshift_ = (zmax_-zmin_)/siz[2];
  outgridz_.SetSize(3,siz);
  outgridzN_.SetSize(3,siz);
  iozX0_ = getCentralIndex(outgridz_.Size(0));
  iozY0_ = getCentralIndex(outgridz_.Size(1));
  iozZ0_ = getCentralIndex(outgridz_.Size(2));
}

void ReProjectGrid::FillOutputGridRedshift() {
  outgridz_=0.;
  outgridzN_=0;
  //cout<<mgal_grid_.SizeX()<<" "<<mgal_grid_.SizeY()<<" "<<mgal_grid_.SizeZ()<<" "<<odx_<<" "<<ody_<<" "<<odz_<<" "<<dx_<<" "<<dy_<<" "<<dz_<<" "<<iozX0_<<" "<<iozY0_<<" "<<iozZ0_<<" "<<" "<<odzredshift_<<" "<<zmin_<<endl;
  size_t okcnt=0;
  for(size_t i=0; i<mgal_grid_.SizeX(); i++) {
    double dX = (i-iX0_)*dx_;
    for(size_t j=0; j<mgal_grid_.SizeY(); j++) {
      double dY = (j-iY0_)*dy_;
      for(size_t k=0; k<mgal_grid_.SizeZ(); k++) {
	double dZ = LOScenter_+(k-iZ0_)*dz_;
	double theta, phi, los, redshift;
	csu_.getAngleRedshift(dX, dY, dZ, theta, phi, los, redshift);
	//sa_size_t io=int((theta/odtheta_)+ioX0AZ_);
	sa_size_t io=int(los*theta*cos(phi)/odx_+iozX0_);
	if ((io<0)||(io>=outgridz_.SizeX())) continue;
	//sa_size_t jo=int((phi/odphi_)+ioY0AZ_);
	sa_size_t jo=int(los*theta*sin(phi)/ody_+iozY0_);
	if ((jo<0)||(jo>=outgridz_.SizeY())) continue;
	sa_size_t ko=int((redshift-zmin_)/odzredshift_);
	if ((ko<0)||(ko>=outgridz_.SizeZ()))  continue;
	outgridz_(io,jo,ko) += mgal_grid_(i,j,k);
	outgridzN_(io,jo,ko) += 1.;
	okcnt++;
      }
    }
  }
  
  cout << "ReProjectGrid::FillOutputGrid() okcount="<<okcnt<<" / Total="<<mgal_grid_.Size()<<endl;
  POutPersist  podbg("outgridzN.ppf");  podbg<<outgridzN_;
  outgridz_.DivElt(outgridzN_, outgridz_, false, true);  // division par zero protegee
}

static int myprtcnt=0;

TArray< TF > ReProjectGrid::AddPhotoZ(TArray <TF> ingrid, double phoz_err){
  sa_size_t siz[3];  siz[0]=ingrid.SizeX(); siz[1]=ingrid.SizeY();
  siz[2]=ingrid.SizeZ();
  TArray< TF > outgridphoz;
  outgridphoz.SetSize(3,siz);
  outgridphoz=0.;
  for(size_t i=0; i<ingrid.SizeX(); i++) {
    for(size_t j=0; j<ingrid.SizeY(); j++) {
      for(size_t k=0; k<ingrid.SizeZ(); k++) {
	double z0 = zmin_+k*odzredshift_;  // redshift of central pixel
	double sigmaz = phoz_err*(1+z0);   // photo-z error of this pixel
	double A = odzredshift_/(sigmaz*sqrt(2*Pi)); // normalize probability to unity
	WGauss w(z0,sigmaz,A);
	//WTopHat w(z0-2*odzredshift_,z0+odzredshift_*2,0.25);
	//cout<<z0<<" "<<sigmaz<<" "<<A<<" "<<odzredshift_<<endl;
	//for(size_t kk=0; kk<outgridphoz.SizeZ(); kk++) {
	//if (myprtcnt<5) cout<<"DBG sigmaz="<<sigmaz<< " z0="<<z0<<endl;
	for(size_t kk=max(0,int((z0-6*sigmaz)/odzredshift_)); kk<min(int(outgridphoz.SizeZ()),int((z0+6*sigmaz)/odzredshift_)); kk++) {
	  //for(size_t kk=max(0,int(k-10)); kk<min(int(outgridphoz.SizeZ()),int(k+10)); kk++) {
	  double z = zmin_+kk*odzredshift_;
	  outgridphoz(i,j,kk) += w(z)*ingrid(i,j,k);
	  //if (myprtcnt<5 && k>10) cout <<"DBG["<<myprtcnt<<"]- k="<<k<<" kk="<<kk<<" z0="<<z0<<" z="<<z<<" w(z)="<<w(z)<<endl;
	}
	myprtcnt++;
      }
    }
  }
  return(outgridphoz);
}

