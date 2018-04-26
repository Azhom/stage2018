==========
PZBAOUtils
==========
R. Ansari ,  Univ. Paris-Sud,  LAL, IN2P3/CNRS

small utility programs for BAO/PhotoZ catalog analysis -
Creation du module en Fevrier 2015
MAJ en septembre 2015, Février 2016
Modifications majeures : Novembre/Décembre 2016 - Janvier 2017 
Ajout repertoire LFCode ( tirage magnitude / type de galaxies selon
des fonctions de luminosite en Juillet 2017  

-------------------------------------------------------------------------
LFCode
-------------------------------------------------------------------------
lfparam.txt		lfparamRamosAll.txt	makefile		mlf.cc
multitypzlf.cc		multitypzlf.h
myschechter.cc		myschechter.h

Classe  MultiType_Z_LF , MySchechter , fichier de parametres LF et
programme de test mlf.cc 

-------------------------------------------------------------------------
mergedtpknosc.cc
-------------------------------------------------------------------------
mergedtpknosc.cc : a small utilit program to merge two P(k) DataTables computed by tpkprj.cc
program, for BAO and noBAO - P(k) and P(k)-NoOsc input power spectra. Computes also the
ratios P(k)/P(k)-noOsc and expected errors 

 ---- mergedtpknosc: merging P(k) with and without oscillation datatables ---- 
 Usage: mergedtpknosc DtPkPPFFile DtPkNoOscPPFFile VolSurvey OutPPFFile 


-------------------------------------------------------------------------
terrpk.cc 
-------------------------------------------------------------------------
terrpk.cc : Computes a DataTable with P(k) and P(k)-NoOsc, the ratio and associated errors 
 Usage: terrpk [options] InputPkFile Out_PPFFile 
   options: [-grid Nx,Ny,Nz,dx,dy,dz] [-k kmin,kmax,nbin] [-p density] 
   InputPkFile : text file with Pk as output by SimLSS

-------------------------------------------------------------------------
tpkln.cc 
-------------------------------------------------------------------------
tpkln.cc : Computes some of the distortions affecting the power spectrum, once a density
inhomogeneity cube is generated and altered. Can be used to compute the following effects: 
  - Gaussian smearing along the z direction in k-space (applied on Fourier coefficients)
  - d rho/rho generation through log-normal or suppression of <-1 values
  - d rho/rho  sampling by galaxies (shot noise)
  
   Usage: Usage: tpkln [options] InputPkFile Out_PPFFile  
   InputPkFile : text file with pairs of values k  Pk on each line 
 options: [-N Nx,Ny,Nz] [-d dx,dy,dz] [-k kmin,kmax,nbin] [-k2d nbin2d,kmax2d] 
            [-neg] [-sneg] [-lognormal] [-thr lowval] [-r rfac] [-rgi] [-prt lev]

tpkln -h provides a more complete description of the arguments and options. 

-------------------------------------------------------------------------
tpkprj.cc 
-------------------------------------------------------------------------
tpkprj.cc : Computes distortions of the recomputed power spectrum, after generating an
input d rho/rho cube, and reprojecting this cube into one or several other cubes,
which might be rotated with respect to the original cube. Includes also the effect of
photoZ smearing (smearing along the radial direction.

 Usage: tpkprj [options] InputPkFile Out_PPFFile 
   InputPkFile : text file with pairs of values k  Pk on each line
 options: [-grid Nx,Ny,Nz,dx,dy,dz] [-in Nx,Ny,Nz,dx,dy,dz] [-out Nx,Ny,Nz,dx,dy,dz] 
          [-center center_dist,tet0,phi0] [-ic center_dist,tet0,phi0]  
          [-oc center_dist,tet0,phi0] [-moc center_dist,tet0,phi0] 
          [-k kmin,kmax,nbin] [-k2d nbin,kmax] 
          [-neg] [-sneg] [-lognormal] [-ltn lowval] [-r rfac] 
          [-s sigmaz] [-ngal ngal_per_cell] [-poiss]  
          [-nosc] [-N nloop] [-nthr NThreads] [-rgi] [-prt lev] [-dbg lev]  
   -grid Nx,Ny,Nz,dx,dy,dz : define 3D grid  default=400x400x400,6,6,6  (6 MPc)
   -in Nx,Ny,Nz,dx,dy,dz : define input 3D grid  default=401x401x401,6,6,6  (6 MPc)
   -out Nx,Ny,Nz,dx,dy,dz : define out 3D grid  default=251x251x251,8,8,8 (8 Mpc)
   -center center_dist,theta0,phi0: define grid center distance(in Mpc) and direction default=3000,0.,0.
   -ic center-dist,theta0,phi0: define input grid center distance(in Mpc) and direction default=3000,0.,0.
   -oc center-dist,theta0,phi0: define output grid center distance(in Mpc) and direction default=3000,0.,0.
   -moc center-dist,theta0,phi0: definition of multiple output grids 
   -k nbin,kmin,kmax : define number of bins and k-range for P(k) computation (def=200,0.,1.) 
   -k2d nbin,kmax : define number of bins and k-max for 2D computation (def=100,1.) 
   -ltn lowval : set input map cells with val<lowval to lowval before computing P(k) (def=-1)
   -neg: set input map cells with val<lowval to lowval (def=-1) before computing P(k) - default behaviour 
   -sneg: set input map cells with val<-0.9 smoothly to -1 (exp(-3.5 x^4)) before computing P(k) 
   -lognormal: replace input map cells with exp(cell_val) 
   -r rfac : P(k) renormalisation factor 
   -s sigmaz : gaussian smearing along z (in Mpc units) default=0.
   -ngal density : Specify galaxy density (number of galaxies / cell) 
   -poiss : Generate Poisson randoms using the specified per cell galaxy density 
   -nosc : use Pk-no-oscillation from the SimLSS output  
   -N nloop : number of different sky cubes generated (default=1) 
   -nthr NThreads: number of parallel threads used for computation 
   -rgi : Auto_initialize Random number generator 
   -dbg lev : activate debug flag and defines debug level 
   -prt lev : define print level (0,1,2..) 


-------------------------------------------------------------------------
tpk2d.cc 
-------------------------------------------------------------------------
tpk2d.cc : small utility program to compute correlation function from power spectrum,
generate grid according to a power spectrum, compute back P(k), 2D power spectrum P(k)-2D
and 1D and 2D correlation function

 Usage: tpk2d [options] InputPkFile Out_PPFFile 
   InputPkFile : text file with pairs of values k  Pk on each line 
       Specify  InputPkFile = '.' or '-' for default P(k) defined in corfunc.h 
 options: [-N Nx,Ny,Nz] [-d dx,dy,dz] [-k kmin,kmax,nbin] [-k2d nbin2d,kmax2d] 
          [-s sigmaz] [-t lowval] [-r rfac] [-p lev] 
   -N Nx,Ny,Nz : define input 3D grid size default=400x400x400 
   -d dx,dy,dz : define input 3D grid cell size (in Mpc) default=5x5x5 Mpc
   -k nbin,kmin,kmax : define number of bins and k-range for P(k) computation (def=200,0.,1.) 
   -k2d nbin2d,kmax2d : define number of bins and k-max for 2D-P(k) computation (def=100,1.) 
   -s sigmaz : gaussian smearing along z (in Mpc units) default=0.
   -p lev : define print level (0,1,2..) 
   -t lowval : set input map cells with val<lowval to lowval before computing P(k) **NOT USED** 
   -r rfac : P(k) renormalisation factor **NOT USED**  

##  Compute reconstructed P(k) and xsi(r) starting from Xsi1 test auto-correlation function
# with 7.5 Mpc smearing along z
./Objs/tpk2d -s 7.5 - totos.ppf

## Check and display some of the results using piapp
spiapp -term

setaxesatt 'font=helvetica,bold,20 fixedfontsize minorticks'
openppf totos.ppf
newwin 2 1 1000 500
disp vpkinterp "arra_xlim=0.,${vpkinterp.info.kmax}  red"
disp recPk 'same black'
disp vxsipki "arra_xlim=0.,${vxsipki.info.rmax}  red"
disp recXsi  'same black'

newwin 2 1 1000 500
disp recPk2D 'xylimits=-0.4,0.4,-0.4,0.4 h2disp=img colbr128'
disp recXsi2D  'xylimits=-200,200,-200,200 h2disp=img colbr128'



-------------------------------------------------------------------------
galcatext.cc : to extract rows from a fits bin table (Galaxy catalog) and write 
 the subset to a new fits file
csh> Objs/galcatext -h
 ---- galcatext: fits catalog extraction program ---- 
 Usage: galcatext InFitsFile OutFitsFile Range [HDU=2] [SegSize=8192] 
    Range: start,end,step  (starting from zero) 

##  To extract the first 10000 rows (galaxies)
csh> rm smallcat.fits  
csh> Objs/galcatext bigcat.fits smallcat.fits 0,10000,1 
##  To extract one every 50 rows (galaxies) - 2%
csh> rm smallcat.fits  
csh> Objs/galcatext bigcat.fits smallcat.fits 0,0,50

 
-------------------------------------------------------------------------
grid2pk.cc : compute power spectrum (P(k)) from grids of de rho/ro or ngals

 Objs/grid2pk -h
 SophyaInitiator::SophyaInitiator() BaseTools Init
 PIOPersist::Initialize() Starting Sophya Persistence management service 
 Usage: grid2pk [options] In3DMapFitsName OutPkTextFile [OutPk_PPFFile] 
 options: [-d dx,dy,dz] [-k kmin,kmax,nbin] [-t lowval] [-r rfac] [-p lev] 
   -d dx,dy,dz : define input 3D map cell size (in Mpc) 
   -k nbin,kmin,kmax : define number of bins and k-range for P(k) computation 
   -t lowval : set input map cells with val<lowval to lowval before computing P(k) 
   -r rfac : P(k) renormalisation factor 
   -p lev : define print level (0,1,2..) 



-------------------------------------------------------------------------
Files :
gpkspec.h  gpkspec.cc : class **GFour3DPk**   grid to P(k) and P(k) to grid computations
corfunc.h : correlation function to P(k) and reverse - 1D and 2D computation
myinteg2d.h myinteg2d.cc : Integrator for P(k)-2D <> xsi-2D(r)
hsplfit.h : utility class to interpolate histograms to be represented as functions
Pk_z_1_0.txt : set of points (k,P(k) representing cosmological power spectrum
