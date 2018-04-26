/*  ------------------------ Projet BAO/PhotoZ/LSST -------------------- 
  Programme de calcul du spectre de puissance (3D) a partir d'un 
  cube de donnees (delta rho/rho ou NGal )
    R. Ansari (pour Adline Choyer) - Feb 2015 
  Usage: grid2pkgrid2pk [options] In3DMap_FitsName OutPk_TextFile [OutPk_PPFFile]
         Options: [-d dx,dy,dz] [-k kmin,kmax,nbin] [-t lowval] [-r rfac] [-p lev]
---------------------------------------------------------------  */

#include "jregrid.h"


HProf Grid2Pk(TArray<TF> grid, double odx, double ody, double odz, int nkbin, double kmin, double kmax, SLinInterp1D Pk, bool gen=false, int prtlev=0, double lowval=-9.e19, bool fgtrunclowval=false) {

  GFour3DPk  gpkc(grid);
  gpkc.SetPrtLevel(prtlev);
  gpkc.SetGridCellSize(odx,ody,odz);
  if(gen) {
    cout << "Generating Fourier amp from Pk " << endl;
    gpkc.generateFourierAmp(Pk);
    gpkc.doInverseFFT();
  }
  if (fgtrunclowval) {  // truncating grid value below threshold 
    double mean, sigma;
    //cout << "jgrid2pk[2] : calling CleanNegatives() ... " << endl;    
    gpkc.CleanNegatives(lowval);
    MeanSigma(grid, mean, sigma);
    //cout << "... After CleanNegatives grid Mean=" << mean << " Sigma=" << sigma << endl;
  }
  //cout << "jgrid2pk[2.c] : Computing Fourier coefficients ... " << endl;    
  gpkc.doFFT();
  
  //cout << "calcpk[3] : computing power spectrum ... " << endl;
  HProf hpk = gpkc.ComputePk(nkbin,kmin,kmax,true);
  return(hpk);
}


//-------------------------------------------------------------------------
//      ------------------ MAIN PROGRAM ------------------------------
//-------------------------------------------------------------------------
int main(int narg, const char* arg[])
{
  if ((narg<3)||((narg>1)&&(strcmp(arg[1],"-h")==0))) {
    cout << " Usage: jgrid2pk [options] InPk_textfile OutPk_TextFile [OutPk_PPFFile] \n" 
	 << " options: [-redshift z,dz] [-d dx,dy,dz] [-k kmin,kmax,nbin] [-t lowval] [-r rfac] [-p lev] \n"
	 << "   -redshift z,dz : input/output 3D grid center redshift, redshift step for output cube \n"
      	 << "   -i nx,ny,nz,dx,dy,dz : define input cube number of cells and cell size (in Mpc) \n"
	 << "   -o nx,ny,nz,dx,dy,dz : define output cube number of cells and cell size \n"
	 << "   -k nbin,kmin,kmax : define number of bins and k-range for P(k) computation \n"
	 << "   -t lowval : set input map cells with val<lowval to lowval before computing P(k) \n"
	 << "   -r rfac : P(k) renormalisation factor \n"
	 << "   -p lev : define print level (0,1,2..) \n" << endl;
    return 1;
  }
  Timer tm("jgrid2pk");
  int rc = 0;
  try { 
    string inpkname;
    string outtextname;
    string outppfname;

    bool fgsetcellsize=false;
    int nx=100,ny=100,nz=100;
    double dx=1.,dy=1.,dz=1.;
    int onx=100,ony=100,onz=100;
    double odx=1.,ody=1.,odz=1.;
    double zcenter=0.5;  // central redshift of input/output cube
    double phoz_err=0.02;  // photo-z uncertainty
    int nkbin=200;
    double kmin=0.001, kmax=1;
    bool fgtrunclowval=false;
    double lowval=-9.e19;
    bool fgrfac=false;
    double rfac=1.;
    int prtlev=0;
    double mean, sigma;
    //----------------------------------------------------------
    // decodage arguments optionnel 
    bool fgoptarg=true;
    while (fgoptarg&&(narg>1)) {
      string fbo = arg[1];
      if (fbo=="-i")  {   // specification taille de cellules (Mpc) drho/rho ou ngal 
	if (narg<3) { cout << " jgrid2pk/missing/bad argument, jgrid2pk -h for help " << endl;  return 2; }
	sscanf(arg[2],"%d,%d,%d,%lf,%lf,%lf",&nx,&ny,&nz,&dx,&dy,&dz); 
	arg+=2; narg-=2; 
      }
      else if (fbo=="-o")  {   // specification taille de cellules (Mpc) drho/rho ou ngal 
	if (narg<3) { cout << " jgrid2pk/missing/bad argument, jgrid2pk -h for help " << endl;  return 2; }
	sscanf(arg[2],"%d,%d,%d,%lf,%lf,%lf",&onx,&ony,&onz,&odx,&ody,&odz); 
	arg+=2; narg-=2; 
      }
      else if (fbo=="-redshift")  {   // redshift central du cube output et erreur photoz
	if (narg<3) { cout << " jgrid2pk/missing/bad argument, jgrid2pk -h for help " << endl;  return 2; }
	sscanf(arg[2],"%lf,%lf",&zcenter,&phoz_err);  arg+=2; narg-=2; 
      }
      else if (fbo=="-k")  {   // specification taille de cellules (Mpc) drho/rho ou ngal 
	if (narg<3) { cout << " jgrid2pk/missing/bad argument, jgrid2pk -h for help " << endl;  return 2; }
	sscanf(arg[2],"%d,%lf,%lf",&nkbin,&kmin,&kmax);  arg+=2; narg-=2; 
      }
      else if (fbo=="-t")  { 
	if (narg<3) { cout << " jgrid2pk/missing/bad argument, jgrid2pk -h for help " << endl;  return 2; }
	lowval=atof(arg[2]);  fgtrunclowval=true; arg+=2; narg-=2; 
      }
      else if (fbo=="-r")  { 
	if (narg<3) { cout << " jgrid2pk/missing/bad argument, jgrid2pk -h for help " << endl;  return 2; }
	rfac=atof(arg[2]);  fgrfac=true; arg+=2; narg-=2; 
      }
      else if (fbo=="-p")  { 
	if (narg<3) { cout << " jgrid2pk/missing/bad argument, jgrid2pk -h for help " << endl;  return 2; }
	prtlev=atoi(arg[2]);  fgrfac=true; arg+=2; narg-=2; 
      }
      else fgoptarg=false;
    }
    //----------------------------------------------------------
    if (narg<3) { cout << " jgrid2pk/missing/bad argument, jgrid2pk -h for help " << endl;  return 2; }
    inpkname=arg[1];
    outtextname=arg[2];
    if (narg>3) outppfname=arg[3];

    //------ Test de la classe InterpCosmoClacl
    SimpleUniverse su;
    InterpCosmoCalc isu(su);
    CosmoCoord csu(su);
    /*
    // on construit un NTuple pour tester
    const char * names[2]={"z","dlos"};
    NTuple nt(2,names);
    double xnt[10];
    for(int i=-10; i<10; i++) {
      double z = zcenter + i*dzredshift;
      double Dc = csu.getLOS(z);
      double zcalc = csu.getRedshift(Dc);
      cout << "z,Dc: "<<z<<" "<< Dc << " "<< zcalc << endl;
      xnt[0]=z;  xnt[1]=Dc;
      nt.Fill(xnt);
    }
    cout << nt;
    FitsInOutFile fos("!jnt.fits",FitsInOutFile::Fits_Create);
    fos<<nt;
    cout << "  nt ecrit ds fichier jnt.fits "<<endl;
    */


    //------ Initialisation du P(k)
    /*
    cout << "jgrid2pk[1] : reading input Pk from file " << inpkname << endl;
    TMatrix<r_8> pkin;
    sa_size_t nr, nc;
    ifstream is(inpkname.c_str());

    pkin.ReadASCII (is, nr, nc);
    Vector ks = pkin.Column(0);
    Vector Pks = pkin.Column(1);
    vector<double> vx = ks.ConvertTostdvec();
    vector<double> vy = Pks.ConvertTostdvec();
    SLinInterp1D Pk(vx,vy);
    */
    
    //PowerSpec Pk(inpkname,zcenter,kmin,kmax);
    PowerSpecCalc Pk(su,true,zcenter,kmin,kmax,nkbin);
    PowerSpecCalc Pknosc(su,false,zcenter,kmin,kmax,nkbin);

    /*
    GFour3DPk  gpkci(ingrid);
    gpkci.SetPrtLevel(prtlev);
    gpkci.SetGridCellSize(dx,dy,dz);
    cout << "jgrid2pk[1.a] : Generating Fourier amp from Pk "<< endl;
    gpkci.generateFourierAmp(Pk);
    gpkci.doInverseFFT();
    double mean, sigma;
    MeanSigma(ingrid, mean, sigma);
    cout << "jgrid2pk[1.b] Input grid sizes " << ingrid.InfoString() << endl;
    ingrid.Show(); 
    cout << "... Input grid Mean=" << mean << " Sigma=" << sigma << endl;
    //tm.Split(" After generating input grid ");
    cout << "jgrid2pk[1.c] : Computing Fourier coefficients ... " << endl;    
    gpkci.doFFT();
    cout << "calcpk[1.d] : computing power spectrum ... " << endl;
    HProf hpki = gpkci.ComputePk(nkbin,kmin,kmax,true);
    */

    // --------- Compute initial grids --------
    sa_size_t nosc = 10;
    sa_size_t siz[3]; siz[0]=nx; siz[1]=ny; siz[2]=nz;
    TArray<TF> ingrid(nx,ny,nz);
    std::vector< TArray<TF> > ingridnosc(nosc);
    ingridnosc[0].SetSize(3,siz);
    cout << "jgrid2pk[2.b] : GFour3DPk  gpkc(outgrid) ... " << endl;
    HProf hpki = Grid2Pk(ingrid,dx,dy,dz,nkbin,kmin,kmax,Pk,true,prtlev,lowval,fgtrunclowval);
    cout << "jgrid2pk[2.c] : GFour3DPk first NoOsc" << endl;
    HProf hpkinosc = Grid2Pk(ingridnosc[0],dx,dy,dz,nkbin,kmin,kmax,Pknosc,true,prtlev,lowval,fgtrunclowval);
    for(int n=1; n<nosc; n++) {
      cout << "jgrid2pk[2.c] : GFour3DPk NoOsc No" << n+1<<endl;
      ingridnosc[n].SetSize(3,siz);
      hpkinosc += Grid2Pk(ingridnosc[n],dx,dy,dz,nkbin,kmin,kmax,Pknosc,true,prtlev,lowval,fgtrunclowval);
    }
    HProf hpkiratio = hpki; hpkiratio /= hpkinosc;

    // --------- Reproject cartesian grid into another cartesian grid --------
    cout << "jgrid2pk[3.a] : ReProjectGrid reproj(ingrid, csu, zcenter) ... " << endl;    
    ReProjectGrid reproj(ingrid, csu, zcenter);
    reproj.SetInputGridCellSize(dx,dy,dz);
    reproj.DefineOutputGrid(onx,ony,onz,odx,ody,odz);
    reproj.FillOutputGrid();
    TArray<TF> outgrid=reproj.getOutputGrid();
    cout << "jgrid2pk[3.b] : GFour3DPk  gpkc(outgrid) ... " << endl;
    HProf hpk = Grid2Pk(outgrid,odx,ody,odz,nkbin,kmin,kmax,Pk,false,prtlev,lowval,fgtrunclowval);
    
    ReProjectGrid reprojnosc(ingridnosc[0], csu, zcenter);
    reprojnosc.SetInputGridCellSize(dx,dy,dz);
    reprojnosc.DefineOutputGrid(onx,ony,onz,odx,ody,odz);
    reprojnosc.FillOutputGrid();
    std::vector< TArray<TF> > outgridnosc(nosc);
    outgridnosc[0] = reprojnosc.getOutputGrid();
    cout << "jgrid2pk[3.c] : GFour3DPk first NoOsc"<<endl;
    HProf hpknosc = Grid2Pk(outgridnosc[0],odx,ody,odz,nkbin,kmin,kmax,Pknosc,false,prtlev,lowval,fgtrunclowval);
    for(int n=1; n<nosc; n++) {
      reprojnosc.UpdateInputGrid(ingridnosc[n]);
      reprojnosc.SetInputGridCellSize(dx,dy,dz);
      reprojnosc.DefineOutputGrid(onx,ony,onz,odx,ody,odz);
      reprojnosc.FillOutputGrid();
      outgridnosc[n] = reprojnosc.getOutputGrid();
      cout << "jgrid2pk[3.c] : GFour3DPk NoOsc No" << n+1<<endl;
      hpknosc += Grid2Pk(outgridnosc[n],odx,ody,odz,nkbin,kmin,kmax,Pknosc,false,prtlev,lowval,fgtrunclowval);
    }  
    HProf hpkratio = hpk; hpkratio /= hpknosc;
    
    /*
    GFour3DPk  gpkc(outgrid);
    gpkc.SetPrtLevel(prtlev);
    gpkc.SetGridCellSize(odx,ody,odz);
  
    if (fgtrunclowval) {  // truncating grid value below threshold 
      cout << "jgrid2pk[2] : calling CleanNegatives() ... " << endl;    
      gpkc.CleanNegatives(lowval);
      MeanSigma(outgrid, mean, sigma);
      cout << "... After CleanNegatives grid Mean=" << mean << " Sigma=" << sigma << endl;
      tm.Split(" After CleanNegatives ");
    }
    cout << "jgrid2pk[2.c] : Computing Fourier coefficients ... " << endl;    
    gpkc.doFFT();
    //tm.Split(" After doFFT ");

    cout << "calcpk[3] : computing power spectrum ... " << endl;
    HProf hpk = gpkc.ComputePk(nkbin,kmin,kmax,true);
    DataTable dtpk; 
    Histo hrap = gpkc.FillPkDataTable(dtpk, rfac);
    tm.Split(" After ComputePk ");
    if (prtlev>0)  { 
      dtpk.SetShowMinMaxFlag(true);
      cout << dtpk; 
    }
    */

    // --------- Reproject cartesian grid into redshift grid --------
    cout << "jgrid2pk[4.a] : ReProjectGrid into redshift grids ... " << endl;   
    reproj.DefineOutputGridRedshift();
    reproj.FillOutputGridRedshift();
    TArray<TF> outgridz=reproj.getOutputGridRedshift();
    tm.Split(" After FillOutputGridAngleRedshift ");
    cout << "jgrid2pk[4.b] : GFour3DPk  gpkc(outgridz) ... " << endl;
    HProf hpkz = Grid2Pk(outgridz,odx,ody,odz,nkbin,kmin,kmax,Pk,false,prtlev,lowval,fgtrunclowval);
    
    reprojnosc.UpdateInputGrid(ingridnosc[0]);
    reprojnosc.SetInputGridCellSize(dx,dy,dz);
    reprojnosc.DefineOutputGrid(onx,ony,onz,odx,ody,odz);
    reprojnosc.DefineOutputGridRedshift();
    reprojnosc.FillOutputGridRedshift();
    std::vector< TArray<TF> > outgridznosc(nosc);
    outgridznosc[0] = reprojnosc.getOutputGridRedshift();
    cout << "jgrid2pk[4.c] : GFour3DPk first NoOsc"<<endl;
    HProf hpkznosc = Grid2Pk(outgridznosc[0],odx,ody,odz,nkbin,kmin,kmax,Pknosc,false,prtlev,lowval,fgtrunclowval);
    for(int n=1; n<nosc; n++) {
      reprojnosc.UpdateInputGrid(ingridnosc[n]);
      reprojnosc.SetInputGridCellSize(dx,dy,dz);
      reprojnosc.DefineOutputGrid(onx,ony,onz,odx,ody,odz);
      reprojnosc.DefineOutputGridRedshift();
      reprojnosc.FillOutputGridRedshift();
      outgridznosc[n] = reprojnosc.getOutputGrid();
      cout << "jgrid2pk[4.c] : GFour3DPk NoOsc No" << n+1<<endl;
      hpkznosc += Grid2Pk(outgridznosc[n],odx,ody,odz,nkbin,kmin,kmax,Pknosc,false,prtlev,lowval,fgtrunclowval);
    }
    hpkznosc *= nosc;
    HProf hpkzratio = hpkz; hpkzratio /= hpkznosc;
    /*
    GFour3DPk  gpkcz(outgridz);
    gpkcz.SetPrtLevel(prtlev);
    gpkcz.SetGridCellSize(odx,ody,odz); 
    if (fgtrunclowval) {  // truncating grid value below threshold 
      cout << "jgrid2pk[4] : calling CleanNegatives() ... " << endl;    
      gpkcz.CleanNegatives(lowval);
      MeanSigma(outgridz, mean, sigma);
      cout << "... After CleanNegatives grid Mean=" << mean << " Sigma=" << sigma << endl;
      tm.Split(" After CleanNegatives ");
    }
    cout << "jgrid2pk[4.b] : Computing Fourier coefficients ... " << endl;    
    gpkcz.doFFT();
    //tm.Split(" After doFFT ");
    cout << "calcpk[4.c] : computing power spectrum ... " << endl;
    HProf hpkz = gpkcz.ComputePk(nkbin,kmin,kmax,true);
    //tm.Split(" After ComputePk ");
    */
   
    // --------- Add photo-z uncertainty --------
    cout << "jgrid2pk[5] : Adding photo-z uncertainty ... " << endl;    
    TArray< TF > outgridphoz = reproj.AddPhotoZ(outgridz,phoz_err);
    cout << "jgrid2pk[5] : GFour3DPk  gpkc(outgridz) ... " << endl;
    HProf hpkphoz = Grid2Pk(outgridphoz,odx,ody,odz,nkbin,kmin,kmax,Pk,false,prtlev,lowval,fgtrunclowval);
    
    std::vector< TArray< TF > > outgridphoznosc(nosc);
    cout << "jgrid2pk[5.c] : AddPhotoZ first NoOsc"<<endl;
    reprojnosc.UpdateInputGrid(ingridnosc[0]);
    reprojnosc.SetInputGridCellSize(dx,dy,dz);
    reprojnosc.DefineOutputGrid(onx,ony,onz,odx,ody,odz);
    outgridphoznosc[0] = reprojnosc.AddPhotoZ(outgridznosc[0],phoz_err);
    HProf hpkphoznosc = Grid2Pk(outgridphoznosc[0],odx,ody,odz,nkbin,kmin,kmax,Pknosc,false,prtlev,lowval,fgtrunclowval);
    for(int n=1; n<nosc; n++) {
      cout << "jgrid2pk[5.c] : AddPhotoZ NoOsc No "<<n+1<<endl;
      reprojnosc.UpdateInputGrid(ingridnosc[n]);
      reprojnosc.SetInputGridCellSize(dx,dy,dz);
      reprojnosc.DefineOutputGrid(onx,ony,onz,odx,ody,odz);
      outgridphoznosc[n] = reprojnosc.AddPhotoZ(outgridznosc[n],phoz_err);
      hpkphoznosc += Grid2Pk(outgridphoznosc[n],odx,ody,odz,nkbin,kmin,kmax,Pknosc,false,prtlev,lowval,fgtrunclowval);
    }
    hpkphoznosc *= nosc;
    HProf hpkphozratio = hpkphoz; hpkphozratio /= hpkphoznosc;
    
    /*
    GFour3DPk  gpkcphoz(outgridphoz);
    gpkcphoz.SetPrtLevel(prtlev);
    gpkcphoz.SetGridCellSize(odx,ody,odz); 
    if (fgtrunclowval) {  // truncating grid value below threshold 
      cout << "jgrid2pk[5] : calling CleanNegatives() ... " << endl;    
      gpkcphoz.CleanNegatives(lowval);
      MeanSigma(outgridphoz, mean, sigma);
      cout << "... After CleanNegatives grid Mean=" << mean << " Sigma=" << sigma << endl;
      tm.Split(" After CleanNegatives ");
    }
    cout << "jgrid2pk[5] : Computing Fourier coefficients ... " << endl;    
    gpkcphoz.doFFT();
    //tm.Split(" After doFFT ");
    cout << "calcpk[5] : computing power spectrum ... " << endl;
    HProf hpkphoz = gpkcphoz.ComputePk(nkbin,kmin,kmax,true);
    // dernier argument a true -> on peut recuperer l'histo count
    Histo hcntmodphoz = gpkcphoz.getModCount();
    */
    
    //tm.Split(" After ComputePk ");
   
    // Saving computed P(k) to text (ascii) file
    /*
    cout << "calcpk[4] : wrting P(k) to text file " << outtextname << endl;
    ofstream tof(outtextname.c_str());
    dtpk.WriteASCII(tof);
    */
    
    // Saving computed P(k) to PPF file
    if (outppfname.length()>0)  {  
      cout << "calcpk[5] : writing profile histo P(k) (hpk) and DataTable P(k) dtpk to PPF file  " << outppfname << endl;
      POutPersist po(outppfname);
      po<<PPFNameTag("nt")<<Pk.GetNt();
      //po<<PPFNameTag("pkin")<<pkin;
      po<<PPFNameTag("hpki")<<hpki;
      po<<PPFNameTag("hpkinosc")<<hpkinosc;
      po<<PPFNameTag("hpkiratio")<<hpkiratio;
      po<<PPFNameTag("hpk")<<hpk;
      po<<PPFNameTag("hpknosc")<<hpknosc;
      po<<PPFNameTag("hpkratio")<<hpkratio;
      po<<PPFNameTag("hpkz")<<hpkz;
      po<<PPFNameTag("hpkznosc")<<hpkznosc;
      po<<PPFNameTag("hpkzratio")<<hpkzratio;
      po<<PPFNameTag("hpkphoz")<<hpkphoz;
      po<<PPFNameTag("hpkphoznosc")<<hpkphoznosc;
      po<<PPFNameTag("hpkphozratio")<<hpkphozratio;
      //po<<PPFNameTag("hcntmodphoz")<<hcntmodphoz;
      //po<<PPFNameTag("dtpk")<<dtpk;
    }
  }  // End of try bloc 
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " jgrid2pk.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " jgrid2pk.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " jgrid2pk.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
  }
  cout << " ==== End of jgrid2pk.cc program  Rc= " << rc << endl;
  return rc;    
}
