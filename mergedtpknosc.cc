/*--- Petit programme utilitaire pour fusionner deux DataTable de P(k) 
   calcules par tpkprj.cc , pour P(k) avec et sans oscillation --- 
       R. Ansari, Univ. Paris-Sud, LAL/IN2P3-CNRS      Jan 2017        */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <iostream>
#include <cmath>
#include <complex>

#include "array.h"
#include "histats.h"

// #include "tarrinit.h"
// #include "histinit.h"

#include "swfitsdtable.h"
#include "fiosinit.h"

using namespace std; 
using namespace SOPHYA; 

DataTable MergeDtPks(DataTable& dtpk, DataTable& dtpknos, double V, double sn);

int Usage()
{
  cout << " ---- mergedtpknosc: merging P(k) with and without oscillation datatables ---- \n"
       << " Usage: mergedtpknosc DtPk_PPFFile DtPkNoOsc_PPFFile VolSurvey ShotNoise OutPPFFile \n" <<endl;
  return 0;
}
/* -- MAIN -- */ 
int main (int narg, char* arg[])
{
  if ((narg < 6)||(strcmp(arg[1],"-h")==0)) {
    return Usage();
  }
  cout << " ------ mergedtpknosc : fits catalog extraction program ------- " << endl;

  int rc = 0;
  try {
    SophyaInit();
    Timer tm("mergedtpknosc.cc");
    string dtpkfile = arg[1];
    string dtpknosfile = arg[2];
    double VSurvey=atof(arg[3]);
    double ShotNoise=atof(arg[4]);
    string outfile=arg[5];

    DataTable dtpk;
    DataTable dtpknos;
    {
      cout << "[1] Opening DtPk file: " << dtpkfile <<endl;
      PInPersist pin(dtpkfile);
      pin>>PPFNameTag("dtpk")>>dtpk;
      dtpk.SetShowMinMaxFlag(true);
      cout<<dtpk;
    }
    {
      cout << "[2] Opening DtPk No-Oscillation file: " << dtpknosfile <<endl;
      PInPersist pin(dtpknosfile);
      pin>>PPFNameTag("dtpk")>>dtpknos;
      dtpknos.SetShowMinMaxFlag(true);
      cout<<dtpknos;      
    }
    if (dtpk.NRows() != dtpknos.NRows())  {
      cout<<"ERROR: Different number of rows in the two DataTables -> exit"<<endl;
      return 91;
    }
    cout << "[3] Merging the two DataTables, computing erros with VSurvey=" <<VSurvey<<" Mpc^3" 
	 << " ShotNoise="<<ShotNoise<< " Mpc^3"<<endl;
    DataTable dt=MergeDtPks(dtpk, dtpknos, VSurvey, ShotNoise);
    cout << "[4] Saving output DataTable to PPF file:"<<outfile<<endl;
    POutPersist po(outfile);
    po<<dt;
  }
  catch (PThrowable & exc) {
    cerr << " mergedtpknosc.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << " - Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {
    cerr << " mergedtpknosc.cc: Catched std::exception "  
         << " - what()= " << e.what() << endl;
    rc = 98;
  }
  catch (...) {
    cerr << " mergedtpknosc.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
  }
  cout << " ------------ END OF mergedtpknosc  (Rc= " << rc << ") ------------ " << endl;
  return rc;
}


/*--Fonction--*/
DataTable MergeDtPks(DataTable& dtpk, DataTable& dtpknos, double V, double sn)
{
  const char* nomcol[25] = {"k","inpk","recpkin","prjpk","sig_recpkin","sig_prjpk",
			    "inpk_nos","recpkin_nos","prjpk_nos","sig_pkin_nos","sig_prjpk_nos",
                            "errinpk", "errinpk_nos","errrecpkin","errrecpkin_nos", "errprjprk", "errprjpk_nos",
                            "rapp_inpk","err_rapp_inpk","rapp_recpkin","err_rapp_recpkin", "rapp_prjpk","err_rapp_prjpk",
                            "sig_rapp_recpkin", "sig_rapp_prjpk" };
  DataTable dt;
  for(int i=0; i<25; i++)   dt.AddDoubleColumn(nomcol[i]);    

  DataTableRow 	dtr = dt.EmptyRow();
  DataTableRow 	dtra = dtpk.EmptyRow();
  DataTableRow 	dtrb = dtpknos.EmptyRow();
  dtpk.GetRow(0,dtra);
  double k0=(double)dtra[0];
  dtpk.GetRow(1,dtra);
  double k1=(double)dtra[0];
  double deltak = (k1-k0);
  cout<<"MergeDtPks: Computing P(k) errors with VSurvey="<<V<<" Mpc^3 and Delta-k="<<deltak<<" Mpc^-1"<<endl;
  double CstSigmaPk=2.*M_PI/sqrt(deltak*V);

  for(size_t i=0; i<dtpk.NRows(); i++) {
    dtpk.GetRow(i,dtra);
    dtpknos.GetRow(i,dtrb);
    double k=dtra[0];   dtr[0]=k;
    for(int j=0; j<6; j++)  dtr[j]=(double)dtra[j];
    for(int j=1; j<6; j++)  {
      dtr[j]=(double)dtra[j];
      dtr[j+5]=(double)dtrb[j];
    }

    // Calcul erreurs P(k) et rapport P(k) / P(k)-noOsc
    for(int j=0; j<3; j++) {
      double pk=dtra[1+j];
      double pknos=dtrb[1+j];
      double sigmaPk=(CstSigmaPk/k)*pk;
      double sigmaPknos=(CstSigmaPk/k)*pknos;
      dtr[11+2*j]=sigmaPk;   dtr[12+2*j]=sigmaPknos;
      dtr[11+6+2*j]=pk/pknos;   dtr[12+6+2*j]=sigmaPk/pknos;
      if (j==2)  dtr[11+6+2*j]=(pk-sn)/pknos; 
      if (j>0)  {
	double siga=dtra[3+j];
	dtr[22+j]=siga/pknos; 
      }
    }
    dt.AddRow(dtr);
  }
  dt.Info()=dtpk.Info();
  dt.SetShowMinMaxFlag(true);
  cout << dt;
  return dt;
}
