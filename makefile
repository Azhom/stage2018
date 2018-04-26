#  Makefile pour des programmes utilitaires BAO/PhotoZ
#      R. Ansari , (C) LAL+UPS , Sep 2014    
include $(SOPHYABASE)/include/sophyamake.inc

####  Decommenter cette ligne et defnir les chemins pour SimLSS si compilation de jgrid2pk et l'ajouter dans les targets 
## include make_simlss.inc

#  Define our target list 
# all :  Objs/galcatext Objs/grid2pk Objs/tpk2d Objs/clcut Objs/pk2cl Objs/jgrid2pk Objs/tpkln Objs/jgrid2pk
all :  Objs/galcatext Objs/grid2pk Objs/tpk2d Objs/clcut Objs/tpkln Objs/tpkshift Objs/tpksph Objs/terrpk Objs/tpkprj \
       Objs/mergedtpknosc Objs/tkbaoerr 
clean :
	rm -f Objs/*


##############
Objs/jgrid2pk : Objs/jgrid2pk.o Objs/jregrid.o Objs/gpkspec.o 
	$(CXXLINK) -o Objs/jgrid2pk Objs/jgrid2pk.o Objs/jregrid.o Objs/gpkspec.o $(SOPHYAEXTSLBLIST) $(SIMLSSLIBS)

Objs/jgrid2pk.o : jgrid2pk.cc jregrid.h gpkspec.h 
	$(CXXCOMPILE) $(SIMLSSINC)  -o Objs/jgrid2pk.o  jgrid2pk.cc
###
Objs/jregrid.o : jregrid.cc jregrid.h gpkspec.h 
	$(CXXCOMPILE) $(SIMLSSINC) -o Objs/jregrid.o  jregrid.cc
##############
Objs/galcatext : Objs/galcatext.o 
	$(CXXLINK) -o Objs/galcatext Objs/galcatext.o  $(SOPHYAEXTSLBLIST) 

Objs/galcatext.o : galcatext.cc 
	$(CXXCOMPILE)  -o Objs/galcatext.o  galcatext.cc

##############
Objs/grid2pk : Objs/grid2pk.o Objs/gpkspec.o
	$(CXXLINK) -o Objs/grid2pk Objs/grid2pk.o Objs/gpkspec.o $(SOPHYAEXTSLBLIST) 

Objs/grid2pk.o : grid2pk.cc gpkspec.h 
	$(CXXCOMPILE)  -o Objs/grid2pk.o  grid2pk.cc
##############
Objs/clcut : Objs/clcut.o Objs/gpkspec.o
	$(CXXLINK) -o Objs/clcut Objs/clcut.o Objs/gpkspec.o $(SOPHYAEXTSLBLIST) 

Objs/clcut.o : clcut.cc gpkspec.h 
	$(CXXCOMPILE)  -o Objs/clcut.o  clcut.cc
##############
Objs/pk2cl : Objs/pk2cl.o Objs/gpkspec.o
	$(CXXLINK) -o Objs/pk2cl Objs/pk2cl.o Objs/gpkspec.o $(SOPHYAEXTSLBLIST)

Objs/pk2cl.o : pk2cl.cc gpkspec.h 
	$(CXXCOMPILE)  -o Objs/pk2cl.o  pk2cl.cc
##############
Objs/tpk2d : Objs/tpk2d.o Objs/gpkspec.o Objs/myinteg2d.o 
	$(CXXLINK) -o Objs/tpk2d Objs/tpk2d.o Objs/gpkspec.o Objs/myinteg2d.o $(SOPHYAEXTSLBLIST) 

Objs/tpk2d.o : tpk2d.cc gpkspec.h corfunc.h myinteg2d.h hsplfit.h 
	$(CXXCOMPILE)  -o Objs/tpk2d.o  tpk2d.cc
##############
Objs/tpkln : Objs/tpkln.o Objs/gpkspec.o Objs/gpkutil.o  Objs/myinteg2d.o 
	$(CXXLINK) -o Objs/tpkln Objs/tpkln.o Objs/gpkspec.o Objs/gpkutil.o Objs/myinteg2d.o $(SOPHYAEXTSLBLIST) 

Objs/tpkln.o : tpkln.cc gpkspec.h corfunc.h myinteg2d.h hsplfit.h 
	$(CXXCOMPILE)  -o Objs/tpkln.o  tpkln.cc
##############
Objs/tpkshift : Objs/tpkshift.o Objs/gpkspec.o Objs/gpkutil.o  Objs/myinteg2d.o 
	$(CXXLINK) -o Objs/tpkshift Objs/tpkshift.o Objs/gpkspec.o Objs/gpkutil.o Objs/myinteg2d.o $(SOPHYAEXTSLBLIST) 

Objs/tpkshift.o : tpkshift.cc gpkspec.h corfunc.h myinteg2d.h hsplfit.h 
	$(CXXCOMPILE)  -o Objs/tpkshift.o  tpkshift.cc
##############
Objs/tpksph : Objs/tpksph.o Objs/gpkspec.o Objs/gpkutil.o 
	$(CXXLINK) -o Objs/tpksph Objs/tpksph.o Objs/gpkspec.o Objs/gpkutil.o $(SOPHYAEXTSLBLIST) 

Objs/tpksph.o : tpksph.cc gpkspec.h gpkutil.h corfunc.h myinteg2d.h hsplfit.h 
	$(CXXCOMPILE)  -o Objs/tpksph.o  tpksph.cc
##############
Objs/tpkprj : Objs/tpkprj.o Objs/gpkspec.o Objs/gpkutil.o 
	$(CXXLINK) -o Objs/tpkprj Objs/tpkprj.o Objs/gpkspec.o Objs/gpkutil.o $(SOPHYAEXTSLBLIST) 

Objs/tpkprj.o : tpkprj.cc gpkspec.h gpkutil.h corfunc.h myinteg2d.h hsplfit.h 
	$(CXXCOMPILE)  -o Objs/tpkprj.o  tpkprj.cc
##############
Objs/tkbaoerr : Objs/tkbaoerr.o Objs/gpkspec.o Objs/gpkutil.o Objs/fitkbao.o 
	$(CXXLINK) -o Objs/tkbaoerr Objs/tkbaoerr.o Objs/gpkspec.o Objs/gpkutil.o Objs/fitkbao.o $(SOPHYAEXTSLBLIST) 

Objs/tkbaoerr.o : tkbaoerr.cc gpkspec.h gpkutil.h fitkbao.h corfunc.h myinteg2d.h hsplfit.h 
	$(CXXCOMPILE)  -o Objs/tkbaoerr.o  tkbaoerr.cc
##############
Objs/terrpk : Objs/terrpk.o Objs/gpkspec.o Objs/gpkutil.o Objs/fitkbao.o 
	$(CXXLINK) -o Objs/terrpk Objs/terrpk.o Objs/gpkspec.o Objs/gpkutil.o Objs/fitkbao.o $(SOPHYAEXTSLBLIST) 

Objs/terrpk.o : terrpk.cc gpkspec.h gpkutil.h fitkbao.h corfunc.h myinteg2d.h hsplfit.h 
	$(CXXCOMPILE)  -o Objs/terrpk.o  terrpk.cc
##############
Objs/mergedtpknosc : Objs/mergedtpknosc.o Objs/mergedtpknosc.o
	$(CXXLINK) -o Objs/mergedtpknosc Objs/mergedtpknosc.o $(SOPHYAEXTSLBLIST) 

Objs/mergedtpknosc.o : mergedtpknosc.cc 
	$(CXXCOMPILE)  -o Objs/mergedtpknosc.o  mergedtpknosc.cc

###
Objs/gpkutil.o : gpkutil.cc  gpkutil.h 
	$(CXXCOMPILE)  -o Objs/gpkutil.o  gpkutil.cc

###
Objs/fitkbao.o : fitkbao.cc  fitkbao.h 
	$(CXXCOMPILE)  -o Objs/fitkbao.o  fitkbao.cc

###
Objs/myinteg2d.o : myinteg2d.cc myinteg2d.h 
	$(CXXCOMPILE)  -o Objs/myinteg2d.o  myinteg2d.cc

###
Objs/gpkspec.o : gpkspec.cc gpkspec.h 
	$(CXXCOMPILE)  -o Objs/gpkspec.o  gpkspec.cc
