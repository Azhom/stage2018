###################################################################
##     Project   LSST/BAO/PhotoZ
##     Makefile to compile multiTypeZlf
## F. Habibi & R.Ansari ,  LAL (CNRS/IN2P3)  & Univ. Paris Sud
##                   April - July 2017      
###################################################################

include $(SOPHYABASE)/include/sophyamake.inc
# contains symbolic links to compile SOPHYA's libraries

MYOBJS = Objs/multitypzlf.o Objs/myschechter.o \
         Objs/kcorr.o Objs/galcnt.o \
         Objs/sinterp.o Objs/reddening.o  Objs/sedfilter.o

MYINCS = myschechter.h multitypzlf.h kcorr.h galcnt.h  sedfilter.h  reddening.h sinterp.h

all :   mlf mapngal

clean :
	rm Objs/* *.ppf 

depend :
	mkdir Objs/


mlf :   Objs/mlf	
	@echo 'Makefile:    mlf made'
	
mapngal : Objs/mapngal
	@echo 'Makefile:	mapngal made'

######    executables made here
### mlf : Basic test of MultiType_Z_LF  
Objs/mlf: Objs/mlf.o Objs/multitypzlf.o Objs/myschechter.o 
	$(CXXLINK) -o Objs/mlf Objs/mlf.o Objs/multitypzlf.o Objs/myschechter.o $(SOPHYAEXTSLBLIST)

Objs/mlf.o: mlf.cc  multitypzlf.h  myschechter.h 
	$(CXXCOMPILE) -I . -o Objs/mlf.o mlf.cc 

#	mapngal : computing map of galaxy number density from dust map
Objs/mapngal: Objs/mapngal.o $(MYOBJS)
	$(CXXLINK) -o Objs/mapngal Objs/mapngal.o $(MYOBJS) $(SOPHYAEXTSLBLIST)
	
Objs/mapngal.o: mapngal.cc  $(MYINCS)
	$(CXXCOMPILE) -I . -o Objs/mapngal.o mapngal.cc 

#   Multi-type, per redshift bin LF's 
Objs/multitypzlf.o: multitypzlf.cc  multitypzlf.h  myschechter.h  
	$(CXXCOMPILE) -I . -o Objs/multitypzlf.o multitypzlf.cc
 
Objs/myschechter.o: myschechter.cc  myschechter.h 
	$(CXXCOMPILE) -I . -o Objs/myschechter.o myschechter.cc

#   K-correction code 
Objs/galcnt.o: galcnt.cc  galcnt.h  multitypzlf.h  myschechter.h  kcorr.h 
	$(CXXCOMPILE) -I . -o Objs/galcnt.o galcnt.cc

#  Galaxy number density count  
Objs/kcorr.o: kcorr.cc  kcorr.h 
	$(CXXCOMPILE) -I . -o Objs/kcorr.o kcorr.cc

#   Files extracted from Alexandra Abate's FZGalSim code
Objs/sinterp.o: sinterp.cc  sinterp.h 
	$(CXXCOMPILE) -I . -o Objs/sinterp.o sinterp.cc

Objs/reddening.o: reddening.cc  reddening.h 
	$(CXXCOMPILE) -I . -o Objs/reddening.o reddening.cc

Objs/sedfilter.o: sedfilter.cc  sedfilter.h 
	$(CXXCOMPILE) -I . -o Objs/sedfilter.o sedfilter.cc
