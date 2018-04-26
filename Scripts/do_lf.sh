#!/bin/bash

##########################################################################
## Script de calcul de nGal pour differentes fct de luminosite
##########################################################################

## Initializing SOPHYA au CC
# source /sps/lsst/Library/Sophya/env.sh new 

PKBAOBASE=/Users/ansari/Work/LSST/PkBAOUtils/
OUTP=/Users/ansari/Work/LSST
## Au CC 
# PKBAOBASE=/sps/lsst/users/ansari/PkBAOUtils
# OUTP=/sps/lsst/users/ansari/Cecile
EXEPATH=${PKBAOBASE}/LFCode/Objs/

### Pour LF-Dahlen 
## Repartition Standard -fEll 1,0,0,0,0,0 -fSp 0,0.4,0.4,0.2,0.,0. -fSB 0,0,0.,0.,0.5,0.5
SEDDIST="-fEll 1,0,0,0,0,0 -fSp 0,0.4,0.4,0.2,0.,0. -fSB 0,0,0.,0.,0.5,0.5"
rm ${OUTP}/dahlen.ppf  ${OUTP}/dahlen.fits
${EXEPATH}/tngal -i ${PKBAOBASE}/LFCode/lfparamDahlenAll.txt -sedp ${PKBAOBASE}/SEDs/ -filtp ${PKBAOBASE}/Filters/ $SEDDIST -limag 25.3,0. -o ${OUTP}/dahlen.ppf -fits ${OUTP}/dahlen.fits

### Pour LF-Ramos 
SEDDIST="-fEll 1,0,0,0,0,0 -fSp 0,0.4,0.4,0.2,0.,0. -fSB 0,0,0.,0.,0.5,0.5"
rm ${OUTP}/ramos.ppf  ${OUTP}/ramos.fits
${EXEPATH}/tngal -i ${PKBAOBASE}/LFCode/lfparamRamosAll.txt -LFiS -sedp ${PKBAOBASE}/SEDs/ -filtp ${PKBAOBASE}/Filters/ $SEDDIST -limag 25.3,0. -o ${OUTP}/ramos.ppf -fits ${OUTP}/ramos.fits

### Pour LF-Zucca
SEDDIST="-fEll 1,0,0,0,0,0 -fSp 0,0.4,0.4,0.2,0.,0. -fSB 0,0,0.,0.,0.5,0.5"
rm ${OUTP}/zucca.ppf  ${OUTP}/zucca.fits
${EXEPATH}/tngal -i ${PKBAOBASE}/LFCode/lfparamZuccaAllFixed.txt -sedp ${PKBAOBASE}/SEDs/ -filtp ${PKBAOBASE}/Filters/ $SEDDIST -limag 25.3,0. -o ${OUTP}/zucca.ppf -fits ${OUTP}/zucca.fits

