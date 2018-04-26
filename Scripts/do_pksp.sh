#!/usr/local/bin/bash

##########################################################################
## Script pour soumission de jobs de simulation/calcul de P(k) au CC-IN2P3
##   Programme PkBAOUtils/tpkprj avec option -spfits 
## Option -spfits pour garder les spectres individuls sous forme d'un 
## tableau 3D au format fits : speca(ix, jy, kz)
## ix : numero de bin en k , jy : numero de simulation
## kz=1 : Spectre P(k) grille avant galaxies/projection
## kz=2 : Spectre P(k) apres projection
## Definir les variables et choisir un des cas ( comment= ... cmd=...)
## avant de lancer l'execution , eventuellement en batch 
##                              R. Ansari  Mars 2017 
##########################################################################


### Displaying all the limits (ulimit is a bash command) 
ulimit -a 
## Initializing SOPHYA 
source /sps/lsst/Library/Sophya/env.sh new 
# source /sps/lsst/Library/Sophya/env.sh V2.410
# source /afs/in2p3.fr/throng/baoradio/Library/Sophya/env.sh new 
# rehash

echo $LD_LIBRARY_PATH 

EXEPATH=/sps/lsst/dev/ansari/PkBAOUtils/Objs/
#  Define here the path for input and output files 
INOUT=/sps/lsst/dev/ansari/Cecile/
## Input P(k) definition file
PKDEFFILE=${INOUT}/simlss.txt
## Nombre de generation 
NLOOP=5
## Nombre de galaxies par cellule
NGAL=0.1
# -- PhotoZ smearing   (en Mpc) 
SIGPZ=50

#### FOR Grid Engine submission 
## set LOGP = /sps/lsst/dev/ansari/Cecile/
# qsub -l ct=3000,os=sl6,vmem=6G,sps=1 -o ${LOGP} -e ${LOGP} -j y  ./do_pkprj.sh 


echo '----------- START of  do_pksp.sh script -----------'
date 

#############################################################
# Projection en passant par les galaxies avec smearing radiale (photoZ)
comment=" (projection through galaxies/Poisson ngal=${NGAL}, Pz-sigma=${SIGPZ}Mpc)"
cmd="${EXEPATH}tpkprj -in 501,501,250,6,6,6 -out 351,351,150 -ic 3000,0.,0. -oc 3000,0.,0. -k 120,0.005,0.605 -neg -N $NLOOP -ngal $NGAL -poiss -rgi -s ${SIGPZ} -spfits ${INOUT}singlegrid.fits ${PKDEFFILE} ${INOUT}singlegrid.ppf" 


#############################################################
# Projection en passant par les galaxies - grilles multiples 
# comment=" (projection galaxies/Poisson ngal=${NGAL}, MultipleGrids)" 
# cmd="${EXEPATH}tpkprj -in 1100,1100,300,6,6,6 -out 300,300,150,6,6,6 -ic 3000,90.,90. -moc 3000,75.,75. -moc 3000,75.,105.  -moc 3000,105.,75. -moc 3000,105.,105. -k 120,0.005,0.605 -neg -N $NLOOP -ngal 10 -poiss -rgi -spfits ${INOUT}multigrid.fits ${PKDEFFILE} ${INOUT}multigrid.ppf"


echo '--------- (1) Executing command for ' $comment
echo $cmd 
echo '---- Starting...'
time $cmd 

echo '------------- END of do_pkprj.sh script ----------------'
qstat -j ${JOB_ID} -nenv


