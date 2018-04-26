#!/usr/local/bin/bash

##########################################################################
## Script pour soumission de jobs de simulation/calcul de P(k) au CC-IN2P3
##   Programme PkBAOUtils/tpkln 
##                              R. Ansari  Janvier 2017 
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
INOUT=/sps/lsst/dev/ansari/Cecile/Res23Jan/
NLOOP=25
# -- PhotoZ smearing 
SIGPZ=20


echo '[1] ---------------------------------------------------'
date 
cmd="${EXEPATH}tpkln -grid 400,400,200,6,6,6 -k 120,0.005,0.605 -neg -N ${NLOOP} -nosc ${INOUT}simlss.txt ${INOUT}Apknos_neg.ppf"
echo '[1] Executing '
echo $cmd
time $cmd 

echo '[2] ---------------------------------------------------'
date 
cmd="${EXEPATH}tpkln -grid 400,400,200,6,6,6 -k 120,0.005,0.605 -s $SIGPZ -N ${NLOOP} -nosc ${INOUT}simlss.txt ${INOUT}Apknos_s${SIGPZ}.ppf"
echo 'Executing '
echo $cmd
time $cmd 

echo '[3] ---------------------------------------------------'
date 
cmd="${EXEPATH}tpkln -grid 400,400,200,6,6,6 -k 120,0.005,0.605 -neg -s $SIGPZ -N ${NLOOP} -nosc ${INOUT}simlss.txt ${INOUT}Apknos_neg_s${SIGPZ}.ppf"
echo 'Executing '
echo $cmd
time $cmd 

echo '------ END of do_tpkln.sh script ------
qstat -j ${JOB_ID} -nenv


