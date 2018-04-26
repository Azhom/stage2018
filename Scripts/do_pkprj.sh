#!/usr/local/bin/bash

##########################################################################
## Script pour soumission de jobs de simulation/calcul de P(k) au CC-IN2P3
##   Programme PkBAOUtils/tpkprj 
## Il faut decommenter une des series de 2 lignes definissant la commande
## et la description comment= ... cmd=...   avant de soumettre le job
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
NGAL=0.1
# -- PhotoZ smearing 
SIGPZ=50

#### FOR Grid Engine submission 
## set LOGP = /sps/lsst/dev/ansari/Cecile/
# qsub -l ct=3000,os=sl6,vmem=6G,sps=1 -o ${LOGP} -e ${LOGP} -j y  ./do_pkprj.sh 


echo '----------- START of  do_pkprj.sh script -----------'
date 
# Projection directe grille a grille 
# comment='(direct grid to grid projection)'
# cmd="${EXEPATH}tpkprj -in 501,501,250,6,6,6 -out 351,351,150 -ic 3000,0.,0. -oc 3000,0.,0. -k 120,0.005,0.605 -neg -N 16 -ngal 0.5 ${INOUT}simlss.txt ${INOUT}pkprjB0.ppf"
# Projection en passant par les galaxies 
# comment=' (projection through galaxies ngal=10)'
# cmd="${EXEPATH}tpkprj -in 501,501,250,6,6,6 -out 351,351,150 -ic 3000,0.,0. -oc 3000,0.,0. -k 120,0.005,0.605 -neg -N 16 -ngal 10 ${INOUT}simlss.txt ${INOUT}pkprjgalB0.ppf"
# Projection en passant par les galaxies - avec Poisson
# comment=' (projection through galaxies/Poisson ngal=10)'
# cmd="${EXEPATH}tpkprj -in 501,501,250,6,6,6 -out 351,351,150 -ic 3000,0.,0. -oc 3000,0.,0. -k 120,0.005,0.605 -neg -N 16 -ngal 10 -poiss ${INOUT}simlss.txt ${INOUT}pkprjgalB0poiss.ppf"
# Projection en passant par les galaxies - avec Poisson, Pk-NoOsc
# comment=' (projection through galaxies/Poisson NoOsc ngal=10)'
# cmd="${EXEPATH}tpkprj -in 501,501,250,6,6,6 -out 351,351,150 -ic 3000,0.,0. -oc 3000,0.,0. -k 120,0.005,0.605 -neg -N 16 -ngal 10 -poiss -nosc ${INOUT}simlss.txt ${INOUT}pknoscprjgalB0poiss.ppf"

# Projection en passant par les galaxies avec smearing radiale (photoZ)
# comment=" (projection through galaxies/Poisson ngal=${NGAL}, Pz-sigma=${SIGPZ}Mpc)"
# cmd="${EXEPATH}tpkprj -in 501,501,250,6,6,6 -out 351,351,150 -ic 3000,0.,0. -oc 3000,0.,0. -k 120,0.005,0.605 -neg -N $NLOOP -ngal $NGAL -poiss -s ${SIGPZ} ${INOUT}simlss.txt ${INOUT}pkprjB0_s${SIGPZ}.ppf"
# Projection en passant par les galaxies avec smearing radiale (photoZ) pour Pk-NoOsc
# comment=" (projection through galaxies/Poisson ngal=${NGAL}, Pz-sigma=${SIGPZ}Mpc Pk-NoOsc)"
# cmd="${EXEPATH}tpkprj -in 501,501,250,6,6,6 -out 351,351,150 -ic 3000,0.,0. -oc 3000,0.,0. -k 120,0.005,0.605 -neg -N $NLOOP -ngal $NGAL -poiss -s ${SIGPZ} -nosc ${INOUT}simlss.txt ${INOUT}pknosprjB0_s${SIGPZ}.ppf"

# Projection en passant par les galaxies - Grande grille pour comparaison avec grilles multiples
# comment=' (projection galaxies/Poisson ngal=10, LargeGrid)' 
# cmd="${EXEPATH}tpkprj -in 1000,1000,250,6,6,6 -out 600,600,150,6,6,6 -ic 3000,0.,0. -oc 3000,0.,0. -k 120,0.005,0.605 -neg -N 16 -ngal 10 -poiss ${INOUT}simlss.txt ${INOUT}pkprjSLG1.ppf"
# comment=' (projection galaxies/Poisson ngal=10, LargeGrid, Pk-NoOsc)' 
# cmd="${EXEPATH}tpkprj -in 1000,1000,250,6,6,6 -out 600,600,150,6,6,6 -ic 3000,0.,0. -oc 3000,0.,0. -k 120,0.005,0.605 -neg -N 16 -ngal 10 -poiss -nosc ${INOUT}simlss.txt ${INOUT}pknosprjSLG1.ppf"
# comment=' (projection galaxies/Poisson ngal=10, LargeGrid- Rotated 10deg)' 
# cmd="${EXEPATH}tpkprj -in 1000,1000,250,6,6,6 -out 600,600,150,6,6,6 -ic 3000,0.,0. -oc 3000,10.,0. -k 120,0.005,0.605 -neg -N 16 -ngal 10 -poiss ${INOUT}simlss.txt ${INOUT}pkprjSLG1_R10deg.ppf"
# comment=' (projection galaxies/Poisson ngal=10, LargeGrid- Rotated 10deg Pk-NoOsc)' 
# cmd="${EXEPATH}tpkprj -in 1000,1000,250,6,6,6 -out 600,600,150,6,6,6 -ic 3000,0.,0. -oc 3000,10.,0. -k 120,0.005,0.605 -neg -N 16 -ngal 10 -poiss -nosc ${INOUT}simlss.txt ${INOUT}pknosprjSLG1_R10deg.ppf"

### Grande grille avec cellules de 8 Mpc 
# comment=' (projection galaxies/Poisson ngal=${NGAL}, LargeGrid-8Mpc)' 
# cmd="${EXEPATH}tpkprj -in 750,750,200,8,8,8 -out 450,450,120,8,8,8 -ic 3000,0.,0. -oc 3000,0.,0. -k 120,0.005,0.605 -neg -N $NLOOP  -ngal $NGAL -poiss ${INOUT}simlss.txt ${INOUT}pkprjSLG2.ppf"
# comment=' (projection galaxies/Poisson ngal=${NGAL}, LargeGrid-8Mpc Rotated 10 deg)' 
# cmd="${EXEPATH}tpkprj -in 750,750,200,8,8,8 -out 450,450,120,8,8,8 -ic 3000,0.,0. -oc 3000,10.,0. -k 120,0.005,0.605 -neg -N $NLOOP -ngal $NGAL -poiss ${INOUT}simlss.txt ${INOUT}pkprjSLG2_R10deg.ppf"
# comment=' (projection galaxies/Poisson ngal=${NGAL}, LargeGrid-8Mpc PkNoOsc)' 
# cmd="${EXEPATH}tpkprj -in 750,750,200,8,8,8 -out 450,450,120,8,8,8 -ic 3000,0.,0. -oc 3000,0.,0. -k 120,0.005,0.605 -neg -N $NLOOP  -ngal $NGAL -poiss -nosc ${INOUT}simlss.txt ${INOUT}pknosprjSLG2.ppf"
# comment=' (projection galaxies/Poisson ngal=${NGAL}, LargeGrid-8Mpc Rotated 10 deg PkNoOsc)' 
# cmd="${EXEPATH}tpkprj -in 750,750,200,8,8,8 -out 450,450,120,8,8,8 -ic 3000,0.,0. -oc 3000,10.,0. -k 120,0.005,0.605 -neg -N $NLOOP -ngal ${NGAL} -poiss -nosc ${INOUT}simlss.txt ${INOUT}pknosprjSLG2_R10deg.ppf"

# comment=" (projection galaxies/Poisson ngal=${NGAL}, LargeGrid+Smearing=${SIGPZ}Mpc)"
# cmd="${EXEPATH}tpkprj -in 1000,1000,250,6,6,6 -out 600,600,150,6,6,6 -ic 3000,0.,0. -oc 3000,0.,0. -k 120,0.005,0.605 -neg -N 16 -ngal $NGAL  -poiss -s ${SIGPZ} ${INOUT}simlss.txt ${INOUT}pkprjSLG1_s${SIGPZ}.ppf"
# comment=" (projection galaxies/Poisson ngal=${NGAL}, LargeGrid-Rotated10deg +Smearing=${SIGPZ}Mpc)"
# cmd="${EXEPATH}tpkprj -in 1000,1000,250,6,6,6 -out 600,600,150,6,6,6 -ic 3000,0.,0. -oc 3000,10.,0. -k 120,0.005,0.605 -neg -N 16 -ngal $NGAL -poiss -s ${SIGPZ} ${INOUT}simlss.txt ${INOUT}pkprjSLG1_R10deg_s${SIGPZ}.ppf"

##########################################
### Grande grille avec cellules de 8 Mpc 
# comment=' (projection galaxies/Poisson ngal=10, LargeGrid-8Mpc)' 
# cmd="${EXEPATH}tpkprj -in 750,750,200,8,8,8 -out 450,450,120,8,8,8 -ic 3000,0.,0. -oc 3000,0.,0. -k 120,0.005,0.605 -neg -N $NLOOP  -ngal 10 -poiss ${INOUT}simlss.txt ${INOUT}pkprjSLG2.ppf"
# comment=' (projection galaxies/Poisson ngal=10, LargeGrid-8Mpc Rotated 10 deg)' 
# cmd="${EXEPATH}tpkprj -in 750,750,200,8,8,8 -out 450,450,120,8,8,8 -ic 3000,0.,0. -oc 3000,10.,0. -k 120,0.005,0.605 -neg -N $NLOOP -ngal 10 -poiss ${INOUT}simlss.txt ${INOUT}pkprjSLG2_R10deg.ppf"
# comment=' (projection galaxies/Poisson ngal=10, LargeGrid-8Mpc PkNoOsc)' 
# cmd="${EXEPATH}tpkprj -in 750,750,200,8,8,8 -out 450,450,120,8,8,8 -ic 3000,0.,0. -oc 3000,0.,0. -k 120,0.005,0.605 -neg -N $NLOOP  -ngal 10 -poiss -nosc ${INOUT}simlss.txt ${INOUT}pknosprjSLG2.ppf"
# comment=' (projection galaxies/Poisson ngal=10, LargeGrid-8Mpc Rotated 10 deg PkNoOsc)' 
# cmd="${EXEPATH}tpkprj -in 750,750,200,8,8,8 -out 450,450,120,8,8,8 -ic 3000,0.,0. -oc 3000,10.,0. -k 120,0.005,0.605 -neg -N $NLOOP -ngal 10 -poiss -nosc ${INOUT}simlss.txt ${INOUT}pknosprjSLG2_R10deg.ppf"
# comment=" (projection galaxies/Poisson ngal=10, LargeGrid-8Mpc  +Smearing=${SIGPZ}Mpc) "
# cmd="${EXEPATH}tpkprj -in 750,750,200,8,8,8 -out 450,450,120,8,8,8 -ic 3000,0.,0. -oc 3000,0.,0. -k 120,0.005,0.605 -neg -N $NLOOP  -ngal 10 -poiss -s ${SIGPZ} ${INOUT}simlss.txt ${INOUT}pkprjSLG2_s${SIGPZ}.ppf"
# comment=" (projection galaxies/Poisson ngal=10, LargeGrid-8Mpc Rotated 10 deg +Smearing=${SIGPZ}Mpc) "
# cmd="${EXEPATH}tpkprj -in 750,750,200,8,8,8 -out 450,450,120,8,8,8 -ic 3000,0.,0. -oc 3000,10.,0. -k 120,0.005,0.605 -neg -N $NLOOP  -ngal 10 -poiss -s ${SIGPZ} ${INOUT}simlss.txt ${INOUT}pkprjSLG2_R10deg_s${SIGPZ}.ppf"
# comment=" (projection galaxies/Poisson ngal=10, LargeGrid-8Mpc  +Smearing=${SIGPZ}Mpc PkNoOsc) "
# cmd="${EXEPATH}tpkprj -in 750,750,200,8,8,8 -out 450,450,120,8,8,8 -ic 3000,0.,0. -oc 3000,0.,0. -k 120,0.005,0.605 -neg -N $NLOOP  -ngal 10 -poiss -s ${SIGPZ} -nosc ${INOUT}simlss.txt ${INOUT}pknosprjSLG2_s${SIGPZ}.ppf"
# comment=" (projection galaxies/Poisson ngal=10, LargeGrid-8Mpc Rotated 10 deg +Smearing=${SIGPZ}Mpc PkNoOsc) "
# cmd="${EXEPATH}tpkprj -in 750,750,200,8,8,8 -out 450,450,120,8,8,8 -ic 3000,0.,0. -oc 3000,10.,0. -k 120,0.005,0.605 -neg -N $NLOOP  -ngal 10 -poiss -s ${SIGPZ} -nosc ${INOUT}simlss.txt ${INOUT}pknosprjSLG2_R10deg_s${SIGPZ}.ppf"


#############################################################
# Projection en passant par les galaxies - grilles multiples 
# comment=" (projection galaxies/Poisson ngal=${NGAL}, MultipleGrids)" 
# cmd="${EXEPATH}tpkprj -in 1100,1100,300,6,6,6 -out 300,300,150,6,6,6 -ic 3000,90.,90. -moc 3000,75.,75. -moc 3000,75.,105.  -moc 3000,105.,75. -moc 3000,105.,105. -k 120,0.005,0.605 -neg -N $NLOOP -ngal 10 -poiss ${INOUT}simlss.txt ${INOUT}pkprj_4G.ppf"
# comment=' (projection galaxies/Poisson ngal=10, MultipleGrids Pk-NoOsc)' 
# cmd="${EXEPATH}tpkprj -in 1100,1100,300,6,6,6 -out 300,300,150,6,6,6 -ic 3000,90.,90. -moc 3000,75.,75. -moc 3000,75.,105.  -moc 3000,105.,75. -moc 3000,105.,105. -k 120,0.005,0.605 -neg -N 16 -ngal 10 -poiss -nosc ${INOUT}simlss.txt ${INOUT}pknosprj_4G.ppf"
# comment=" (projection galaxies/Poisson ngal=10, MultipleGrids +Smearing=${SIGPZ}Mpc)"
# cmd="${EXEPATH}tpkprj -in 1100,1100,300,6,6,6 -out 300,300,150,6,6,6 -ic 3000,90.,90. -moc 3000,75.,75. -moc 3000,75.,105.  -moc 3000,105.,75. -moc 3000,105.,105. -k 120,0.005,0.605 -neg -N 16 -ngal 10 -poiss -s ${SIGPZ} ${INOUT}simlss.txt ${INOUT}pkprj_4G_s${SIGPZ}.ppf"
# comment=" (projection galaxies/Poisson ngal=10, MultipleGrids +Smearing=${SIGPZ}Mpc PkNoOsc)"
# cmd="${EXEPATH}tpkprj -in 1100,1100,300,6,6,6 -out 300,300,150,6,6,6 -ic 3000,90.,90. -moc 3000,75.,75. -moc 3000,75.,105.  -moc 3000,105.,75. -moc 3000,105.,105. -k 120,0.005,0.605 -neg -N 16 -ngal 10 -poiss -s ${SIGPZ} -nosc ${INOUT}simlss.txt ${INOUT}pknosprj_4G_s${SIGPZ}.ppf"

#############################################################
# Projection en passant par les galaxies - grilles multiples - cellules de 8 Mpc 
# comment=" Projection galaxies/Poisson ngal=${NGAL}, MultipleGrids, 8Mpc cells " 
# cmd="${EXEPATH}tpkprj -in 800,800,250,8,8,8 -out 200,200,100,8,8,8 -ic 3000,90.,90. -moc 3200,75,75 -moc 3200,75,105 -moc 3200,105,75 -moc 3200,105,105 -k 120,0.005,0.605 -neg -N $NLOOP -ngal $NGAL -poiss -rgi ${INOUT}simlss.txt ${INOUT}pkprj_4Gb.ppf"
# comment=" Projection galaxies/Poisson ngal=${NGAL}, MultipleGrids, 8Mpc cellsPkNoOsc " 
# cmd="${EXEPATH}tpkprj -in 800,800,250,8,8,8 -out 200,200,100,8,8,8 -ic 3000,90.,90. -moc 3200,75,75 -moc 3200,75,105 -moc 3200,105,75 -moc 3200,105,105 -k 120,0.005,0.605 -neg -N $NLOOP -ngal $NGAL -poiss -nosc -rgi ${INOUT}simlss.txt ${INOUT}pknosprj_4Gb.ppf"
# comment=" Projection galaxies/Poisson ngal=${NGAL}, MultipleGrids, 8Mpc cells  +Smearing=${SIGPZ}Mpc " 
# cmd="${EXEPATH}tpkprj -in 800,800,250,8,8,8 -out 200,200,100,8,8,8 -ic 3000,90.,90. -moc 3200,75,75 -moc 3200,75,105 -moc 3200,105,75 -moc 3200,105,105 -k 120,0.005,0.605 -neg -N $NLOOP -ngal $NGAL -poiss -s ${SIGPZ} -rgi ${INOUT}simlss.txt ${INOUT}pkprj_4Gb_s${SIGPZ}.ppf"
# comment=" Projection galaxies/Poisson ngal=${NGAL}, MultipleGrids, 8Mpc cells  +Smearing=${SIGPZ}Mpc PkNoOsc" 
# cmd="${EXEPATH}tpkprj -in 800,800,250,8,8,8 -out 200,200,100,8,8,8 -ic 3000,90.,90. -moc 3200,75,75 -moc 3200,75,105 -moc 3200,105,75 -moc 3200,105,105 -k 120,0.005,0.605 -neg -N $NLOOP -ngal ${NGAL} -poiss -s ${SIGPZ} -nosc -rgi ${INOUT}simlss.txt ${INOUT}pknosprj_4Gb_s${SIGPZ}.ppf"


#############################################################
# Projection en passant par les galaxies - PETITE grille - cellules de 8 Mpc 
# comment=" Projection galaxies/Poisson ngal=${NGAL}, Small Grid, 8Mpc cells " 
# cmd="${EXEPATH}tpkprj -in 800,800,250,8,8,8 -out 200,200,100,8,8,8 -ic 3000,0.,0. -oc 3200,0.,0. -k 120,0.005,0.605 -neg -N $NLOOP -ngal $NGAL -poiss -rgi ${INOUT}simlss.txt ${INOUT}pkprj_SSG.ppf"
# comment=" Projection galaxies/Poisson ngal=${NGAL}, Small Grid, 8Mpc cellsPkNoOsc " 
# cmd="${EXEPATH}tpkprj -in 800,800,250,8,8,8 -out 200,200,100,8,8,8 -ic 3000,0.,0. -oc 3200,0,0 -k 120,0.005,0.605 -neg -N $NLOOP -ngal $NGAL -poiss -nosc -rgi ${INOUT}simlss.txt ${INOUT}pknosprj_SSG.ppf"
# comment=" Projection galaxies/Poisson ngal=${NGAL}, Small Grid, 8Mpc cells  +Smearing=${SIGPZ}Mpc " 
# cmd="${EXEPATH}tpkprj -in 800,800,250,8,8,8 -out 200,200,100,8,8,8 -ic 3000,0.,0. -oc 3200,0,0 -k 120,0.005,0.605 -neg -N $NLOOP -ngal $NGAL -poiss -s ${SIGPZ} -rgi ${INOUT}simlss.txt ${INOUT}pkprj_SSG_s${SIGPZ}.ppf"
# comment=" Projection galaxies/Poisson ngal=${NGAL}, Small Grid, 8Mpc cells  +Smearing=${SIGPZ}Mpc PkNoOsc" 
# cmd="${EXEPATH}tpkprj -in 800,800,250,8,8,8 -out 200,200,100,8,8,8 -ic 3000,0.,0. -oc 3200,0,0 -k 120,0.005,0.605 -neg -N $NLOOP -ngal ${NGAL} -poiss -s ${SIGPZ} -nosc -rgi ${INOUT}simlss.txt ${INOUT}pknosprj_SSG_s${SIGPZ}.ppf"



##### Pour debugging #### 
# Projection en passant par les galaxies - grille petite dans grande 1
# comment=' (projection galaxies/Poisson ngal=20, BigGrid to SmallGrid G1f)' 
# cmd="${EXEPATH}tpkprj -in 1200,1200,300,6,6,6 -out 450,450,225,4,4,4 -ic 3000,90.,90. -oc 3000,90.,90. -k 120,0.005,0.605 -neg -N 16 -ngal 20 -poiss ${INOUT}simlss.txt ${INOUT}pkprjgal_G1f.ppf"
# Projection en passant par les galaxies - grille petite dans grande 2
# comment=' (projection galaxies/Poisson ngal=20, BigGrid to SmallGrid G2f)' 
# cmd="${EXEPATH}tpkprj -in 1100,1100,280,6,6,6 -out 450,450,225,4,4,4 -ic 3000,90.,90. -oc 3000,75.,75. -k 120,0.005,0.605 -neg -N 16 -ngal 20 -poiss ${INOUT}simlss.txt ${INOUT}pkprjgal_G2f.ppf"
# ---- Projection directe de grille a grille 
# comment=' (projection grid to grid Big to Small 2' 
# cmd="${EXEPATH}tpkprj -in 900,900,250,6,6,6 -out 250,250,120,8,8,8 -ic 3000,90.,90. -oc 3000,75.,75. -k 120,0.005,0.605 -neg -N 25 -ngal 0. ${INOUT}simlss.txt ${INOUT}pkprj_G2.ppf"

echo '--------- (1) Executing command for ' $comment
echo $cmd 
echo '---- Starting...'
time $cmd 

echo '------------- END of do_pkprj.sh script ----------------'
qstat -j ${JOB_ID} -nenv


