#!/bin/bash

DATE=`date +%Y%m%d`
echo "$DATE"

CH=$1
FUNC=$2
MASS=$3
NTOYS=$4
NSIGFAC=$5
JOBID=$6

echo '------------------------------------------------------'
echo 'Arg: CH                           '$CH
echo 'Arg: FUNC                         '$FUNC
echo 'Arg: MASS                         '$MASS
echo 'Arg: NTOYS                        '$NTOYS
echo 'Arg: NSIGFAC                      '$NSIGFAC
echo 'Arg: JOBID                        '$JOBID
echo '------------------------------------------------------'
echo 'PBS: qsub is running on           '$PBS_O_HOST
echo 'PBS: originating queue is         '$PBS_O_QUEUE
echo 'PBS: executing queue is           '$PBS_QUEUE
echo 'PBS: working directory is         '$PBS_O_WORKDIR
echo 'PBS: execution mode is            '$PBS_ENVIRONMENT
echo 'PBS: job identifier is            '$PBS_JOBID
echo 'PBS: job name is                  '$PBS_JOBNAME
echo 'PBS: current home directory is    '$PBS_O_HOME
echo '------------------------------------------------------'

#-- PBS parameters
#PBS -l walltime=48:00:00,cput=48:00:00
#PBS -V
#PBS -k o
#PBS -j oe

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/work/j/jschulte/bias/CMSSW_10_2_13/src/
eval `scramv1 runtime -sh`
which cmsRun
which root
echo ""

cd $PBS_O_WORKDIR

echo "python runBiasTestSigBkg.py --channel $CH -f $FUNC --mass $MASS --nToys $NTOYS --fix 1 --nSigFac $NSIGFAC --jobId $JOBID"
python runBiasTestSigBkg.py --channel $CH -f $FUNC --mass $MASS --nToys $NTOYS --fix 1 --nSigFac $NSIGFAC --jobId $JOBID

echo 'script done'
