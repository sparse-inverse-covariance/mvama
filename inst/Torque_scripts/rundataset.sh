#!/bin/bash
FNAME=butyrate_residuals1.Rdata
VMEM=800GB   #should be 32*NODES
WALLTIME=24:00:00
CPCYCLES=800  #CPCYCLES*NP*JOBS (number of jobs) should be greater than the p dimension of the data
NODES=25
NP=100  #NODES*4 (we ask for 8 processors but only use 4 per node so that we have more available memory
QUEUE=NORMAL  #or  NAMD or extended
JOB01=`qsub -q $QUEUE -l nodes=$NODES:ppn=8,walltime=$WALLTIME,vmem=$VMEM -v NP=$NP,FILENAME=$FNAME,RESTART=FALSE,CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB01
JOB02=`qsub -q $QUEUE -l nodes=$NODES:ppn=8,walltime=$WALLTIME,vmem=$VMEM -W depend=afterany:$JOB01 -v NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB01"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB02
JOB03=`qsub -q $QUEUE -l nodes=$NODES:ppn=8,walltime=$WALLTIME,vmem=$VMEM -W depend=afterany:$JOB02 -v NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB02"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB03
JOB04=`qsub -q $QUEUE -l nodes=$NODES:ppn=8,walltime=$WALLTIME,vmem=$VMEM -W depend=afterany:$JOB03 -v NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB03"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB04
JOB05=`qsub -q $QUEUE -l nodes=$NODES:ppn=8,walltime=$WALLTIME,vmem=$VMEM -W depend=afterany:$JOB04 -v NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB04"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB05
JOB06=`qsub -q $QUEUE -l nodes=$NODES:ppn=8,walltime=$WALLTIME,vmem=$VMEM -W depend=afterany:$JOB05 -v NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB05"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB06
JOB07=`qsub -q $QUEUE -l nodes=$NODES:ppn=8,walltime=$WALLTIME,vmem=$VMEM -W depend=afterany:$JOB06 -v NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB06"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB07
JOB08=`qsub -q $QUEUE -l nodes=$NODES:ppn=8,walltime=$WALLTIME,vmem=$VMEM -W depend=afterany:$JOB07 -v NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB07"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB08
JOB09=`qsub -q $QUEUE -l nodes=$NODES:ppn=8,walltime=$WALLTIME,vmem=$VMEM -W depend=afterany:$JOB08 -v NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB08"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB09
JOB10=`qsub -q $QUEUE -l nodes=$NODES:ppn=8,walltime=$WALLTIME,vmem=$VMEM -W depend=afterany:$JOB09 -v NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB09"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB10
JOB11=`qsub -q $QUEUE -l nodes=$NODES:ppn=8,walltime=$WALLTIME,vmem=$VMEM -W depend=afterany:$JOB10 -v NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB10"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB11
JOB12=`qsub -q $QUEUE -l nodes=$NODES:ppn=8,walltime=$WALLTIME,vmem=$VMEM -W depend=afterany:$JOB11 -v NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB11"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB12
JOB13=`qsub -q $QUEUE -l nodes=$NODES:ppn=8,walltime=$WALLTIME,vmem=$VMEM -W depend=afterany:$JOB12 -v NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB12"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB13
JOB14=`qsub -q $QUEUE -l nodes=$NODES:ppn=8,walltime=$WALLTIME,vmem=$VMEM -W depend=afterany:$JOB13 -v NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB13"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB14
JOB15=`qsub -q $QUEUE -l nodes=$NODES:ppn=8,walltime=$WALLTIME,vmem=$VMEM -W depend=afterany:$JOB14 -v NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB14"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB15
JOB16=`qsub -q $QUEUE -l nodes=$NODES:ppn=8,walltime=$WALLTIME,vmem=$VMEM -W depend=afterany:$JOB15 -v NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB15"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB16
JOB17=`qsub -q $QUEUE -l nodes=$NODES:ppn=8,walltime=$WALLTIME,vmem=$VMEM -W depend=afterany:$JOB16 -v NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB16"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB17
JOB18=`qsub -q $QUEUE -l nodes=$NODES:ppn=8,walltime=$WALLTIME,vmem=$VMEM -W depend=afterany:$JOB17 -v NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB17"",CPEXIT=FALSE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB18

exit 0
