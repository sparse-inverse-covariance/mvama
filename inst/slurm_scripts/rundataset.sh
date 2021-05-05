#!/bin/bash
FNAME=/home/dun280/butyrate/butyrate_residuals1.Rdata
VMEM=32768   # 32GB 
WALLTIME=24:00:00
CPCYCLES=800  #CPCYCLES*NP*JOBS (number of jobs) should be greater than the p dimension of the data
NODES=25
NP=100  #NODES*4 (we ask for 8 processors but only use 4 per node so that we have more available memory
#QUEUE=defq  #  Don't specify partition at all.  Jobs are automatically allocated to the right one.
JOB01=`sbatch --parsable  --nodes=$NODES --time=$WALLTIME --ntasks-per-node=8 --mem-per-cpu=$VMEM  --export=NP=$NP,FILENAME=$FNAME,RESTART=FALSE,CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB01
JOB02=`sbatch --parsable  --nodes=$NODES --time=$WALLTIME --ntasks-per-node=8   --mem-per-cpu=$VMEM --dependency=afterany:$JOB01 --export=NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB01"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB02
JOB03=`sbatch --parsable  --nodes=$NODES --time=$WALLTIME --ntasks-per-node=8   --mem-per-cpu=$VMEM --dependency=afterany:$JOB02 --export=NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB02"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB03
JOB04=`sbatch --parsable  --nodes=$NODES --time=$WALLTIME --ntasks-per-node=8   --mem-per-cpu=$VMEM --dependency=afterany:$JOB03 --export=NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB03"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB04
JOB05=`sbatch --parsable  --nodes=$NODES --time=$WALLTIME --ntasks-per-node=8   --mem-per-cpu=$VMEM --dependency=afterany:$JOB04 --export=NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB04"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB05
JOB06=`sbatch --parsable  --nodes=$NODES --time=$WALLTIME --ntasks-per-node=8   --mem-per-cpu=$VMEM --dependency=afterany:$JOB05 --export=NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB05"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB06
JOB07=`sbatch --parsable  --nodes=$NODES --time=$WALLTIME --ntasks-per-node=8   --mem-per-cpu=$VMEM --dependency=afterany:$JOB06 --export=NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB06"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB07
JOB08=`sbatch --parsable  --nodes=$NODES --time=$WALLTIME --ntasks-per-node=8   --mem-per-cpu=$VMEM --dependency=afterany:$JOB07 --export=NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB07"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB08
JOB09=`sbatch --parsable  --nodes=$NODES --time=$WALLTIME --ntasks-per-node=8   --mem-per-cpu=$VMEM --dependency=afterany:$JOB08 --export=NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB08"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB09
JOB10=`sbatch --parsable  --nodes=$NODES --time=$WALLTIME --ntasks-per-node=8   --mem-per-cpu=$VMEM --dependency=afterany:$JOB09 --export=NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB09"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB10
JOB11=`sbatch --parsable  --nodes=$NODES --time=$WALLTIME --ntasks-per-node=8   --mem-per-cpu=$VMEM --dependency=afterany:$JOB10 --export=NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB10"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB11
JOB12=`sbatch --parsable  --nodes=$NODES --time=$WALLTIME --ntasks-per-node=8   --mem-per-cpu=$VMEM --dependency=afterany:$JOB11 --export=NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB11"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB12
JOB13=`sbatch --parsable  --nodes=$NODES --time=$WALLTIME --ntasks-per-node=8   --mem-per-cpu=$VMEM --dependency=afterany:$JOB12 --export=NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB12"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB13
JOB14=`sbatch --parsable  --nodes=$NODES --time=$WALLTIME --ntasks-per-node=8   --mem-per-cpu=$VMEM --dependency=afterany:$JOB13 --export=NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB13"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB14
JOB15=`sbatch --parsable  --nodes=$NODES --time=$WALLTIME --ntasks-per-node=8   --mem-per-cpu=$VMEM --dependency=afterany:$JOB14 --export=NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB14"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB15
JOB16=`sbatch --parsable  --nodes=$NODES --time=$WALLTIME --ntasks-per-node=8   --mem-per-cpu=$VMEM --dependency=afterany:$JOB15 --export=NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB15"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB16
JOB17=`sbatch --parsable  --nodes=$NODES --time=$WALLTIME --ntasks-per-node=8   --mem-per-cpu=$VMEM --dependency=afterany:$JOB16 --export=NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB16"",CPEXIT=TRUE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB17
JOB18=`sbatch --parsable  --nodes=$NODES --time=$WALLTIME --ntasks-per-node=8   --mem-per-cpu=$VMEM --dependency=afterany:$JOB17 --export=NP=$NP,FILENAME=$FNAME,RESTART=TRUE,RESTARTDIR="./progress_"$FNAME"_"$JOB17"",CPEXIT=FALSE,CPCYCLES=$CPCYCLES regression.q`
echo $JOB18

exit 0
