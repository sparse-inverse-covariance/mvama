#!/bin/sh
#PBS -j oe 
cd ${PBS_O_WORKDIR:-.}
#cd $PBS_O_WORKDIR
module load openmpi
module load R/2.11.1
PROCS=8
CPE=""
CPC=""
if [ "$NP" != "" ]
then
	PROCS=$NP
fi
if [ "$CPEXIT" == "TRUE" ]
then
	CPE="exitoncheckpoint="$CPEXIT
fi
if [ "$CPCYCLES" != "" ]
then
	CPC="checkpointcycles="$CPCYCLES
fi
echo "Processes = "$PROCS
if [ "$FILENAME" != "" ]
then

	mkdir "./progress_"$FILENAME"_"$PBS_JOBID

	if [ "$RESTART" == "TRUE" ]
	then
		if [ "$RESTARTDIR" != "" ]
		then
			mpirun --bynode -n $PROCS R  CMD BATCH --no-restore --no-save --slave "--args filename=\""$FILENAME"\" outdir=\"./progress_"$FILENAME"_"$PBS_JOBID"\" "$CPE" "$CPC" restart="$RESTART" restartdir=\""$RESTARTDIR"\"" regression.s "./progress_"$FILENAME"_"$PBS_JOBID"/output"$PBS_JOBID".Rout"
		else
			echo "RESTART=TRUE must be accompanied by RESTARTDIR=...."
		fi
	else
		mpirun --bynode -n $PROCS R  CMD BATCH --no-restore --no-save --slave "--args filename=\""$FILENAME"\" outdir=\"./progress_"$FILENAME"_"$PBS_JOBID"\" "$CPE" "$CPC"" regression.s "./progress_"$FILENAME"_"$PBS_JOBID"/output"$PBS_JOBID".Rout"
	fi
else
	echo "Usage regression.q FILENAME=.... [RESTART=TRUE|FALSE RESTARTDIR=.....] [CPCYCLES=N] [CPEXIT=TRUE|FALSE]"
fi

