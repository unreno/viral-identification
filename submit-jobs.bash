#!/usr/bin/env bash


#set -x
set -e  #       exit if any command fails
set -u  #       Error on usage of unset variables
set -o pipefail


for inurl in $( tail list ) ; do


	jobname=${inurl//\//_}
	jobname=${jobname//./_}
	jobname=${jobname//:/_}
	echo $jobname

	aws batch submit-job --job-name ${jobname} --job-definition myJobDefinition --job-queue myJobQueue --container-overrides '{ "command": ["viral_identification.bash","'${inurl}'"]}'

	#	Test
	#viral_identification.bash ${inurl}

done

