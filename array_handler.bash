#!/usr/bin/env bash

list=${1:-1000genomes}
page=${2:-1}

page_length=10000

echo "0:command:$0"
echo "1:list:$1"
echo "2:page:$2"
echo "list:$list"
echo "page:$page"
echo "AWS_BATCH_JOB_ARRAY_INDEX:$AWS_BATCH_JOB_ARRAY_INDEX"

item=$(( ($page - 1) * $page_length + $AWS_BATCH_JOB_ARRAY_INDEX + 1 ))

echo "item#:$item"

if [ -f ${list} ] ; then

	echo "List file :${list}: found"
	echo "Picking out the :${item}th: line"

	url=$( tail -n +${item} ${list} | head -1 )

	if [ -n ${url} ] ; then
		echo "Found :${url}:"
		#viral_identification.bash ${url}
	else
		echo "No line match these specs"
	fi

else

	echo "List file :${list}: not found"

fi

