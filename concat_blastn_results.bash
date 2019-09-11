#!/usr/bin/env bash



set -x
set -e  #       exit if any command fails
set -u  #       Error on usage of unset variables
set -o pipefail


#	python keeps everything in memory, whereas bash/awk is just a stream so faster?
#

rm merged.masked.tsv.gz

for p in $( find . -name \*.viral.masked.tsv.gz ); do

	echo $f

	gzcat $f | awk -v p=${p} '\
		BEGIN{
			f=p
			split(f,fparts,"/")
			fname=fparts[length(fparts)]

			s=fparts[length(fparts)]
			sub(/.viral.masked.tsv.gz/,"",s)
			sub(/.gz/,"",s)
			split(s,sparts,".")
			format=sparts[length(sparts])


		}
		{
			print fname"\t"source"\t"$0
		}' | gzip --best >> merged.masked.tsv.gz

done

