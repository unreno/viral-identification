#!/usr/bin/env bash


set -x
set -e  #       exit if any command fails
set -u  #       Error on usage of unset variables
set -o pipefail


echo "$0 $*"
echo "env:"
env
echo

outbucket="s3://viral-identification"


#	Expecting $1 to include the COMPLETE PATH WITH PROTOCOL
#	of gzipped fastq file. Otherwise, this script will need modified.
#	Testing bam files.


infile=$1

outbase=$1

#	Check if file exists? Or just let crash when doesn't?

if [ ${infile:0:2} == "s3" ]; then
	echo "Infile is S3"
	outbase=${outbase#s3://}
	#fasta_input_stream="aws s3 cp --debug ${infile} -"
	fasta_input_stream="aws s3 cp ${infile} -"
elif [ ${infile:0:3} == "ftp" ]; then
	echo "Infile is FTP"
	outbase=${outbase#ftp://}
	fasta_input_stream="curl -s ${infile}"
elif [ ${infile:0:7} == "http://" ]; then
	echo "Infile is HTTP"
	outbase=${outbase#http://}
	fasta_input_stream="curl -s ${infile}"
elif [ ${infile:0:8} == "https://" ]; then
	echo "Infile is HTTPS"
	outbase=${outbase#https://}
	fasta_input_stream="curl -s ${infile}"
else
	echo "Infile has an unknown protocol"
	exit 1
fi



if [ ${infile:(-9)} == '.fastq.gz' ] ; then
	#	Keep extensions so know source when aggregating
	#outbase=${outbase%.fastq.gz}   # 1000genomes and geuvadis

	#	Sadly, not all read names are of the same format
	#	Some like ...
	#		@SRR003916.100412 HWI-EAS109_303KG:5:60:511:1416 length=36			NO LANE (SINGLE, not PAIRED)
	#		>SRR003916.100412
	#	or
	#		@ERR188406.1 HWI-ST537:99:D0ACLACXX:8:1101:1207:2198/2
	#		>ERR188406.1/2
	#	or
	#		@SRR449493.136 136/1
	#		>SRR449493.136/1

	fasta_input_stream="${fasta_input_stream} | zcat | paste - - - - | cut -f 1,2 | tr '\t' '\n' | sed -E -e 's/^@/>/' -e 's/ .+\/(.)$/\/\1/' -e 's/ .*$//'"

elif [ ${infile:(-4)} == '.bam' ] ; then
	#	Keep extensions so know source when aggregating
	#outbase=${outbase%.bam}
	#fasta_input_stream="${fasta_input_stream} | samtools fasta -N -"
	# Most bam files are aligned to human. Only grab those that are unaligned?
	fasta_input_stream="${fasta_input_stream} | samtools fasta -f 4 -N -"
else
	echo "Infile is of unknown file type"
	exit 1
fi  

#outbase=${outbase%.filt}   # 1000genomes


#	There are no duplicates in the unmapped bams from 1000genomes.
#awk -F/ '{print $6}' 1000genomes.unmapped  | sort | uniq -d
#
#	s3://1000genomes/phase3/data/HG00096/alignment/HG00096.unmapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam
#	So most of the name is unnecesary, but I'm keeping it




outbase=${outbucket}/${outbase}		#.viral.csv.gz
#	s3://viral-identification/ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188040/ERR188040_1
#	s3://viral-identification/1000genomes/phase3/data/NA19238/sequence_read/SRR002989
#	s3://1000genomes/phase3/data/HG00096/alignment/HG00096.unmapped.ILLUMINA.bwa.GBR.low_coverage.20120522

#	It is always possible that this will not finish, so leave a marker and remove it at the end.
echo "Leaving marker"
date | aws s3 cp - ${outbase}.Running

outbase_dir=$( dirname ${outbase} )
#	s3://viral-identification/ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188040
#	s3://viral-identification/1000genomes/phase3/data/NA19238/sequence_read
#	s3://1000genomes/phase3/data/HG00096/alignment

outbase_file=$( basename ${outbase} )
#	ERR188040_1
#	SRR002989
#	HG00096.unmapped.ILLUMINA.bwa.GBR.low_coverage.20120522

#ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188022/ERR188022_1.fastq.gz

echo "Infile:${infile}:"
echo "Outbase:${outbase}:"
echo "Outbase dir:${outbase_dir}:"
echo "Outbase file:${outbase_file}:"

date
echo "Piping infile through several utilities and then blastn to outfile"



#	sadly, diamond will not output daa files (outfmt 100) to STDOUT??

#eval ${fasta_input_stream} | \
#	diamond blastx --block-size 1.5 --threads 1 --outfmt 6 --db $DIAMOND/viral | gzip | 
#		aws s3 cp - ${outbase_dir}/${outbase_file}.diamond.viral.csv.gz;

#	blastn only takes fasta (no fastq, no gzipped, no sam or bam or cram)
#
#	blastn's default evalue is 10
#	diamond's is 0.001
#	Setting's blastn's to 0.001
#
#eval ${fasta_input_stream} | \
#	blastn -evalue 0.001 -num_threads 1 -outfmt 6 -db viral.masked | gzip | \
#		aws s3 cp - ${outbase_dir}/${outbase_file}.blastn.viral.masked.csv.gz;

#	We found that up to 50% of these unmapped read will map to hg38 with bowtie2.
#	This is unfortunate, so I will try to add a hg38 alignment to the mix.
#	Be sure to use different output filename

#bowtie2 -f -U <( eval ${fasta_input_stream} ) -x hg38 --very-sensitive \
#	| samtools fasta -f 4 - | \
#	blastn -evalue 0.001 -num_threads 1 -outfmt 6 -db viral.masked | gzip | \
#		aws s3 cp - ${outbase_dir}/${outbase_file}.unmapped.blastn.viral.masked.csv.gz;

#	Will need the numbers to normalize later.
#	No reason I can't do multiple steps in this script.


#	hg38 reference is kinda too large to include in docker image.
#	So, copy in the hg38 reference

#	date

#	it appears that these Batch instances somehow share the Docker instance and this causes a collision
#	when run on multiple instances.
#	Try again to include in docker image
#aws s3 sync --no-sign-request --exclude \*fa s3://ngi-igenomes/igenomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/ /bowtie2/
#rename genome hg38 /bowtie2/genome.*.bt2


date
bowtie2 --xeq -f -U <( eval ${fasta_input_stream} ) -x hg38 --very-sensitive \
	| samtools fasta -f 4 - | gzip \
	| aws s3 cp - ${outbase_dir}/${outbase_file}.bowtie2-e2e.hg38.unmapped.fasta.gz

date
blastn -query <( aws s3 cp ${outbase_dir}/${outbase_file}.bowtie2-e2e.hg38.unmapped.fasta.gz - | zcat ) \
	-evalue 0.001 -outfmt 6 -db viral.masked | gzip \
	| aws s3 cp - ${outbase_dir}/${outbase_file}.bowtie2-e2e.hg38.unmapped.blastn.viral.masked.csv.gz

date
aws s3 cp ${outbase_dir}/${outbase_file}.bowtie2-e2e.hg38.unmapped.fasta.gz - | zcat | paste - - | wc -l \
	| aws s3 cp - ${outbase_dir}/${outbase_file}.bowtie2-e2e.hg38.unmapped.fasta.gz.read_count.txt

sample=${outbase_file%%.*}

for index in /bowtie2/NC_00*.rev.2.bt2 ; do
	date
	index=$( basename ${index} .rev.2.bt2 )
	bowtie2 --xeq -x ${index} --very-sensitive-local --no-unal \
		--rg-id ${sample}.loc --rg "SM:${sample}" \
		-U <( aws s3 cp ${outbase_dir}/${outbase_file}.bowtie2-e2e.hg38.unmapped.fasta.gz - ) \
		| samtools sort - \
		| aws s3 cp - ${outbase_dir}/${outbase_file}.bowtie2-e2e.hg38.unmapped.bowtie2-loc.${index}.bam
	aws s3 cp ${outbase_dir}/${outbase_file}.bowtie2-e2e.hg38.unmapped.bowtie2-loc.${index}.bam - \
		| samtools view -c - \
		| aws s3 cp - ${outbase_dir}/${outbase_file}.bowtie2-e2e.hg38.unmapped.bowtie2-loc.${index}.bam.read_count.txt
done



#	eval ${fasta_input_stream} | tee \
#		>( blastn -num_threads 1 -outfmt 6 -db viral.masked | gzip --best | \
#			aws s3 cp - ${outbase_dir}/${outbase_file}.viral.masked.tsv.gz;
#			touch viral.masked.done ) \
#		>( blastn -num_threads 1 -outfmt 6 -db viral.raw | gzip --best | \
#			aws s3 cp - ${outbase_dir}/${outbase_file}.viral.raw.tsv.gz; \
#			touch viral.raw.done ) > /dev/null
#	
#	until [ -e viral.masked.done ] && [ -e viral.raw.done ] ; do
#		echo "Waiting"
#		sleep 1
#	done
#	
#	echo "Complete"
#	date
#	
#	echo "Removing markers"
#	rm viral.masked.done viral.raw.done
#	

aws s3 rm ${outbase}.Running
echo "Complete"

