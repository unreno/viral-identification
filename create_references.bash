#!/usr/bin/env bash

set -x
set -e  #       exit if any command fails
set -u  #       Error on usage of unset variables
set -o pipefail




#	Update Dfam in RepeatMasker/Libraries
#	wget http://www.dfam.org/releases/Dfam_3.1/families/Dfam.hmm.gz 
#	wget http://www.dfam.org/releases/Dfam_3.1/families/Dfam.embl.gz
#	gunzip ...
#	./configure


rsync -avz --progress --include="*fna.gz" --exclude="*" rsync://ftp.ncbi.nih.gov/refseq/release/viral/ ./

zcat viral*.fna.gz > viral.genomic.fa

sed -e 's/<[^>]*>/_/g' viral.genomic.fa > viral.raw.fa

grep -c "^>" viral.raw.fa > viral.raw.count.txt

RepeatMasker -s -pa 40 viral.raw.fa > RepeatMasker.viral.raw.out

ln -s viral.raw.fa.masked viral.masked.fa


for f in viral.raw.fa viral.masked.fa ; do

	b=$( basename $f .fa )

	faSplit size -extra=50 ${f} 50 ${b}-100bp -oneFile > ${b}.faSplit.out 2> ${b}.faSplit.err
	
	bowtie2 --threads 40 -x hg38 -f -U ${b}-100bp.fa --very-sensitive --no-unal -S ${b}-100bp.hg38.e2e.sam > ${b}.bowtie2.e2e.out 2> ${b}.bowtie2.e2e.err
	
	bowtie2 --threads 40 -x hg38 -f -U ${b}-100bp.fa --very-sensitive-local --no-unal -S ${b}-100bp.hg38.loc.sam > ${b}.bowtie2.loc.out 2> ${b}.bowtie2.loc.err
	
	samtools view -c ${b}-100bp.hg38.e2e.sam > ${b}-100bp.hg38.e2e.count.txt
	
	samtools view -c ${b}-100bp.hg38.loc.sam > ${b}-100bp.hg38.loc.count.txt
	
	samtools view ${b}-100bp.hg38.e2e.sam | awk '{print $10}' | sort | uniq -c | sort -n | tail -n 20 > ${b}-100bp.hg38.e2e.common.txt
	
	samtools view ${b}-100bp.hg38.loc.sam | awk '{print $10}' | sort | uniq -c | sort -n | tail -n 20 > ${b}-100bp.hg38.loc.common.txt

	makeblastdb -in ${f} -dbtype nucl -out ${b} -title ${b} -parse_seqids

	ls -1s ${b}.n* > ${b}.blast_files.txt

	tar cvf - ${b}.n* | gzip --best > ${b}.tar.gz

done
