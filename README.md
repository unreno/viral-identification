#	Viral Identification


Thorough, unique and unambiguous identification of viral content in a genomic sample.




#	Reference Preparation



##	Viral Sequence Acquisition From NCBI RefSeq



```BASH
rsync -avz --progress --include="*fna.gz" --exclude="*" rsync://ftp.ncbi.nih.gov/refseq/release/viral/ ./

zcat viral*.fna.gz > viral.genomic.fa
```




##	Viral Sequence Extraction From NCBI's nt Database


The complete nt database is quite large. ~100GB. Extracting just the viruses.


```BASH
blastn -version
blastn: 2.9.0+
 Package: blast 2.9.0, build Mar 11 2019 15:20:05

update_blastdb.pl -version
/home/jake/.local/bin/update_blastdb.pl version 581818

blastdbcmd -version
blastdbcmd: 2.9.0+
 Package: blast 2.9.0, build Mar 11 2019 15:20:05

update_blastdb.pl nt

blastdbcmd -db nt -entry all -outfmt "%K" | sort | uniq -c


blastdbcmd -db nt -entry all -outfmt "%K,%i" | awk -F, '( $1 == 'Viruses' ){ print $2 }' > viruses.seqidlist

blastdbcmd -db nt -entry_batch viruses.seqidlist > nt.viruses.fa
```




##	Reference Cleanup




















#	Script Preparation




When using AWS Batch, the job command cannot include the use of pipes as they will be converted to strings.

This basically means that your job will probably need a script.





#	Docker Image Creation and Storage









#	AWS Preparation






#	Job Submission






#	Job Debugging






