#	Viral Identification


Thorough, unique and unambiguous identification of viral content in a genomic sample.


#	AWS FYI

This analysis is designed to use AWS resources from the command line.
`aws` usually requires your credentials in `~/.aws/config` and `~/.aws/credentials`.
These can be set up using `aws configure` or simply creating and editing those files.

If you have multiple sets of credentials, as I do, you will need to add `--profile` and/or `--region` options to each of these `aws` commands, otherwise they will use the default.

To setup multiple credentials, you could use `aws configure` and then edit the files or simply just edit them.

```BASH
cat ~/.aws/config

[default]
output = json
region = us-east-1
cli_history = enable
[profile oneofmyprofiles]
output = json
region = us-east-1
[profile anotherprofile]
output = json
region = us-west-2

cat ~/.aws/credentials

[default]
aws_access_key_id = ABCDEF
aws_secret_access_key = GHIJKL
[oneofmyprofiles]
aws_access_key_id = ABCDEF
aws_secret_access_key = GHIJKL
[anotherprofile]
aws_access_key_id = MNOPQR
aws_secret_access_key = STUVWX
```







#	Reference Preparation


What is the most appropriate viral reference?

RefSeq only has about 12,000 sequences.



##	Viral Sequence Acquisition From NCBI RefSeq


On August 28, 2019 ...


```BASH
rsync -avz --progress --include="*fna.gz" --exclude="*" rsync://ftp.ncbi.nih.gov/refseq/release/viral/ ./

ls -1s viral*genomic.fna.gz
63744 viral.1.1.genomic.fna.gz
30208 viral.2.1.genomic.fna.gz

zcat viral*.fna.gz > viral.genomic.fa
```



##	Reference Cleanup


Some of these entries have sequence names which are invalid to some software.


The most common issue is including a greater-than symbol after the first character.


```BASH
sed -e 's/<[^>]*>/_/g' viral.genomic.fa > viral.raw.fa

diff viral.genomic.fa viral.raw.fa
3440738c3440738
< >NC_042059.1 Halobacterium phage phiH T4, T4', and T<down>LX1</down> genes, complete sequence; and orf75 (T<down>LX3</down>) gene, complete cds
---
> >NC_042059.1 Halobacterium phage phiH T4, T4', and T_LX1_ genes, complete sequence; and orf75 (T_LX3_) gene, complete cds

```


##	Human Homology


Many viruses have some sequence similarity to parts of the human genome.
In order to minimize any mis-alignment of human sequence to the viral reference, I will mask it.
But first, I'll check the homology by chopping the viral reference into smaller sequences and attempt to align them to a human reference.

For this example, let's chop the reference into 100bp reads.
`faSplit` gives us the opportunity take a 50bp sequence and then add the next 50bp as "extra".
This overlap will result in the majority (not the first or last 50bp of each sequence) of being checked twice.





```BASH

faSplit size -extra=50 viral.raw.fa 50 viral.raw-100bp -oneFile
6256732 pieces of 6257355 written

grep -c "^>" viral.raw-100bp.fa
6256732

bowtie2 --version

/usr/local/bin/bowtie2-align-s version 2.3.4.1
64-bit
Built on system76-server
Tue Apr 17 15:46:04 MDT 2018
Compiler: gcc version 5.4.0 20160609 (Ubuntu 5.4.0-6ubuntu1~16.04.9)
Options: -O3 -m64 -msse2 -funroll-loops -g3 -std=c++98 -DPOPCNT_CAPABILITY -DWITH_TBB -DNO_SPINLOCK -DWITH_QUEUELOCK=1
Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}


bowtie2 -x hg38 -f -U viral.raw-100bp.fa --very-sensitive --no-unal -S viral.raw-100bp.hg38.e2e.sam
6256732 reads; of these:
  6256732 (100.00%) were unpaired; of these:
    6256188 (99.99%) aligned 0 times
    73 (0.00%) aligned exactly 1 time
    471 (0.01%) aligned >1 times
0.01% overall alignment rate


bowtie2 -x hg38 -f -U viral.raw-100bp.fa --very-sensitive-local --no-unal -S viral.raw-100bp.hg38.loc.sam
6256732 reads; of these:
  6256732 (100.00%) were unpaired; of these:
    6246279 (99.83%) aligned 0 times
    3311 (0.05%) aligned exactly 1 time
    7142 (0.11%) aligned >1 times
0.17% overall alignment rate


samtools --version
samtools 1.8
Using htslib 1.8
Copyright (C) 2018 Genome Research Ltd.

samtools view -c viral.raw-100bp.hg38.e2e.sam
544

samtools view -c viral.raw-100bp.hg38.loc.sam
10453

samtools view viral.raw-100bp.hg38.e2e.sam | awk '{print $10}' | sort | uniq -c | sort -n | tail
      2 TGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTG
      3 CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT
      3 CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA
      3 GTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA
      4 GGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTT
      6 GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
      6 GTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT
      7 TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC
     10 ACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC
     10 CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA

samtools view viral.raw-100bp.hg38.loc.sam | awk '{print $10}' | sort | uniq -c | sort -n | tail
      3 CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT
      3 CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAA
      3 GTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTA
      3 TTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTAT
      4 GGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTT
      6 GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
      6 GTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGTGT
      7 TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC
     10 ACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC
     10 CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTA
```






How to properly analyze these masked sequences?
Does masking matter in post alignment analysis?



###	RepeatMasking































https://github.com/yjx1217/RMRB
RepBaseRepeatMaskerEdition-20170127.tar.gz
Can't use. Only includes RMRBSeqs.embl and RepeatMasker seems to be expecting RMRBMeta.embl


Trying to get the latest RepBase

wget http://www.dfam.org/releases/Dfam_3.1/families/Dfam.embl.gz
wget http://www.dfam.org/releases/Dfam_3.1/families/Dfam.hmm.gz
cd RepeatMasker/Libraries; gunzip Dfam.embl.gz; gunzip Dfam.hmm.gz ; cd .. ; ./configure


Install RepeatMasker Libraries
RepeatMasker now comes with the open Dfam database and will work out-of-the box with this library. However it is advised that you also obtain a license for the RepBase RepeatMasker Edition to supplement these sequences. To obtain a license and download the library go to http://www.girinst.org. Once you have obtaianed the library ( current version is RepBaseRepeatMaskerEdition-20181026.tar.gz ) file from GIRI unpack it in the RepeatMasker directory and it will automatically place the contents in the correct subdirectories. Then complete the installation by running ( or rerunning ) the configure program to prepare the libraries for RepeatMasker use.
cp RepBaseRepeatMaskerEdition-########.tar.gz /usr/local/RepeatMasker/
cd /usr/local/RepeatMasker
gunzip RepBaseRepeatMaskerEdition-########.tar.gz
tar xvf RepBaseRepeatMaskerEdition-########.tar
rm RepBaseRepeatMaskerEdition-########.tar


Check for Dfam Updates ( optional )
Download the Dfam.hmm.gz library from http://www.dfam.org and save it to the RepeatMasker/Libraries directory. Uncompress the file before using RepeatMasker. The RepeatMasker distribution contains the Dfam 2.0 library.




```BASH



RepeatMasker -s -pa 40 viral.raw.fa > RepeatMasker.out





head -3 RepeatMasker.out
RepeatMasker version open-4.0.9
Search Engine: NCBI/RMBLAST [ 2.9.0+ ]
Master RepeatMasker Database: /home/jake/.local/RepeatMasker-open-4-0-9-p2/Libraries/RepeatMaskerLib.embl ( Complete Database: CONS-Dfam_3.0 )






cat RepeatMasker.out | sed -E 's/ in batch [[:digit:]]+ of [[:digit:]]+//' | sort | uniq
cat RepeatMasker.out | sed 's/ in batch [[:digit:]]\+ of [[:digit:]]\+//' | sort | uniq

Checking for E. coli insertion elements
identifying ancient repeats
identifying full-length ALUs
identifying full-length interspersed repeats
identifying long interspersed repeats
identifying most interspersed repeats
identifying remaining ALUs
identifying retrovirus-like sequences
identifying Simple Repeats

The following E coli IS elements could not be confidently clipped out:
  IS1#ARTEFACT in NC_005856.1frag-1: 22650 - 22848
  IS1#ARTEFACT in NC_005856.1frag-1: 23125 - 23314
  IS1#ARTEFACT in NC_022749.1: 35265 - 35605
  IS1#ARTEFACT in NC_022749.1: 35731 - 35929
  IS1#ARTEFACT in NC_042128.1frag-2: 50923 - 51121
  IS1#ARTEFACT in NC_042128.1frag-2: 51247 - 51360
  IS1#ARTEFACT in NC_042128.1frag-2: 51375 - 51570
  IS5#ARTEFACT in NC_005856.1frag-1: 7667 - 8702
  IS5#ARTEFACT in NC_042128.1frag-2: 50115 - 50731
  IS5#ARTEFACT in NC_042128.1frag-2: 57463 - 58079



faSplit size -extra=50 viral.raw.fa.masked 50 viral.masked-100bp -oneFile
6224051 pieces of 6257355 written

grep -c "^>" viral.masked-100bp.fa
6224051


bowtie2 -x hg38 -f -U viral.masked-100bp.fa --very-sensitive --no-unal -S viral.masked-100bp.hg38.e2e.sam
6224051 reads; of these:
  6224051 (100.00%) were unpaired; of these:
    6223929 (100.00%) aligned 0 times
    62 (0.00%) aligned exactly 1 time
    60 (0.00%) aligned >1 times
0.00% overall alignment rate

bowtie2 -x hg38 -f -U viral.masked-100bp.fa --very-sensitive-local --no-unal -S viral.masked-100bp.hg38.loc.sam
6224051 reads; of these:
  6224051 (100.00%) were unpaired; of these:
    6221488 (99.96%) aligned 0 times
    1494 (0.02%) aligned exactly 1 time
    1069 (0.02%) aligned >1 times
0.04% overall alignment rate


samtools view -c viral.masked-100bp.hg38.e2e.sam
122

samtools view -c viral.masked-100bp.hg38.loc.sam
2563


samtools view viral.masked-100bp.hg38.e2e.sam | awk '{print $10}' | sort | uniq -c | sort -n | tail
      1 TGGATCGGGAAACTGGTTCTATCAAGGTTGTAGTGTCCAAATAGTGTATTTTGTAGAATTCCAGTAAAGATGGCAAAGAATCAAACACCTGGTCACCTAT
      1 TGGCACAACCATCTGAATCCAGAAGTGAAGAAAACCTCCTGGACAGAAGAGGAAGATAGAATTATTTACCAGGCACACAAGAGACTGGGAAACAGATGGG
      1 TGGCGGAAGGACCCTGAGGAGCGGCCCACTTTTGAGTACCTGCAGGCCTTCCTGGAGGACTACTTCACCTCGACAGAGCCCCCAGTACCAGCCTGGAGAG
      1 TGGCTCGCTGCTCCTGGGAAGTCCCCGGGCTTCGGGTCACAGCCCGTGCAGCTGCCACTATCTCACACTTGCATGCCAGGTGGTCCTCCAGCGTCACCGT
      1 TTCCATCCATCAGCGCAAAGTAGGTGATTTTGAGGCCCAACATGCTTGACTCTGCCCAAAAGTCACCTTCCTCAGCAGGATGCAGCCTATTACACTCAGC
      1 TTGCAACAGCCGGAGCAGCGCTGCACCTCCACGCAGGGCGGCCACACCAGGAAGTTGGCATTGGTGCGGTCGATGAGGCGCCGGGAGATCTCGAACACCT
      1 TTTTAACCAGTGAAATTGACCTGCCCGTGAAGAGGCGGGCATGACACAGCAAGACGAGAAGACCCTATGGAGCTTTAATTTATTAATGCAAACAGTACCT
      2 CCAACCCTAACCCTAACCCTAGCTCTAAGCCTAACCCCAACCCTAACCCTAACCCTAGCTCTAAGCCTAACCCCAACCCTAACCCTAACCCTAGCTCTAA
      2 TAACCCTAACCCTAACCCTAAGTCTAACCCTAACCCTAAGTCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTGACCCTAACCCTAGCTCTAAC
      2 TAGGGTTAGACTTAGGGTTAGGGTTAGACTTAGGGTTAGGGTTAGGGTTAGACCTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGG

samtools view viral.masked-100bp.hg38.loc.sam | awk '{print $10}' | sort | uniq -c | sort -n | tail
      2 GTTAGAGCTAGGGTTAGGGTCAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGACTTAGGGTTAGGGTTAGACTTAGGGTTAGGGTTAGGGTTA
      2 GTTTAAAATTTTGTCCAACGTTCCTTTAGGTACACCAGACTTCTCGGATAATTGTTTTGAAGTAAATCCTTTTTCTTTTTTTAATTTTTCTATTATTTCC
      2 TAAGTCCTCATCCAATTTATAAAAAATATTTGAATTGTTTAATTTAATTATTAATATACTTTTTTCTTTTTCTGTTAAATTAGAATTTTCTAAAATTGTA
      2 TAATATCTCTTAATAATTCTTTAGAGCCAATATCTACAATATTTTTATTATTTATTTTTAAATATACTCTTTCTTCACCTATAAATATACTATTTCTATT
      2 TAGCTCTAACCCTAGCCCTAACCCTAACCCTAGCTCTAACCCTAGCCCTAACCCTAACCCTAGCTCTAACCCTAGCCCTAACCCTAACCCTAGCTCTAAC
      2 TAGGGTTAGACTTAGGGTTAGGGTTAGACTTAGGGTTAGGGTTAGGGTTAGACCTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGG
      2 TCATCAGACTTCTCCGTTTTTTCTTCTTTGCTGTCTTCTTCTTCTTTCTTGTCTTTCTTGTCTTTCTTACCTTTTTTGTTCTCGTCTTCTTCACTGTCTT
      2 TCTCCTCCCCTCTCTCCTCCCCTCTCTCCTCCCCTCTCTCCTCCCCTCTCTCCTCCCCTCTCTCCTCCCCTCTCTCCTCCCCTCTCTCCTCCCCTCTCTC
      2 TTATTAGATATAGCTTCGTGGATATCTATTGAATTTTATAATAAATGTTATAATATAATAGAAAAATATTATATAGACATATTCACTAAAGAGTATAAAG
      2 TTTTATTATACAAAGCTAGTAAAAAATATAAATAAGGGGGTATTGGTTTAAAAATTTTTTTTTTATATTTTTTTAAACTCATTCATTTCTTCAAACACCT
```










updated to Dfam 3.1








```BASH
RepeatMasker -s -pa 40 viral.raw.fa > RepeatMasker.out

head -3 RepeatMasker.out
RepeatMasker version open-4.0.9
Search Engine: NCBI/RMBLAST [ 2.9.0+ ]
Master RepeatMasker Database: /Users/jakewendt/.local/RepeatMasker-open-4-0-9-p2/Libraries/RepeatMaskerLib.embl ( Complete Database: CONS-Dfam_3.1 )

cat RepeatMasker.out | sed -E 's/ in batch [[:digit:]]+ of [[:digit:]]+//' | sort | uniq ... moderate

Checking for E. coli insertion elements
identifying ancient repeats
identifying full-length ALUs
identifying full-length interspersed repeats
identifying long interspersed repeats
identifying most interspersed repeats
identifying remaining ALUs
identifying retrovirus-like sequences
identifying Simple Repeats

The following E coli IS elements could not be confidently clipped out:
  IS1#ARTEFACT in NC_005856.1frag-1: 22650 - 22848
  IS1#ARTEFACT in NC_005856.1frag-1: 23125 - 23314
  IS1#ARTEFACT in NC_022749.1: 35265 - 35605
  IS1#ARTEFACT in NC_022749.1: 35731 - 35929
  IS1#ARTEFACT in NC_042128.1frag-2: 50923 - 51121
  IS1#ARTEFACT in NC_042128.1frag-2: 51247 - 51360
  IS1#ARTEFACT in NC_042128.1frag-2: 51375 - 51570
  IS5#ARTEFACT in NC_005856.1frag-1: 7667 - 8702
  IS5#ARTEFACT in NC_042128.1frag-2: 50115 - 50731
  IS5#ARTEFACT in NC_042128.1frag-2: 57463 - 58079

faSplit size -extra=50 viral.raw.fa.masked 50 viral.masked-100bp -oneFile
6224051 pieces of 6257355 written

grep -c "^>" viral.masked-100bp.fa
6224051

bowtie2 -x hg38 -f -U viral.masked-100bp.fa --very-sensitive --no-unal -S viral.masked-100bp.hg38.e2e.sam
6224051 reads; of these:
  6224051 (100.00%) were unpaired; of these:
    6223929 (100.00%) aligned 0 times
    62 (0.00%) aligned exactly 1 time
    60 (0.00%) aligned >1 times
0.00% overall alignment rate

bowtie2 -x hg38 -f -U viral.masked-100bp.fa --very-sensitive-local --no-unal -S viral.masked-100bp.hg38.loc.sam
6224051 reads; of these:
  6224051 (100.00%) were unpaired; of these:
    6221488 (99.96%) aligned 0 times
    1494 (0.02%) aligned exactly 1 time
    1069 (0.02%) aligned >1 times
0.04% overall alignment rate

samtools view -c viral.masked-100bp.hg38.e2e.sam
122

samtools view -c viral.masked-100bp.hg38.loc.sam
2563
```

No difference in the numbers comparing Dfam 3.0 and Dfam 3.1.


Perhaps the RepBaseRepeatMaskerEdition-20181026.tar.gz would make a difference?
Cost thousands of dollars!



I've run RepeatModeler on hg38 to make my own mask.
I've run RepeatMasker multiple times.
Both mask a few thousand more base pairs, but neither reduce the number of aligned 100bp snippets.













##	Creating blastn Reference


```BASH

makeblastdb -in viral.raw.fa -dbtype nucl -out viral.raw -title viral.raw -parse_seqids

ls -1s viral.raw.n*
 1444 viral.raw.nhr
  144 viral.raw.nin
   48 viral.raw.nog
  384 viral.raw.nsd
   12 viral.raw.nsi
82176 viral.raw.nsq

tar cvf - viral.raw.n* | gzip > viral.raw.tar.gz





Masked or MultiMasked?





makeblastdb -in viral.raw.fa.masked -dbtype nucl -out viral.masked -title viral.masked -parse_seqids

ls -1s viral.masked.n*
 1444 viral.masked.nhr
  144 viral.masked.nin
   48 viral.masked.nog
  384 viral.masked.nsd
   12 viral.masked.nsi
76944 viral.masked.nsq

tar cvf - viral.masked.n* | gzip > viral.masked.tar.gz
```








#	Script Preparation


When using AWS Batch, the job command cannot include the use of pipes as they will be converted to strings.

This basically means that your job will probably need a script.

But don't need a complex script. Processing is a series of pipes from S3(or ftp) to S3.





#	Docker Image Creation and Storage


Create a Dockerfile defining the contents of your environment.

Create a Docker image containing anything needed by any of your job scripts.



```BASH
docker build -t viral_identification .
```


Manually inspect your new creation.
Given that this docker container is not meant to be interactive, you MAY need to redefine the entrypoint.

```BASH
docker run --rm -it viral_identification
```

If you added an ENTRYPOINT to your Dockerfile, you'll need to override it with something like ...

```BASH
docker run --rm --entrypoint "/bin/bash" -it viral_identification
```

I do not use entrypoints as I feel that they limit what the container can do.
It can't be overridden in the AWS Batch world.
At least not as far as I have found anyway.

The fetch_and_run example is built on this which is fine, but it effectively ignores the first part of the passed command.
This would work well if you can't or don't want to edit the docker image.
It is stuck to getting its job script from S3.





##	Create an ECR Repository

Strongly advised to use the Amazon Elastic Container Registry (ECR) to store your docker image.

I imagine that you can store it elsewhere but that is outside the focus of this project.


```BASH
aws ecr create-repository --repository-name viral_identification
{
    "repository": {
        "repositoryArn": "arn:aws:ecr:us-east-1:XXXXXXXXXX:repository/viral_identification",
        "registryId": "XXXXXXXXXX",
        "repositoryName": "viral_identification",
        "repositoryUri": "XXXXXXXXXX.dkr.ecr.us-east-1.amazonaws.com/viral_identification",
        "createdAt": 1564519088.0,
        "imageTagMutability": "MUTABLE"
    }
}
```


##	Push Docker Image


You'll need the repository URI.
If you weren't paying attention when you created it, you can get it with the following command.

```BASH
aws ecr describe-repositories | jq -r '.repositories | map(select( .repositoryName == "viral_identification" ).repositoryUri )[]'

XXXXXXXXXX.dkr.ecr.us-east-1.amazonaws.com/viral_identification
```


Login, tag the image and push it to the repository.
Be sure to use the added --no-include-email option as without it the output will include deprecated "-e none"

For privacy in this README, I used a command to set my Account ID / ECR repository.


```BASH
ecr_repo=$( aws ecr describe-repositories | jq -r '.repositories | map(select( .repositoryName == "viral_identification" ).repositoryUri )[]' )

aws ecr get-login --no-include-email --region us-east-1 | bash

docker tag viral_identification:latest ${ecr_repo}:latest

docker push ${ecr_repo}:latest
```






#	AWS Preparation



Your AWS Account ID ...
```BASH
aws sts get-caller-identity | jq -r '.Account'
```







##	Create an EC2 KeyPair

This is kinda optional, but again strongly advised.

This keypair will allow you connect to instances, predominantly for debugging.

This MUST be done BEFORE the cloudformation, as the keypair name is referenced in the template.


```BASH
private_key=$( aws ec2 create-key-pair --key-name batchKeyPair )

echo $private_key | jq -r '.KeyMaterial' > ~/.ssh/batchKeyPair

chmod 400 ~/.ssh/batchKeyPair
```









##	Use CloudFormation to Build All Your Resources



Previously, IAM Roles, Compute Environments, Queues, etc. were created with bash scripts.
They were very busy, but they checked things, were rerunnable and dealt with existing resources and errors.
CloudFormation doesn't seem to do that. If a resource already exists, it crashes.

I recently learned of an AWS resource called CloudFormation.
Previously, I've used bash scripts to create and destroy AWS resources.
CloudFormation is not perfect, but it is VERY effective for what is needed here.
Is it better? Worse? Cleaner? Or just different? I do love scripts.
Could, but not using it for creating the container repository though.



```BASH
aws cloudformation create-stack --template-body file://batch-setup.template.yaml --stack-name batch --capabilities CAPABILITY_NAMED_IAM
```


To destroy everything ...


```BASH
aws cloudformation delete-stack --stack-name batch
```













##	Limits



I was recently awarded AWS research credits and I am currently developing the pipeline to do the processing. Later this year, I intend on beginning this processing with over 100,000 jobs using SPOT instances and Batch. Some of these jobs take several minutes while some take several days. My jobs will stream data mostly from S3 but also from the internet and then push the results back to S3. They seem optimized for the C5 family having each job use 2 cores and about 3.5GB of memory.

I'd like to up all of my limits for the US East 1 (N Virginia) region to enable the fastest processing.
I would like to run about 1000 jobs at a time so that I can complete this processing in about a week or so.

Increase my instance limits to 100 (or higher) for each in the C5 family.
Increase my General Purpose (SSD) volume storage (TiB) from 300	to 500 ( or higher).
Increase my SPOT instances from 20 to 100 ( or higher ).

Are there any other obvious resources that I may need to increase the limit on to optimize this processing?


22GB gp2 per container. 8GB gp2 per instance.
if running 1000 jobs at once, and if 1 container on 1 instance, 30TB needed. Have 300TB

IOPS is 100 per volume so 100 per instance plus 100 per container. Seems enough.
Again, 200 IOPS per job, 1000 jobs. 200,000 IOPS. Have 300,000 "Provisioned" IOPS. Same?

I'm not sure that I need to up my instance limits. I've started c5.2xlarge, which my limit for is 0.


Many of AWS' resources have limits attached to them.
I'm not certain how Batch abides by these limits as I know some instances have been created that I have a 0 limit on.
How many instances will Batch create?
How many jobs will it try to run at once?
How much disk space will each use?

Nevertheless, one should be aware.




https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/using-spot-limits.html?icmpid=docs_ec2_console#spot-limits-general

Spot Instance limits are dynamic. When your account is new, your limit might be lower than 20 to start, but can increase over time. In addition, your account might have limits on specific Spot Instance types. If you submit a Spot Instance request and you receive the error Max spot instance count exceeded, you can complete the AWS Support Center Create case form to request a Spot Instance limit increase. For Limit type, choose EC2 Spot Instances. For more information, see Amazon EC2 Service Limits.

8/19/2019, 4:07:13 PM	information	launchSpecTemporarilyBlacklisted	Repeated errors have occurred processing the launch specification "c5.18xlarge, ami-0d09143c6fc181fe3, Linux/UNIX, us-east-1a while launching spot instance". It will not be retried for at least 13 minutes. Error message: Spot Max Instance Count Exceeded


Stuck in a "loop".

```
8/19/2019, 7:18:47 PM	error	spotInstanceCountLimitExceeded
8/19/2019, 7:08:59 PM	information	launchSpecTemporarilyBlacklisted	Repeated errors have occurred processing the launch specification "c5.18xlarge, ami-0d09143c6fc181fe3, Linux/UNIX, us-east-1a while launching spot instance". It will not be retried for at least 13 minutes. Error message: Spot Max Instance Count Exceeded
8/19/2019, 7:08:59 PM	error	allLaunchSpecsTemporarilyBlacklisted	Several attempts to launch instances have failed. Either the request could not be satisfied or the configuration is not valid. We will retry the request again later. For more information, see the description of the event.
8/19/2019, 7:07:58 PM	information	launchSpecTemporarilyBlacklisted	Repeated errors have occurred processing the launch specification "c5.18xlarge, ami-0d09143c6fc181fe3, Linux/UNIX, us-east-1c while launching spot instance". It will not be retried for at least 13 minutes. Error message: Spot Max Instance Count Exceeded
8/19/2019, 7:06:58 PM	information	launchSpecTemporarilyBlacklisted	Repeated errors have occurred processing the launch specification "c5.18xlarge, ami-0d09143c6fc181fe3, Linux/UNIX, us-east-1b while launching spot instance". It will not be retried for at least 13 minutes. Error message: Spot Max Instance Count Exceeded
8/19/2019, 7:05:38 PM	information	launchSpecTemporarilyBlacklisted	Repeated errors have occurred processing the launch specification "c5.18xlarge, ami-0d09143c6fc181fe3, Linux/UNIX, us-east-1f while launching spot instance". It will not be retried for at least 13 minutes. Error message: Spot Max Instance Count Exceeded
8/19/2019, 7:04:37 PM	information	launchSpecTemporarilyBlacklisted	Repeated errors have occurred processing the launch specification "c5.18xlarge, ami-0d09143c6fc181fe3, Linux/UNIX, us-east-1d while launching spot instance". It will not be retried for at least 13 minutes. Error message: Spot Max Instance Count Exceeded
8/19/2019, 7:03:37 PM	error	spotInstanceCountLimitExceeded
8/19/2019, 7:03:26 PM	information	launchSpecUnusable	c5.18xlarge, ami-0d09143c6fc181fe3, Linux/UNIX, us-east-1e, Spot bid price is less than Spot market price $-1.0000
8/19/2019, 6:48:26 PM	error	spotInstanceCountLimitExceeded
```







#	Job Submission



Something like.

```BASH
aws batch submit-job --job-name testing4 --job-definition myJobDefinition --job-queue myJobQueue --container-overrides '{ "command": ["echo","testing","command","line"]}'

aws batch submit-job --job-name TESTING1 --job-definition myJobDefinition --job-queue myJobQueue --container-overrides command="echo testing command line"   NOPE

aws batch submit-job --job-name s3___1000genomes_phase3_data_HG01863_sequence_read_SRR395998_1_filt_fastq_gz --job-definition myJobDefinition --job-queue myJobQueue --container-overrides { "command": ["viral_identification.bash","s3://1000genomes/phase3/data/HG01863/sequence_read/SRR395998_1.filt.fastq.gz"]}


aws batch submit-job --job-name ftp___ftp_sra_ebi_ac_uk_vol1_fastq_ERR188_ERR188022_ERR188022_1_fastq_gz --job-definition myJobDefinition --job-queue myJobQueue --container-overrides { "command": ["viral_identification.bash","ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188022/ERR188022_1.fastq.gz"]}
```









This is effective for several, perhaps hundreds of jobs.
However, when submitting ~100,000 jobs, it just isn't very efficient.
If it is the intention to submit massive numbers of jobs, the "array job" should be considered.

https://docs.aws.amazon.com/batch/latest/userguide/array_index_example.html


An array job will submit the number of jobs defined by `size` each setting the environment variable `AWS_BATCH_JOB_ARRAY_INDEX` from 0 to `size-1`.
The submitted script just needs to deal with this variable.
`size` needs to be between 2 and 10,000.
I have the need for creating an array job of about 109,000 plus or minus a few so I've added a `page` argument.





```BASH
aws batch submit-job --job-name array_test_1 --job-definition myJobDefinition --job-queue myJobQueue --array-properties size=5 --container-overrides '{ "command": ["array_handler.bash","1000genomes","1"]}'

aws batch submit-job --job-name array_test_2 --job-definition myJobDefinition --job-queue myJobQueue --array-properties size=5 --container-overrides '{ "command": ["array_handler.bash","Other","10"]}'

aws batch submit-job --job-name array_test_2 --job-definition myJobDefinition --job-queue myJobQueue --array-properties size=5 --container-overrides '{ "command": ["array_handler.bash","Other","10"], "vcpus": 1, "memory": 100, "instanceType": "optimal" }'

aws batch submit-job --job-name array_test_1 --job-definition myJobDefinition --job-queue myJobQueue --array-properties size=100 --container-overrides '{ "command": ["array_handler.bash","1000genomes","1"], "vcpus": 1, "memory": 100, "instanceType": "optimal"}'
```


can also override the docker image



SPOT Instance Requests Limit was increased from "default" to 180.
Attempt at 100 item array job started and executed nicely this morning.


```
However, the limits for Spot Instances are not raised depending upon the type of instance but the total Spot Instance limits in general. For example, in this case, I've raised the limits for c5.18xlarges Spot Instances upto 100 but you can also launch other c5 instances depending upon the size. For your convenience, I've listed the number of Spot Instances you can launch depending upon its size
> upto 74 c5.24xlarges Spot Instances with the current Spot Instance limits
> upto 100 c5.18xlarges Spot Instances with the current Spot Instance limits
> upto 112 c5.16xlarges Spot Instances with the current Spot Instance limits
> upto 150 c5.12xlarges Spot Instances with the current Spot Instance limits
Not sure how this translates into the 180.
Also, given that Batch is starting these instances, the numbers are beyond my control.
It is still unclear to me as to why it wouldn't start my SPOT instance even when I had none running.

8/20/2019, 12:51:31 PM	instanceChange	launched	{"instanceType":"c5.18xlarge","image":"ami-0d09143c6fc181fe3","productDescription":"Linux/UNIX","availabilityZone":"us-east-1d"}	i-0593e7f9bf14cfee1
8/20/2019, 12:51:31 PM	instanceChange	launched	{"instanceType":"c5.18xlarge","image":"ami-0d09143c6fc181fe3","productDescription":"Linux/UNIX","availabilityZone":"us-east-1d"}	i-0ed25aae31baf8885
8/20/2019, 12:51:31 PM	bidChange	active	BidId sir-2k384zik, PreviousState: active
8/20/2019, 12:51:31 PM	bidChange	active	BidId sir-49zr6mfk, PreviousState: active
8/20/2019, 12:51:29 PM	fleetRequestChange	progress	c5.18xlarge, ami-0d09143c6fc181fe3, Linux/UNIX, us-east-1d, capacityUnitsRequested: 2.0, totalCapacityUnitsRequested: 2.0, totalCapacityUnitsFulfilled: 0.0, targetCapacity: 2
8/20/2019, 12:51:28 PM	fleetRequestChange	active
8/20/2019, 12:51:28 PM	information	launchSpecUnusable	c5.18xlarge, ami-0d09143c6fc181fe3, Linux/UNIX, us-east-1e, Spot bid price is less than Spot market price $-1.0000
8/20/2019, 12:51:18 PM	fleetRequestChange	submitted
```









#	Job Debugging


If you created a KeyPair and want to connect to your running instance and container,
get the ip address from the web console or command line and `ssh` to it.

Once connected, list the docker containers. There should be at least 2.
One for the "agent" and one for each of your running containers.
Connect to it with a `docker exec` command.

```BASH

aws ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress'

[
    "3.91.106.157",
    "18.207.137.221"
]

ssh -i ~/.ssh/batchKeyPair -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ec2-user@PUBLIC_IP_ADDRESS


docker ps -a
CONTAINER ID        IMAGE                  COMMAND             CREATED             STATUS              PORTS               NAMES
235487af1717        viral_identification   "/bin/bash"         8 seconds ago       Up 7 seconds                            wonderful_fermi

docker exec -it 235487af1717 /bin/bash
```

If you used a different Linux base for your docker container, the username will likely be different.




#	Total time

I processed the 7th largest file in 13 hours.

```BASH
31468335171 data/HG00251/sequence_read/ERR018427_2.filt.fastq.gz

31468335171 bases / 13 hours
2420641167 bases / hour

awk '{s+=$1}END{print s}' 1kg.files
71033560137387 total bases
71033560137387 bases / 2420641167(bases/hour)
29344.94 hours
```

Approximately 30,000 hours of total processing time.




#	AWS Questions

##	What happens if spot instance going to be killed? how to react?


##	Why don't Job queues iterate over multiple spot instances compute environments?
Not sure if I worded that clearly enough.


##	Will cancelled spot instances rerun the job?

Not on the default MANAGED Compute Environment. Possibly if self managed.


##	How to change the bid percentage of a Batch Compute Environment?

Seems like bidPercentage isn't valid in the update-compute-environment operation.
Can't?


##	How to delete Batch Job Definitions?

I can only deregister them and then they don't go away?
Can't?


##	How to clear "Deletion failed" for AWSServiceRoleForTrustedAdvisor?

I was cleaning up existing roles and tried to delete them all. This one failed, produced a message and highlighted in red. This message and style will not go away.

I ended up turning off the feature, deleting the role, then recreating.


##	Why can't delete invalid Compute Envirionments?

Cloud Formation or just by accident, if roles or job queue deleted first, can't delete compute environment


##	Batch Job Queue doesn't cascade over different Compute Environments

I created 3 identical compute environments each with higher SPOT bids. 35, 40 and 45%.
Recently, the first was too low and it just sat there for 30 minutes.
Investigation showed that the SPOT request errored.
The bid price was too low.
Batch didn't even try the second and third compute environments.
Does this not work with SPOT?
Do I just need to wait longer?

I disabled the "stuck" compute environment and the spot request was canceled.
The jobs stayed RUNNABLE for about 10 minutes and then a new spot request
was issued with the next compute environment. Then they ran.






#	References


https://aws.amazon.com/research-credits/


https://aws.amazon.com/blogs/compute/creating-a-simple-fetch-and-run-aws-batch-job/

https://stackify.com/aws-batch-guide/

https://docs.aws.amazon.com/batch/latest/userguide/Batch_GetStarted.html

https://docs.aws.amazon.com/batch/latest/userguide/troubleshooting.html


https://aws.amazon.com/blogs/compute/creating-a-simple-fetch-and-run-aws-batch-job/

https://stackify.com/aws-batch-guide/

https://docs.aws.amazon.com/batch/latest/userguide/Batch_GetStarted.html

https://docs.aws.amazon.com/batch/latest/userguide/troubleshooting.html

https://github.com/aws-samples/aws-batch-genomics/blob/v1.0.0/batch/setup/iam.template.yaml

https://github.com/aws-samples/aws-batch-genomics/blob/v1.0.0/batch/setup/batch_env.template.yaml


https://www.reddit.com/r/aws/comments/a3vhtp/cloudformation_is_not_linking_my_routetable_with/

https://keithmsharp.wordpress.com/2016/11/15/building-a-vpc-with-aws-cloudformation/


https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/aws-properties-batch-computeenvironment-computeresources.html#cfn-batch-computeenvironment-computeresources-instancerole



https://docs.aws.amazon.com/AWSCloudFormation/latest/UserGuide/pseudo-parameter-reference.html




