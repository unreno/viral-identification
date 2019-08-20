#	Viral Identification


Thorough, unique and unambiguous identification of viral content in a genomic sample.




#	Reference Preparation


What is the most appropriate viral reference?
RefSeq only has about 12,000 sequences.
NCBI's nt has nearly 3 million.

Why the huge difference?

What about NCBI's nr?


##	Viral Sequence Acquisition From NCBI RefSeq



```BASH
rsync -avz --progress --include="*fna.gz" --exclude="*" rsync://ftp.ncbi.nih.gov/refseq/release/viral/ ./

ls -1s viral*genomic.fna.gz
63744 viral.1.1.genomic.fna.gz
29952 viral.2.1.genomic.fna.gz

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


update_blastdb.pl --decompress nt
...



blastdbcmd -db nt -entry all -outfmt "%K" | sort | uniq -c
 451708 Archaea
8858629 Bacteria
46649264 Eukaryota
 939111 N/A
2918320 Viruses


blastdbcmd -db nt -entry all -outfmt "%K,%i" | awk -F, '( $1 == "Viruses" ){ print $2 }' > nt.viruses.seqidlist

wc -l nt.viruses.seqidlist
2918320 nt.viruses.seqidlist


blastdbcmd -db nt -entry_batch nt.viruses.seqidlist > nt.viruses.fa

grep -c "^>" nt.viruses.fa
2918320
```



```BASH
update_blastdb.pl --decompress nr
...







blastdbcmd -db nr -entry all -outfmt "%K" | sort | uniq -c
6009090 Archaea
670701365 Bacteria
65740422 Eukaryota
 312376 N/A
5694044 Viruses


blastdbcmd -db nr -entry all -outfmt "%K,%i" | awk -F, '( $1 == "Viruses" ){ print $2 }' > nr.viruses.seqidlist

wc -l nr.viruses.seqidlist


blastdbcmd -db nr -entry_batch nr.viruses.seqidlist > nr.viruses.fa

grep -c "^>" nr.viruses.fa

```




##	Reference Cleanup




















#	Script Preparation


When using AWS Batch, the job command cannot include the use of pipes as they will be converted to strings.

This basically means that your job will probably need a script.

Don't even need a script
Actually, kinda do NEED a script. The command converts pipes to strings.
But don't need a complex script. Processing is a series of pipes from S3(or ftp) to S3.
I don't think that including pipes in the command is acceptable.
They just get converted to strings and then passed as params to the first command.

viral_identification.bash

array_handler.bash





#	Docker Image Creation and Storage


Create a Docker image containing anything needed by any of your job scripts.

Add blastn executable
ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz
Add blastn masked viral reference
viral.masked.gz
Add script


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

XXXXXXXXXX.dkr.ecr.us-east-1.amazonaws.com/viral_identification:latest
```


Login, tag the image and push it to the repository.
Be sure to use the added --no-include-email option as without it the output will include deprecated "-e none"


```BASH
aws ecr get-login --no-include-email --region us-east-1 | bash

docker tag viral_identification:latest XXXXXXXXXX.dkr.ecr.us-east-1.amazonaws.com/viral_identification:latest

docker push XXXXXXXXXX.dkr.ecr.us-east-1.amazonaws.com/viral_identification:latest
```






#	AWS Preparation




##	Create an EC2 KeyPair

This is kinda optional.

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







#	Job Submission



Something like.

```BASH
aws batch submit-job --job-name testing4 --job-definition JobDefinition --job-queue JobQueue --container-overrides '{ "command": ["echo","testing","command","line"]}'

aws batch submit-job --job-name TESTING1 --job-definition JobDefinition --job-queue JobQueue --container-overrides command="echo testing command line"   NOPE

aws batch submit-job --job-name s3___1000genomes_phase3_data_HG01863_sequence_read_SRR395998_1_filt_fastq_gz --job-definition JobDefinition --job-queue JobQueue --container-overrides { "command": ["viral_identification.bash","s3://1000genomes/phase3/data/HG01863/sequence_read/SRR395998_1.filt.fastq.gz"]}


aws batch submit-job --job-name ftp___ftp_sra_ebi_ac_uk_vol1_fastq_ERR188_ERR188022_ERR188022_1_fastq_gz --job-definition JobDefinition --job-queue JobQueue --container-overrides { "command": ["viral_identification.bash","ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR188/ERR188022/ERR188022_1.fastq.gz"]}
```









This is effective for several, perhaps hundreds of jobs.
However, when submitting ~100,000 jobs, it just isn't very efficient.
If it is the intention to submit massive numbers of jobs, the "array job" should be considered.

https://docs.aws.amazon.com/batch/latest/userguide/array_index_example.html

```BASH
aws batch submit-job --job-name array_test_1 --job-definition JobDefinition --job-queue JobQueue --array-properties size=5 --container-overrides '{ "command": ["array_handler.bash","1000genomes","1"]}'

aws batch submit-job --job-name array_test_2 --job-definition JobDefinition --job-queue JobQueue --array-properties size=5 --container-overrides '{ "command": ["array_handler.bash","Other","10"]}'

aws batch submit-job --job-name array_test_2 --job-definition JobDefinition --job-queue JobQueue --array-properties size=5 --container-overrides '{ "command": ["array_handler.bash","Other","10"], "vcpus": 1, "memory": 100, "instanceType": "optimal" }'

aws batch submit-job --job-name array_test_1 --job-definition JobDefinition --job-queue JobQueue --array-properties size=100 --container-overrides '{ "command": ["array_handler.bash","1000genomes","1"], "vcpus": 1, "memory": 100, "instanceType": "optimal"}'
```

AWS_BATCH_JOB_ARRAY_INDEX=0


the array index starts at 0







#	Job Debugging



```BASH

> aws ec2 describe-instances --query 'Reservations[].Instances[].PublicIpAddress' 

[
    "3.91.106.157",
    "18.207.137.221"
]

> ssh -i ~/.ssh/batchKeyPair -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no ec2-user@PUBLIC_IP_ADDRESS


> docker ps -a


> docker exec -it CONTAINER_ID /bin/bash

```



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



