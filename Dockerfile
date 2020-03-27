FROM amazonlinux:latest
RUN yum -y install which unzip aws-cli tar gzip bzip2 gcc g++ zlib-devel bzip2-devel xz-devel make libcurl-devel ncurses-devel openssl-devel wget procps htop python3 boto3

RUN rm -f /usr/bin/python
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN ln -s /usr/bin/pip3 /usr/bin/pip
#	aws-cli/1.16.300 Python/2.7.16 Linux/4.9.125-linuxkit botocore/1.13.36

RUN pip install --upgrade awscli pip numpy scipy boto3
#	aws-cli/1.18.18 Python/3.7.6 Linux/4.9.125-linuxkit botocore/1.15.18


ARG SAMTOOLS_VERSION=1.9
ARG SAMTOOLS_URL=https://github.com/samtools/samtools/releases/download
ARG HTSLIB_URL=https://github.com/samtools/htslib/releases/download
 
RUN cd / \
	&& wget ${HTSLIB_URL}/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xvfj htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd htslib-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& cd ~ \
	&& /bin/rm -rf /htslib-${SAMTOOLS_VERSION}*
 
RUN cd / \
	&& wget ${SAMTOOLS_URL}/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar xvfj samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& cd samtools-${SAMTOOLS_VERSION} \
	&& ./configure \
	&& make \
	&& make install \
	&& cd ~ \
	&& /bin/rm -rf /samtools-${SAMTOOLS_VERSION}*




#	
#	RUN cd / \
#		&& wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz \
#		&& tar xvfz ncbi-blast-2.9.0+-x64-linux.tar.gz \
#		&& mv ncbi-blast-2.9.0+/bin/* /usr/local/bin/ \
#		&& /bin/rm -rf ncbi-blast-2.9.0+-x64-linux.tar.gz ncbi-blast-2.9.0+ 
#	
#	ENV BLASTDB=/blastdb
#	
#	#	ADD will actually untar and gunzip for you.
#	ADD references/viral.masked.tar.gz /blastdb/
#	ADD references/viral.raw.tar.gz /blastdb/
#	
#	USING DIAMOND INSTEAD OF BLASTN
#




RUN cd / \
	&& wget https://github.com/bbuchfink/diamond/releases/download/v0.9.30/diamond-linux64.tar.gz \
	&& tar xvfz diamond-linux64.tar.gz \
	&& mv diamond /usr/local/bin/ \
	&& /bin/rm -rf diamond-linux64.tar.gz

ENV DIAMOND=/diamond

#	ADD will actually untar and gunzip for you.
#	doesn't seem to gunzip
ADD references/viral.dmnd.tar.gz /diamond/







ADD array_handler.bash /usr/local/bin/array_handler.bash
ADD viral_identification.bash /usr/local/bin/viral_identification.bash


#	RUN aws s3 cp --no-sign-request s3://1000genomes/sequence.index - | awk -F"\t" '(NR>1) && ($13 ~ /ILLUMINA/) && ($25 != "not available") && ($26 ~ /coverage/){print $25,"s3://1000genomes/phase3/"$1}' | sort -n -r | awk '{print $2}' > /tmp/1000genomes.urls
ADD urls/1000genomes.urls /tmp/

#	RUN curl https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt | awk -F"\t" '(NR>1){print $35}' > /tmp/geuvadis.urls
ADD urls/geuvadis.urls /tmp/

#	RUN curl https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt | awk -F"\t" '(NR>1){print "https://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/"$27".bam"}' | uniq > geuvadis.bam.urls 
ADD urls/geuvadis.bam.urls /tmp/

#	RUN aws s3 ls --recursive --no-sign-request s3://1000genomes/phase3/ | grep ".unmapped.ILLUMINA." | grep "_coverage." | grep ".bam$" | awk '{print "s3://1000genomes/"$4}' > 1000genomes.unmapped.urls
#	RUN aws s3 cp --no-sign-request s3://1000genomes/alignment.index - | grep ".unmapped.ILLUMINA." | grep ".low_coverage." | awk '{print "s3://1000genomes/phase3/"$1}' > 1000genomes.unmapped.urls
ADD urls/1000genomes.unmapped.urls /tmp/

#wc -l *.urls
# 108873 1000genomes.urls
#    924 geuvadis.urls
#   2535 1000genomes.unmapped.urls
#    462 geuvadis.bam.urls



#	DO THIS ONLY FOR LOCAL TESTING
#	tar xfvz aws.tar.gz
#	chmod 666 .aws/c*
#COPY .aws /.aws/



#
#	docker run --rm gwendt/fetch_and_run will execute the following command and remove itself ...
#	Either include the script to be run in the docker image itself 
#	or place on S3 and pass as environment variable to job.
#
#ADD fetch_and_run.sh /usr/local/bin/fetch_and_run.sh
WORKDIR /tmp
USER nobody
#
#ENTRYPOINT ["/usr/local/bin/fetch_and_run.sh"]
