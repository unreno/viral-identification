#!/usr/bin/env python


import pandas as pd
import numpy as np
import glob
import os

def lane_or_dot(row):
	parts=row['qaccver'].split('/')
	if len(parts) > 1:
		return parts[-1]
	else:
		return '.'


startdir = os.getcwd()

dfs=[]

#	expecting files like ...
#./1000genomes/phase3/data/NA11931/sequence_read/ERR243205.filt.viral.raw.tsv.gz
#./1000genomes/phase3/data/NA11931/sequence_read/ERR243205.filt.viral.masked.tsv.gz
#./1000genomes/phase3/data/HG01384/alignment/HG01384.unmapped.ILLUMINA.bwa.CLM.low_coverage.20120522.viral.masked.tsv.gz
#./1000genomes/phase3/data/HG01384/alignment/HG01384.unmapped.ILLUMINA.bwa.CLM.low_coverage.20120522.viral.raw.tsv.gz
#./1000genomes/phase3/data/HG01377/alignment/HG01377.unmapped.ILLUMINA.bwa.CLM.low_coverage.20120522.viral.raw.tsv.gz
#./1000genomes/phase3/data/HG01377/alignment/HG01377.unmapped.ILLUMINA.bwa.CLM.low_coverage.20120522.viral.masked.tsv.gz
#./1000genomes/phase3/data/NA18968/sequence_read/SRR003916.filt.viral.masked.tsv.gz
#./1000genomes/phase3/data/NA18968/sequence_read/SRR003916.filt.viral.raw.tsv.gz
#./1000genomes/phase3/data/HG02420/sequence_read/SRR449493.filt.viral.raw.tsv.gz
#./1000genomes/phase3/data/HG02420/sequence_read/SRR449493.filt.viral.masked.tsv.gz
#./1000genomes/phase3/data/HG03108/sequence_read/SRR594817.filt.viral.raw.tsv.gz
#./1000genomes/phase3/data/HG03108/sequence_read/SRR594817.filt.viral.masked.tsv.gz



#	gonna need to change some stuff when add other datasets



base="/Users/jakewendt/s3/viral-identification"

os.chdir(base)

for filepath in glob.iglob('**/*.viral.*.tsv.gz', recursive=True):
	print(filepath)

	#	The file may exist even if empty so ...
	try:
		f = pd.read_csv(filepath, header=None, sep="\t", 
			names=['qaccver','saccver','pident','length','mismatch','gapopen',
				'qstart','qend','sstart','send','evalue','bitscore'] )
		filepath_pieces=filepath.split('/')
		f['subject']=filepath_pieces[3]
		f['source']=filepath_pieces[0]

		file_pieces=filepath_pieces[5].split('.')

		if filepath_pieces[4] == 'sequence_read':
			f['format']='fastq'
			#	assign lane (the unlaned file contains "unpaired" reads, but the raw file shows a lane)
			filelane = file_pieces[0].split('_')
			if len(filelane) > 1:
				f['filelane'] = filelane[1]
			else:
				f['filelane'] = '.'
		else:
			f['format']='unmapped_bam'

		f['reference'] = file_pieces[-3]
		dfs.append(f)
	except:
		continue

df = pd.concat(dfs,sort=False)	#sort=True)	# do not sort column names
#df=df.fillna('.')
#df.sort_index(inplace=True)

df['lane'] = df.apply (lambda row: lane_or_dot(row), axis='columns')

df.insert(0, 'subject', df.pop('subject') )

os.chdir(startdir)

df.to_csv('viral.merged.csv',index=False)

