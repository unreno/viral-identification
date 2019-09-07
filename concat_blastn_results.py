#!/usr/bin/env python

import pandas as pd
import numpy as np
import glob
import os


startdir = os.getcwd()

dfs=[]

base="/Users/jakewendt/s3/viral-identification"

os.chdir(base)

#for filepath in glob.iglob('**/*.viral.*.tsv.gz', recursive=True):
for filepath in glob.iglob('**/*.viral.masked.tsv.gz', recursive=True):
	print(filepath)

	f = pd.read_csv(filepath, header=None, sep="\t",
		names=['qaccver','saccver','pident','length','mismatch','gapopen',
			'qstart','qend','sstart','send','evalue','bitscore'] )

	if len(f.index) == 0:
		continue

	filepath_pieces=filepath.split('/')
	f.insert( 0, 'file', filepath_pieces[-1] )
	f.insert( 1, 'source', filepath_pieces[0] )

	# raw or masked ( ...viral.masked.tsv.gz )
	#f['reference'] = f.file.str.split('.',expand=True).iloc[:,-3]

	dfs.append(f)

###########################################################################

df = pd.concat(dfs, sort=False, ignore_index=True)

#df=df.fillna('.')
#df.sort_index(inplace=True)


pd.set_option('display.max_rows', None)

#
#	Is it faster to do all at once after concat, or individually before concat?
#
# "expand" option is required here

df['lane'] = (df.qaccver + "/0").str.split('/',expand=True)[1]
df['sample'] = df.qaccver.str.split('.', expand=True)[0]

#df['reference'] = df.file.str.split('.', expand=True)[-3]	#	# raw or masked ( ...viral.masked.tsv.gz )
#	Creates a dataframe, so negative subscripts are based on the longest array.
#	Easier just to do this individually.
#print(df.file.str.split('.', expand=True) )
#                   0      1      2      3      4       5    6     7     8     9    10    11
#0        ERR188231_2  fastq     gz  viral    raw     tsv   gz  None  None  None  None  None
#1        ERR188231_2  fastq     gz  viral    raw     tsv   gz  None  None  None  None  None
#2        ERR188231_2  fastq     gz  viral    raw     tsv   gz  None  None  None  None  None
#3        ERR188231_2  fastq     gz  viral    raw     tsv   gz  None  None  None  None  None
#4        ERR188231_2  fastq     gz  viral    raw     tsv   gz  None  None  None  None  None
#...              ...    ...    ...    ...    ...     ...  ...   ...   ...   ...   ...   ...
#8868846    SRR584240   filt  fastq     gz  viral  masked  tsv    gz  None  None  None  None
#8868847    SRR584240   filt  fastq     gz  viral  masked  tsv    gz  None  None  None  None
#8868848    SRR584240   filt  fastq     gz  viral  masked  tsv    gz  None  None  None  None
#8868849    SRR584240   filt  fastq     gz  viral  masked  tsv    gz  None  None  None  None
#8868850    SRR584240   filt  fastq     gz  viral  masked  tsv    gz  None  None  None  None



#df['source_format']=''
#df['source_format'][df['file'].str.contains('fastq')] = 'fastq'
#df['source_format'][df['file'].str.contains('bam')] = 'bam'
#df['source_format'] = df.apply(lambda row: 'bam' if row['file'].str.contains('bam') else 'fastq')

#	Would really like to do this in a one liner. Faster?

#df.loc[df['file'].str.contains('fastq'), 'source_format'] = 'fastq'
#df.loc[df['file'].str.contains('bam'), 'source_format'] = 'bam'

df['source_format'] = df.file.apply(lambda f: 'bam' if 'bam.viral' in f else 'fastq')

#df['elderly'] = np.where(df['age']>=50, 'yes', 'no')
#df['color'] = np.where(df['Set']=='Z', 'green', 'red')




os.chdir(startdir)

df.to_csv('viral.merged.csv',index=False)

