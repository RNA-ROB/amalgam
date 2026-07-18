#!/usr/bin/env python3

'''
Author: Robert Wang (Xing Lab)
Date: 2026.07.17
Version: 1.0.0

This is a script designed to filter the GTF file produced by Build_Transcriptome.py to include annotations
for transcripts detected with at least X reads in at least one input sample (by default, X = 1). This script
requires the following three inputs:
    1. A GTF file containing reference transcript annotations
    2. A text file in which each line specifies the path to each input sample's transcript counts file
    3. Minimum read count threshold to decide if a transcript is detected or not (by default, this value is 1)

This script will generate a smaller GTF file containing annotations for transcripts detected with at least X reads in
at least one input sample. 
'''

# =====================================================================================================================
#                                                       PACKAGES 
# =====================================================================================================================

# Load required packages
from datetime import datetime
import argparse, csv
import pandas as pd

# =====================================================================================================================
#                                                   HELPER FUNCTIONS
# =====================================================================================================================

def outlog(myString):
    '''
    This is a function to print some output message string (myString) with the date and time
    '''
    print('[', datetime.now().strftime("%Y-%m-%d %H:%M:%S"), '] ', myString, sep = '', flush = True)

def PullFeature(infoString, feature):
    '''
    This is a function designed to pull out the value of a particular feature in a GTF info field
    '''
    return next(iter([item.split('"')[1] for item in infoString.split(';') if feature in item]), '.')

# =====================================================================================================================
#                                                    MAIN FUNCTIONS
# =====================================================================================================================

def main():
    message = 'Filters input annotations for transcripts detected with at least X reads in at least one input sample'
    parser = argparse.ArgumentParser(description = message)

    # Add arguments
    parser.add_argument('-i', metavar = '/path/to/input/GTF/file', required = True,
        help = 'path to input GTF file')
    parser.add_argument('-s', metavar = '/path/to/list/of/files', required = True,
        help = 'path to file containing paths to transcript counts files for input samples')
    parser.add_argument('-x', metavar = 'threshold', type = float, default = 1.0, 
        help = 'minimum read count threshold for a transcript to be considered detected')
    parser.add_argument('-o', metavar = '/path/to/output/GTF/file', required = True,
        help = 'path to output GTF file')
    
    # Parse command-line arguments
    args = parser.parse_args()
    infile, fplist, threshold, outfile = args.i, args.s, args.x, args.o

    # Read in files in fplist and construct a transcript read count matrix across input samples
    idx, outDF = 0, pd.DataFrame(columns = ['gene_id', 'transcript_id'])
    with open(fplist, 'r') as file:
        for line in file:
            currDF = pd.read_csv(line.strip(), sep = '\t', header = 0)
            currDF.columns = ['gene_id', 'transcript_id', idx]
            outDF = pd.merge(outDF, currDF, how = 'outer', on = ['gene_id', 'transcript_id'])
            idx += 1
            outlog('Reading in ' + line.strip())
    
    outlog('Reading in ' + infile)

    # Only keep transcripts where the maximum read count across all samples is at least threshold
    outDF = outDF[outDF.drop(['gene_id', 'transcript_id'], axis = 1).max(axis = 1) >= threshold]
    keepGenes, keepTranscripts = set(outDF['gene_id']), set(outDF['transcript_id'])

    # Read in infile as a dataframe and extract gene ID and transcript ID fields for each entry
    annoDF = pd.read_csv(infile, sep = '\t', header = None, comment = '#')
    annoDF['gene_id'] = annoDF[8].apply(lambda x: PullFeature(x, 'gene_id'))
    annoDF['transcript_id'] = annoDF[8].apply(lambda x: PullFeature(x, 'transcript_id'))

    # Filter annoDF for entries involving genes and transcripts in keepGenes and keepTranscripts
    annoDF = annoDF[annoDF['gene_id'].isin(keepGenes) & annoDF['transcript_id'].isin(keepTranscripts)]
    annoDF = annoDF.drop(['gene_id', 'transcript_id'], axis = 1)

    # Save annoDF to output file
    annoDF.to_csv(outfile, sep = '\t', index = False, header = False, quoting = csv.QUOTE_NONE)

if __name__ == '__main__':
    main()
