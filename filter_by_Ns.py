#!/usr/bin/env python 
import argparse

parser=argparse.ArgumentParser(description='Script for filtering a barcode fasta with number of Ns')
parser.add_argument('-i','--infile',help='Path to input fasta file',dest="infile",required=True)
parser.add_argument('-n','--nambs',help='Number of ambiguities allowed',type=int,dest="nambs",required=True)

args=parser.parse_args()
seqdict={}
with open(args.infile) as infile:
	with open(args.infile.split(".")[0]+"_Nfilter.fa",'w') as outfile:
		l=infile.readlines()
		for i,j in enumerate(l):
			if ">" in j:
				if l[i+1].upper().count("N")<=int(args.nambs):
					outfile.write(j+l[i+1])

