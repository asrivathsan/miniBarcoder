#!/usr/bin/env python
import argparse,os
parser=argparse.ArgumentParser(description='Script for splitting barcodes to separate fasta for graphmap')
parser.add_argument('-i','--infasta',help='Path to input racon corrected barcode fasta file',dest="infasta",required=True)
parser.add_argument('-o','--outdir',help='Path to output file containing statistics',dest="outdir",required=True)

args=parser.parse_args()
if os.path.isdir(args.outdir):
	print "output directory exists, please delete or name another"
	sys.exit()
else:
	os.system("mkdir "+args.outdir)
with open(args.infasta) as infile:
	l=infile.readlines()
	for i,j in enumerate(l):
		if ">" in j:
			with open(args.outdir+"/"+j[1:].split(";")[0].strip()+".fa",'w') as outfile:
				outfile.write(j+l[i+1])
