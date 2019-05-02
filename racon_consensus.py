#!/usr/bin/env python
import argparse,os
parser=argparse.ArgumentParser(description='Script for obtaining racon barcodes')
parser.add_argument('-i','--indir',help='output folder path of mb_parallel_demultiplex',dest="indir",required=True)
parser.add_argument('-b','--barcodefile',help='input barcode set',dest="barcodefile",required=True)
parser.add_argument('-d','--outdir',help='outputdirectory',dest="outdir",required=True)
parser.add_argument('-o','--outfile',help='output filen name',dest="outfile",required=True)

args=parser.parse_args()
if os.path.isdir(args.outdir):
	print "output directory exists, please delete or name another"
	sys.exit()
else:
	os.system("mkdir "+args.outdir)

with open(args.barcodefile) as infile:
	l=infile.readlines()
	for i,j in enumerate(l):
		if ">" in j:
			with open(args.outdir+"/"+j[1:].split(";")[0].strip()+".fa",'w') as outfile:
				outfile.write(j+l[i+1])	
dirlist=os.listdir(args.outdir)
for f in dirlist:
	os.system("graphmap align --max-error 0.15 -r "+args.outdir+"/"+f+" -d "+ args.indir+"/demultiplexed/"+f+" -o "+args.outdir+"/"+f.split(".")[0]+".sam")
	os.system("racon "+args.indir+"/demultiplexed/"+f+" "+args.outdir+"/"+f.split(".")[0]+".sam " +args.outdir+"/"+f+">"+args.outdir+"/"+f.split(".")[0]+"_racon.fasta")

os.system("cat "+args.outdir+"/*racon.fasta > "+args.outfile)

