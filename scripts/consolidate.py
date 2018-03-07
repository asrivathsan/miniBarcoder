import sys,os,argparse
parser=argparse.ArgumentParser(description='Script for obtaining consolidated barcodes from MAFFT+AA and RACON+AA barcodes')
parser.add_argument('-m','--mafft',help='Path to input mafft corrected barcode fasta file',dest="mafftb",required=True)
parser.add_argument('-r','--racon',help='Path to input racon corrected barcode fasta file',dest="raconb",required=True)
parser.add_argument('-o','--outfile',help='Path to output corrected barcode fasta file',dest="outfile",required=True)
args=parser.parse_args()
def get_consensus(seq1,seq2):
	newseq=''
	for i,j in enumerate(seq1):
		flag=False
		try:
			if j==seq2[i]:
				bp=j	
			elif j=="N":
				bp=seq2[i]
			elif seq2[i]=="N":
				bp=j
			else:
				flag=True
				break
			newseq+=bp
		except:
			flag=True
	if flag==True:
		return seq1
	else:
		return newseq
		
mafftbarcodes={}
raconbarcodes={}
with open(args.mafftb) as mfile:
	l=mfile.readlines()
	for i,line in enumerate(l):
		if ">" in line:
			mafftbarcodes[line.split(";")[0][1:].strip()]=l[i+1].strip()
with open(args.raconb) as rfile:
	l=rfile.readlines()
	for i,line in enumerate(l):
		if ">" in line:
			raconbarcodes[line.split(";")[0][1:].strip()]=l[i+1].strip()
with open(args.outfile,'w') as outfile:	
	for barcode in mafftbarcodes:
		if barcode in raconbarcodes:
			print barcode
			conbarcode=get_consensus(mafftbarcodes[barcode],raconbarcodes[barcode])
			outfile.write(">"+barcode+'\n'+conbarcode+'\n')
		
