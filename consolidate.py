import sys,os,argparse
parser=argparse.ArgumentParser(description='Script for obtaining consolidated barcodes from MAFFT+AA and RACON+AA barcodes')
parser.add_argument('-m','--mafft',help='Path to input mafft corrected barcode fasta file',dest="mafftb",required=True)
parser.add_argument('-r','--racon',help='Path to input racon corrected barcode fasta file',dest="raconb",required=True)
parser.add_argument('-o','--outfile',help='Path to output corrected barcode fasta file',dest="outfile",required=True)
parser.add_argument('-t','--tempdir',help='Path to tempdirectory',dest="temp",required=True)

args=parser.parse_args()
os.system("mkdir "+args.temp)
def change_ext_gaps(sequence):
	bps_base=['A','T','G','C','N']
	start_pos, end_pos = 0, 0
	for i,bp in enumerate(sequence):
		if bp in bps_base:
			start_pos = i - 1
			break
	for i,bp in enumerate(sequence[:: - 1]):
		if bp in bps_base:
			end_pos = len(sequence) - i
			break
	new_sequence = "N" * (start_pos + 1) + sequence[start_pos + 1 : end_pos] + "N" * (len(sequence) - end_pos)
	return new_sequence
def reformat(inputfile, outputfile):
	with open(inputfile) as input:
		with open(outputfile,'w') as output:
			l=input.readlines()
			n=0
			for line in l:
				if n==0 and ">" in line:
					output.write(line)
				elif ">" in line:
					output.write("\n"+line)
				else:
					output.write(line.strip())
				n+=1
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
				bp="N"
	#			break
			newseq+=bp
		except:
			pass
	return newseq

mafftbarcodes={}
raconbarcodes={}
with open(args.mafftb) as mfile:
	l=mfile.readlines()
	for i,line in enumerate(l):
		if ">" in line:
			mafftbarcodes[line.split(";")[0][1:].strip()]=l[i+1].strip()
print mafftbarcodes.keys()
with open(args.raconb) as rfile:
	l=rfile.readlines()
	for i,line in enumerate(l):
		if ">" in line:
			raconbarcodes[line.split(";")[0][1:].strip()]=l[i+1].strip()
print raconbarcodes.keys()
with open(args.outfile,'w') as outfile:	
	for barcode in mafftbarcodes:
		if barcode in raconbarcodes:
			with open(args.temp+"/"+barcode,'w') as ofile:
				ofile.write(">"+barcode+'\n'+mafftbarcodes[barcode]+'\n'+">"+barcode+'\n'+raconbarcodes[barcode]+'\n')
			os.system("mafft "+args.temp+"/"+barcode+">"+args.temp+"/"+barcode+"_mafft")
			reformat(args.temp+"/"+barcode+"_mafft",args.temp+"/"+barcode+"_mafftreformat")
			with open(args.temp+"/"+barcode+"_mafftreformat") as infile:
				l=infile.readlines()
				m=change_ext_gaps(l[1].strip().upper())
				r=change_ext_gaps(l[3].strip().upper())
				print m
				print r
				if "-" not in m:
					if "-" not in r:
						conbarcode=get_consensus(m,r)
						outfile.write(">"+barcode+'\n'+conbarcode+'\n')
				else:
					print barcode
		
