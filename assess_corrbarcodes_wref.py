import sys,os,fileinput,re,argparse
parser=argparse.ArgumentParser(description='Script for assessing correct barcodes using Sanger/illumina references')
parser.add_argument('-m','--minionfasta',help='Path to input mafft corrected barcode fasta file',dest="minionfasta",required=True)
parser.add_argument('-r','--reffasta',help='Path to input racon corrected barcode fasta file',dest="reffasta",required=True)
parser.add_argument('-t','--tempdir',help='Path to input racon corrected barcode fasta file',dest="outdir",required=True)
parser.add_argument('-o','--outfile',help='Path to output file containing statistics',dest="outfile",required=True)

args=parser.parse_args()
if os.path.isdir(args.outdir):
	print "output directory exists, please delete or name another"
	sys.exit()
else:
	os.system("mkdir "+args.outdir)
def measurepair_distancem(seq1, seq2):
	num_d = 0
	num_s = 0
	diff=[]
	num_g=0
	for ele_x, ele_y in zip(seq1,seq2):
		if ele_x!=ele_y:
			if ele_x in ["A","T","G","C"] and ele_y in ["A","T","G","C"]:
					num_d += 1
					num_s+=1
			if ele_x =="-" and ele_y in ["A","T","G","C","N"]:
					num_g += 1
			if ele_y =="-" and ele_x in ["A","T","G","C","N"]:
					num_g += 1
	return num_d,num_g
def reformat(inputfile, outputfile):
	with open(outputfile,'w') as output:
		n=0
		for line in fileinput.input([inputfile]):
			if n==0 and ">" in line:	
				output.write(line)
			elif ">" in line:
				output.write("\n"+line)
			else:
				output.write(line.strip())
			n+=1
referencebarcodes={}

with open(args.reffasta) as referencefile:
	l=referencefile.readlines()
	for i,each in enumerate(l):
		if ">" in each:
			referencebarcodes[each.strip()[1:].split("_all")[0]]=l[i+1]

minionbarcodes={}
with open(args.minionfasta) as minionfile:
	l=minionfile.readlines()
	for i,each in enumerate(l):
		if ">" in each:
			minionbarcodes[each.strip()[1:].split("_all")[0]]=l[i+1]
print referencebarcodes
with open(args.outfile,'w') as statsfile:
	statsfile.write('ID\tn_diff\tn_gaps\tlen_aln\n')
	for k in minionbarcodes.keys():
		try:
			with open(args.outdir+"/"+k+"_withref",'w') as outfile:
				outfile.write(">ref_"+k+'\n'+referencebarcodes[k]+">minion_"+k+'\n'+minionbarcodes[k])
			cmd="mafft --adjustdirection "+args.outdir+"/"+k+"_withref > "+args.outdir+"/"+k+"_withrefaln"
			os.system(cmd)
			reformat(args.outdir+'/'+k+"_withrefaln",args.outdir+"/"+k+"_withrefaln_reformat")
			with open(args.outdir+"/"+k+"_withrefaln_reformat") as infile:
				l=infile.readlines()
				minionbarcode=l[3].strip().upper()
				refbarcode=l[1].strip().upper()
				startpos1=0
				endpos1=0
				startpos2=0
				endpos2=0
				for n,bp in enumerate(minionbarcode):
					if bp!="-":
						startpos1=n
						break
				for n,bp in enumerate(minionbarcode[:: - 1]):
					if bp !="-":
						endpos1 = len(minionbarcode) - n 
						break
				for n,bp in enumerate(refbarcode):
					if bp!="-":
						startpos2=n
						break
				for n,bp in enumerate(refbarcode[:: - 1]):
					if bp !="-":
						endpos2 = len(refbarcode) - n 
						break
				startpos=max([startpos1,startpos2])
				endpos=min([endpos1,endpos2])
				if startpos<0:
					startpos=0
				with open(args.outdir+"/"+k+"_withrefaln_reformat_trimmednoext",'w') as outfile:
					for i,j in enumerate(l):
						if ">" in j:
							outfile.write(j+l[i+1].strip()[startpos:endpos]+'\n')
				pdist,gaps=measurepair_distancem(minionbarcode[startpos:endpos],refbarcode[startpos:endpos])
				statsfile.write(k+'\t'+str(pdist)+'\t'+str(gaps)+'\t'+str(len(minionbarcode[startpos:endpos]))+'\n')
		except KeyError:
			pass
