import sys, os,argparse
parser=argparse.ArgumentParser(description='Script for getting fastq after demultiplexing')
parser.add_argument('-fq','--infastq',help='Path to input fastq file',dest="infastq",required=True)
parser.add_argument('-dr','--demreads',help='Path to folder containing demultiplexed files',dest="demreads",required=True)
parser.add_argument('-se','--startend',help='file containing start and end position of primers, \"COIpred file\"',dest="sefile",required=True)
parser.add_argument('-o','--outdir',help='output directory',dest="outdir",required=True)
args=parser.parse_args()

dirlist=os.listdir(args.demreads)

if os.path.isdir(args.outdir):
	print "output directory exists, please delete or name another"
	sys.exit()
else:
	os.system("mkdir "+args.outdir)
def retrieve_partseq_fn(seq,range,direction):
	if direction==1:
		try:
			return seq[range[0]-1:range[1]]
		except IndexError:
			pass
	
def retrieve_partseq(inputseqs,inputfile):
	newdict={}
	with open(inputfile) as startendfile:
		startendlines=startendfile.readlines()
		for each in startendlines:
			m=each.strip().split(' ')
			outseq=retrieve_partseq_fn(inputseqs[int(m[0])*2][1],[int(m[2]),int(m[3])],int(m[1]))
			outquals=retrieve_partseq_fn(inputseqs[int(m[0])*2][3],[int(m[2]),int(m[3])],int(m[1]))
			newdict[int(m[0])*2]=inputseqs[int(m[0])*2][0]+outseq+'\n'+inputseqs[int(m[0])*2][2]+outquals+'\n'
	return newdict

fqlines={}
with open(args.infastq) as infile:
	l=infile.readlines()
	for i,j in enumerate(l):
		if i%4 == 0:
			fqlines[i]=[j,l[i+1],l[i+2],l[i+3]]
			
newfqlines=retrieve_partseq(fqlines,args.sefile)

print newfqlines.keys()[0:100]
for fname in dirlist:
	with open(args.demreads+"/"+fname) as infile:
		with open(args.outdir+"/"+fname+"stq",'w') as outfile:
			l=infile.readlines()
			for i,j in enumerate(l):
				if ">" in j:
					id=int(j.strip().split(";")[0].replace(">",""))*2
					outfile.write(newfqlines[id])
