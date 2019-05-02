import sys,os,fileinput,re,argparse
from distutils.spawn import find_executable

parser=argparse.ArgumentParser(description='Script for assessing uncorrected barcodes using sanger illumina references')
parser.add_argument('-m','--minionfasta',help='Path to input mafft corrected barcode fasta file',dest="minionfasta",required=True)
parser.add_argument('-r','--reffasta',help='Path to input racon corrected barcode fasta file',dest="reffasta",required=True)
parser.add_argument('-t','--tempdir',help='Path to input racon corrected barcode fasta file',dest="outdir",required=True)
parser.add_argument('-o','--outfile',help='Path to output file containing statistics',dest="outfile",required=True)

args=parser.parse_args()
def is_inpath(toolname):
	return find_executable(toolname) is not None
if is_inpath("dnadiff")==False:
	print "dnadiff is not in path"
	sys.exit()

if is_inpath("mafft")==False:
	print "mafft is not in path"
	sys.exit()
if os.path.isdir(args.outdir):
	print "temp directory exists, please delete or name another"
	sys.exit()
else:
	os.system("mkdir "+args.outdir)

def cleanup_dnadiff(inline):
	return filter(None,inline.split(" "))[2]
def get_cor_lines(infilename,outfile):
	print infilename
	with open(infilename) as infile:
		l=infile.readlines()
		sample=os.path.basename(infilename)
		outfile.write(sample)
		for each in l:
			if "AlignedBases" in each:
				len_aln=cleanup_dnadiff(each.strip())
			if "TotalIndels" in each:
				n_gaps=cleanup_dnadiff(each.strip())
			if "TotalSNPs" in each:
				n_diff=cleanup_dnadiff(each.strip())
		outfile.write('\t'+n_diff+'\t'+n_gaps+'\t'+len_aln+'\n')

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

reformat(args.reffasta,args.reffasta+"_reformat")
reformat(args.minionfasta,args.minionfasta+"_reformat")

referencebarcodes={}

with open(args.reffasta+"_reformat") as referencefile:
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
			
with open(args.outfile,'w') as statsfile:
	statsfile.write('ID\tn_diff\tn_gaps\tlen_aln\n')
	for k in minionbarcodes.keys():
		try:
			with open(args.outdir+"/"+k+"_minion",'w') as outfile:
				outfile.write(">minion_"+k+'\n'+minionbarcodes[k])
			with open(args.outdir+"/"+k+"_reference",'w') as outfile:
				outfile.write(">ref_"+k+'\n'+referencebarcodes[k])
			cmd="dnadiff "+args.outdir+"/"+k+"_reference "+args.outdir+"/"+k+"_minion -p "+args.outdir+"/"+k+"_dnadiff "
			os.system(cmd)
			get_cor_lines(args.outdir+"/"+k+"_dnadiff.report",statsfile)
		#	statsfile.write(k+'\t'+str(pdist)+'\t'+str(gaps)+'\t'+str(ncount)+'\n')
		except KeyError:
			pass
