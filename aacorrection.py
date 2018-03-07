import sys,os,fileinput,argparse
import correction_funcs as corr
from distutils.spawn import find_executable
from Bio import Entrez

parser=argparse.ArgumentParser(description='Script for correcting barcodes using conserved amino acids')
parser.add_argument('-b','--barcodes',help='Path to input barcode fasta file',dest="infasta",required=True)
parser.add_argument('-o','--outfile',help='outfile file name',dest="outfile",required=True)
parser.add_argument('-p','--threads',help='number of threads for BLAST,default=4',dest="threads",default="4")
blastinputgroup = parser.add_mutually_exclusive_group(required=True)
blastinputgroup.add_argument('-bo','--blastout',help='Path to blast output file, outputformat 6',dest="blastoutfile")
parser.add_argument('-bf','--blastfasta',help='Path to fasta file containing sequences of BLAST hits, required if -bo or --blastout is given',dest="blastaccfile")
blastinputgroup.add_argument('-d','--db',help='Path to nucleotide database with database prefix, if local copy is unavailable you can try typing \'nt -remote\'. note that remote has not been extensively tested and is slower',dest="path_to_db")
parser.add_argument('-a','--amb',help='proportion of ambiguities allowed per barcode, default=0.01',dest="nambs",default="0.01")
parser.add_argument('-l','--minlen',help='exclude barcodes shorter than this length, default=640',dest="minlen",default="640")
parser.add_argument('-L','--maxlen',help='exclude barcodes longer than this length, default=670',dest="maxlen",default="670")
parser.add_argument('-c','--congaps',help='exclude sequences containing any gap of length >= value, default=5',dest="congaps",default="5")
parser.add_argument('-n','--namino',help='number of flanking amino acids around the gap used for correction, default=3',dest="namino",default="3")
parser.add_argument('-g','--gencode',help='genetic code https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi, default=5, invertebrate mitochondrial',dest="gencode",default="5")
parser.add_argument('-e','--evalue',help='e-value for BLAST search, default=1e-5',dest="evalue",default="1e-5")
parser.add_argument('-H','--hplen',help='minimum homopolymer length, default=2',dest="hplen",default="2")
parser.add_argument('-s','--support',help='minimum support for indel in references. Reducing this increases chances of errors and improper detection of reading frames',dest="support",type=int,default=2)

args=parser.parse_args()
if args.blastoutfile is not None:
	if args.blastaccfile is None:
		parser.error("--blastout/-bo requires --blastfasta/-bf.")
		
def is_inpath(toolname):
	return find_executable(toolname) is not None

if is_inpath("blastn")==False:
	print "blastn is not in path"
	sys.exit()
if is_inpath("blastdbcmd")==False:
	print "blastdbcmd is not in path"
	sys.exit()
if is_inpath("mafft")==False:
	print "mafft is not in path"
	sys.exit()
if os.path.exists(args.infasta+"_w10blasthits"):
	print "this script needs to create a folder named input fasta filename+_w10blasthits, this folder exists, please delete existing or rename old folder"
	sys.exit()
## This script is hard-coded to work with 10 BLAST hits where it terminates for any sequence if <5 are available. 

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

def reverse_comp(inseq):
	comp_seq=''
	comps={'A':'T','T':'A','G':'C','C':'G','N':'N','M':'K','R':'Y','W':'W','S':'S','Y':'R','K':'M','V':'B','H':'D','D':'H','B':'V'}
	for nucl in inseq:
		comp_seq=comp_seq+comps[nucl]
	return comp_seq[::-1]

def get_acc_list(infilename):
	accdict={}
	for line in fileinput.input([infilename]):
		m=line.split('\t')[1]
		accdict[m]=''
	fileinput.close()
	with open(infilename+".list",'w') as acclistfile:
		acclistfile.write('\n'.join(accdict.keys()))

def build_accdict(infastaname):
	# can be improved.
	reformat(infastaname,infastaname+"reformat.fasta")
	accdict={}
	with open(infastaname+"reformat.fasta") as accfile:
		l=accfile.readlines()
		for i,line in enumerate(l):
			if ">" in line:
				accdict[line.split(" ")[0][1:]]=l[i+1].strip()
	return accdict

def build_blastoutdict(blastoutfile):
	iddict,lendict={},{}
	for line in fileinput.input([blastoutfile]):
		m=line.split('\t')
		try:
			iddict[m[0].split(";")[0]].append(m)
			lendict[m[0].split(";")[0]].append(int(m[7])-int(m[6])-1)
		except KeyError:
			iddict[m[0].split(";")[0]]=[m]
			lendict[m[0].split(";")[0]]=[m]
	iddict_sorted={}
	lendict_sorted={}
	for id in iddict.keys():
		newlist=sorted(iddict[id], key=lambda x:float(x[2]))
		newlist2=sorted(newlist, key=lambda x:float(x[7])-float(x[6]))
		iddict_sorted[id]=newlist[::-1]
		lendict_sorted[id]=newlist2[::-1]
	return iddict_sorted,lendict_sorted

def get_top5(iddict_sorted,lendict_sorted):
	diffseqset={}	
	startend={}
	for k in iddict_sorted.keys():
		minstart,maxend =700,0
		startend[k]=[minstart,maxend]
		diffseqset[k]=[]
		for m in iddict_sorted[k]:
			if len(diffseqset[m[0].split(";")[0]])<5:
				try:
					seq=accdict[m[1]]
					if int(m[8])<int(m[9]):
						subseq=seq[int(m[8])-1:int(m[9])]
						if "N" not in subseq:
							if subseq not in diffseqset[m[0].split(";")[0]]:
								diffseqset[m[0].split(";")[0]].append(subseq)
								if int(m[6])<minstart:
									minstart=int(m[6])
								if int(m[7])>maxend:
									maxend=int(m[7])
					else:
						subseq=reverse_comp(seq[int(m[9])-1:int(m[8])])
						if "N" not in subseq:
							if subseq not in diffseqset[m[0].split(";")[0]]:
								diffseqset[m[0].split(";")[0]].append(subseq)
								if int(m[6])<minstart:
									minstart=int(m[6])
								if int(m[7])>maxend:
									maxend=int(m[7])
				except KeyError:
					pass
		for m in lendict_sorted[k]:
			if len(diffseqset[m[0].split(";")[0]])<10:
				try:
					seq=accdict[m[1]]
					if int(m[8])<int(m[9]):
						subseq=seq[int(m[8])-1:int(m[9])]
						if "N" not in subseq:
							if subseq not in diffseqset[m[0].split(";")[0]]:
								diffseqset[m[0].split(";")[0]].append(subseq)
								if int(m[6])<minstart:
									minstart=int(m[6])
								if int(m[7])>maxend:
									maxend=int(m[7])
					else:
						subseq=reverse_comp(seq[int(m[9])-1:int(m[8])])
						if "N" not in subseq:
							if subseq not in diffseqset[m[0].split(";")[0]]:
								diffseqset[m[0].split(";")[0]].append(subseq)
								if int(m[6])<minstart:
									minstart=int(m[6])
								if int(m[7])>maxend:
									maxend=int(m[7])
				except KeyError:
					pass
		startend[k]=[minstart,maxend]
	return diffseqset,startend

def get_top5_blast(infilename,outdir,diffseqset,startend):
	mseqdict={}
	with open(args.infasta) as minionfile:
		minionlines=minionfile.readlines()
		for i,line in enumerate(minionlines): 
				if ">" in line:
					mseqdict[line.strip().replace(">","").split(";")[0]]=corr.remove_ext_Ns(minionlines[i+1].strip().upper())
				
		for each in diffseqset.keys():
			try:
				with open(outdir+"/"+each.split(';')[0],'w') as outfile:
					outfile.write(">"+each+'\n' +mseqdict[each]+'\n')
					for i,seq in enumerate(diffseqset[each]):
						outfile.write(">" +str(i)+'\n'+seq+'\n')
			except KeyError:
				pass

def revert_hp(infile):
	o=open(infile)
	clustset={}
	seqclust={}
	seqset=[]
	l=o.readlines()
	for i,j in enumerate(l):
		if ">" in j:
			if i!=0:
				seqset.append(l[i+1].strip().upper().replace("-","N"))
	nseq=corr.change_ext_gaps(l[1].strip().upper()).replace("-","")
	n=0
	while n<len(nseq)-1:
		if nseq[n]=="N" and nseq[n+1]=="N":
			strict_consensus,k1,k2=corr.retrieve_multns(nseq,seqset,n,args.hplen)
			nseq=nseq[:k1]+strict_consensus+nseq[k2:]
			n+=len(strict_consensus)
		else:
			n+=1
	return l[0]+nseq.replace("?","").replace("-","")+'\n'
	


#def trim_aln(infilename,outfilename):
#	seqset=[]
#	trimmed_seqset={}
#	l=infile.readlines()
#	seqset = [
#	for i,j in  enumerate(l):
#		filled_sequence,start, end = corr.get_start_end(l[i+1].strip())
#		seqset.append(corr.get_start_end(l[i+1].strip()))
#		if i>1:
#			if len(l[i+1].strip().replace("-",""))>maxlen:
#				maxlen=len(l[i+1].strip().replace("-",""))
#				refseq=i/2
#	for i,j in  enumerate(l):
#		if i
#		if ">" in j:
#			trimmed_seqset[
#	startpos=0
#	endpos=0
#	for n,bp in enumerate(l[1].strip()):
#		if bp!="-":
#			startpos=n
#			break
#	for n,bp in enumerate(l[1].strip()[:: - 1]):
#		if bp !="-":
#			endpos = len(l[1].strip()) - n 
#			break
#	print startpos,endpos
#	with open(sys.argv[1]+"/"+fname+"_reformat_trimmed",'w') as outfile:
#		for i,j in enumerate(l):
#			if ">" in j:
#				outfile.write(j+l[i+1].strip()[startpos:endpos]+'\n')
		

def get_fasta_from_acclist(listfile,outfile):
	with open(outfile,'w') as outfile:
		with open(listfile) as inputfile:
			l=inputfile.readlines()
			for line in l:
				m=Entrez.efetch(db="nucleotide", id=line.strip(), rettype="fasta", retmode="text")
				for ind in m.readlines():
					outfile.write(ind)				


if __name__ == '__main__':
# BLAST of uncorrected barcodes to nt
	if args.path_to_db is not None:
		if args.path_to_db!="nt -remote":
			os.system('blastn -outfmt 6 -query '+ args.infasta+ ' -db ' + args.path_to_db + ' -evalue ' + args.evalue + ' -num_threads ' + args.threads + ' -out '+args.infasta+'_megablast')
		else:
			os.system('blastn -outfmt 6 -query '+ args.infasta+ ' -db ' + args.path_to_db + ' -evalue ' + args.evalue + ' -out '+args.infasta+'_megablast')			
# Getting list of accession of subject sequences
		get_acc_list(args.infasta+"_megablast")

# Getting nucleotide sequences corresponding to hits
		if args.path_to_db!="nt -remote":
			os.system('blastdbcmd -entry_batch '+ args.infasta+'_megablast.list -db ' + args.path_to_db + ' -out ' + args.infasta+"_megablast.fasta")
		else:
			get_fasta_from_acclist(args.infasta+'_megablast.list',args.infasta+"_megablast.fasta")
# Building sequence dictionary of the best hits
		accdict=build_accdict(args.infasta+'_megablast.fasta')

# Reading BLAST output file to store sorted information based on identities and hit length
		iddict,lendict=build_blastoutdict(args.infasta+"_megablast")
	else:
		accdict=build_accdict(args.blastaccfile)
		iddict,lendict=build_blastoutdict(args.blastoutfile)

# Parsing BLAST hits to retain only unique sequences and five best identity and hit lengths
	diffseqset,startend=get_top5(iddict,lendict)

# Make directory for storing files needed for correction
	if not os.path.exists(args.infasta+"_w10blasthits"):
		os.mkdir(args.infasta+"_w10blasthits")

# Generating files per barcode and corresponding 10 hits
	get_top5_blast(args.infasta,args.infasta+"_w10blasthits",diffseqset,startend)

# Executing correction
	dirlist=os.listdir(args.infasta+"_w10blasthits")
	outfile=open(args.outfile,'w')
	for fname in dirlist:
		try:
			print fname
			os.system("mafft --globalpair --adjustdirection " + args.infasta+"_w10blasthits/"+ fname +"> " + args.infasta+"_w10blasthits/"+ fname + "_mafft")
			reformat(args.infasta+"_w10blasthits/"+ fname + "_mafft",args.infasta+"_w10blasthits/"+ fname + "_mafft_reformat")
		#	trim_aln(args.infasta+"_w10blasthits/"+ fname + "_mafft_reformat",args.infasta+"_w10blasthits/"+ fname + "_mafft_reformat_aln"
			with open(args.infasta+"_w10blasthits/"+ fname+ "_mafft_reformat") as infile:
				l=infile.readlines()
				seqid=l[0]
				if l[1].upper().count("N")<=len(l[1].strip().replace("-",""))*float(args.nambs) and len(l[1].replace("-","").strip())<int(args.maxlen) and len(l[1].replace("-","").strip())>int(args.minlen) and "-"*int(args.congaps) not in l[1]:
					seqset,retainflag=corr.check_alignment(l,int(args.support))					
					maxlen=0
			#		print seqset
					for i,j in enumerate(seqset):
						if i>0:
			#				print i
							if len(corr.remove_ext_Ns(j.strip().replace("-","")))>maxlen:
								maxlen=len(corr.remove_ext_Ns(j.strip().replace("-","")))
								refseq=i
					mseq,refseqset,orseq=corr.translate_corframe(seqset,refseq,int(args.gencode))
					newseq,orseq,nrefseqset=corr.runcorrection(mseq,refseqset,orseq,int(args.gencode),int(args.namino))
					with open(args.infasta+"_w10blasthits/" + fname+"_wambs",'w') as toutfile:
						toutfile.write(l[0].strip()+"_corr"+'\n'+newseq+'\n'+l[0]+orseq+'\n')
						for i,each in enumerate(nrefseqset):
							toutfile.write(">"+str(i)+'\n'+each+'\n')

					fbarcode=revert_hp(args.infasta+"_w10blasthits/" + fname+"_wambs")
					outfile.write(fbarcode)


#					seqset,retainflag=corr.check_alignment(l,int(args.support))					
#					if retainflag==False:
#						with open(args.infasta+"_w10blasthits/"+ fname+ "2",'w') as outfile:
#							outfile.write(seqid+seqset[0]+'\n')
#							for i,j in enumerate(seqset[1:]):
#								outfile.write(">"+str(i)+"\n"+j+'\n')
#						os.system("mafft --globalpair --adjustdirection " + args.infasta+"_w10blasthits/"+ fname +"2> " + args.infasta+"_w10blasthits/"+ fname + "2_mafft")
#						reformat(args.infasta+"_w10blasthits/"+ fname + "2_mafft",args.infasta+"_w10blasthits/"+ fname + "2_mafft_reformat")
#						with open(args.infasta+"_w10blasthits/"+ fname+ "2_mafft_reformat") as infile:
#							l=infile.readlines()
#							seqset,retainflag=corr.check_alignment(l,int(args.support))
#					if len(seqset)<6:
#						break
#					mseq,refseqset,orseq,flag=corr.translate_corframe(seqset,int(args.gencode))
#					if flag==True:
#						newseq,orseq,nrefseqset=corr.runcorrection(mseq,refseqset,orseq,int(args.gencode),int(args.namino),int(args.congaps))
#						with open(args.infasta+"_w10blasthits/" + fname+"_wambs",'w') as toutfile:
#							toutfile.write(l[0].strip()+"_corr"+'\n'+newseq+'\n'+l[0]+orseq+'\n')
#							for i,each in enumerate(nrefseqset):
#								toutfile.write(">"+str(i)+'\n'+each+'\n')
#						fbarcode=revert_hp(args.infasta+"_w10blasthits/" + fname+"_wambs")
#						outfile.write(fbarcode)
		except IndexError:
			pass
