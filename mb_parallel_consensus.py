#!/usr/bin/env python 
import sys,os,fileinput,re,time,numpy,argparse,random,fnmatch,multiprocessing
from distutils.spawn import find_executable
from multiprocessing import Pool
from collections import Counter 


def is_inpath(toolname):
	return find_executable(toolname) is not None

if is_inpath("glsearch36")==False:
	print "glsearch36 is not in path"
	sys.exit()

def reverse_comp(inseq):
	comp_seq=''
	flag=True
	ambiguity_codes=['N','K','M','R','Y','S','W','B','V','H','D','N','X']
	for base in ambiguity_codes:
		if base in inseq:
			flag==False
	if flag==True:	
		for nucl in inseq:
			if nucl=='A':
				comp_seq=comp_seq+'T'
			elif nucl=='T':
				comp_seq=comp_seq+'A'
			elif nucl=='G':
				comp_seq=comp_seq+'C'
			elif nucl=='C':
				comp_seq=comp_seq+'G'
	elif flag==False:
		comp_seq=0
	return comp_seq[::-1]

def allunique(i,o):
	with open(i,'r') as infile:
		l=infile.readlines()
		seqdict={}
		poslist=[]
		for i,j in enumerate(l):
			if ">" in j:
				seqdict[j.strip().split('>')[1]]=l[i+1].strip()
		seqlist=seqdict.values()
		seqlistuniques={}
		for seq in seqlist:
			if seq not in seqlistuniques.keys():
				if reverse_comp(seq) not in seqlistuniques.keys():
					seqlistuniques[seq]=0
		for each in seqlist:
			try:
				seqlistuniques[each]=seqlistuniques[each]+1
			except KeyError:
				seqlistuniques[reverse_comp(each)]=seqlistuniques[reverse_comp(each)]+1
		with open(o,'w') as outfile:
			for i,each in enumerate(seqlistuniques.keys()):
				outfile.write(">Unique_"+str(i)+"; "+str(seqlistuniques[each])+'\n'+each+'\n')
			outfile.close()
def consensus(indict,perc_thresh,abs_thresh):
	if len(indict.keys())>=abs_thresh:
		poslist=[]
		n=0
		while n<len(indict.values()[0]):
			newlist=[]
			for each in indict.values():
				newlist.append(each[n])
			poslist.append(newlist)
			n+=1
		sequence=[]
		countpos=1
		for character in poslist:
			charcounter=Counter(character)
			baseset={}
			for k,v in charcounter.iteritems():
				percbp=float(v)/float(len(character))
				if percbp>perc_thresh:
					baseset[k]=v
			if len(baseset)==0:
				bp='N'
			if len(baseset)==1:
				bp=baseset.keys()[0]
			if len(baseset)>1:
				bp='N'
			sequence.append(bp)
			countpos+=1
		return ''.join(sequence)
	else:
		return ''
							
def callconsensus(i,perc_thresh,abs_thresh,o,name):
	with open(i) as infile:
		l=infile.readlines()
		seqdict,poslist={},[]
		for i,j in enumerate(l):
			if ">" in j:
				poslist.append(i)
		for i,j in enumerate(poslist):
			ambcounts=0
			k01=l[j].strip().split('>')[1]
			if i!=len(poslist)-1:
				k3=l[j+1:poslist[i+1]]
			if i==len(poslist)-1:
				k3=l[j+1:]
			k4=''.join(k3).replace('\n','')
			seqdict[k01]=k4
		conseq=consensus(seqdict,perc_thresh,abs_thresh)
		conseq=conseq.replace("-",'')
		if len(conseq)!=0:
			with open(o,'w') as outfile:
				outfile.write(">"+name.split(".")[0]+";"+str(len(conseq))+";"+str(len(seqdict.keys()))+'\n'+conseq+'\n')
def subset_randomly (infile,outfile,n):
	with open(infile) as fulldata:
		with open(outfile,'w') as subsetdata:
			l=fulldata.readlines()
			if len(l)/2<=n:
				subsetdata.write("".join(l))
			else:
				subsetlist=[x*2 for x in random.sample(range(len(l)/2),n)]
				for i in subsetlist:
					subsetdata.write(l[i]+l[i+1])


def runparts(name):
#merging sequences to uniques
	if args.maxdepth!=0:
		subset_randomly(args.indir+"/demultiplexed/"+name,args.indir+"/demultiplexed_"+str(args.maxdepth)+"/"+name,int(args.maxdepth))
	else:
		os.system("cp "+args.indir+"/demultiplexed/"+name +" "+args.indir+"/demultiplexed_"+str(args.maxdepth)+"/"+name)
	allunique(args.indir+"/demultiplexed_"+str(args.maxdepth)+"/"+name, args.indir+"/demultiplexed_uniqs/"+name+".uniq")
	os.system('mafft --adjustdirection --op 0 '+args.indir+"/demultiplexed_"+str(args.maxdepth)+"/"+ name + ' > ' +args.indir+'/demultiplexed_uniqs_mafft/'+name + '_aln.fasta')

parser=argparse.ArgumentParser(description='Script for obtaining barcodes')
parser.add_argument('-i','--indir',help='output folder path of mb_parallel_demultiplex',dest="indir",required=True)
parser.add_argument('-D','--maxdepth',type=int,help='set max depth per coverage to improve speed, default 100, must be >2',dest="maxdepth",default=100)
parser.add_argument('-t','--threads',help='number of threads',dest="threads",default="4")
args=parser.parse_args()

nthreads=int(args.threads)


#	#aligning and consensus calling
os.system("mkdir "+args.indir+"/"+"demultiplexed_"+str(args.maxdepth))
os.system("mkdir "+args.indir+"/demultiplexed_uniqs")
os.system("mkdir "+args.indir+"/demultiplexed_uniqs_mafft")
os.system("mkdir "+args.indir+"/demultiplexed_uniqs_mafft_consensus")
dirlist=os.listdir(args.indir+"/demultiplexed")
partlist=fnmatch.filter(os.listdir(args.indir+"/demultiplexed"),"*")
lens=range(0,len(partlist),nthreads)
n=0
for i,j in enumerate(lens[:-1]):
	print "Processing parts "+str(n) +" to "+str(n+nthreads-1)
	parts=partlist[j:lens[i+1]]
	print parts
	p=Pool(nthreads)
	p.map(runparts,parts)
	p.close()
	n+=nthreads
parts=partlist[lens[-1]:]
p=Pool(nthreads)
p.map(runparts,parts)
p.close()
for name in dirlist:
	callconsensus(args.indir+"/demultiplexed_uniqs_mafft/"+name+"_aln.fasta",0.5,5,args.indir+"/demultiplexed_uniqs_mafft_consensus/"+name.split(".")[0]+"_consensus.fa",name)


os.system("cat "+args.indir+"/demultiplexed_uniqs_mafft_consensus/* > " + args.indir+"/all_barcodes.fa")
