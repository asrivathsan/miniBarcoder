#!/usr/bin/env python 
# minibarcoder: a pipeline for BLAST based read identifications

# Copyright 2017 Amrita Srivathsan

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

# Contact: asrivathsan@gmail.com


import sys,os,fileinput,re,time,numpy,argparse,random,fnmatch,multiprocessing,itertools
from distutils.spawn import find_executable
from multiprocessing import Pool
from collections import Counter
parser=argparse.ArgumentParser(description='Script for obtaining barcodes')
parser.add_argument('-f','--infasta',help='Path to input fasta file',dest="infasta",required=True)
parser.add_argument('-d','--demfile',help='Path to demultiplexing file',dest="demfile",required=True)
parser.add_argument('-o','--outdir',help='set output directory path',dest="outdir",required=True)
parser.add_argument('-l','--minlen',help='exclude barcode sequences identified that are shorter than specified length',dest="minlen",required=True)
parser.add_argument('-L','--maxlen',help='exclude barcode sequences identified that are shorter than specified length',dest="maxlen",default="2000")
parser.add_argument('-mm','--mismatch',type=int,choices=range(0,6),help='number of mismatches allowed in tags, must be <=5, default 2',dest="mismatch",default=2)
parser.add_argument('-e','--evalue',help='evalue for primer search using glsearch36,default 1e+6',dest="evalue",default="1000000")
parser.add_argument('-g','--gaps',help='number of gaps allowed for primer identification, default 5',dest="gaps",default="5")
parser.add_argument('-D','--maxdepth',type=int,help='set max depth per coverage to improve speed, default 100, must be >2',dest="maxdepth",default=100)
parser.add_argument('-t','--threads',help='number of threads for glsearch and mafft',dest="threads",default="4")
parser.add_argument('-s','--splitlength',type=int,help='length to split reads for ligation of products',dest="slen",default="650")

args=parser.parse_args()




def is_inpath(toolname):
	return find_executable(toolname) is not None

if is_inpath("glsearch36")==False:
	print "glsearch36 is not in path"
	sys.exit()


if os.path.isdir(args.outdir):
	print "output directory exists, please delete or name another"
 	sys.exit()
else:
	os.system("mkdir "+args.outdir)
def rev_comp(inseq,comps):
	comp_seq=''
	for nucl in inseq:
		comp_seq=comp_seq+comps[nucl]
	return comp_seq[::-1]
def retrieve_partseq_fn(seq,range,direction):
	if direction==1:
		try:
			return seq[range[0]-1:range[1]]
		except IndexError:
			pass
	elif direction==2:
		try:
			return rev_comp(seq[range[0]-1:range[1]],comps)
		except IndexError:
			pass
	
def parseglsearch_primermatch (infilename,seqdict,start,end):
	with open(infilename) as infile:
		l=infile.readlines()
		hitseqs,hitseqs2,inseqs,inseqs2={},{},{},{}

		for line in l:
			m=line.split('\t')
			if len(m)==12:
				hitseqs[m[1]],hitseqs2[m[1]]=["",""],["",""]
				inseqs[m[1]],inseqs2[m[1]]="",""

		for each in l:
			m=each.split('\t')
			if len(m)==12:
				qse, se = [int(m[6]),int(m[7])] , [int(m[8]),int(m[9])]
				if qse[0]<qse[1]:
					if se[0]-taglen>0:
						hitseqs[m[1]][0]=se[0]-taglen-7
						hitseqs[m[1]][1]=se[0]-1
						inseqs[m[1]]=se[1]+1
				else:
					hitseqs2[m[1]][0]=len(seqdict[m[1]])-end+ se[1]+1
					hitseqs2[m[1]][1]=len(seqdict[m[1]])-end+ se[1]+taglen+7
					inseqs2[m[1]]=len(seqdict[m[1]])-end+ se[0]-1


	with open(infilename+"parsed_f",'w') as outfile1:
		with open(infilename+"start_f","w") as outfile2:
			for each in hitseqs.keys():
				f=[]
				for v in hitseqs[each]:
					if isinstance( v, ( int, long ) ) == True:
						f.append(1)
				if f==[1,1]:
					outfile1.write(">"+each+'\n'+retrieve_partseq_fn(seqdict[each],hitseqs[each],1)+'\n')
					outfile2.write(">"+each+'\n'+seqdict[each][inseqs[each]:]+'\n')
			for each in hitseqs2.keys():
				f=[]
				for v in hitseqs2[each]:
					if isinstance( v, ( int, long ) ) == True:
						f.append(1)
				if f==[1,1]:
					outfile1.write(">"+each+'\n'+retrieve_partseq_fn(seqdict[each],hitseqs2[each],2)+'\n')
					outfile2.write(">"+each+'\n'+rev_comp(seqdict[each][:inseqs2[each]],comps)+'\n')

def parseglsearch_primermatchR(infilename,seqdict,start,end):
	with open(infilename) as infile:
		l=infile.readlines()
		hitseqs,hitseqs2,inseqs,inseqs2={},{},{},{}

		for line in l:
			m=line.split('\t')
			if len(m)==12:
				hitseqs[m[1]],hitseqs2[m[1]]=["",""],["",""]
				inseqs[m[1]],inseqs2[m[1]]="",""

		for each in l:
			m=each.split('\t')
			if len(m)==12:
				qse, se = [int(m[6]),int(m[7])] , [int(m[8]),int(m[9])]
				if qse[0]<qse[1]:
					if se[0]-taglen>0:
						hitseqs[m[1]][0]=se[0]-taglen-7
						hitseqs[m[1]][1]=se[0]-1
						inseqs[m[1]]=se[1]+1
				else:
					hitseqs2[m[1]][0]=len(seqdict[m[1]])-end+ se[1]+1
					hitseqs2[m[1]][1]=len(seqdict[m[1]])-end+ se[1]+taglen+7
					inseqs2[m[1]]=len(seqdict[m[1]])-end+ se[0]-1

	with open(infilename+"parsed_r",'w') as outfile1:
		with open(infilename+"endr","w") as outfile2:
			for each in hitseqs2.keys():
				f=[]
				for v in hitseqs2[each]:
					if isinstance( v, ( int, long ) ) == True:
						f.append(1)
				if f==[1,1]:
					outfile1.write(">"+each+'\n'+retrieve_partseq_fn(seqdict[each],hitseqs2[each],2)+'\n')
					outfile2.write(">"+each+'\n'+seqdict[each][:inseqs2[each]]+'\n')

def parse_by_id(InFileName,gapCutoff):
	PassedHits=[] 

	# Hits that meet the overlap requirement
	
	for line in fileinput.input([InFileName]):
		if "\t" in line:
			HitValues=line.split('\t')
			if int(HitValues[5])<=gapCutoff:
				PassedHits.append(line)
			
	fileinput.close()
	IdsWithHits={} 

	# Dictionary of unique ids that meet the overlap requirement
	for hit in PassedHits:
		HitValues=hit.split('\t')
		IdsWithHits[HitValues[1]]=[]

	for hit in PassedHits:
		HitValues=hit.split('\t')
		IdsWithHits[HitValues[1]].append(hit)
	
	FilteredHits={} 
	# Dictionary of unique ids such that the hit with best ID is retained.

	for ID in IdsWithHits.keys():
		FilteredHits[ID]=[]

	for ID in IdsWithHits.keys():
		IDsublist=[]
		Gaplist=[]
		for hit in IdsWithHits[ID]:
			HitValues=hit.split('\t')
			IdentityToSubject=HitValues[4]
			Identity=IdentityToSubject.strip()
			IDsublist.append(float(IdentityToSubject))
		BestIdentity=min(IDsublist)
		for hit in IdsWithHits[ID]:
			HitValues=hit.split('\t')
			IdentityToSubject=HitValues[4]
			IdentityToSubject=IdentityToSubject.replace('\n','')
			if BestIdentity==float(IdentityToSubject):
				Gaplist.append(float(HitValues[5]))
		Bestgaps=min(Gaplist)
		for hit in IdsWithHits[ID]:
			HitValues=hit.split('\t')
			IdentityToSubject=HitValues[4]
			IdentityToSubject=IdentityToSubject.replace('\n','')
			if BestIdentity==float(IdentityToSubject) and Bestgaps==float(HitValues[5]):
				FilteredHits[ID].append(hit) 
	outfile=open(InFileName+".parsed.lencutoff"+str(5),'w')
	for ID in FilteredHits.keys():
		for hit in FilteredHits[ID]:
			outfile.write(hit)

	outfile.close()

def crmutant_m2(tfile,nbp):
	tagdict={}
	sampledict={}
	with open(tfile) as tagfile:
		t=tagfile.readlines()
		for each in t:
			tagdict[each.split(',')[1]]=''
			tagdict[each.split(',')[2]]=''
	counter=1
	typedict={}
	for each in tagdict.keys():
		tagdict[each]="t"+str(counter)
		typedict[each]=0
		counter+=1
	with open(tfile) as tagfile:
		t=tagfile.readlines()
		for each in t:
			sampledict[(tagdict[each.split(',')[1]],tagdict[each.split(',')[2]])]=each.split(',')[0]
	n=1
	while n<=nbp:
		muttags_fr,newtags_fr=create_all_mutants(tagdict)
#		newtags_fr_counts=Counter(newtags_fr)	

		with open("conflicts",'w') as conflictfile:
			for k in newtags_fr.keys():
				if len(newtags_fr[k])>1:
					conflictfile.write(k+'\t'+",".join(newtags_fr[k])+'\n')
					del muttags_fr[k]
		ntagset=list(set(muttags_fr.keys())-set(tagdict.keys()))
		print len(tagdict),len(muttags_fr),len(ntagset)
		for k in ntagset:
			tagdict[k]=muttags_fr[k]
			typedict[k]=n
		n+=1
	with open(args.outdir+"/temp1.fas",'w') as outfile:
		for k in tagdict.keys():
			outfile.write(k+'\t'+tagdict[k]+'\n')
	return tagdict,muttags_fr,sampledict,typedict

def demultiplex(pf,pr,seqdict,tagdict,muttags_fr,sampledict,typedict):
	indict1,indict2,indict=readprimertagfasta(pf,pr)

	for each in indict.keys():
		try: 
			indict[each][0]=indict1[each]
			indict[each][1]=indict2[each]
		except KeyError:
			del indict[each]
			
	try:
		with open(args.outdir+"/dmpfile_tagfr",'w') as tagfrdmp:	
			for tag in muttags_fr.keys():
				tagfrdmp.write(tag+'\t'+muttags_fr[tag]+'\n')
		findmatch_m2(tagdict, muttags_fr,indict,seqdict,sampledict,typedict)
	except UnboundLocalError:
		findmatch_m2(tagdict, {},indict,seqdict,sampledict,typedict)
def readprimertagfasta(filef,filer):
	indict2,indict1,indict={},{},{}
	with open(filef) as infile1:
		with open(filer) as infile2:
			l1=infile1.readlines()
			l2=infile2.readlines()
			for i,j in enumerate(l1):
				if ">" in j:
					indict1[j.strip().replace(">","")]=l1[i+1].strip()
					indict[j.strip().replace(">","")]=[0,0]
			for i,j in enumerate(l2):
				if ">" in j:
					indict2[j.strip().replace(">","")]=l2[i+1].strip()
					indict[j.strip().replace(">","")]=[0,0]
	return indict1,indict2,indict

def create_all_mutants(tagdict):
	muttags,newtags={},{}
	for tag in tagdict.keys():
		submutant=createsubmutant(tag)
		delmutant=createdelmutant(tag)
		insmutant=createinsmutant(tag)
		mutantset=list(set(submutant)|set(delmutant)|set(insmutant))
		mutantset[:]=[x for x in mutantset if x != tag]
		for mutant in mutantset:
			try:
				if tagdict[tag] not in newtags[mutant]:
					newtags[mutant].append(tagdict[tag])
			except:
				newtags[mutant]=[tagdict[tag]]
			muttags[mutant]= tagdict[tag]
	return muttags,newtags
def createsubmutant(tag):
	newlist=[]
	bpset=["A","T","G","C"]
	for i,v in enumerate(tag):
		for bp in [x for x in bpset if x != v]:
			if i==0:
				newlist.append(bp+tag[i+1:])
			if i>0 and i<len(tag)-1:
				newlist.append(tag[0:i]+bp+tag[i+1:])
			if i==len(tag)-1:
				newlist.append(tag[0:i]+bp)
	return newlist

def createdelmutant(tag):
	newlist=[]
	bpset=["A","T","G","C"]
	for i,v in enumerate(tag):
		for bp in bpset:
			if i>0 and i<len(tag)-1:
				newlist.append(bp+tag[0:i]+tag[i+1:])
			if i==len(tag)-1:
				newlist.append(bp+tag[0:i])
	return newlist
	
def createinsmutant(tag):
	newlist=[]
	bpset=["A","T","G","C"]
	for i,v in enumerate(tag):
		for bp in bpset:
			if i>0 and i<len(tag)-1:
				newlist.append(tag[1:i]+bp+tag[i:])
			if i==len(tag)-1:
				newlist.append(tag[1:]+bp)
	return newlist

def findmatch_m2(taglist,muttags_fr,indict,seqdict,sampledict,typedict):
	dmpfile=open(args.outdir+"/dmpfile",'a')
	idcombs={}
	cumscores={}
	for each in indict.keys():
		try:
			idcombs[each]=(taglist[indict[each][0]], taglist[indict[each][1]])
			cumscore=typedict[indict[each][0]]+typedict[indict[each][1]]
			cumscores[each]=[str(typedict[indict[each][0]]),str(typedict[indict[each][1]]),str(cumscore)]
		except KeyError:
			pass

	matched_idcombs={}
	for each in idcombs.keys():
		try:
			with open(args.outdir+"/"+sampledict[idcombs[each]]+"_all.fa",'a') as outfile:
				outfile.write(">"+each+"_"+"_".join(cumscores[each])+"\n"+seqdict[each]+"\n")
		except KeyError:
			pass
				
	dmpfile.close()	
# checks length of file to ensure the file is not empty
def chklen(inputfile):
	with open(inputfile) as input:
		l=input.readlines()
		return len(l)
def createprimers(seq,ambcodes):
	seqset=[seq[0]]
	for bp in seq[1:]:
		newseqset=[]
		for item in seqset:
			for amb in ambcodes[bp]:
				newseq=item+amb
				newseqset.append(newseq)
		seqset[:]=newseqset
	return seqset

def builddict_sequences(infile):
	seqdict={}
	with open(infile) as inseqs:
		l=inseqs.readlines()
		for i,j in enumerate(l):
			if ">" in j:
				seqdict[j.strip().replace(">","")]=l[i+1].strip()
	return seqdict

def cleantagfile(infile):
	with open(infile) as primertagfile:
		with open(infile+"cleaned",'w') as primertagfile_cleaned:
			taglines=primertagfile.readlines()
			ids=[]
			seqs={}
			for n,line in enumerate(taglines):
				if ">" in line:
					ids.append(line)
					seqs[line]=modseq(taglines[n+1])[-(taglen+1):]

			idcounts=Counter(ids)
			nseqs=0
			for id in ids:
				if idcounts[id]==1:
					nseqs+=1
					primertagfile_cleaned.write(id+seqs[id])

def modseq(seq):
	for bp in ["A","T","G","C"]:
		n=20
		while n>2:
			seq=seq.replace(bp*n,bp*2)
			n-=1
	return seq

def runglsearch2(infile):
	inputseqs=builddict_sequences(args.outdir+"/"+infile)
	with open(args.outdir+"/"+infile+"1",'w') as outfile1:
			with open(args.outdir+"/"+infile+"2",'w') as outfile2:
				for each in inputseqs.keys():
					outfile1.write(">"+each+'p1\n'+inputseqs[each][0:-650]+'\n')
					outfile2.write(">"+each+'p2\n'+inputseqs[each][650:]+'\n')
	runglsearchperfile(infile+"1",100,500)
	runglsearchperfile(infile+"2",100,500)


def runglsearch1(infile):
	runglsearchperfile(infile,100,100)

def runglsearchperfile(infile, start,end):
	inputseqs=builddict_sequences(args.outdir+"/"+infile)
	with open(args.outdir+"/"+infile+"_first100",'w') as outfile1:
			with open(args.outdir+"/"+infile+"_last100",'w') as outfile2:
				for each in inputseqs.keys():
					outfile1.write(">"+each+'\n'+inputseqs[each][0:start]+'\n')
					outfile2.write(">"+each+'\n'+inputseqs[each][-end:]+'\n')
	seqlens=[len(i.strip()) for i in inputseqs.values()]
	searchfilefirst100=infile+"_first100"
	searchfilelast100=infile+"_last100"
	for pf in primerfset:
		with open(args.outdir+"/"+infile+"_primerf.fa",'w') as outfile:
			outfile.write(">primerf\n"+pf+'\n')
		cmd="glsearch36 -n -3 -m 8 -E "+args.evalue+" -T 1 "+args.outdir+"/"+infile+"_primerf.fa"+" "+args.outdir+"/"+searchfilefirst100+ " > "+args.outdir+"/"+searchfilefirst100+"_glsearchF"
		os.system(cmd)
		cmd="glsearch36 -n -i -m 8 -E "+args.evalue+" -T 1 "+args.outdir+"/"+infile+"_primerf.fa"+" "+args.outdir+"/"+searchfilelast100 + " > "+args.outdir+"/"+searchfilelast100+"_glsearchR"
		os.system(cmd)
		os.system("cat " +args.outdir+"/"+searchfilefirst100+"_glsearchF "+args.outdir+"/"+searchfilelast100+"_glsearchR>" +args.outdir+"/"+infile+"_"+pf+"_glsearch1")
		parse_by_id(args.outdir+"/"+infile+"_"+pf+"_glsearch1",int(args.gaps))
		parseglsearch_primermatch(args.outdir+"/"+infile+"_"+pf+"_glsearch1.parsed.lencutoff5",inputseqs,start,end)
		doneseqs=builddict_sequences(args.outdir+"/"+infile+"_"+pf+"_glsearch1.parsed.lencutoff5start_f")
		for id in doneseqs:
			del inputseqs[id]
		if len(inputseqs)==0:
			break
		with open(args.outdir+"/"+infile+"_first100_post"+pf,'w') as tfile:
			for each in inputseqs.keys():
				tfile.write(">" +each+'\n'+inputseqs[each]+'\n')
		with open(args.outdir+"/"+infile+"_last100_post"+pf,'w') as tfile:
			for each in inputseqs.keys():
				tfile.write(">" +each+'\n'+inputseqs[each]+'\n')
		searchfilefirst100=infile+"_first100_post"+pf
		searchfilelast100=infile+"_last100_post"+pf
	cmd="cat "+args.outdir+"/"+infile+"*_glsearch1.parsed.lencutoff5start_f > "+args.outdir+"/"+ infile+"_all_glsearch1.parsed.lencutoff5start_f"
	os.system(cmd)
	cmd="cat "+args.outdir+"/"+infile+"*_glsearch1.parsed.lencutoff5parsed_f > "+args.outdir+"/"+ infile+"_all_glsearch1.parsed.lencutoff5parsed_f"
	os.system(cmd)
	inputseqs=builddict_sequences(args.outdir+"/"+infile+"_all_glsearch1.parsed.lencutoff5start_f")
	with open(args.outdir+"/"+infile+"_glsearch1.parsed.lencutoff5start_f_last100",'w') as outfile2:
		for each in inputseqs.keys():
			outfile2.write(">"+each+'\n'+inputseqs[each][-end:]+'\n')
	searchfilelast100=infile+"_glsearch1.parsed.lencutoff5start_f_last100"
	for pr in primerrset:
		with open(args.outdir+"/"+infile+"_primerr.fa",'w') as outfile:
			outfile.write(">primerr\n"+pr+'\n')
		cmd="glsearch36 -n -i -m 8 -E "+args.evalue+" -T 1 "+args.outdir+"/"+infile+"_primerr.fa"+" "+args.outdir+"/"+searchfilelast100+" > "+args.outdir+"/"+infile+"_"+pr+"_glsearchR"
		os.system(cmd)
		parse_by_id(args.outdir+"/"+infile+"_"+pr+"_glsearchR",int(args.gaps))
		parseglsearch_primermatchR(args.outdir+"/"+infile+"_"+pr+"_glsearchR.parsed.lencutoff5",inputseqs,start,end)	
		doneseqs=builddict_sequences(args.outdir+"/"+infile+"_"+pr+"_glsearchR.parsed.lencutoff5endr")
		for id in doneseqs:
			del inputseqs[id]
		if len(inputseqs)==0:
			break
		with open(args.outdir+"/"+infile+"_glsearch1.parsed.lencutoff5start_f_last100"+"_post"+pr,'w') as tfile:
			for each in inputseqs.keys():
				tfile.write(">" +each+'\n'+inputseqs[each]+'\n')
		searchfilelast100=infile+"_glsearch1.parsed.lencutoff5start_f_last100"+"_post"+pr

		
	cmd="cat "+args.outdir+"/"+infile+"*_glsearchR.parsed.lencutoff5endr"+" > "+args.outdir+"/"+ infile+"_all_glsearchR.parsed.lencutoff5endr"
	os.system(cmd)
	cmd="cat "+args.outdir+"/"+infile+"*_glsearchR.parsed.lencutoff5parsed_r"+" > "+args.outdir+"/"+ infile+"_all_glsearchR.parsed.lencutoff5parsed_r"
	os.system(cmd)
	cleantagfile(args.outdir+"/"+infile+"_all_glsearch1.parsed.lencutoff5parsed_f")
	cleantagfile(args.outdir+"/"+ infile+"_all_glsearchR.parsed.lencutoff5parsed_r")
	inputseqs=builddict_sequences(args.outdir+"/"+ infile+"_all_glsearchR.parsed.lencutoff5endr")
	demultiplex(args.outdir+"/"+infile+"_all_glsearch1.parsed.lencutoff5parsed_fcleaned",args.outdir+"/"+infile+"_all_glsearchR.parsed.lencutoff5parsed_rcleaned", inputseqs,tagdict,muttags_fr,sampledict,typedict)
	cmd="rm "+args.outdir+"/"+infile+"*"
	os.system(cmd)


	#identifies primerF in all reads
#	search_parse_primerseq(primerfset,"primerf",inputseqs,infile)
#	print "Identified primerF in sequences"

	#identifies primerR in all reads
#	search_parse_primerseq(primerrset,"primerr",inputseqs,infile)
#	print "Identified primerR in sequences"


if __name__ == '__main__':
	plen=735
	n=0

	minlen=int(args.minlen)
	nthreads=int(args.threads)
	#defines nucleotides
	ambcodes={'N':['A','G','C','T'],'K':['G','T'],'M':['A','C'],'R':['A','G'],'Y':['C','T'],'S':['C','G'],'W':['A','T'],'B':['C','G','T'],'V':['A','C','G'],'H':['A','C','T'],'D':['A','G','T'],'G':['G'],'C':['C'],'T':['T'],'A':['A']}
	comps={'A':'T','T':'A','C':'G','G':'C'}
	basename=os.path.basename(args.infasta)
	start_time=time.time()
	os.system("mkdir "+args.outdir+"/demultiplexed")
	with open(args.demfile) as demulfile:
		demullines=demulfile.readlines()
		primerf=demullines[0].split(',')[3]
		primerr=demullines[0].split(',')[4].strip()
		taglen=len(demullines[0].split(',')[1])
		primerfset=createprimers(primerf,ambcodes)
		primerrset=createprimers(primerr,ambcodes)
	print taglen
	with open(args.infasta) as infile:
		with open(args.outdir+"/"+basename+"_reformat_out",'w') as outfile:
			n=1
			for line1,line2 in itertools.izip_longest(*[infile]*2):
				if len(line1.strip())!=0:
					seqid=">"+str(n)+'\n'
			 		sequence=line2.strip()
			 		if len(sequence)>=minlen:
			 			outfile.write(seqid+line2)
					n+=1


	#reformats MinION input fasta to non-interleaved fasta and fix headers

	with open(args.outdir+"/"+basename+"_reformat_out_1pdt",'w') as file1:
		with open(args.outdir+"/"+basename+"_reformat_out_2pdt",'w') as file2:
			with open(args.outdir+"/"+basename+"_reformat_out") as infile:
			    for line1,line2 in itertools.izip_longest(*[infile]*2):
			    	if len(line1.strip())!=0:
						seqid=line1
				 		sequence=line2.strip()
						if len(sequence)<int(args.slen)*2:
							file1.write(seqid+line2)
						elif len(sequence)>int(args.slen)*2 and len(sequence)<int(args.maxlen):
							file2.write(seqid+line2)

	tagdict,muttags_fr,sampledict,typedict=crmutant_m2(args.demfile,int(args.mismatch))
	print "Non-interleaved file generated"

	os.system("split -l 40000 "+args.outdir+"/"+basename+"_reformat_out_1pdt " +args.outdir+"/"+basename+"_reformat_out_1pdt_p")
	os.system("split -l 40000 "+args.outdir+"/"+basename+"_reformat_out_2pdt " +args.outdir+"/"+basename+"_reformat_out_2pdt_p")

	#Builds dictionary of sequences from reformated file
	prefix1=basename+"_reformat_out_1pdt_p"
	prefix2=basename+"_reformat_out_2pdt_p"
	print prefix1,prefix2
	partlist1=fnmatch.filter(os.listdir(args.outdir), prefix1+"*")
	n=0
#	runglsearch1(partlist1[0])
	lens=range(0,len(partlist1),nthreads)
#	if len(partlist)%nthreads!=0:
	for i,j in enumerate(lens[:-1]):
		print "Processing parts "+str(n) +" to "+str(n+nthreads-1)+ " of " + str(len(partlist1))
		parts=partlist1[j:lens[i+1]]
		p=Pool(nthreads)
		p.map(runglsearch1,parts)
		n+=nthreads
		p.close()
	parts=partlist1[lens[-1]:]
	p=Pool(nthreads)
	p.map(runglsearch1,parts)
	p.close()	

	try:
		partlist2=fnmatch.filter(os.listdir(args.outdir), prefix2+"*")
		print partlist2
		n=0
	#	runglsearch2(partlist2[0])
		lens=range(0,len(partlist2),nthreads)
	#	if len(partlist)%nthreads!=0:
		for i,j in enumerate(lens[:-1]):
			print "Processing parts "+str(n) +" to "+str(n+nthreads-1) + " of " + str(len(partlist2))
			parts=partlist2[j:lens[i+1]]
			p=Pool(nthreads)
			p.map(runglsearch2,parts)
			p.close()
			n+=nthreads
		parts=partlist2[lens[-1]:]
		p=Pool(nthreads)
		p.map(runglsearch2,parts)
		p.close()
	except IndexError:
		pass


	os.system("mv "+args.outdir+"/*all.fa "+ args.outdir+"/demultiplexed")

	with open("runtime",'w') as runtimefile:
		runtimefile.write(str(round(time.time() - start_time, 2)))
