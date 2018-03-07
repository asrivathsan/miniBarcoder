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


import sys,os,fileinput,re,time,numpy,argparse,random
from distutils.spawn import find_executable
from collections import Counter
parser=argparse.ArgumentParser(description='Script for obtaining barcodes')
parser.add_argument('-f','--infasta',help='Path to input fasta file',dest="infasta",required=True)
parser.add_argument('-d','--demfile',help='Path to demultiplexing file',dest="demfile",required=True)
parser.add_argument('-o','--outdir',help='set output directory path',dest="outdir",required=True)
parser.add_argument('-l','--minlen',help='exclude barcode sequences identified that are shorter than specified length',dest="minlen",required=True)
parser.add_argument('-m','--mode',help='run with unique tag mode (1) or dual tag mode (2)',dest="mode",default="2")
parser.add_argument('-mm','--mismatch',type=int,choices=range(0,6),help='number of mismatches allowed in tags, must be <=5, default 2',dest="mismatch",default=2)
parser.add_argument('-e','--evalue',help='evalue for primer search using glsearch36,default 1e+6',dest="evalue",default="1000000")
parser.add_argument('-g','--gaps',help='number of gaps allowed for primer identification, default 5',dest="gaps",default="5")
parser.add_argument('-D','--maxdepth',type=int,help='set max depth per coverage to improve speed, default 100, must be >2',dest="maxdepth",default=100)
parser.add_argument('-t','--threads',help='number of threads for glsearch and mafft',dest="threads",default="4")
parser.add_argument('-bl','--barcodelength',type=int,help='estimated barcode length, used for unique tag mode only, please keep it slightly shorter than actual barcode length, due to indel errors',dest="blen",default="300")

args=parser.parse_args()

def is_inpath(toolname):
	return find_executable(toolname) is not None

if is_inpath("glsearch36")==False:
	print "glsearch36 is not in path"
	sys.exit()
if is_inpath("mafft")==False:
	print "mafft is not in path"
	sys.exit()

if os.path.isdir(args.outdir):
	print "output directory exists, please delete or name another"
	sys.exit()
else:
	os.system("mkdir "+args.outdir)
	
# checks length of file to ensure the file is not empty
def chklen(inputfile):
	with open(inputfile) as input:
		l=input.readlines()
		return len(l)

#reads input file and filters all sequences shorter than specified length threshold and writes output file
def parselength(i,o,length):
	with open(i) as infile:
		with open(o,'w') as outfile:
			l=infile.readlines()
			for i,j in enumerate(l):
				if ">" in j:
					if len(l[i+1].strip())>length:
						outfile.write(j+l[i+1])

#collates information on glsearch of primer f and r and predicts start and end of barcode sequence. NEEDS RENAMING STUFF
def collate(f1,f2,f3,f4,o):

	# builds dictionaries of primer positions
	primerf_startdict=builddict_position(f1)
	primerf_enddict=builddict_position(f2)
	primerr_startdict=builddict_position(f3)
	primerr_enddict=builddict_position(f4)

	# identifies sequences with both primerf and primerr recovered
	intdict={}
	for id in list(set(primerf_startdict.keys())&set(primerr_enddict.keys())):
		intdict[id]=[primerf_startdict[id],primerr_enddict[id]]
	

	for id in list(set(primerr_startdict.keys())&set(primerf_enddict.keys())):
		intdict[id]=[primerr_startdict[id],primerf_enddict[id]]

	dmp=[]
	for id in list(set(primerf_startdict.keys())&set(primerr_startdict.keys())):
		dmp.append(id)
	for id in list(set(primerf_enddict.keys())&set(primerr_enddict.keys())):
		dmp.append(id)

	# identifies sequences with primerf or primerr only
	if args.mode=="1":
		nonintdict={}
		for id in list(set(primerf_startdict.keys())-set(primerr_enddict.keys())):
			nonintdict[id]=[primerf_startdict[id],str(int(primerf_startdict[id])+args.blen)]
		for id in list(set(primerr_startdict.keys())-set(primerf_enddict.keys())):
			nonintdict[id]=[primerr_startdict[id],str(int(primerr_startdict[id])+args.blen)]
		for id in list(set(primerf_enddict.keys())-set(primerr_startdict.keys())):
			if int(primerf_enddict[id])-args.blen >=0:
				nonintdict[id]=[str(int(primerf_enddict[id])-args.blen),primerf_enddict[id]]
		for id in list(set(primerr_enddict.keys())-set(primerf_startdict.keys())):
			if int(primerr_enddict[id])-args.blen >=0:
				nonintdict[id]=[str(int(primerr_enddict[id])-args.blen),primerr_enddict[id]]
		for id in dmp:
			nonintdict.pop(id, None)
			
	# writes output with start and end positions for retrieving sequences. 
	with open(o,'w') as outfile:
		for each in intdict.keys():
			outfile.write(each+" 1 "+ intdict[each][0]+' '+intdict[each][1]+'\n')
		if args.mode=="1":
			for each in nonintdict.keys():
				outfile.write(each+" 1 "+ nonintdict[each][0]+' '+nonintdict[each][1]+'\n')
		
#parses glsearch output for each primer to get start end positions of alignments
def parseglsearch_primermatch (infilename):
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
						hitseqs[m[1]][0]=se[0]-taglen
						hitseqs[m[1]][1]=se[0]-1
						inseqs[m[1]]=se[1]+1
				else:
					hitseqs2[m[1]][0]=se[1]+1
					hitseqs2[m[1]][1]=se[1]+taglen
					inseqs2[m[1]]=se[0]-1

	with open(infilename+"parsed_f",'w') as outfile1:
		with open(infilename+"start_f","w") as outfile3:
			for each in hitseqs.keys():
				f=[]
				for v in hitseqs[each]:
					if isinstance( v, ( int, long ) ) == True:
						f.append(1)
				if f==[1,1]:
					if inseqs[each]<avglen/2  :
						outfile1.write(each+' 1 '+str(hitseqs[each][0])+' '+str(hitseqs[each][1])+'\n')
						outfile3.write(each+' ' + str(inseqs[each])+'\n')

	with open(infilename+"parsed_r",'w') as outfile2:
		with open(infilename+"start_r","w") as outfile4:
			for each in hitseqs2.keys():
				f=[]
				for v in hitseqs2[each]:
					if isinstance( v, ( int, long ) ) == True:
						f.append(1)
				if f==[1,1]:
					if inseqs2[each]>avglen/2:
						outfile2.write(each+' 2 '+str(hitseqs2[each][0])+' '+str(hitseqs2[each][1])+'\n')
						outfile4.write(each+' ' + str(inseqs2[each])+'\n')
			

#reformats input fasta to a noninterleaved file
def reformat_minion(inputfile, outputfile):
	with open(outputfile,'w') as output:
		n=0
		for line in fileinput.input([inputfile]):
			if n==0 and ">" in line:
				output.write(">"+str(n)+'\n')
			elif ">" in line:
				output.write("\n>"+str(n)+'\n')
			else:
				output.write(line.strip())
			n+=1

# reversecomplments input sequence
def rev_comp(inseq,comps):
	comp_seq=''
	for nucl in inseq:
		comp_seq=comp_seq+comps[nucl]
	return comp_seq[::-1]

# given sequence, length, and orientation, function retrieves a subsection in the given orientation	
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
	
# uses retrieve_partseq_fn for getting sections of sequences from an entire file
def retrieve_partseq(seqdict,inputfile,outfile):
	with open(outfile,'w') as retrievedfile:
		with open(inputfile) as startendfile:
			startendlines=startendfile.readlines()
			for each in startendlines:
				m=each.strip().split(' ')
				outseq=retrieve_partseq_fn(inputseqs[m[0]],[int(m[2]),int(m[3])],int(m[1]))
				retrievedfile.write(">"+m[0]+'\n'+outseq+'\n')

#create the various possible primer sequences if a primer has an ambiguity code
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

# runs glsearch36 for a given primer sequence fasta to the fasta file
def runglsearch(primerfile,infile,i,iter):
	cmd="glsearch36 -m 8 -E "+args.evalue+" -T "+args.threads+" "+primerfile+" "+infile+" > "+args.outdir+"/temp"+str(i)+"_"+str(iter)+"_glsearch"
	os.system(cmd)


# builds dictionary of start/end position of primer match
def builddict_position(infile):
	outdict={}
	with open(infile) as file1:
		for each in file1.readlines():
			each=each.strip()
			m=each.split(" ")
			outdict[m[0]]=m[1]
	return outdict

# filter glsearch output to retain matches with best identity
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

def search_parse_primerseq(primerset,name,inputseqs):
	print inputseqs.keys()[0:10]
	findprimer(primerset,name,inputseqs)
	parse_by_id(args.outdir+"/"+"all_"+name+"_primerglsearch",int(args.gaps))
	os.system("rm "+args.outdir+"/temp*")
	parseglsearch_primermatch(args.outdir+"/"+"all_"+name+"_primerglsearch.parsed.lencutoff5")
	retrieve_partseq(inputseqs, args.outdir+"/"+"all_"+name+"_primerglsearch.parsed.lencutoff5parsed_f",args.outdir+"/"+"all_"+name+"_primerglsearch.parsed.lencutoff5parsed_f.retrieved")
	retrieve_partseq(inputseqs, args.outdir+"/"+"all_"+name+"_primerglsearch.parsed.lencutoff5parsed_r",args.outdir+"/"+"all_"+name+"_primerglsearch.parsed.lencutoff5parsed_r.retrieved")
	os.system("cat "+ args.outdir+"/all_"+name+"_primerglsearch.parsed.lencutoff5parsed_f.retrieved " +  args.outdir+"/all_"+name+"_primerglsearch.parsed.lencutoff5parsed_r.retrieved > "+ args.outdir+"/all_"+name+"_primerglsearch.parsed.lencutoff5parsed_fr.retrieved")
	with open( args.outdir+"/"+"all_"+name+"_primerglsearch.parsed.lencutoff5parsed_fr.retrieved") as primertagfile:
		with open( args.outdir+"/"+"all_"+name+"_primerglsearch.parsed.lencutoff5parsed_fr.retrievedcleaned",'w') as primertagfile_cleaned:
			taglines=primertagfile.readlines()
			ids=[]
			seqs={}
			for n,line in enumerate(taglines):
				if ">" in line:
					ids.append(line)
					seqs[line]=taglines[n+1]
			idcounts=Counter(ids)
			nseqs=0
			for id in ids:
				if idcounts[id]==1:
					nseqs+=1
					primertagfile_cleaned.write(id+seqs[id])


def findprimer(primerset,name,indict2):
	indict=dict(indict2)
	infasta=args.outdir+"/"+basename+"_reformat_out"
	for i,each in enumerate(primerset):
		iternum=0
		while True:
			with open(args.outdir+"/"+"tempfile.fa",'w') as temp:
				temp.write(">primerf\n"+each+'\n')
			runglsearch(args.outdir+"/"+"tempfile.fa",infasta,i,iternum)
			parse_by_id(args.outdir+"/"+"temp"+str(i)+"_"+str(iternum)+"_glsearch",int(args.gaps))		
			with open(args.outdir+"/"+"temp"+str(i)+"_"+str(iternum)+"_glsearch.parsed.lencutoff5") as compfile:
				t=compfile.readlines()
				if len(t)==0:
					break
				toremove={}
				for line in t:
					if "\t" in line:
						toremove[line.split('\t')[1]]=''
				for id in toremove.keys():
					indict.pop(id)
			with open(args.outdir+"/"+"tempfastasubset"+str(i)+"_"+str(iternum)+".fa",'w') as newfastafile:
				for id in indict.keys():
					newfastafile.write(">"+id+'\n'+indict[id]+'\n')
			if iternum!=0:
				os.system("rm "+args.outdir+"/"+"tempfastasubset"+str(i)+"_"+str(iternum-1)+".fa")
			infasta=args.outdir+"/"+"tempfastasubset"+str(i)+"_"+str(iternum)+".fa"
			if os.path.getsize(infasta)==0:
				break
			iternum+=1
		if os.path.getsize(infasta)==0:
			break
	os.system("cat "+args.outdir+"/temp*_glsearch > "+args.outdir+"/all_"+name+"_primerglsearch")

# this function applies a hierarchical approach to demultiplexing, used only for unique tag mode. Highest demultiplexing priority is given to perfect matches on both tags, to perfect match to one and imperfect match to other, to imperfect match to both, to perfect match to only one and no match to other and finally to imperfect match to one and not the other.
def findmatch_m1(taglistf,taglistr,muttags_f,muttags_r,indict,seqdict):
	dmpfile=open(args.outdir+"/dmpfile",'a')

	for each in indict.keys():
	
		#search perfect tag matches
		if indict[each][0] in taglistf.keys() and indict[each][1] in taglistr.keys():
			if taglistf[indict[each][0]] != taglistr[indict[each][1]]:
				dmpfile.write(each+'\t'+indict[each][0]+'\t'+indict[each][1]+'\t'+taglistf[indict[each][0]]+'\t'+taglistr[indict[each][1]]+'\n')
			elif taglistf[indict[each][0]] == taglistr[indict[each][1]]:
				with open(args.outdir+"/"+taglistf[indict[each][0]]+"_all.fa",'a') as outfile:
					outfile.write(">"+each+"; "+taglistf[indict[each][0]]+"; "+indict[each][0]+"\n"+seqdict[each]+"\n")
					
		# search perfect forward, mutant reverse
		elif indict[each][0] in taglistf.keys() and indict[each][1] in muttags_r.keys():
			if taglistf[indict[each][0]] != muttags_r[indict[each][1]]:
				dmpfile.write(each+'\t'+indict[each][0]+'\t'+indict[each][1]+'\t'+taglistf[indict[each][0]]+'\t'+muttags_r[indict[each][1]]+'\n')
			elif taglistf[indict[each][0]] == muttags_r[indict[each][1]]:
				with open(args.outdir+"/"+taglistf[indict[each][0]]+"_all.fa",'a') as outfile:
					outfile.write(">"+each+"; "+taglistf[indict[each][0]]+"; "+indict[each][0]+"\n"+seqdict[each]+"\n")
		
		# search perfect reverse, mutant forwards
		elif indict[each][0] in muttags_f.keys() and indict[each][1] in taglistr.keys():
			if muttags_f[indict[each][0]] != taglistr[indict[each][1]]:
				dmpfile.write(each+'\t'+indict[each][0]+'\t'+indict[each][1]+'\t'+muttags_f[indict[each][0]]+'\t'+taglistr[indict[each][1]]+'\n')
			elif muttags_f[indict[each][0]] == taglistr[indict[each][1]]:
				with open(args.outdir+"/"+muttags_f[indict[each][0]]+"_all.fa",'a') as outfile:
					outfile.write(">"+each+"; "+muttags_f[indict[each][0]]+"; "+indict[each][0]+"\n"+seqdict[each]+"\n")
					
		# search mutant forward and mutant reverse
		elif indict[each][0] in muttags_f.keys() and indict[each][1] in muttags_r.keys():
			if muttags_f[indict[each][0]] != muttags_r[indict[each][1]]:
				dmpfile.write(each+'\t'+indict[each][0]+'\t'+indict[each][1]+'\t'+muttags_f[indict[each][0]]+'\t'+muttags_r[indict[each][1]]+'\n')
			elif muttags_f[indict[each][0]] == muttags_r[indict[each][1]]:
				with open(args.outdir+"/"+muttags_f[indict[each][0]]+"_all.fa",'a') as outfile:
					outfile.write(">"+each+"; "+muttags_f[indict[each][0]]+"; "+indict[each][0]+"\n"+seqdict[each]+"\n")

		#search single end matches
		#perfect forward only
		elif indict[each][0] in taglistf.keys():
			with open(args.outdir+"/"+taglistf[indict[each][0]]+"_all.fa",'a') as outfile:
				outfile.write(">"+each+"; "+taglistf[indict[each][0]]+"; " +indict[each][0]+"\n"+seqdict[each]+"\n")
		
		# perfect reverse only
		elif indict[each][1] in taglistr.keys():
			with open(args.outdir+"/"+taglistr[indict[each][1]]+"_all.fa",'a') as outfile:
				outfile.write(">"+each+"; "+taglistr[indict[each][1]]+"; " +indict[each][1]+"\n"+seqdict[each]+"\n")
				
		# mutant f only
		elif indict[each][0] in muttags_f.keys():
			with open(args.outdir+"/"+muttags_f[indict[each][0]]+"_all.fa",'a') as outfile:
				outfile.write(">"+each+"; "+muttags_f[indict[each][0]]+"; " +indict[each][0]+"\n"+seqdict[each]+"\n")
				
		# mutant r only
		elif indict[each][1] in muttags_r.keys():
			with open(args.outdir+"/"+muttags_r[indict[each][1]]+"_all.fa",'a') as outfile:
				outfile.write(">"+each+"; "+muttags_r[indict[each][1]]+"; " +indict[each][1]+"\n"+seqdict[each]+"\n")
	dmpfile.close()

def findmatch_m2(taglist,muttags_fr,indict,seqdict,sampledict):
	dmpfile=open(args.outdir+"/dmpfile",'a')
	idcombs={}
	for each in indict.keys():
		try:
			idcombs[each]=(taglist[indict[each][0]], taglist[indict[each][1]])
		except KeyError:
			pass
	print len(idcombs)
	matched_idcombs={}
	for each in idcombs.keys():
		try:
			with open(args.outdir+"/"+sampledict[idcombs[each]]+"_all.fa",'a') as outfile:
				outfile.write(">"+each+"\n"+seqdict[each]+"\n")
		except KeyError:
			pass
				
	dmpfile.close()

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

def builddict_sequences(infile):
	seqdict={}
	with open(infile) as inseqs:
		l=inseqs.readlines()
		for i,j in enumerate(l):
			if ">" in j:
				seqdict[j.strip().replace(">","")]=l[i+1].strip()
	return seqdict


def create_all_mutants(tagdict):
	muttags,newtags={},[]
	for tag in tagdict.keys():
		submutant=createsubmutant(tag)
		delmutant=createdelmutant(tag)
		insmutant=createinsmutant(tag)
		mutantset=list(set(submutant)|set(delmutant)|set(insmutant))
		mutantset[:]=[x for x in mutantset if x != tag]
		for mutant in mutantset:
			newtags.append(mutant)
			muttags[mutant]= tagdict[tag]
	return muttags,newtags

def crmutant_m1(tfile,pf,pr,seqdict):
	tagfdict,tagrdict={},{}
	with open(tfile) as tagfile:
		t=tagfile.readlines()
		for each in t:
			tagfdict[each.split(',')[1]]=each.split(',')[0].strip()
			tagrdict[each.split(',')[2]]=each.split(',')[0].strip()

	indict1,indict2,indict=readprimertagfasta(pf,pr)
	for each in indict.keys():
		try: 
			indict[each][0]=indict1[each]
		except KeyError:
			pass
		try:
			indict[each][1]=indict2[each]
		except KeyError:
			pass
			
	for k in indict.keys():
		if k not in seqdict.keys():
			del indict[k]
	muttags_f,newtags_f=create_all_mutants(tagfdict)
	muttags_r,newtags_r=create_all_mutants(tagrdict)

	newtags_f_counts=Counter(newtags_f)		
	newtags_r_counts=Counter(newtags_r)	

	for k in newtags_f_counts.keys():
		if newtags_f_counts[k]>1:
			del muttags_f[k]

	for k in newtags_r_counts.keys():
		if newtags_r_counts[k]>1:
			del muttags_r[k]

	with open("dmpfile_tagf",'w') as tagfdmp:	
		for tag in muttags_f.keys():
			tagfdmp.write(tag+'\t'+muttags_f[tag]+'\n')

	with open("dmpfile_tagr",'w') as tagrdmp:	
		for tag in muttags_r.keys():
			tagrdmp.write(tag+'\t'+muttags_r[tag]+'\n')
			
	findmatch_m1(tagfdict, tagrdict,muttags_f,muttags_r,indict,seqdict)	
def crmutant_m2(tfile,pf,pr,seqdict,nbp):
	tagdict={}
	sampledict={}
	with open(tfile) as tagfile:
		t=tagfile.readlines()
		for each in t:
			tagdict[each.split(',')[1]]=''
			tagdict[each.split(',')[2]]=''
	counter=1
	for each in tagdict.keys():
		tagdict[each]="t"+str(counter)
		counter+=1
	with open(tfile) as tagfile:
		t=tagfile.readlines()
		for each in t:
			sampledict[(tagdict[each.split(',')[1]],tagdict[each.split(',')[2]])]=each.split(',')[0]
	indict1,indict2,indict=readprimertagfasta(pf,pr)
	for each in indict.keys():
		try: 
			indict[each][0]=indict1[each]
			indict[each][1]=indict2[each]
		except KeyError:
			del indict[each]
			
	n=1
	while n<=nbp:
		muttags_fr,newtags_fr=create_all_mutants(tagdict)
		newtags_fr_counts=Counter(newtags_fr)	

#		for k in newtags_fr_counts.keys():
#			if newtags_fr_counts[k]>1:
#				del muttags_fr[k]
		tagdict.update(muttags_fr)
		n+=1
	with open(args.outdir+"/temp1.fas",'w') as outfile:
		for k in tagdict.keys():
			outfile.write(k+'\t'+tagdict[k]+'\n')
	with open(args.outdir+"/dmpfile_tagfr",'w') as tagfrdmp:	
		for tag in muttags_fr.keys():
			tagfrdmp.write(tag+'\t'+muttags_fr[tag]+'\n')

	findmatch_m2(tagdict, muttags_fr,indict,seqdict,sampledict)

#def cr_tag_mutants(tag)
	
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
	if sum(indict.values())>=abs_thresh:
		poslist=[]
		n=0
		while n<len(indict.keys()[0]):
			newlist=[]
			for each in indict.keys():
				m=0
				while m<indict[each]:
					newlist.append(each[n])
					m+=1
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
			seqdict[k4]=int(k01.split(";")[-1].strip())
		conseq=consensus(seqdict,perc_thresh,abs_thresh)
		conseq=conseq.replace("-",'')
		if len(conseq)!=0:
			with open(o,'w') as outfile:
				outfile.write(">"+name.split(".")[0]+";"+str(len(conseq))+";"+str(sum(seqdict.values()))+'\n'+conseq+'\n')
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
if __name__ == '__main__':
	minlen=int(args.minlen)

	#defines nucleotides
	ambcodes={'N':['A','G','C','T'],'K':['G','T'],'M':['A','C'],'R':['A','G'],'Y':['C','T'],'S':['C','G'],'W':['A','T'],'B':['C','G','T'],'V':['A','C','G'],'H':['A','C','T'],'D':['A','G','T'],'G':['G'],'C':['C'],'T':['T'],'A':['A']}
	comps={'A':'T','T':'A','C':'G','G':'C'}
	basename=os.path.basename(args.infasta)
	start_time=time.time()
	#reformats MinION input fasta to non-interleaved fasta and fix headers
	reformat_minion(args.infasta,args.outdir+"/"+basename+"_reformat_out")

	print "Non-interleaved file generated"

	#Builds dictionary of sequences from reformated file
	inputseqs=builddict_sequences(args.outdir+"/"+basename+"_reformat_out")
	print "Input sequence dictionary built"

	seqlens=[len(i.strip()) for i in inputseqs.values()]
	avglen=numpy.mean(seqlens)
	print avglen
#	#Reads demultiplexing file for generating primer sets
	with open(args.demfile) as demulfile:
		demullines=demulfile.readlines()
		primerf=demullines[0].split(',')[3]
		primerr=demullines[0].split(',')[4].strip()
		taglen=len(demullines[0].split(',')[1])
		primerfset=createprimers(primerf,ambcodes)
		primerrset=createprimers(primerr,ambcodes)
	print taglen
	print "Generated the primersets, predicting primer start end"

	#identifies primerF in all reads
	search_parse_primerseq(primerfset,"primerf",inputseqs)
	print "Identified primerF in sequences"

	#identifies primerR in all reads
	search_parse_primerseq(primerrset,"primerr",inputseqs)
	print "Identified primerR in sequences"

	#parses sequences to get tag sequences and COIbarcode sequences
	collate(args.outdir+"/"+"all_primerf_primerglsearch.parsed.lencutoff5start_f",args.outdir+"/"+"all_primerf_primerglsearch.parsed.lencutoff5start_r",args.outdir+"/"+"all_primerr_primerglsearch.parsed.lencutoff5start_f",args.outdir+"/"+"all_primerr_primerglsearch.parsed.lencutoff5start_r", args.outdir+"/"+basename+"_reformat_out"+"_COIpred")
	print "Identified tags and COI barcodes"
	#retrieves COI and tag sequences
	retrieve_partseq(inputseqs,args.outdir+"/"+basename+"_reformat_out_COIpred",args.outdir+"/"+basename+"_reformat_out_COIpred.retrieved")
	print "COI barcodes retrieved"

   #excludes short sequences
	parselength(args.outdir+"/"+basename+"_reformat_out_COIpred.retrieved", args.outdir+"/"+basename+"_reformat_out_COIpred.retrieved_len"+str(minlen), minlen)
	print "Sequences < (minlen) bp excluded"
	inputseqs=builddict_sequences(args.outdir+"/"+basename+"_reformat_out_COIpred.retrieved_len"+str(minlen))
	#demultiplexes data
	if int(args.mode)==1:
		crmutant_m1(args.demfile, args.outdir+"/"+"all_primerf_primerglsearch.parsed.lencutoff5parsed_fr.retrievedcleaned",args.outdir+"/"+"all_primerr_primerglsearch.parsed.lencutoff5parsed_fr.retrievedcleaned", inputseqs)
	elif int(args.mode)==2:
		crmutant_m2(args.demfile, args.outdir+"/"+"all_primerf_primerglsearch.parsed.lencutoff5parsed_fr.retrievedcleaned",args.outdir+"/"+"all_primerr_primerglsearch.parsed.lencutoff5parsed_fr.retrievedcleaned", inputseqs,int(args.mismatch))
	
	print "Sequences have been demultiplexed"
	os.system("mkdir "+args.outdir+"/temp_alls")
	os.system("mv "+args.outdir+"/*all.fa "+ args.outdir+"/temp_alls")

#	#aligning and consensus calling
	if args.maxdepth!=0:
		os.system("mkdir "+args.outdir+"/"+"temp_alls_"+str(args.maxdepth))
	os.system("mkdir "+args.outdir+"/temp_alls_uniqs")
	os.system("mkdir "+args.outdir+"/temp_alls_uniqs_mafft")
	os.system("mkdir "+args.outdir+"/temp_alls_uniqs_mafft_consensus")
	dirlist=os.listdir(args.outdir+"/temp_alls")

	#merging sequences to uniques
	for name in dirlist:
		if args.maxdepth!=0:
			subset_randomly(args.outdir+"/temp_alls/"+name,args.outdir+"/temp_alls_"+str(args.maxdepth)+"/"+name,int(args.maxdepth))
			allunique(args.outdir+"/temp_alls_"+str(args.maxdepth)+"/"+name, args.outdir+"/temp_alls_uniqs/"+name+".uniq")
		else:
			allunique(args.outdir+"/temp_alls/"+name, args.outdir+"/temp_alls_uniqs/"+name+".uniq")
		
	#align sequences
	dirlist=os.listdir(args.outdir+"/temp_alls_uniqs")
	for name in dirlist:
		os.system('mafft --adjustdirection --op 0 --thread '+args.threads+' '+args.outdir+'/temp_alls_uniqs/'+ name + ' > ' +args.outdir+'/temp_alls_uniqs_mafft/'+name + '_aln.fasta')

    #majority rule consensus
	dirlist=os.listdir(args.outdir+"/temp_alls_uniqs_mafft")
	for name in dirlist:
		callconsensus(args.outdir+"/temp_alls_uniqs_mafft/"+name,0.5,10,args.outdir+"/temp_alls_uniqs_mafft_consensus/"+name.split(".")[0]+"_consensus.fa",name)
	print "Alignment and consensus calling completed"
	os.system("cat "+args.outdir+"/temp_alls_uniqs_mafft_consensus/* > " + args.outdir+"/all_barcodes.fa")
	with open("runtime",'w') as runtimefile:
		runtimefile.write(str(round(time.time() - start_time, 2)))
