import sys,os,itertools
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC,generic_dna
bps = ['A','T','G','C','U','N','M','R','W','S','Y','K','V','H','D','P','B']

#move gaps in all positions of a window. return a non-duplicated set of kmers (nodups)
def create_gap_combs(seq,ngaps):
	n=1
	seqset=[seq.replace("-",'')]
	while n<=ngaps:
		nseqset=[]
		for k,each in enumerate(seqset):
			pos=0
			while pos <= len(each):
				nseqset.append(each[:pos]+"-"+each[pos:])
				pos+=1
		seqset=nseqset
		n+=1
	nodups={}
	for each in seqset:
		nodups[each]=''
	return nodups.keys()

#get the reference codons corresponding to window around each gap. If reference codons contain "N" in the window, it is excluded.
def get_ref_codons(refset,startpos,frame,n):
	endposset=[]
	for seq in refset:
		counter,endposseq=0,startpos
		for i,j in enumerate(seq[frame:][startpos:]):
			if counter <=n:
				if j!="-":
					counter+=1
				endposseq+=1
			else:
				break
		endposset.append(endposseq)
	endpos=min(endposset)-1
	refcodonset=[]
	gapcount=[]
	for each in refset:
		if "N" not in each[frame:][startpos:endpos]:
			gapcount.append(each[frame:][startpos:endpos].count('-'))
	mingapcount=0
	try:
		mingapcount=min(gapcount)
		for each in refset:
			if each[frame:][startpos:endpos].count('-')==mingapcount:
				refcodonset.append(each[frame:][startpos:endpos])
	except ValueError:
		pass
	nrefcodonset=[]
	for each in refcodonset:
		if "N" not in each:
			nrefcodonset.append(each)
	return endpos, nrefcodonset,mingapcount

#for insertions, create a set of kmers that change the position to be deleted
def create_del_combs(seq,ngaps):
	n=1
	seqset=[seq]
	while n<=ngaps:
		nseqset=[]
		for k,each in enumerate(seqset):
			pos=0
			while pos < len(each):
				nseqset.append(each[:pos]+each[pos+1:])
				pos+=1
		seqset=nseqset
		n+=1
	nodups={}
	for each in seqset:
		nodups[each]=''
	return nodups.keys()

def parsealn(inset):
	newset=[]
	cols_to_remove=[]
	n=0
	while n<len(inset[0]):
		tlist=[each[n] for each in inset]
		nored=[]
		for each in tlist:
			if each not in nored:
				nored.append(each)
		if nored==['-']:
			cols_to_remove.append(n)
		n+=1
	for each in inset:
		for n in cols_to_remove:
			each=''.join(list(each).pop(n))
		newset.append(each)
	return newset
		
# compares a set of amino acids to a reference sequence set and calculates the best score.
def find_amino_similarity(seqset, compset):
	best_score=0
	comp_by_pos=[]
	maxlen=max([len(each) for each in compset])
	#excludes reference codons with partial overlaps
	newcompset=[each for each in compset if len(each)==maxlen]
	n=0
	while n<len(newcompset[0]):
		comp_by_pos.append([each[n] for each in newcompset])
		n+=1

	for seq in seqset:
		score=0
		for i,j in enumerate(seq):
			try:
				if j in comp_by_pos[i]:
					score+=1
			except IndexError:
				pass
		if score>best_score:
			best_score=score

	return best_score
#creates strict consensus of a sequence set
def strict_con(seqset):
	con=''
	n=0
	while n<len(seqset[0]):
		tempset={}
		for each in seqset:
			each=each.replace("-","n")
			tempset[each[n]]=''
		if len(tempset.keys())==1:
			con+=tempset.keys()[0]
		else:
			con+="n"
		n+=1
	return con
		
#get the correct reading frame.
def get_cor_frame(seq,gencode):
	seqset=[Seq(seq,generic_dna),Seq(seq[1:],generic_dna),Seq(seq[2:],generic_dna),Seq(seq,generic_dna).reverse_complement(),Seq(seq[:-1],generic_dna).reverse_complement(),Seq(seq[:-2],generic_dna).reverse_complement()]
	aminoset=[]
	for each in seqset:
		a=each.translate(table=gencode,to_stop=True)
		aminoset.append(a.__str__())
	maxlen,corframe=0,0
	for i,each in enumerate(aminoset):
		if len(each)>maxlen:
			maxlen=len(each)
			corframe=i+1
	return corframe
def remove_ext_gaps(seq):
	startpos=0
	endpos=0
	for n,bp in enumerate(seq):
		if bp!="-":
			startpos=n
			break
	for n,bp in enumerate(seq[:: - 1]):
		if bp !="-":
			endpos = len(seq) - n 
			break
	newseq=seq[startpos:endpos]
	return newseq,startpos,endpos	
	
def remove_ext_Ns(seq):
	startpos=0
	endpos=0
	for n,bp in enumerate(seq):
		if bp!="N":
			startpos=n
			break
	for n,bp in enumerate(seq[:: - 1]):
		if bp !="N":
			endpos = len(seq) - n 
			break
	newseq=seq[startpos:endpos]
	return newseq
#locate the start and end position of a sequence in a multiple sequence alignment, replace terminal gaps with Ns
def get_start_end2(seq):
	startpos=0
	endpos=0
	for n,bp in enumerate(seq):
		if bp!="-":
			startpos=n
			break
	for n,bp in enumerate(seq[:: - 1]):
		if bp !="-":
			endpos = len(seq) - n 
			break
	newseq=startpos*'N'+seq[startpos:endpos]+(len(seq)-endpos)*'N'
	return newseq,startpos, endpos
	
#locate the start and end position of a sequence in a multiple sequence alignment, replace terminal gaps with Ns
def get_start_end(seq):
	startpos=0
	endpos=0
	for n,bp in enumerate(seq):
		if bp!="-":
			startpos=n
			break
	for n,bp in enumerate(seq[:: - 1]):
		if bp !="-":
			endpos = len(seq) - n 
			break
	newseq=startpos*'N'+seq[startpos:endpos]+(len(seq)-endpos)*'N'
	return newseq

# identify positions of gaps in uncorrected sequence window: This step ensures that there are at least 5 references that support deletion
def identify_uncor_gaps(window,flag,winlen):
	gappos=[]
	for n,bp in enumerate(window):
		if bp =="-":
			gappos.append(n)
	for pos in gappos:
		if pos in range(0,winlen*3+3):
			flag=1
	return gappos,flag

# identify positions of gaps in reference sequences
def identify_ref_gaps(refset,flag,winlen):
	refgappos=[]
	for each in refset:
		for n,bp in enumerate(each):
			if bp=="-":
				if n not in refgappos:
					refgappos.append(n)	
	if len(refgappos)>0 and len(refset)>=5:
		for pos in refgappos:
			if pos in range(1,winlen*3+3):
				if flag==False:
					flag=2
			elif flag==1:
				flag=3
	return refgappos,flag

# correct insertions in a given window. if corrected, return flag=4, so that the length of the sequence can be adjusted
def correct_ins(c3set,refgappos,gencode,aminosetrefs):
	con=c3set
	nseqset=create_del_combs(c3set,len(refgappos))
	combaminoset={}
	allcombaminolist=[]
	for var in nseqset:
		varset=[]
		for v in var:
			if v=="-":
				varset.append(["A","T","G","C"])
			else:
				varset.append([v])
		aminoset=[Seq(''.join(list(each)),generic_dna).translate(table=gencode,to_stop=True).__str__() for each in list(itertools.product(*varset))]
		combaminoset[var]= aminoset
		allcombaminolist+=aminoset
	corr_sets=[]
	best_score=find_amino_similarity(allcombaminolist,aminosetrefs)
	for k in combaminoset.keys():
		tscore=find_amino_similarity(combaminoset[k],aminosetrefs)
		if tscore==best_score:
			corr_sets.append(k)
	if len(corr_sets)!=0:
		con=strict_con(corr_sets)
		flag=4
	else:
		con=c3set
	return con,flag

# correct deletions in a given window. 
def correct_del(c3set,gencode,aminosetrefs):
	nseqset=create_gap_combs(c3set,c3set.count("-"))
	combaminoset={}
	allcombaminolist=[]
	for var in nseqset:
		varset=[]
		for v in var:
			if v=="-" or v=="n":
				varset.append(["A","T","G","C"])
			else:
				varset.append([v])

		aminoset=[Seq(''.join(list(each)),generic_dna).translate(table=gencode,to_stop=True).__str__() for each in list(itertools.product(*varset))]
		combaminoset[var]= aminoset
		allcombaminolist+=aminoset

	corr_sets=[]
	best_score=find_amino_similarity(allcombaminolist,aminosetrefs)
	for k in combaminoset.keys():
		tscore=find_amino_similarity(combaminoset[k],aminosetrefs)
		if tscore==best_score:
			corr_sets.append(k)
	if len(corr_sets)!=0:
		con=strict_con(corr_sets)
	else:
		con=c3set
	return con
def runcorrection(mseq,refseqset,orseq,gencode,winlen):
	nrefseqset=refseqset
	i=0
	newseq=''
	while i<len(mseq)-winlen*3:
		refend,refcodonset,gapcount=get_ref_codons(refseqset,i,0,winlen*6+3)
		if i!=len(mseq)-winlen*3-3:
			c3set=mseq[i:refend]
		else:
			c3set=mseq[i:]
		if c3set.replace("-","")=="":
			break

		if len(refcodonset)!=0:
			#define type of sliding window: flag: 1: deletion only 2: insertion only 3: insertion and deletion
			flag=False
			#identify positions in the sliding window that are gaps
			gappos,flag=identify_uncor_gaps(c3set,flag,winlen)
			refgappos,flag=identify_ref_gaps(refcodonset,flag,winlen)
			# modify sliding window based on gap position: if 
			temprefcodonset=[]
			if flag==2:
				if len(refgappos)>4:
					break
			if flag==1:
				for num in range((winlen+2)*3,(winlen*2+1)*3+3):
					if num in gappos:
						c3set=c3set[:num/3*3]
						for each in refcodonset:
							temprefcodonset.append(each[:num/3*3])
					elif (winlen+1)*3 in gappos:
						c3set=mseq[i:refend+3]
						refend,refcodonset,gapcount=get_ref_codons(refseqset,i,0,winlen*6+6)
				if temprefcodonset!=[]:
					refcodonset=temprefcodonset

			aminosetrefs=[Seq(each.replace("-",""),generic_dna).translate(table=gencode,to_stop=True).__str__() for each in refcodonset]

			if flag==1 or flag==3:
				con=correct_del(c3set,gencode,aminosetrefs)
			elif flag==2 or flag==3:
				con,flag=correct_ins(c3set,refgappos,gencode,aminosetrefs)
			elif flag==False:
				con=c3set
			if flag==1:
				newseq+=con
				i+=len(con)

			elif flag==4:
				nrefseqset=[]
				for each in refseqset:
					nrefseqset.append(each[:i+min(refgappos)]+each[i+max(refgappos)+1:])
				refseqset=nrefseqset
				orseq=orseq[:i+min(refgappos)]+each[i+max(refgappos)+1:]
				if i==len(mseq)-(winlen+1)*3:
					newseq+=con
				else:
					newseq+=con[0:3]
				i+=3

			elif flag==2 or flag==3:
				if i==len(mseq)-(winlen+1)*3:
					newseq+=con
				else:
					newseq+=con[0:3]
				i+=3+gapcount

			else:
				if i==len(mseq)-(winlen+1)*3:
					newseq+=con
				else:
					newseq+=con[0:3]
				i+=3
			mseq=mseq.replace(c3set,con)
		else:
			i+=3
	return newseq,orseq,nrefseqset
def check_alignment(inlist,support):
	aligned_refs=[]
	uncor_minion_seq=''
	for i,j in enumerate(inlist):
		if i>0 and ">" in j:
			aligned_refs.append(inlist[i+1].strip().upper())
		if i==0 and ">" in j:
			uncor_minion_seq=inlist[i+1].strip().upper()
	#Trim alignment to region covered by BLAST hits.
	start,end=[],[]
	for seq in aligned_refs:
		newseq,startseq,endseq=remove_ext_gaps(seq)
		start.append(startseq)
		end.append(endseq)
	minstart,maxend=min(start), max(end)
	seqset=[get_start_end(uncor_minion_seq[minstart:maxend])]
	for seq in aligned_refs:
		seqset.append(get_start_end(seq[minstart:maxend]))
	#exclude positions in alignments if and only if 
	#1) only sequences <support threshold create an insertion and 
	#2) those positions are also gaps in the minion uncorrected barcodes. 
	#I.E. insertions exclusively created by few reference sequences
	to_remove=[]
	flag=True
	for bp in range(0,len(seqset[0])):
		bpset={"-":[],"N":[]}
		for seq in seqset[1:]:
			if seq[bp]=="-":
				bpset["-"].append(seq)
			else:
				bpset["N"].append(seq)
	#	print bp, len(bpset["-"]), len(bpset["N"])
	#	if len(bpset["-"])>0 and len(bpset["-"])<=support:
	#		for each in bpset["-"]:
	#		if seqset[0][bp]!="-":
	#			to_remove.append(bp)
			#	seqset.remove(each)
			#flag=False
		if len(bpset["N"])>0 and len(bpset["N"])<=support:
		#	print bpset["N"]
		#	for each in bpset["N"]:
			if seqset[0][bp]=="-":
				to_remove.append(bp)
		#	flag=False
#	print to_remove
	trimmed_seqset=[]
	for seq in seqset:
		newseq=''
		for bp in range(0,len(seq)):
			if bp not in to_remove:
				newseq+=seq[bp]
		trimmed_seqset.append(newseq)
#	print trimmed_seqset
	return trimmed_seqset,flag
	
def translate_corframe(seqlist,refseq,gencode):
	corframe=get_cor_frame(seqlist[refseq].replace("-",""),gencode)
#	print corframe,seqlist[refseq]
	translatedset=[]
	for each in seqlist[1:]:
		if corframe==1:
			translatedset.append(Seq(each.replace('-',''),generic_dna).translate(table=gencode,to_stop=True).__str__())
		if corframe==2:
			translatedset.append(Seq(each[1:].replace('-',''),generic_dna).translate(table=gencode,to_stop=True).__str__())
		if corframe==3:
			translatedset.append(Seq(each[2:].replace('-',''),generic_dna).translate(table=gencode,to_stop=True).__str__())
		if corframe==4:
			translatedset.append(Seq(each.replace('-',''),generic_dna).reverse_complement().translate(table=gencode,to_stop=True).__str__())
		if corframe==5:
			translatedset.append(Seq(each[:-1].replace('-',''),generic_dna).reverse_complement().translate(table=gencode,to_stop=True).__str__())
		if corframe==6:
			translatedset.append(Seq(each[:-2].replace('-',''),generic_dna).reverse_complement().translate(table=gencode,to_stop=True).__str__())
	refseqset=[]
	
	if corframe in [1,2,3]:
		mseq=seqlist[0][corframe-1:]
		orseq=seqlist[0][corframe-1:]
		for each in seqlist[1:]:
			refseqset.append(each[corframe-1:])
	else:
		corframe=corframe-3
		orseq=Seq(seqlist[0],generic_dna).reverse_complement().__str__()[corframe-1:]
		mseq=Seq(seqlist[0],generic_dna).reverse_complement().__str__()[corframe-1:]
		for each in seqlist[1:]:
			refseqset.append(Seq(each,generic_dna).reverse_complement().__str__()[corframe-1:])
	return mseq,refseqset,orseq
def change_ext_gaps(sequence):
	start_pos, end_pos = 0, 0
	for i,bp in enumerate(sequence):
		if bp in bps:
			start_pos = i - 1
			break
	for i,bp in enumerate(sequence[:: - 1]):
		if bp in bps:
			end_pos = len(sequence) - i
			break
	new_sequence = "?" * (start_pos + 1) + sequence[start_pos + 1 : end_pos] + "?" * (len(sequence) - end_pos)
	return new_sequence

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
def strict_con2(seqset):
	con=''
	n=0
	while n<len(seqset[0]):
		tempset={}
		for each in seqset:
			each=each.replace("-","N")
			tempset[each[n]]=''
		if len(tempset.keys())==1:
			con+=tempset.keys()[0]
		else:
			con+="N"
		n+=1
	return con

def find_hp(seq1,hplen):
	counter={}
	n=0
	newseq=''
	while n<len(seq1):
		i=0
		seq=seq1[n+i].replace("-","N")
		curbp=seq1[n+i]
		try:
			while seq1[n+i+1]==seq1[n+i] or seq1[n+i+1]=="-" or seq1[n+i]=="-":
				if n+i+1<len(seq1):
					if seq1[n+i+1]=="-":
						seq+="N"
					else:
						seq+=seq1[n+i+1]
					i+=1
				else:
					break
		except IndexError:
			pass
		if len(seq.replace("N",""))>=int(hplen):
			newseq+=seq
		else:
			newseq+="N"*len(seq)
		n=n+i+1
	return newseq

def retrieve_multns(nseq,seqset1,n,hplen):
	k1,k2=n,n
	while nseq[k1-1]=="N":
		k1=k1-1
	while nseq[k2+1]=="N":
		if k2!=len(nseq)-2:
			k2+=1
		else:
			break
	newseqset=[]
	for i,each in enumerate(seqset1):
		if i==0:
			subseq=each[k1:k2+1]
			newseqset.append(find_hp(subseq,hplen))
		else:
			if "N" not in each[k1:k2+1]:
				subseq=each[k1:k2+1]
				newseqset.append(find_hp(subseq,hplen))
	if len(newseqset)>=5:
		test=strict_con2(newseqset).split('N')
		newtest=[]
		flag=False
		for i,item in enumerate(test):
			if len(item)<int(hplen) and len(item)>0:
				newtest.append('N')
			else:
				newtest.append(item)
		strict_consensus='N'.join(newtest)

	else:
		strict_consensus=nseq[k1:k2+1]
	return strict_consensus,k1,k2+1

