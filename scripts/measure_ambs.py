import sys,os
with open(sys.argv[1]+"ambstats",'w') as outfile:
	with open(sys.argv[1]) as infile:
		l=infile.readlines()
		for i,j in enumerate(l):
			if ">" in j:
				outfile.write(j.strip()[1:]+'\t'+str(l[i+1].upper().count("N"))+'\t'+str(len(l[i+1].strip()))+'\n')
