#!/usr/bin/env python 
import sys,csv,xlwt,math,argparse
#style = xlwt.easyxf('pattern: pattern solid, fore_colour 0x7')

parser=argparse.ArgumentParser(description='Script for repooling low coverage data')
parser.add_argument('-i','--infile',help='input tab delimited file see example',dest="infile",required=True)
parser.add_argument('-o','--output',help='prefix for excel outputfile',dest="outfile",required=True)
parser.add_argument('-m','--min',help='minimum coverage',dest="min",required=True)
parser.add_argument('-M','--max',help='minimum coverage',dest="max",required=True)
parser.add_argument('-c','--colour',help='colour to highlight cells ',dest="colour",default="0x22")
args=parser.parse_args()

style = xlwt.easyxf('pattern: pattern solid, fore_colour '+args.colour)
infile=open(args.infile)
l=infile.readlines()

p={}
for each in l:
	m=each.split('\t')
	p[m[3]]=[]

for each in l:
	m=each.strip().split('\t')
	p[m[3]].append(m[0]+':'+m[4])



k={}

def platebuilder(platelist):
	n=0
	outlist=[]
	while n<8:
		t=[]
		for each in range(n,96,8):
			try:
				t.append(platelist[each])
			except IndexError:
				print "not generated"
				break
		n+=1
		outlist.append(t)
	return outlist

wb=xlwt.Workbook()
for each in p.keys():
	o=open(each,'w')
	sheet=wb.add_sheet(each)
#	for n in range(0,12):
#		sheet.col(n).width=1079
#	for n in range(0,8):
#		sheet.row(n).height_mismatch=True
#		sheet.row(n).height=510
	print each
	eachlist=platebuilder(p[each])
	for i,j in enumerate(eachlist):
		for u,v in enumerate(j):
			value=int(v.split(":")[1])
			if value <=int(args.max) and value >=int(args.min):
				sheet.write(i,u,v,style)
			else:
				sheet.write(i,u,v)
wb.save(args.outfile+".xls")
#	for i in eachlist:
#		o.write('\t'.join(i)+'\n')
#	o.close()


