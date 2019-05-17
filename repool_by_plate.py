#!/usr/bin/env python 
# usage is repool_by_plate.py inputrepoolingsheet outputexcelfilename coverage_minimum coverage_maximum
import sys,csv,xlwt,math
#style = xlwt.easyxf('pattern: pattern solid, fore_colour 0x7')
style = xlwt.easyxf('pattern: pattern solid, fore_colour 0x22')
o=open(sys.argv[1])
l=o.readlines()

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
		for each in range(n,95,8):
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
			if value <=int(sys.argv[4]) and value >=int(sys.argv[3]):
				sheet.write(i,u,v,style)
			else:
				sheet.write(i,u,v)
wb.save(sys.argv[2]+".xls")
#	for i in eachlist:
#		o.write('\t'.join(i)+'\n')
#	o.close()


