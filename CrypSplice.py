#!/usr/bin/env python
# Last update 1/31/2016 #
# Overlapping Gains/Loss/Diff #
# Added Parser #
# Check file existence #
# Added default match values #
# Isoform Duplicate junction removal # 

########################### CrypSplice Version 1 ##############################
#    Identification of cryptic splice sites from RNA-Seq data 				  #
#    Copyright (C) 2015  Hari Krishna Yalamanchili							  #
#    Results: JN-Gain/Loss/Diff.txt     Log file : Log_CrypSplice.txt         #
#    Requires bedtools V 2.17 or higher to be accessed from cmd line		  #
#    Requires ibb and MASS R libraries, see README for more information 	  #
#    Requires sorted bam files named <sample.sorted.bam> in working directory #
#    i.e. for a sample XX.bed corresponding bam file should be XX.sorted.bam  #
#    Under the terms of the GNU General Public License as published by		  #
#    the Free Software Foundation (version 3 or later) and No Warranty		  #
###############################################################################


# **************** Adjust BED Files **************** #
def adjustbed(sample,covFilter):
	fileread=open(sample,'r')
	filewrite=open(sample+'.Ad','w')	
	lines=fileread.readlines()
	for line in lines:
		line.strip()
		data1=line.split('\t')
		try:
			if int(data1[4]) >= covFilter and data1[5]!="?":
				data2=data1[10].split(',')
				data1[1]=int(data1[1])+int(data2[0])
				data1[2]=int(data1[2])-int(data2[1])
				s='\t'.join(str(v) for v in data1)
				filewrite.write(s)
		except:
			pass
	fileread.close()
	filewrite.close()

# ************* Remove Known Junctions ************* #
def removeknownjunctions(sample,jnFile,junMatch):
	sample = sample.strip()
	cmd='bedtools intersect -a '+sample+'.Ad -b '+jnFile+' -v -f '+str(junMatch)+' -r -s > Temp_'+sample+'.bed'
	logfile.write(sample+'\t'+cmd+'\n')
	print(sample+'\t'+cmd+'\n')
	os.system(cmd)
	os.system('mv Temp_'+sample+'.bed '+sample+'.Ad.Fl')
	
# ************ Collapse and sort ******************** #
def collapsesort(sample):
	sample = sample.strip()
	fileread=open(sample+'.Ad.Fl','r')
	filewrite=open(sample+'.Ad.Fl'+'.Cl','w')
	junctions=fileread.readlines()
	i=0
	k=0 # Container
	if i == (len(junctions)-1):
		junctions="\n".join(junctions)
		filewrite.write(junctions)
		#print str(i)+"\t"+str(len(junctions))
		
	while i < (len(junctions)-1):
		current=junctions[i].split('\t')
		j=i+1
		while j < len(junctions):
			check=junctions[j].split('\t')
			try:
				if((check[0] != current[0]) or (((int(check[1])-int(current[1])) > 500) and ((int(check[2])-int(current[2])) > 500))):
					break
				if((check[0] == current[0]) and (abs(int(check[1]) - int(current[1])) < 100) and (abs(int(check[2]) - int(current[2])) < 100) and (check[5] == current[5])):
					new_cov=int(check[4]) + int(current[4])
					if(int(check[4]) > int(current[4])):
						current=check
						current[4]=str(new_cov)	
					if(int(check[4]) <= int(current[4])):
						current[4]=str(new_cov)
					junctions.pop(j)	
				else:
					j+=1;
			except:
				break	
		s = '\t'.join(current)
		filewrite.write(s)
		i+=1;
	filewrite.close()
	fileread.close()
	
	# **********Sort the collapsed junction file********** #
	os.system('sort -V -k1,1 -k2,2n '+sample+'.Ad.Fl'+'.Cl > Temp_'+sample+'.temp')
	os.system('mv Temp_'+sample+'.temp '+sample+'.Ad.Fl'+'.Cl')

# ********* Jn Single co-ordinate system  ************* #
# this command is just comparing the sample jn against the annotated jn first
def processjunctionbed(sample,i,jn_list,ex_list):
	fileread1=open(sample+'.Ad.Fl.Cl','r')
	fileread2=open('All_Junctions.bed.So.Cl','r')
	filewrite=open(sample+'.Ex.bed','w')
	jn_dict={}
	ex_dict={}
	sample_junctions=fileread1.readlines() #readlines() outputs a list containing each line as an item of the list
	all_junctions=fileread2.readlines()
	
	#with open(sample+'.txt', 'w') as f:
	#	for item in sample_junctions:
	#		f.write("%s\n" % item)
	#with open(sample+'alljunctions.txt', 'w') as f:
	#	for item in all_junctions:
	#		f.write("%s\n" % item)

	for line in sample_junctions:
		line.strip() #strip removes any leading and trailing spaces of the string
		data=line.split('\t') #data (sample) now contains a list of the items in each line that is seperated by a tab
		ii=0
		while ii < len(all_junctions): #len is the number of lines in already annotated junctions
			current=all_junctions[ii].split('\t')
			ii=ii+1
			# If not same chr continue #
			if (data[0] != current[0]):
				continue
			# If passed that coordinate quit loop #
			if((data[0] == current[0]) and (((int(current[1])-int(data[1])) > 100) and ((int(current[2])-int(data[2])) > 100))):
				break
			# If not check for overlap #
			if((data[0] == current[0]) and (abs(int(data[1]) - int(current[1])) < 100) and (abs(int(data[2]) - int(current[2])) < 100) and (data[5] == current[5])):
				key=data[0]+'\t'+current[1]+'\t'+current[2]
				jn_dict[key]=data[4]
				#********** Generate Exon coverage Bed template *************# 
				x=0;
				y=0;
				if data[5] == '-':
					x=int(data[2])+2 # smallest exon 4nt
					y=int(data[2])+2 # smallest exon 4nt
				if data[5] == '+':
					x=int(data[1])-2 # smallest exon 4nt
					y=int(data[1])-2 # smallest exon 4nt
				val=data[0]+'\t'+str(x)+'\t'+str(y)+'\t'+data[3]+'\t'+'0'+'\t'+data[5]
				ex_dict[key]=val
				filewrite.write(data[0]+'\t'+str(x)+'\t'+str(y)+'\t'+data[3]+'\t'+'0'+'\t'+data[5]+'\n')
				break
	
	jn_list[i]=jn_dict
	ex_list[i]=ex_dict
	fileread1.close()
	fileread2.close()
	filewrite.close()
# ********* Remove extra Columns **********#
def removeextracols(sample):
	#fileread=open(sample+'.Ex.bed2','r')
	#line=fileread.readlines()
	sample = sample.strip()
	cutcol='cut -f1,2,3,4,5,6,13,14,15,16- '+sample+'.Ex.bed2 > Temp_'+sample+'.Ex.bed2'
	#logfile.write(sample+'\t'+cutcol+'\n')
	#print(sample+'\t'+cutcol+'\n')
	os.system(cutcol)
	os.system('mv Temp_'+sample+'.Ex.bed '+sample+'.Ex.bed2')

# ********* Ex Single co-ordinate system  ************* #	
def processexonbed(sample,i,ex_list):
	temp_dict={}
	sample = sample.strip()
	fileread=open(sample+'.Ex.bed2','r')
	temp_dict=ex_list[i]
	exon_junctions=fileread.readlines()

	#with open(sample+'.2.txt', 'w') as f:
	#	for item in exon_junctions:
	#		f.write("%s\n" % item)

	for line in exon_junctions:
		line.strip()
		data=line.split('\t')
		if len(data) == 10:
			for key, value in temp_dict.iteritems():
				#v=value.strip()
				#v1=v.split('\t')
				#if((data[0]==v1[0]) and ())
				if value in line and len(str(value)) > 20:
					value=data[6]
					#print(data[6])
					temp_dict[key]=value	
		else:
			pass
	ex_list[i]=temp_dict
	fileread.close()

# **************** Process BBoutput files **************** #
def processbboutput(sample,ns):
	fileread=open(sample,'r')
	filewrite=open(sample+'_temp.bb.txt','w')	
	lines=fileread.readlines()
	for line in lines:
		line.strip()
		try:
			data1=line.split()
			data2=data1[0].split(':')
			data3=data2[1].split('-')
			data4=data3[1].split('_')
			jd=data1[1:(ns+1)]
			jds=' '.join(str(v) for v in jd)
			ed=data1[(ns+1):((2*ns)+1)]
			eds=' '.join(str(v) for v in ed)
			pd=data1[((2*ns)+1):len(data1)]
			pds='\t'.join(str(v) for v in pd)
			if data4[1] == '+':
				pass
			else:
				data4[1] = '-'
			filewrite.write(data2[0]+'\t'+data3[0]+'\t'+data4[0]+'\t0\t0\t'+data4[1]+'\t'+jds+'\t'+eds+'\t'+pds+'\n')
		except:
			pass
	fileread.close()
	filewrite.close()
	os.system('mv '+sample+'_temp.bb.txt '+sample)
	
# ********************************** MAIN ********************************** #
# import  modules #
try:
	import glob,os,sys, re
	import time
	import subprocess
	import math
	from multiprocessing import Process, Manager
	from time import sleep
	from collections import defaultdict
	from sys import argv, exit
	import argparse 
except ImportError:

    print '\n\nCrypSplice requires glob, os, sys, time and collections modules\n\n'
    exit()

# Argument Parser #
parser = argparse.ArgumentParser(description='Process Arguments')
parser.add_argument('-C','--control', help="Comma separated list of control bed files. ", required=True)
parser.add_argument('-T','--treatment', help="Comma separated list of treatment bed files. ", required=True)
parser.add_argument('-G','--genome', help="Reference Genome Directory with gene model, known junctions and known alternatice splicing events.", required=True)
parser.add_argument('-F','--filter', help="Read filter cutoff. Default 10.",type=int, default=10)
parser.add_argument('-M','--match', help="Junction Overlap cutoff. Default 0.90.",type=float, default=0.90)
parser.add_argument('-P','--ncores', help="No. of cores to use. Default 1.", type=int, default=1)
#parser.add_argument('-R','--results', help="Results directory.",default=".")
args = parser.parse_args()


# Check parsed arguments and filters #
crFile = args.control
cr_samples=crFile.split(',')
for temp in cr_samples:
	if(os.path.isfile(temp)):
		temp2=temp.split(".")
		temp2=temp2[0]+".sorted.bam"
		if(os.path.isfile(temp2)):
			pass
		else:
			print "Error message: File "+temp2+" not found."
	else:
		print "Error message: File "+temp+" not found."
		exit()
ncr=len(cr_samples)

trFile = args.treatment
tr_samples=trFile.split(',')
for temp in tr_samples:
	if(os.path.isfile(temp)):
		temp2=temp.split(".")
		temp2=temp2[0]+".sorted.bam"
		if(os.path.isfile(temp2)):
			pass
		else:
			print "Error message: File "+temp2+" not found."
	else:
		print "Error message: File "+temp+" not found."
		exit()
ntr=len(tr_samples)

samples=cr_samples+tr_samples
ns=len(samples)

genome = args.genome
jnFile = genome+'/'+genome+'-Known-SpliceJunctions.txt'
gmFile = genome+'/'+genome+'.genes.bed'
aeFile = genome+'/'+genome+'-Known-AE.txt'
covFilter = int(args.filter)
junMatch = float(args.match)
np = int(args.ncores)
if np ==0:
	np = 2 
#rdir=args.results

# check for bedtools #
try:
    os.system('bedtools --version')
except:
    print "\n\nbedtools cannot be accessed from cmd line\n\n"
    exit()
    
# ********* Logging the script **************** #
logfile = open('Log_CrypSplic.txt','w')

# ********* Print Start of Script ************* #
localtime = time.strftime("%a:%H:%M:%S")
logfile.write('Started CrypSplice at: ' + localtime + '\n')


#  ********* Define an output Manager ************* #
manager = Manager()

# ****** Data structures to store jn coverage ****** #
jn_list = manager.list()
ex_list = manager.list()
for i in range(0,ns):
	jn_list.insert(i,0)
for i in range(0,ns):
	ex_list.insert(i,0)
	
#************* Adjust BED Files and Filter Reads less than Filter -F *************** #
Pros = []
nloops=math.ceil(ns/float(np))
for x in range(0,int(nloops)):
	if (x+1)*np <= ns:
		lsamples=samples[(x*np):(x+1)*np]
	else:
		if ns < np:
			lsamples=samples
		if ns > np:
			lsamples=samples[(x*np):ns]
	for i in range(0,len(lsamples)):
		sample = lsamples[i]
		p = Process(target=adjustbed, args=(sample,covFilter))
		p.start()
		logfile.write('Adjusting bed coordinates for: '+ sample+ '\n')
		Pros.append(p)
	for p in Pros:
		p.join()
	
localtime = time.strftime("%a:%H:%M:%S")
logfile.write('Finished adjusting and coverage filtering of junction bed files at: ' + localtime + '\n')

# ******* Remove known Splice Junctions ******** # # Corrected for strandedness #
Pros = []
nloops=math.ceil(ns/float(np))
for x in range(0,int(nloops)):
	if (x+1)*np <= ns:
		lsamples=samples[(x*np):(x+1)*np]
	else:
		if ns < np:
			lsamples=samples
		if ns > np:
			lsamples=samples[(x*np):ns]
	for i in range(0,len(lsamples)):
		sample = lsamples[i]
		p = Process(target=removeknownjunctions, args=(sample,jnFile,junMatch))
		p.start()
		logfile.write('Filtering known junctions from : '+ sample+ '\n')
		Pros.append(p)
	for p in Pros:
		p.join()
os.system('rm *.Ad')
	
localtime = time.strftime("%a:%H:%M:%S")
logfile.write('Finished filtering out of known Splice Junctions at: ' + localtime + '\n')

#********** Collapse overlapping JNs and take Sum *************# # Corrected for strandedness #
Pros = []
nloops=math.ceil(ns/float(np))
for x in range(0,int(nloops)):
	if (x+1)*np <= ns:
		lsamples=samples[(x*np):(x+1)*np]
	else:
		if ns < np:
			lsamples=samples
		if ns > np:
			lsamples=samples[(x*np):ns]
	for i in range(0,len(lsamples)):
		sample = lsamples[i]
		p = Process(target=collapsesort, args=(sample,))
		p.start()
		logfile.write('Collapsing overlapping for: '+ sample+ '\n')
		Pros.append(p)
	for p in Pros:
		p.join()
os.system('rm *.Fl')
		
localtime = time.strftime("%a:%H:%M:%S")
logfile.write('Finished collapsing overlapping junctions at: ' + localtime + '\n')

# ********** Concatenate each sample bed file to a big file ********** #
junctionfile=open('All_Junctions.bed','w') 
junctionfile.close()
logfile.write('Concatenating all bed files\n')
for sample in samples:
	sample = sample.strip()
	os.system('cat '+sample+'.Ad.Fl.Cl >> All_Junctions.bed')


# **********Sort the pooled junction file********** #
logfile.write('Sorting pooled junctions\n')
os.system('sort -V -k1,1 -k2,2n All_Junctions.bed > All_Junctions.bed.So')
	
# ********** Collapse the sorted junction file ********** # 
bedfile = open('All_Junctions.bed.So','r')
bedfileout = open('All_Junctions.bed.So.Cl','w')
logfile.write('Collapsing pooled junctions\n')
d = defaultdict(list)
junctions=bedfile.readlines()
i=0
k=0 # Container
while i < (len(junctions)-1):
	current=junctions[i].split('\t')
	j=i+1
	while j < len(junctions):
		check=junctions[j].split('\t')
		try:
			if((check[0] != current[0]) or (((int(check[1])-int(current[1])) > 500) and ((int(check[2])-int(current[2])) > 500))):
				break
			if((check[0] == current[0]) and (abs(int(check[1]) - int(current[1])) < 100) and (abs(int(check[2]) - int(current[2])) < 100)and (check[5] == current[5])):
				junctions.pop(j)	
			else:
				j+=1;
		except:
			break	
	if len(current) >5:
		d[current[0]+'_'+current[1]+'_'+current[2]+'_'+current[5]].append('')
	s = '\t'.join(current)
	bedfileout.write(s)
	i+=1;
bedfileout.close()
bedfile.close()

localtime = time.strftime("%a:%H:%M:%S")
logfile.write('Finished processing pooles junctions: ' + localtime + '\n')



# ********** Get junction beds to a single co-ordinate system ***************** #
Pros = []
nloops=math.ceil(ns/float(np))
for x in range(0,int(nloops)):
	if (x+1)*np <= ns:
		lsamples=samples[(x*np):(x+1)*np]
	else:
		if ns < np:
			lsamples=samples
		if ns > np:
			lsamples=samples[(x*np):ns]
	for i in range(0,len(lsamples)):
		sample = lsamples[i]
		p = Process(target=processjunctionbed, args=(sample,i,jn_list,ex_list))
		p.start()
		logfile.write('Extracting junction coverage from: ' + sample + '\n')
		Pros.append(p)
	for p in Pros:
		p.join()
os.system('rm *.Fl.Cl')

localtime = time.strftime("%a:%H:%M:%S")
logfile.write('Finished transforming junction beds to a single co-ordinate system at: ' + localtime + '\n')


# aFter this command, jnlist contains 6 lists, each list is for each sample, MT1-3,WT1-3.  
# within each list, it contains strings of 'chr#\start pos\end pos' : 'score'
# basically all the scores of already known junctions found within the list itself  
#
#exlist contians also 6 lists
#but within each list, it contains the same start as each entry for jnlist 'chr#\start pos\end pos' : 
#what comes after differs, supposed to contain the score but it contains first 6 columns of the exon bed files, but surprisingly the score column is just 't0'


#********** Get coverage of 3' Exon from BAM *************# 
pids = []
nloops=math.ceil(ns/float(np))
for x in range(0,int(nloops)):
	if (x+1)*np <= ns:
		lsamples=samples[(x*np):(x+1)*np]
	else:
		if ns < np:
			lsamples=samples
		if ns > np:
			lsamples=samples[(x*np):ns]
	for sample in lsamples:
		sample = sample.strip()
		bamsample=sample
		bamsample=bamsample.split('.')
		cmd='bedtools coverage -a '+sample+'.Ex.bed -b '+bamsample[0]+'.sorted.bam   -hist > '+sample+'.Ex.bed2'
		#cmd='bedtools coverage -abam '+bamsample[0]+'.sorted.bam -b '+sample+'.Ex.bed -hist > '+sample+'.t.Ex.bed2'
		logfile.write(cmd+'\n')
		p = subprocess.Popen(cmd, shell=True)
		logfile.write('Extracting exon coverage from: ' + sample + '\n')
		pids.insert(0,p) #at index 0 of pids, insert p
	# Wait for the sub-processes #
	exitcode = []
	for i in xrange(len(lsamples)):
		exitcode.insert(len(exitcode),pids[i].wait())
os.system('rm *.Ex.bed')

localtime = time.strftime("%a:%H:%M:%S")
logfile.write('Finished get coverage of 3\' exons from BAMs  at: ' + localtime + '\n')	

#******* removing the extra cols *********#
'''Pros = []
nloops=math.ceil(ns/float(np))
for x in range(0,int(nloops)):
	if (x+1)*np <= ns:
		lsamples=samples[(x*np):(x+1)*np]
	else:
		if ns < np:
			lsamples=samples
		if ns > np:
			lsamples=samples[(x*np):ns]
	for i in range(0,len(lsamples)):
		sample = lsamples[i]
		p = Process(target=removeextracols, args=(sample))
		p.start()
		logfile.write('Removing extra columns from: '+ sample+ '\n')
		Pros.append(p)
	for p in Pros:
		p.join()

pids = []
nloops=math.ceil(ns/float(np))
for x in range(0,int(nloops)):
	if (x+1)*np <= ns:
		lsamples=samples[(x*np):(x+1)*np]
	else:
		if ns < np:
			lsamples=samples
		if ns > np:
			lsamples=samples[(x*np):ns]
	for sample in lsamples:
		sample = sample.strip()
		cutcol='cut -f1,2,3,4,5,6,13,14,15,16- '+sample+'.t.Ex.bed2 > '+sample+'.Ex.bed2'
		p = subprocess.Popen(cutcol, shell=True)
		logfile.write('Removing extra columns from: ' + sample + '\n')
		pids.insert(0,p) #at index 0 of pids, insert p
	# Wait for the sub-processes #
	exitcode = []
	for i in xrange(len(lsamples)):
		exitcode.insert(len(exitcode),pids[i].wait())
#os.system('rm *.t.Ex.bed2')
localtime = time.strftime("%a:%H:%M:%S")
logfile.write('Finished removing extra columns from Ex.bed2 files at: ' + localtime + '\n')
'''

# ********** Get exon coverage to a single co-ordinate system ***************** #
Pros = []
nloops=math.ceil(ns/float(np))
for x in range(0,int(nloops)):
	if (x+1)*np <= ns:
		lsamples=samples[(x*np):(x+1)*np]
	else:
		if ns < np:
			lsamples=samples
		if ns > np:
			lsamples=samples[(x*np):ns]
	for i in range(0,len(lsamples)):
		sample = lsamples[i]
		p = Process(target=processexonbed, args=(sample,i,ex_list))
		p.start()
		logfile.write('Processing exon coverage for: ' + sample + '\n')
		Pros.append(p)
	for p in Pros:
		p.join()
os.system('rm *.bed2')
		
localtime = time.strftime("%a:%H:%M:%S")
logfile.write('iFinished transforming exon coverage beds to a single co-ordinate system at: ' + localtime + '\n')

#********check what is within the two lists*********#
'''
with open('jnlist.txt', 'w') as f:
	for item in jn_list:
		print('jnlist')
		f.write("%s\n" % item)
with open('exlist.txt', 'w') as f:
	for item in ex_list:
		print('exlist')
		f.write("%s\n" % item)

localtime = time.strftime("%a:%H:%M:%S")
logfile.write('did you make the list??? ' + localtime + '\n')

'''

# ********** Match to individual sample and extract count ratios ********** # 
filewrite=open('BBInput.txt','w')
fileread=open('All_Junctions.bed.So.Cl','r')
all_junctions=fileread.readlines()
for line in all_junctions:
	line.strip()
	data=line.split('\t')
	key=data[0]+'\t'+data[1]+'\t'+data[2]
	jn_line=''
	ex_line=''
	temp_jn_dict={}
	temp_ex_dict={}
	for i in range(0,ns):
		temp_jn_dict=jn_list[i]
		temp_ex_dict=ex_list[i]
		try:
			jn_line=jn_line+temp_jn_dict[key]+'\t'
		except:
			jn_line=jn_line+'0\t'
			pass
		try:
			ex_line=ex_line+temp_ex_dict[key]+'\t'
		except:
			ex_line=ex_line+'0\t'
			pass
	
	#print "jn_line: "
	#print jn_line+"|"
	#print "ex_line: "
	#print ex_line+"|"
         
	# Check for consistency of Jn expression # 
	jn_array=jn_line.split('\t')
	ec_array=ex_line.split('\t')
	printstatus1 = 1 # Status of Cr
	printstatus2 = 1 # Status of Tr
	
	# Count number of zeros in Cr samples#
	count1 = 0
	for check in  range(0,ncr): # Note #
		if int(jn_array[check]) >= covFilter:  # Check for 0s in exon expression due to read coverage cutoff
			count1 = count1 +1
	
	# Check for 50% >= filter cutoff 		
	if count1 >= math.ceil(ncr/float(2)):
		printstatus1 = 1 
	else:
		printstatus1 = 0
	
	# Count number of zeros in Tr samples#
	count2 = 0
	for check in  range(ncr,ns): # Note #
		if int(jn_array[check]) >= covFilter:  # Check for 0s in exon expression due to read coverage cutoff
			count2 = count2 +1
	
	# Check for 50% >= filter cutoff 		
	if count2 >= math.ceil(ntr/float(2)):
		printstatus2 = 1
	else:
		printstatus2 = 0
	
	# Do correction to EX coverage if less than JN coverage #
	jc=jn_array
	ec=ec_array
	for p in range(0,ns):
		if int(jc[p]) > int(ec[p]):
			ec[p]=jc[p]
	#with open('jn.txt', 'w') as f:
	#	for item in jc:
	#		f.write("%s\n" % item)
	#with open('ex.txt', 'w') as f:
	#	for item in ec:
	#		f.write("%s\n" % item)
	
	# ****************** For gain of junctions ****************** #
	if printstatus1 == 0 and printstatus2 ==1:
		
		# Get average of non zero Ex counts for Tr samples#
		avgcov=0
		for p in range(ncr,ns):
			if int(ec[p]) >= covFilter:
				avgcov=avgcov+int(ec[p])
		
		try:
			avgcov=avgcov/count2
			# Do correction to EX coverage if 0 #
			for p in range(0,ns):
				if int(ec[p]) < covFilter:
					ec[p]=str(avgcov)
			# Prepare the printing line #
			jn_line='\t'.join(jc)
			ex_line='\t'.join(ec)		
			filewrite.write(data[0]+':'+data[1]+'-'+data[2]+'_'+data[5]+'\t'+jn_line+'\t'+ex_line+'\tJN-Gain\n')	
		except:
			logfile.write("ERROR check required at\t"+line+'\n'+jn_line+'\n'+ex_line+'\t'+str(avgcov)+'\t'+str(count2)+'\t'+str(ntr)+'\n')
			pass
	
	# ****************** For loss of junctions ****************** #
	if printstatus1 == 1 and printstatus2 ==0:
		
		# Get average of non zero Ex counts for Cr samples#
		avgcov=0
		for p in range(0,ncr):
			if int(ec[p]) >= covFilter:
				avgcov=avgcov+int(ec[p])
		
		try:
			avgcov=avgcov/count1
			# Do correction to EX coverage if 0 #
			for p in range(0,ns):
				if int(ec[p]) < covFilter:
					ec[p]=str(avgcov)
			# Prepare the printing line #
			jn_line='\t'.join(jc)
			ex_line='\t'.join(ec)		
			filewrite.write(data[0]+':'+data[1]+'-'+data[2]+'_'+data[5]+'\t'+jn_line+'\t'+ex_line+'\tJN-Loss\n')
		except:
			logfile.write("ERROR check required at\t"+line+'\n'+jn_line+'\n'+ex_line+'\t'+str(avgcov)+'\t'+str(count1)+'\t'+str(ncr)+'\n')
			pass
	
	# ****************** For differential junctions ****************** #
	if printstatus1 == 1 and printstatus2 == 1:
		
		# Get average of non zero Ex counts for Cr samples#
		avgcov=0
		for p in range(0,ncr):
			if int(ec[p]) >= covFilter:
				avgcov=avgcov+int(ec[p])
		try:
			avgcov=avgcov/count1
			# Do correction to EX coverage if 0 #
			for p in range(0,ncr):
				if int(ec[p]) < covFilter:
					ec[p]=str(avgcov)
		except:
			logfile.write("ERROR check required at\t"+line+'\n'+jn_line+'\n'+ex_line+'\t'+str(avgcov)+'\t'+str(count1)+'\t'+str(ncr)+'\n')
			pass
		
		# Get average of non zero Ex counts for Tr samples#
		avgcov=0
		for p in range(ncr,ns):
			if int(ec[p]) >= covFilter:
				avgcov=avgcov+int(ec[p])
		
		try:
			avgcov=avgcov/count2
			# Do correction to EX coverage if 0 #
			for p in range(ncr,ns):
				if int(ec[p]) < covFilter:
					ec[p]=str(avgcov)
		except:
			logfile.write("ERROR check required at\t"+line+'\n'+jn_line+'\n'+ex_line+'\t'+str(avgcov)+'\t'+str(count2)+'\t'+str(ntr)+'\n')
			pass
		
		# Prepare the printing line #
		jn_line='\t'.join(jc)
		ex_line='\t'.join(ec)
		filewrite.write(data[0]+':'+data[1]+'-'+data[2]+'_'+data[5]+'\t'+jn_line+'\t'+ex_line+'\tJN-Diff\n')	
filewrite.close()
fileread.close()
os.system('rm All_Junctions.bed *.So *.So.Cl')

localtime = time.strftime("%a:%H:%M:%S")
logfile.write('Finished estracting junction and exon coverage counts at: ' + localtime + '\n')

# ********** Execute Beta binomial test and do Multiple testing correction  ********** # 
pids = []
#ncr=len(cr_samples) 3
#ntr=len(tr_samples) 3
cmd='Rscript BBTestP.r '+str(ncr)+' '+str(ntr)+' BBInput.txt BBOutput.txt > Log_BBTest.txt'
p = subprocess.Popen(cmd, shell=True)
pids.insert(0,p)
logfile.write('Started '+cmd+'\n')

# Wait for the sub-processes #
for i in xrange(1):
	pids[i].wait()

localtime = time.strftime("%a:%H:%M:%S")
logfile.write('Finished Beta Binomial test and multiple testing correction at: ' + localtime + '\n')

#************** Processing BBoutput *******************#
Pros = []
sample = 'BBOutput.txt'
p = Process(target=processbboutput, args=(sample,ns))
p.start()
logfile.write('Processing : ' + sample + '\n')
Pros.append(p)
for p in Pros:
	p.join()

localtime = time.strftime("%a:%H:%M:%S")
logfile.write('Finished processing BBOutput at: ' + localtime + '\n')

#************** Mapping junctions to Genes *******************#
pids = []

cmd='bedtools intersect -a BBOutput.txt -b '+gmFile+'  -wa -wb -s > CrypticJunctions.txt'
p = subprocess.Popen(cmd, shell=True)
pids.insert(0,p)
logfile.write('Started '+cmd+'\n')
	

# Wait for the sub-processes #
for i in xrange(1):
	pids[i].wait()

localtime = time.strftime("%a:%H:%M:%S")
logfile.write('Finished mapping junctions to genes at: ' + localtime + '\n')

#os.system('rm BBInput* BBOutput* Log_BBTest*')


#************** Screen for multi gene Jns *******************#
fileread=open('CrypticJunctions.txt','r')
#<edit>#
os.system("cp CrypticJunctions.txt CrypticJunctions1.txt")
#<edit>#
filewrite=open('CrypTemp1.txt','w')
# Get Counts #
dict_count={}
lines=fileread.readlines()
fileread.close()
for line in lines:
	line.strip('\n')
	data=line.split("\t")
	key=data[0]+'\t'+data[1]+'\t'+data[2]
	if key in dict_count.keys():
		dict_count[key]=dict_count[key]+1
	else:
		dict_count[key]=1
# Change the class and prit #
for line in lines:
	line.strip('\n')
	data=line.split("\t")
	key=data[0]+'\t'+data[1]+'\t'+data[2]
	if int(dict_count[key]) >1:
		if data[12]=="JN-Gain":
			data[12]="JN-GaOv"
		if data[12]=="JN-Loss":
			data[12]="JN-LoOv"
		if data[12]=="JN-Diff":
			data[12]="JN-DfOv"  
	line="\t".join(data)
	filewrite.write(line)
filewrite.close()
#<edit>#
os.system("mv CrypTemp1.txt CrypticJunctions.txt")


localtime = time.strftime("%a:%H:%M:%S")
logfile.write('Finished screening gene fusion events at: ' + localtime + '\n')

#************** Screen and Classify *******************#
fileread1=open('CrypticJunctions.txt','r')
fileread2=open(aeFile,'r')
filewrite1=open('JN-Gain.txt','w')
filewrite2=open('JN-Loss.txt','w')
filewrite3=open('JN-Diff.txt','w')
# <edit> #
#filewrite4=open('JN-Ovlp.txt','w') # Not to Change Gene Fusion Events #

printl=[]
printl.append("JnChr\tJnStart\tJnStop\tJnStrand\t")
printl.append('')
printl.append("JnExp%\tP-value\tP-value(B)\tP-value(BH)\tP-value(BY)\tJnType\tJnClass\tMapGeneChr\tGStart\tGStop\tMapGene\tGStrand\n")
for pi in range(0,ncr):
	pii=pi+1
	printl[1]=printl[1]+"CJn"+str(pii)+"\t"
for pi in range(0,ntr):
	pii=pi+1
	printl[1]=printl[1]+"TJn"+str(pii)+"\t"
for pi in range(0,ncr):
	pii=pi+1
	printl[1]=printl[1]+"CEx"+str(pii)+"\t"
for pi in range(0,ntr):
	pii=pi+1
	printl[1]=printl[1]+"TEx"+str(pii)+"\t"

printline="".join(printl)
filewrite1.write(printline)
filewrite2.write(printline)

printl[2]="P-value\tP-value(B)\tP-value(BH)\tP-value(BY)\tJnType\tJnClass\tMapGeneChr\tGStart\tGStop\tMapGene\tGStrand\n";
printline="".join(printl)
filewrite3.write(printline)
#filewrite4.write(printline)

lines1=fileread1.readlines()
lines2=fileread2.readlines()
fileread1.close()
fileread2.close()

for line1 in lines1:
	data1=line1.split("\t")
	status1=0
	status2=0
	for line2 in lines2:
		line2.strip('\n')
		data2=line2.split("\t")
		if data1[0] == data2[0] and (int(data1[1]) <=int(data2[1]) and int(data1[2]) >= int(data2[1])):
			status1=1
		if data1[0] == data2[0] and (int(data1[1]) <=int(data2[2]) and int(data1[2]) >= int(data2[2])):
			status2=1
		if status1 ==1 and status2 ==1:
			break
	if (status1 + status2) == 2:
		data1.insert(13,'CD-CD')
	if (status1 + status2) == 1:
		data1.insert(13,'CD-CY')
	if (status1 + status2) == 0:
		data1.insert(13,'CY-CY')

	# Remove padded zeros #
	#print data1
	data1.pop(3)
	#print data1
	data1.pop(3)
	#print data1
	data1=data1[0:-6]
	data1.pop(-2)
	#data1[4]=data1[4].replace(" ","\t")
	#data1[5]=data1[5].replace(" ","\t")
	#print data1
	
		
	if data1[10] =="JN-Gain" or data1[10] =="JN-GaOv":
		tdata=data1[4].split(" ")
		data1[4]="\t".join(tdata)
		tdata=map(int,tdata)
		sumtdata=0
		for pi in range(0,ntr):
			sumtdata=sumtdata+tdata[ncr+pi]
		tavg=float(sumtdata)/float(ntr)
		
		edata=data1[5].split(" ")
		data1[5]="\t".join(edata)
		edata=map(int,edata)
		sumedata=0
		for pi in range(0,ntr):
			sumedata=sumedata+edata[ncr+pi]
		eavg=float(sumedata)/float(ntr)
		
		tavg=tavg/eavg
		data1.insert(6,str(tavg))
		
		line1="\t".join(data1)
		filewrite1.write(line1+"\n")
	
	if data1[10] =="JN-Loss" or data1[10] =="JN-LoOv":
		tdata=data1[4].split(" ")
		data1[4]="\t".join(tdata)
		tdata=map(int,tdata)
		sumtdata=0
		for pi in range(0,ncr):
			sumtdata=sumtdata+tdata[pi]
		tavg=float(sumtdata)/float(ncr)
		
		edata=data1[5].split(" ")
		data1[5]="\t".join(edata)
		edata=map(int,edata)
		sumedata=0
		for pi in range(0,ncr):
			sumedata=sumedata+edata[pi]
		eavg=float(sumedata)/float(ncr)
		
		tavg=tavg/eavg
		data1.insert(6,str(tavg))
		
		line1="\t".join(data1)
		filewrite2.write(line1+"\n")
	
	if data1[10] =="JN-Diff" or data1[10] =="JN-DfOv":
		tdata=data1[4].split(" ")
		data1[4]="\t".join(tdata)
		edata=data1[5].split(" ")
		data1[5]="\t".join(edata)
		edata=map(int,edata)
		line1="\t".join(data1)
		filewrite3.write(line1+"\n")
	
	#<edit>#
	# Not to Change Gene Fusion Events #
	#if data1[10] =="Jn-GaOv" or data1[10] =="Jn-LoOv" or data1[10] =="Jn-DfOv" :
	#	line1="\t".join(data1)
	#	filewrite4.write(line1)

os.system("rm CrypticJunctions.txt CrypticJunctions1.txt")	

localtime = time.strftime("%a:%H:%M:%S")
logfile.write('Finished annotating splicing events at: ' + localtime + '\n')

filewrite1.close()
filewrite2.close()
filewrite3.close()

# ***************** Print End of Script ************* #
localtime = time.strftime("%a:%H:%M:%S")
logfile.write('Finished CrypSplice at : ' + localtime + '\n')
logfile.close()