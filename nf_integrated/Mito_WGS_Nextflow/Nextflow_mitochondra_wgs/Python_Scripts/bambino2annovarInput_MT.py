#!/usr/bin/python 

import sys

def build_annovar_input_file(bambinoFile):
	bFile = open(bambinoFile, 'r')
	newfile = open(bambinoFile + ".annovar_input", 'w')
	headers = bFile.readline()
	headCols = headers.split('\t')
	altFwdInd = headCols.index('alternative_fwd_count')
	altRevInd = headCols.index('alternative_rev_count')
	normalRefInd = headCols.index('reference_normal_count')
	normalAltInd = headCols.index('alternative_normal_count')
	tumorRefInd = headCols.index('reference_tumor_count')
	tumorAltInd = headCols.index('alternative_tumor_count')
	typeInd = headCols.index('Type')
	endPosInd = headCols.index('Pos')
	refAlleleInd = headCols.index('Chr_Allele')
	altAlleleInd = headCols.index('Alternative_Allele')
	sizeInd = headCols.index('Size')
	for line in bFile:
		line = line.rstrip()
		cols = line.split('\t')
		mutT = int(cols[tumorAltInd])
		totT = int(cols[tumorAltInd]) + int(cols[tumorRefInd])
		mutN = int(cols[normalAltInd])
		totN = int(cols[normalRefInd]) + int(cols[normalAltInd])
		type = cols[typeInd]
		size = cols[sizeInd]
		endPos = int(cols[endPosInd])
		ref = cols[refAlleleInd]
		alt = cols[altAlleleInd]
		if type == 'deletion':
			alt = '-'
			endPos = endPos-1+int(cols[4])
		elif type == 'insertion':
			ref = '-'
	
		newline = "MT\t%d\t%d\t%s\t%s\t%s\t%d/%d\t%d/%d\n" % (int(cols[endPosInd]), endPos, ref, alt, type + " " + size, mutT, totT, mutN, totN) 
		if cols[0] == 'chrM.8612':
			print cols
			print newline
		newfile.write(newline)
	newfile.close()
		

if __name__ == "__main__":
	if len(sys.argv) < 2:
		print "Usage: bambino2annovarInput_MT.py <bambino_file_name>"
	else:
		fileName =  sys.argv[1]
		build_annovar_input_file(fileName)

