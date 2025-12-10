
def annotateFromMtDB(fileName, relpath):
	file = open(fileName, 'r')
	newfile = open(fileName + "_mtDB", 'w')
	#mtDBFile = open("/research/groups/kundugrp/projects/Mitochondria/common/FilesForPipeline/annotationFiles/mtDB_database.txt", 'r')
        mtDBFile = open(relpath + "/Mito_WGS_Nextflow/Nextflow_mitochondra_wgs/annotationFiles/mtDB_database.txt", 'r')
	mtDBMap = {}
	headers = mtDBFile.readline()
	for line in mtDBFile:
		line = line.strip()
		cols = line.split('\t')
		pos = cols[0]
		refBase = cols[1]
		A = cols[2]
		if A == ' ' or A == '': A = 0
		G = cols[3]
		if G == ' ' or G == '': G = 0
		C = cols[4]
		if C == ' ' or C == '': C = 0
		T = cols[5]
		if T == ' ' or T == '': T = 0
		baseVals = {}
		baseVals['A'] = int(A)
		baseVals['G'] = int(G)
		baseVals['C'] = int(C)
		baseVals['T'] = int(T)
		mtDBMap[pos] = [refBase, baseVals]
	mtDBFile.close()

	for line in file:
		line = line.rstrip()
		cols = line.split('\t')
		pos = cols[5]
		ref = cols[7]
		altBase = cols[8]
		found = False
		if pos in mtDBMap:
			refBase, baseVals = mtDBMap[pos]
                        print(mtDBMap[pos])
			#if ref == refBase and baseVals[altBase] != '0':
                        if ref==refBase and baseVals in ['A','T','C','G'] and baseVals[altBase] != '0': ### TCC
				total = 0
				for key in baseVals:
					total += baseVals[key]
				altNum = baseVals[altBase]
				found = True

		if not found:
			phrase1 = 'notinDB'
			phrase2 = 'notinDB'
		else:
			phrase1 = pos + ref + '>' + altBase
			phrase2 = '"' + str(altNum) + '/' + str(total) + '"' 

		newline = line
		while len(cols) < 19:
			cols.append('')
			newline += '\t'
		newline += '\t' + phrase1 + '\t' + phrase2 + '\n'

		newfile.write(newline)

	newfile.close()
	file.close()

def annotateFromMitoMap(fileName, relpath):
	file = open(fileName, 'r')
	newfile = open(fileName + "_annotatedMitoMapDS", 'w')
        #diseaseFile = open('/research/groups/kundugrp/projects/Mitochondria/common/FilesForPipeline/annotationFiles/MitoMapDiseaseAssociationsAll.txt', 'r')
	#diseaseFile = open('/research/groups/cab/projects/automapper/common/wzhang42/Mitochondra_Pipeline/annotationFiles/MitoMapDiseaseAssociationsAll.txt', 'r')
	diseaseFile = open(relpath + "/Mito_WGS_Nextflow/Nextflow_mitochondra_wgs/annotationFiles/MitoMapDiseaseAssociationsAll.txt", 'r')
        hash = {}
	headers = diseaseFile.readline()
	for line in diseaseFile:
		cols = line.split('\t')
		hash[cols[3]] = (cols[2], cols[8])

	diseaseFile.close()
	for line in file:
		line = line.strip()
		cols = line.split('\t')
		newline = line
		if len(cols) < 17:
			newline += '\t'
		key = cols[7] + cols[6] + cols[8]
		info = ''
		if key in hash:
			info = hash[key]
		if info != '':
			newline += '\t' + info[0] + '\t' + info[1]
		newfile.write(newline + '\n')

	file.close()

def annotateNonSynSnpsSummary(sumFileName, relpath):
        sumFile = open(sumFileName, 'r')
        sumFileFinal = open(sumFileName + "_nonSyn", 'w')
	sampleHash = {}
	samples = []
        for line in sumFile:
                line = line.strip()
                cols = line.split('\t')
                sampleName = cols[0]
                pos = cols[5]
		type = cols[8]
		if (sampleName, pos, type) not in sampleHash:
			if sampleName not in samples:
                		#fileName =  "/research/groups/kundugrp/projects/Mitochondria/common/FilesForPipeline/annotationFiles/MockData.txt.exonic_variant_function"
				fileName =  relpath + "/Mito_WGS_Nextflow/Nextflow_mitochondra_wgs/annotationFiles/MockData.txt.exonic_variant_function"
                                samples.append(sampleName)
                		file = open(fileName, 'r')
                		for line2 in file:
					line2 = line2.rstrip()
                        		cols2 = line2.split('\t')
					type2 = cols2[7]
					sampleHash[(sampleName, cols2[4], type2)] = [cols2[1], cols2[2]]
				file.close()
	
		if (sampleName, pos, type) in sampleHash:
			nonSyn, exon = sampleHash[(sampleName, pos, type)]
               		newline = nonSyn + '\t' + exon + '\t'
                else:
                        newline = '\t\t'
                newline += line
                sumFileFinal.write(newline + '\n')
        sumFileFinal.close()
        sumFile.close()


def addMutationAssessorAnnotations(sumFileName, relpath):
	#fileName = "/research/groups/kundugrp/projects/Mitochondria/common/FilesForPipeline/annotationFiles/MA_scores_rel3_hg19_chrM_full.txt"
        fileName = relpath + "/Mito_WGS_Nextflow/Nextflow_mitochondra_wgs/annotationFiles/MA_scores_rel3_hg19_chrM_full.txt"
	file = open(fileName, 'r')
	mutHash = {}
	headers = file.readline()
	for line in file:
		line = line.rstrip()
		cols = line.split('\t')
		vals = cols[0].split(',')
		mutHash[(vals[2], vals[4][0])] = cols
	file.close()

	#This file was downloaded from the mutation assessor site and modified by me to exclude some columns
	#fileName = "/research/groups/kundugrp/projects/Mitochondria/common/FilesForPipeline/annotationFiles/MutationAssessor.output_fromGenes_PCGPnonSynOnly.txt"
        fileName = relpath + "/Mito_WGS_Nextflow/Nextflow_mitochondra_wgs/annotationFiles/MutationAssessor.output_fromGenes_PCGPnonSynOnly.txt"
	file = open(fileName, 'r')
	mutHashByGene = {}
	headers = file.readline()
	for line in file:
		line = line.rstrip()
		cols = line.split('\t')
		aaChange = cols[0]
		gene = cols[1]
		impact = cols[2]
		score = cols[3]
		mutHashByGene[(gene, aaChange)] = (impact, score)
	file.close()

        sumFile = open(sumFileName, 'r')
        sumFileFinal = open(sumFileName + "_MutationAssessor", 'w')
        for line in sumFile:
                cols = line.split('\t')
                pos = cols[7]
                mutAllele = cols[10]
                newline = ''
		gene = ''
		aaChange = ''
		impact = 'NA'
		score = 'NA'
		if cols[1].find('.') != -1:
			aaChange = cols[1][cols[1].rindex('.')+1:-1]
			gene = cols[1][:cols[1].index(':')]
		if cols[0] == 'nonsynonymous SNV':
			if (pos, mutAllele) in mutHash:
				cols2 = mutHash[(pos, mutAllele)]
				#if not 8 fields, then it doesn't have an impact & score or probably even an AA change value
				if len(cols2) == 8:
					impact = cols2[6]
					score = cols2[7]
					fromMutAssAAchange = cols2[5]
				#use alternate file instead because values aren't available in normal file
				else:
					if (gene, aaChange) in mutHashByGene:
						impact, score = mutHashByGene[(gene, aaChange)]
						fromMutAssAAchange = aaChange
					else:
						#no impact value found in either file, so put NA
						impact = 'NA'	
						score = 'NA'
						fromMutAssAAchange = ''
				if aaChange != fromMutAssAAchange:
					if (gene, aaChange) in mutHashByGene:
						impact2, score2 = mutHashByGene[(gene, aaChange)]
						if impact2 != '':
							#If alt file has value, use alt value = more reliable since amino acid change correct
							impact = impact2
							score = score2
						#if both blank, put NA
						if impact2 == '' and impact == '':
							impact = 'NA'
							score = 'NA'
						#if impact2 != impact and impact != '' and impact2 != '':
						#	print aaChange, gene, impact2, score2, impact, score
					
                        if len(cols) == 18:
                        	line = line[:-1] + '\t'
                        else:
                        	line = line[:-1]
                        newline = line + '\t' + aaChange + '\t' + gene + '\t' + impact + '\t' + score + '\n'
		else:
			#stopgain or synonymous - then list aaChange and gene, no impact or score should be available
			if cols[4] == 'exonic':
                        	if len(cols) == 18:
                        		line = line[:-1] + '\t'
                        	else:
                        		line = line[:-1]
                        	newline = line + '\t' + aaChange + '\t' + gene + '\n'
			#all non-exonic
			else:
                        	newline = line
               	sumFileFinal.write(newline)
        sumFileFinal.close()
        sumFile.close()

def annotateFromLevinEtAl(sumFileName, relpath):
        sumFile = open(sumFileName, 'r')
        sumFileFinal = open(sumFileName + "_LevinEtAl", 'w')
        #fileName =  "/research/groups/kundugrp/projects/Mitochondria/common/FilesForPipeline/annotationFiles/LevinEtAlTable.txt"
        fileName = relpath + "/Mito_WGS_Nextflow/Nextflow_mitochondra_wgs/annotationFiles/LevinEtAlTable.txt"
        file = open(fileName, 'r')
	levinHash = {}
	for line in file:
		cols = line.split('\t')
		pos = cols[1]
		alleles = cols[4]
		refAllele = alleles[0]
		altAllele = alleles[-1]
		levinHash[pos] = refAllele + altAllele
        file.close()
        for line in sumFile:
		line = line.rstrip()
                cols = line.split('\t')
                pos = cols[7]
                refAllele = cols[9]
                mutAllele = cols[10]
		inLevin = 'no'
		if pos in levinHash and levinHash[pos] == refAllele + mutAllele:
			inLevin = 'yes'
                newline = ''
		for col in cols[:23]:
			newline += col + '\t'
		newline += inLevin + '\t'
		for col in cols[23:]:
			newline += col + '\t'
		newline = newline[:-1] + '\n'
		sumFileFinal.write(newline)

        sumFileFinal.close()
        sumFile.close()


def finalizeFile(filename):
	file = open(filename, 'r')
	newfile = open(filename + "_finalized.txt", 'w')

	for line in file:
		line = line.rstrip()
		cols = line.split('\t')
		if cols[9] == '-' or cols[10] == '-':
			continue
		newline = cols[2] + '\t' + cols[0] + '\t' + cols[1] + '\t' + cols[3] + '\t' + cols[4] + '\t' + cols[5] + '\t'	
		for ind in xrange(8, len(cols)):
			newline += cols[ind] + '\t'
		newline = newline[:-1] + '\n'
		newfile.write(newline)
	newfile.close()
	file.close()

def addRnaInfo(filename, relpath):
	#rnaFile = open("/research/groups/kundugrp/projects/Mitochondria/common/FilesForPipeline/annotationFiles/rRNA_mutationPotentialDisruption.txt", 'r')
        rnaFile = open(relpath + "/Mito_WGS_Nextflow/Nextflow_mitochondra_wgs/annotationFiles/rRNA_mutationPotentialDisruption.txt", 'r')
	rnaHash = {}
	for line in rnaFile:
		line = line.rstrip()
		cols = line.split('\t')
		potential = cols[0]
		pos = cols[1]
		ref = cols[2]
		alt = cols[3]
		
		key = ref + pos + alt
		rnaHash[key] = potential
	rnaFile.close()
	#trnaFile = open("/research/groups/kundugrp/projects/Mitochondria/common/FilesForPipeline/annotationFiles/tRNA_mutationBenignDeleterious.txt", 'r')
        trnaFile = open(relpath + "/Mito_WGS_Nextflow/Nextflow_mitochondra_wgs/annotationFiles/tRNA_mutationBenignDeleterious.txt", 'r')
	trnaHash = {}
	for line in trnaFile:
		line = line.rstrip()
		cols = line.split('\t')
		pos = cols[3]
		alt = cols[4]
		potential = cols[5]
		trnaHash[pos+alt] = potential
	trnaFile.close()

	bothP = 0
	bothN = 0
	mtPotherN = 0
	mtNotherP = 0
        #trnaFile = open("/research/groups/kundugrp/projects/Mitochondria/common/FilesForPipeline/annotationFiles/MitoTipRawScores_Aug2017.txt", 'r')
	#trnaFile = open("/research/groups/cab/projects/automapper/common/wzhang42/Mitochondra_Pipeline/annotationFiles/MitoTipRawScores_Oct2024.txt", 'r')
        trnaFile = open(relpath + "/Mito_WGS_Nextflow/Nextflow_mitochondra_wgs/annotationFiles/MitoTipRawScores_Oct2024.txt", 'r')	

	mitoTipHash = {}
	for line in trnaFile:
		line = line.rstrip()
		cols = line.split('\t')
		pos = cols[0]
		alt = cols[2]
		potential = cols[4]
		if alt == 'del':
			continue
		mitoTipHash[pos+alt] = potential
		'''
		if pos+alt in trnaHash:
			if potential == 'pathogenic' and trnaHash[pos+alt] == 'D':
				bothP += 1
			elif potential == 'benign' and trnaHash[pos+alt] == 'B':
				bothN += 1
			elif potential == 'pathogenic' and trnaHash[pos+alt] != 'D':
				mtPotherN += 1
			elif potential == 'benign' and trnaHash[pos+alt] != 'B':
				mtNotherP += 1
		else:
			print 'Not in original file', cols
		'''
	trnaFile.close()
	#print bothP, bothN, mtPotherN, mtNotherP

	file = open(filename, 'r')
	newFile = open(filename[:-4] + "_rnaInfo.txt", 'w')
	for line in file:
		newline = line
		cols = line.split('\t')
		pos = cols[6]
		ref = cols[7]
		alt = cols[8]
		if cols[3] == 'rRNA':
			key = ref + pos + alt
			call = 'NotFound'
			if key in rnaHash:
				call = rnaHash[key]
			newline = line.rstrip() + '\t\t\t' + call + '\n'
			
				
		elif cols[3] == 'tRNA':
			key = pos + alt
			call = 'NotFound'
			if key in trnaHash:
				call = trnaHash[key]
				if call == 'B':
					call = 'benign'
				elif call == 'D':
					call = 'deleterious'
			mitoTipCall = 'NotFound'
			if key in mitoTipHash:
				mitoTipCall = mitoTipHash[key]
			newline = line.rstrip() + '\t\t' + call + '\t' + mitoTipCall + '\n'
		newFile.write(newline)
	newFile.close()


def addMMpopFreq(filename, relpath):
	#mmFile = open('/research/groups/kundugrp/projects/Mitochondria/common/FilesForPipeline/annotationFiles/MitoMapFrequencies_updatedMay_31_2018_45494totalFullSamples.txt', 'r')
        mmFile = open(relpath + "/Mito_WGS_Nextflow/Nextflow_mitochondra_wgs/annotationFiles/MitoMapFrequencies_updatedMay_31_2018_45494totalFullSamples.txt", 'r')
	headers = mmFile.readline()
	popFreqMap = {}
	for line in mmFile:
		line = line.rstrip()
		pos, ref, alt, cnt = line.split('\t')
		freq = cnt + '/45494'
		popFreqMap[ref + pos + alt] = freq

	mmFile.close()
	newFile = open(filename[:-4] + "_MMpopFreq.txt", 'w')
	file = open(filename, 'r')
	#headers = file.readline()
	for line in file:
		line = line.rstrip()
		cols = line.split('\t')
		newline = ''
		for col in cols[:22]:
			newline += col + '\t'
		pos = cols[6]
		code = cols[7] + pos + cols[8]
		freq = ''
		if code in popFreqMap:
			freq = popFreqMap[code]
		newline += freq + '\t'
		for col in cols[22:]:
			newline += col + '\t'
		newFile.write(newline[:-1] + '\n')
	newFile.close()

def removeSnp1Col(filename):
	file = open(filename, 'r')

	newfilename = filename[:filename.index('annotated')] + 'annotated.txt'
	newfile = open(newfilename, 'w')
	for line in file:
		cols = line.split('\t')
		newline = ''
		for col in cols:
			if col != 'SNP 1':
				newline += col + '\t'
		newfile.write(newline[:-1])
	newfile.close()
	file.close()
		
def main():
	# Find the path at which this python script is executing
	import os
	script_dir = os.path.dirname(os.path.abspath(__file__))
	# Now get three directories up from this path
	relpath = os.path.dirname(script_dir)
	relpath = os.path.dirname(relpath)
	relpath = os.path.dirname(relpath)
	filename = "./all.maf.fisher"
	annotateFromMitoMap(filename, relpath)
	annotateFromMtDB(filename + "_annotatedMitoMapDS", relpath)
	annotateNonSynSnpsSummary(filename + "_annotatedMitoMapDS_mtDB", relpath)
	addMutationAssessorAnnotations(filename + "_annotatedMitoMapDS_mtDB_nonSyn", relpath)
	annotateFromLevinEtAl(filename + "_annotatedMitoMapDS_mtDB_nonSyn_MutationAssessor", relpath)
	finalizeFile(filename + "_annotatedMitoMapDS_mtDB_nonSyn_MutationAssessor_LevinEtAl")
	addRnaInfo(filename + "_annotatedMitoMapDS_mtDB_nonSyn_MutationAssessor_LevinEtAl_finalized.txt", relpath)
	addMMpopFreq(filename + "_annotatedMitoMapDS_mtDB_nonSyn_MutationAssessor_LevinEtAl_finalized_rnaInfo.txt", relpath)
	removeSnp1Col(filename + "_annotatedMitoMapDS_mtDB_nonSyn_MutationAssessor_LevinEtAl_finalized_rnaInfo_MMpopFreq.txt")


main()
