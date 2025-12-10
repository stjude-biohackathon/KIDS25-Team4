import scipy
import sys
from scipy import stats

ballPrefix = ['SJBALL', 'SJPHALL', 'SJINF', 'SJHYPO', 'SJHYPER', 'SJETV', 'SJERG', 'SJE2A']

def addSubtypesColumnToSNVs(filename):
 	tumorType = {'SJACT':'Solid', 'SJCBF': 'Leukemia', 'SJEPD':'Brain', 'SJEWS':'Solid', 'SJHYPO':'Leukemia', 'SJMB':'Brain', 'SJOS':'Solid',
                'SJRHB': 'Solid', 'SJAMLM7': 'Leukemia', 'SJCPC':'Brain', 'SJERG': 'Leukemia', 'SJHGG': 'Brain', 'SJINF':'Leukemia', 'SJMEL': 'Solid',
                'SJPHALL': 'Leukemia', 'SJTALL': 'Leukemia', 'SJBALL': 'Leukemia', 'SJE2A': 'Leukemia', 'SJETV': 'Leukemia', 'SJHYPER':'Leukemia',
                'SJLGG': 'Brain', 'SJNBL': 'Solid', 'SJRB': 'Solid'}
        subClass = {'SJBALL011': 'Ph-like', 'SJBALL012': 'Ph-like', 'SJBALL020013': 'Ph-like', 'SJBALL020340': 'SJERG', 'SJBALL020374': 'iAMP21',
        'SJBALL020422': 'Ph-like', 'SJBALL020469': 'SJERG', 'SJBALL020516': 'iAMP21', 'SJBALL020570': 'iAMP21', 'SJBALL020579':'Ph-like', 'SJBALL020589':'Ph-like',
        'SJBALL020595': 'SJERG', 'SJBALL020609': 'SJERG', 'SJBALL020625': 'Ph-like', 'SJBALL020635': 'SJERG', 'SJBALL020649':'Ph-like', 'SJBALL020704': 'Ph-like',
        'SJBALL020828': 'SJERG', 'SJBALL020877': 'Ph-like', 'SJBALL020882': 'Ph-like', 'SJBALL020984':'Ph-like', 'SJBALL021031': 'iAMP21', 'SJBALL021058':'Ph-like',
        'SJBALL021083':'Ph-like', 'SJBALL021108': 'iAMP21', 'SJBALL021130': 'SJERG','SJBALL021170': 'TCF3-HLF', 'SJBALL021305': 'Ph-like', 'SJBALL021358': 'iAMP21',
        'SJBALL021373': 'iAMP21', 'SJBALL021491': 'iAMP21', 'SJBALL021516': 'SJERG', 'SJBALL021893':'TCF3-HLF', 'SJBALL021894': 'TCF3-HLF', 'SJBALL021895': 'TCF3-HLF',
        'SJBALL021896': 'TCF3-HLF','SJBALL021897': 'TCF3-HLF', 'SJBALL021898': 'TCF3-HLF', 'SJBALL021900': 'iAMP21', 'SJBALL021901': 'iAMP21', 'SJBALL063': 'Ph-like',
        'SJBALL101': 'Ph-like', 'SJBALL102': 'iAMP21', 'SJBALL153': 'Ph-like', 'SJBALL231': 'Ph-like', 'SJBALL239': 'Ph-like', 'SJBALL247': 'Ph-like',
        'SJBALL255': 'Ph-like', 'SJBALL263': 'Ph-like', 'SJBALL264': 'Ph-like', 'SJBALL267': 'Ph-like', 'SJBALL271': 'iAMP21', 'SJHYPER003': 'Ph-likeHyp',
        'SJHYPER013': 'Ph-likeHyp', 'SJHYPER021': 'Ph-likeHyp', 'SJHYPER120': 'Ph-likeHyp', 'SJHYPER146': 'Ph-likeHyp', 'SJHYPO020': 'Ph-likeHyp', 'SJHYPO109': 'Ph-likeHyp',
        'SJHYPO110': 'Ph-likeHyp', 'SJHYPO123': 'Ph-likeHyp', 'SJHYPO147': 'Ph-likeHyp'}

        sampleTypes = tumorType.keys()
        sampleTypes.append('iAMP21')
        sampleTypes.append('Ph-like')
        sampleTypes.append('Ph-likeHyp')
        sampleTypes.append('TCF3-HLF')
        sampleTypes.sort()

        allSamples = []
        file = open("/home/cwelsh1/projects/PCGP/Mito/allSampleNames_withNewTall.txt", 'r')
        for line in file:
                line = line.strip()
                allSamples.append(line)
        sampleTypeMap = {}
        for sample in allSamples:
                if sample in subClass:
                        key = subClass[sample]
                        sampleTypeMap[sample] = key
                else:
                        for key in sampleTypes:
				if sample.startswith(key):
					sampleTypeMap[sample] = key
	file.close()

	file = open(filename, 'r')
	newfile = open(filename[:-4] + '_withSubtypes.txt', 'w')
	for line in file:
		line = line.rstrip()
		cols = line.split('\t')
		sample = cols[0]
		if sample == 'Sample':
			newfile.write('Sample\tSubtype\tmt_copy_number_avg\n')
			continue
		else:
			if '_' in sample:
				sample2 = sample[:sample.index('_')]
			else:
				sample2 = sample
			if '-' in sample:
				sample3 = sample[:sample.index('-')]
			else:
				sample3 = sample
		subtype = sampleTypeMap[sample2]
		newline = sample3 + '\t' + subtype + '\t'
		for col in cols[1:]:
			newline += col + '\t'
		newline = newline[:-1] + '\n'
		newfile.write(newline)
	newfile.close()
	file.close()
#This function adds 3 additional columns to the spreadsheet
#Adds Complex if SNV is exonic, tRNA, rRNA, D-loop, or intergenic otherwise
#Adds predicted benign, deleterious, D-loop or unknown column
#Adds flag column for probable sequencing errors
def addComplexAndPredictedImpactAndMarkComSeqErrors(filename):
	file = open(filename, 'r')
	newfile = open(filename[:-4] + '_allColumns.txt', 'w')

	refSeqFile = open('/home/cwelsh1/projects/PCGP/Mito/double-alignment/MT.fa', 'r')
	refSeq = ''
	for line in refSeqFile:
		if line.startswith('>MT'):
			continue
		refSeq += line.rstrip()

	refSeqFile.close()
		
	for line in file:
		line = line.rstrip()
		cols = line.split('\t')
		sample = cols[0]
		print cols
		#pos = int(cols[13])
		#pos = int(cols[8])
		pos = int(cols[7])
                if (pos >= 302 and pos <= 315) or (pos >= 513 and pos <= 525) or (pos >= 3105 and pos <= 3109) or pos == 567:
			flag = 'x'
		elif pos == 16192 or pos == 16183 or pos == 16189:
			flag = 'x'
		else:
			flag = ''
		#gene = cols[9]
		#loc = cols[10]
		gene = cols[3]
		loc = cols[4]
		#gene = cols[4]
		#loc = cols[5]
		if gene == '' and loc == '':
			complex = 'intergenic'
		elif gene == '' and loc != '':
			complex = loc
		else:
			#gene = cols[7]
			gene = cols[6]
			if gene == 'ATP6' or gene == 'ATP8':
				complex = 'Complex 5'
			elif gene == 'CYTB':
				complex = 'Complex 3'
			elif gene == 'COX1' or gene == 'COX2' or gene == 'COX3':
				complex = 'Complex 4'
			elif gene == 'ATP6,COX3':
				complex = 'Complex 4/5'
			else:
				complex = 'Complex 1'

		#synNon = cols[8]
		#if len(cols) > 31:
		#	mutAss = cols[31]
		#synNon = cols[3]
		synNon = cols[2]
		#if len(cols) > 26:
		#	mutAss = cols[26]
		if len(cols) > 25:
			mutAss = cols[25]
		else:
			mutAss = ''
		print loc, mutAss, synNon
		if loc == 'D-loop':
                        typeOfSNV = 'D-loop'
                elif (synNon == 'nonsynonymous SNV' and (mutAss != 'neutral' and mutAss != 'NA')) or synNon == 'stoploss SNV' or synNon == 'stopgain SNV' or (loc == 'tRNA' and mutAss == 'pathogenic') or (loc == 'rRNA' and (mutAss == 'Expectedly' or mutAss == 'Proven')) :
                        typeOfSNV = 'Predicted Pathogenic'
                elif synNon == 'synonymous SNV' or (synNon == 'nonsynonymous SNV' and (mutAss == 'neutral' or mutAss == 'NA')) or (loc == 'tRNA' and mutAss == 'benign') or (loc == 'rRNA' and mutAss == 'NEE'):
                        typeOfSNV = 'Predicted Benign'
                else:
                        typeOfSNV = 'Unknown'

		print typeOfSNV
		#refAllele = cols[9]
		#altAllele = cols[10]
		refAllele = cols[8]
		altAllele = cols[9]
		sequence = refSeq[pos-26:pos-1] + '[' + refSeq[pos-1] + '/' + altAllele + ']' + refSeq[pos:pos+25]
		
		line = flag + '\t' + complex + '\t' + typeOfSNV + '\t' + sequence + '\t' + line + '\n'
		newfile.write(line)
	
	file.close()
	newfile.close()
	

def addModifiedNewGroups(filename, sampleListFile="/home/cwelsh1/projects/PCGP/Mito/allSampleNames_withNewTall.txt"):
	file = open(filename, 'r')
	newfile = open(filename[:-4] + '_withModifiedNewGroups.txt', 'w')
 	groupStats = {}
        groupStats['1'] = [0, {}]
        groupStats['2A'] = [0, {}]
        groupStats['2B'] = [0, {}]
        groupStats['3A'] = [0, {}]
        groupStats['3B'] = [0, {}]
        groupStats['4'] = [0, {}]
        groupStats['5'] = [0, {}]
        groupStats['6'] = [0, {}]

	cvalues = {}

	sampleStats = {}
	sampleFile = open(sampleListFile, 'r')
	for line in sampleFile:
		sample = line.strip()
		sampleStats[sample] = {}
		sampleStats[sample]['1'] = [0,0,0]
               	sampleStats[sample]['2A'] = [0,0,0]
               	sampleStats[sample]['2B'] = [0,0,0]
               	sampleStats[sample]['3A'] = [0,0,0]
               	sampleStats[sample]['3B'] = [0,0,0]
               	sampleStats[sample]['4'] = [0,0,0]
               	sampleStats[sample]['5'] = [0,0,0]
               	sampleStats[sample]['6'] = [0,0,0]
               	sampleStats[sample]['backMutation'] = [0, 0, 0]
	sampleFile.close()

	numDiff = 0
	for line in file:
		line = line.rstrip()
		cols = line.split('\t')
		sample = cols[0]
		germlineCvalue  = cols[3]
		tumorCvalue = cols[4]
		cvalues[sample] = [tumorCvalue, germlineCvalue]
		gC = float(cols[3]) * 5
		tC = float(cols[4]) * 5
		gMaf = float(cols[18])
		tMaf = float(cols[17])
		'''
		#use the following for indels
		gMaf = float(cols[19])
		tMaf = float(cols[18])
		'''
		if gMaf >= .97 and tMaf < .97:
                        sampleStats[sample]['backMutation'][0] += 1
                        sampleStats[sample]['backMutation'][1] += tMaf
                        sampleStats[sample]['backMutation'][2] += gMaf

		#print gC, tC, gMaf, tMaf
		if tMaf > .50:
			tmaf = tMaf + tC
		else:
			tmaf = tMaf - tC
		if gMaf > .50:
			nmaf = gMaf + gC
		else:
			nmaf = gMaf - gC

		if gMaf < 0.03 and tMaf < 0.03:
                        group = '6'
		elif nmaf >= 0.97 and tmaf >= 0.97:
                        group = '1'
                elif tmaf >= 0.03 and nmaf < .01:
                        group = '2A'
                elif tmaf >= 0.03 and nmaf >= 0.01 and tmaf/nmaf >= 3:
                        group = '2B'
                elif tmaf >= .97 and nmaf < .97 and (nmaf == 0 or tmaf/nmaf < 3):
                        group = '3B'
                elif nmaf >= 0.03 and tmaf <= 0.03:
                        group = '3A'
                elif nmaf < 0.03 and tmaf < 0.03:
                        group = '5'
                else:
                        group = '4'

		if group != cols[2]:
			numDiff += 1
		groupStats[group][0] += 1
                groupStats[group][1][sample]=''

		if sample in sampleStats:
			sampleStats[sample][group][0] += 1
			sampleStats[sample][group][1] += tMaf
			sampleStats[sample][group][2] += gMaf
		newline = cols[0] + '\t' + cols[1] + '\t' + cols[2] + '\t' + group + '\t' + str(tmaf) + '\t'  + str(nmaf) + '\t'
		for col in cols[3:]:
			newline += col + '\t'
		newline = newline[:-1] + '\n'
		newfile.write(newline)
	newfile.close()
	file.close()
	print "Number of differences -", numDiff
        print "Group 1 -", groupStats['1'][0], 'in', len(groupStats['1'][1]), 'samples'
        print "Group 2A -", groupStats['2A'][0], 'in', len(groupStats['2A'][1]), 'samples'
        print "Group 2B -", groupStats['2B'][0], 'in', len(groupStats['2B'][1]), 'samples'
        print "Group 3A -", groupStats['3A'][0], 'in', len(groupStats['3A'][1]), 'samples'
        print "Group 3B -", groupStats['3B'][0], 'in', len(groupStats['3B'][1]), 'samples'
        print "Group 4 -", groupStats['4'][0], 'in', len(groupStats['4'][1]), 'samples'
        print "Group 5 -", groupStats['5'][0], 'in', len(groupStats['5'][1]), 'samples'
        print "Group 6 -", groupStats['6'][0], 'in', len(groupStats['6'][1]), 'samples'

	sampleStatsFileName = filename[:filename.rindex("/")+1] + "groupStatsBySample_newGroups.txt"
	sampleStatsFile = open(sampleStatsFileName, 'w')
	sampleStatsFile.write('Sample\tTumor C-Value\tGermline C-Value\tGroup 1\tAvg tMaf\tAvg gMaf\tGroup 2A\tAvg tMaf\tAvg gMaf\tGroup 2B\tAvg tMaf\tAvg gMaf\tGroup 3A\tAvg tMaf\tAvg gMaf\tGroup 3B\tAvg tMaf\tAvg gMaf\tGroup 4\tAvg tMaf\tAvg gMaf\tGroup 5\tAvg tMaf\tAvg gMaf\tBack Mutations\tAvg tMaf\tAvg gMaf\n')
	for sample in sampleStats:
		if sampleStats[sample]['1'][0] == 0:
			continue
		line = sample + '\t' + cvalues[sample][0] + '\t' +  cvalues[sample][1] + '\t'
		for group in ['1', '2A', '2B', '3A', '3B', '4', '5', 'backMutation']:
			line += str(sampleStats[sample][group][0]) + '\t'
			if sampleStats[sample][group][0] != 0:
				line += str(sampleStats[sample][group][1]/sampleStats[sample][group][0]) + '\t'
				line += str(sampleStats[sample][group][2]/sampleStats[sample][group][0]) + '\t'
			else:
				line += '\t\t'
		sampleStatsFile.write(line[:-1] + '\n')
	sampleStatsFile.close()

#This function removes already marked Common Sequencing Errors
def removeCSEs(filename):
	file = open(filename, 'r')
	newfile = open(filename[:-4] + '_CSEsRemoved.txt', 'w')
	cnt = 0
	for line in file:
		if line.startswith('x'):
			continue
		newfile.write(line)
		cnt += 1
	newfile.close()
	file.close()
	print cnt

#This function removes a list of samples that were pre-determined to have too many back mutations
def removeBackMutationSamples(filename):

	removeSamples = ['SJMB031', 'SJBALL021170', 'SJHYPER084', 'SJPHALL020040', 'SJERG020051', 'SJBALL020625', 'SJHYPER119', 'SJE2A006', 'SJHYPER010', 'SJBALL021058', 'SJBALL021516', 'SJBALL021893', 'SJHYPER123', 'SJHYPER095']

	file = open(filename, 'r')
	newfile = open(filename[:-4] + '_BackMutationSamplesRemoved.txt', 'w')
	cnt = 0
	for line in file:
		cols = line.split('\t')
		if cols[4] in removeSamples:
			continue
		newfile.write(line)
		cnt += 1
	newfile.close()
	file.close()
	print cnt

def finalizeGroupsAndMarkTRNAdiffs(filename):
	file = open(filename, 'r')
	newfile = open(filename[:-4] + '_FinalGroups.txt', 'w')
	cnt = 0
	for line in file:
		cols = line.rstrip().split('\t')
		group = cols[6]
		newgroup = cols[7]
		if group != newgroup and group == '4':
			if newgroup.startswith('2') or newgroup.startswith('3'):
				line = ''
				for col in cols[:7]:
					line += col + '\t'
				line += '4\t'
				for col in cols[8:]:
					line += col + '\t'
				line = line[:-1] + '\n'
		#if cols[1] == 'tRNA' and (cols[-2] == 'benign' and cols[-1] == 'pathogenic') or (cols[-2] == 'deleterious' and cols[-1] == 'benign'):
		#	line = 'x' + line
		newfile.write(line)
		cnt += 1
	newfile.close()
	file.close()
	print cnt


def main():
	#addComplexAndPredictedImpactAndMarkComSeqErrors("/home/cwelsh1/projects/PCGP/Mito/Mulligan43More/combinedFiles/all_snvs.maf.fisher_annotated_gMafOrtMaf_greaterThanOrEqual0.03_withNewGroups.txt")
        addComplexAndPredictedImpactAndMarkComSeqErrors("./all.maf.fisher_annotated_gMafOrtMaf_greaterThanOrEqual0.03_withNewGroups.txt")
    

	#UNCOMMENT THESE 6 lines for PCGP runs - left these calls in so you could see some of the code used for these functions
        #Most of this is covered in parseFilesByThreshold.py now, so shouldn't need this here, but leaving as reference
	'''
	addSubtypesColumnToSNVs("/home/cwelsh1/projects/PCGP/Mito/rerunAll/combinedFiles/all_snvs.maf.fisher_annotated_gMafOrtMaf_greaterThanOrEqual0.03_withNewGroups_removedContaminatedSamples.txt")
	addModifiedNewGroups("/home/cwelsh1/projects/PCGP/Mito/rerunAll/combinedFiles/all_snvs.maf.fisher_annotated_gMafOrtMaf_greaterThanOrEqual0.03_withNewGroups_removedContaminatedSamples_withSubtypes.txt")
	addComplexAndPredictedImpactAndMarkComSeqErrors("/home/cwelsh1/projects/PCGP/Mito/rerunAll/combinedFiles/all_snvs.maf.fisher_annotated_gMafOrtMaf_greaterThanOrEqual0.03_withNewGroups_removedContaminatedSamples_withSubtypes_withModifiedNewGroups.txt")
	removeBackMutationSamples("/home/cwelsh1/projects/PCGP/Mito/rerunAll/combinedFiles/all_snvs.maf.fisher_annotated_gMafOrtMaf_greaterThanOrEqual0.03_withNewGroups_removedContaminatedSamples_withSubtypes_withModifiedNewGroups_allColumns.txt")
	finalizeGroupsAndMarkTRNAdiffs("/home/cwelsh1/projects/PCGP/Mito/rerunAll/combinedFiles/all_snvs.maf.fisher_annotated_gMafOrtMaf_greaterThanOrEqual0.03_withNewGroups_removedContaminatedSamples_withSubtypes_withModifiedNewGroups_allColumns_BackMutationSamplesRemoved.txt")
	'''

main()
