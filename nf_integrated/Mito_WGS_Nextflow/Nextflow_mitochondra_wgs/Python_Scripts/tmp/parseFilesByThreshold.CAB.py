import glob


def parseByMafs(filename):
	file = open(filename, 'r')
	outputFileName = filename[:-4] + "_gMafOrtMaf_greaterThanOrEqual0.03.txt"
	outputFile = open(outputFileName, 'w')
	cnt = 0
	for line in file:
		cols = line.split('\t')
		pos = int(cols[6])
		#nmaf = float(cols[15])
		#tmaf = float(cols[14])
		#pvalue = float(cols[16])
		nmaf = float(cols[14])
		tmaf = float(cols[13])
		pvalue = float(cols[15])
		'''
		#use these for indels
		nmaf = float(cols[15])
		tmaf = float(cols[14])
		pvalue = float(cols[16])
		'''
		#remove known common sequencing error 'SNVs'
                #if (pos >= 302 and pos <= 315) or (pos >= 513 and pos <= 525) or (pos >= 3105 and pos <= 3109) or pos == 567:
                #        continue
		#if nmaf < 0.001 and tmaf < 0.001:
			#continue
		if nmaf >= 0.03 or tmaf >= 0.03:
		#if nmaf >= 0.01 or tmaf >= 0.01:
			outputFile.write(line)
			cnt += 1
		'''
		if nmaf >= 0.01 or tmaf >= 0.01:
                	if (pos >= 302 and pos <= 315) or (pos >= 513 and pos <= 525) or (pos >= 3105 and pos <= 3109) or pos == 567:
				outputFile.write("CommonSeqError\t" + line)
			else:
				outputFile.write("\t" + line)
			cnt += 1
		'''
	print cnt
	outputFile.close()
	file.close()

def assignNewGroups(filename):
	file = open(filename, 'r')
 	outputFileName = filename[:-4] + "_withNewGroups.txt"
        outputFile = open(outputFileName, 'w')
	group1cnt = 0
	group2Acnt = 0
	group2Bcnt = 0
	group3Acnt = 0
	group3Bcnt = 0
	group4cnt = 0
	group6cnt = 0

        ### column headers:
        headers=['Sample', 'Group', 'Syn vs NonSyn','Gene Annotation','Location Information.mtANN','Location Information.exonic.func','Location Information.gene','mtDNA_pos','Ref. Allele','Mutant Allele', 'Mut_Allele_Read_Tumor','Total_Read_Tumor', 'Mut_Allele_Read_Normal','Total_Read_Normal','MAF_Tumor','MAF_Normal','P-value','Reported_disease_association', 'Disease_status','is_in_mtDB','MAF_in_mtDB','Reported_by_Levin_et_al(GBE,2013,doi"10.1093/gbe/evt058)', 'MitoMap Pop Freq','AA Variant','Gene/tRNA old','Func.Impact','MutationAssessorScore']

        outputFile.write('\t'.join(headers)+'\n')

        for line in file:
                cols = line.split('\t')
                nmaf = float(cols[14])
                tmaf = float(cols[13])
		'''
		#use these for indels
		nmaf = float(cols[15])
		tmaf = float(cols[14])
		'''
		if nmaf < 0.03 and tmaf < 0.03:
			group = '6'
			group6cnt += 1
                elif nmaf >= 0.97 and tmaf >= 0.97:
			group = '1'
			group1cnt += 1
		elif tmaf >= 0.03 and nmaf < .01:
			group = '2A'
			group2Acnt += 1
		elif tmaf >= 0.03 and nmaf >= 0.01 and tmaf/nmaf >= 3:
			group = '2B'
			group2Bcnt += 1
		elif tmaf >= .97 and nmaf < .97 and (nmaf == 0 or tmaf/nmaf < 3):
			group = '3B'
			group3Bcnt += 1
		elif nmaf >= 0.03 and tmaf <= 0.03:
			group = '3A'
			group3Acnt += 1
		else:
			group = '4'
			group4cnt += 1
		newline = cols[0] + '\t' + group + '\t'
		for col in cols[1:]:
			newline += col + '\t'
		newline = newline[:-1] 
		outputFile.write(newline)

	print "Group 1 -", group1cnt
	print "Group 2A -", group2Acnt
	print "Group 2B -", group2Bcnt
	print "Group 3A -", group3Acnt
	print "Group 3B -", group3Bcnt
	print "Group 4 -", group4cnt
	print "Group 6 -", group6cnt
        outputFile.close()


def main():

	#filename = "/home/cwelsh1/projects/PCGP/Mito/Mulligan43More/combinedFiles/all_snvs.maf.fisher_annotated.txt"
	filename ="all.maf.fisher_annotated.txt"
        parseByMafs(filename)
	assignNewGroups(filename[:-4] + "_gMafOrtMaf_greaterThanOrEqual0.03.txt")

        #once you have a list of C-values for the samples, then add a function to remove samples with C-values > 0.02
        #also check for samples with back mutations - remove any samples with

main()

