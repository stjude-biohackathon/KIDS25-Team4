import sys
from scipy import stats

def calcMAF_fisher(fileName):
	if True:
		file = open(fileName, 'r')
		#if '_chrM' not in fileName:
		#	sampleName = 'sample'
		#else:
		#    sampleName = fileName[:fileName.index("_chrM")]
                sampleName=fileName[:fileName.index("_Combine_filter.tab.annovar_input.variant_function")]
		#outputFileName = "/home/cwelsh1/projects/PCGP/Mito/updatedAnnovar/" + sampleName + "_chrM.maf.fisher"
		#outputFileName = "/home/cwelsh1/projects/PCGP/Mito/SJTALL/combinedFiles/" + sampleName + "_chrM.maf.fisher"
		#outputFileName = "/home/cwelsh1/projects/PCGP/Mito/relapseSamples/combinedFiles/" + sampleName + "_chrM.maf.fisher"
		#outputFileName = "/home/cwelsh1/projects/PCGP/Mito/SCMC_RELAPSE/relapse_tumor/combinedFiles/" + sampleName + "_chrM.maf.fisher"
		#outputFileName = "/home/cwelsh1/projects/PCGP/Mito/RELAPSE_ALL_MULLIGAN/tumor_germline/combinedFiles/" + sampleName + "_chrM.maf.fisher"
		#outputFileName = "/home/cwelsh1/projects/PCGP/Mito/ClinicalGenomics_RTCG/combinedFiles/" + sampleName + "_chrM.maf.fisher"
		#outputFileName = "/home/cwelsh1/projects/PCGP/Mito/ClinicalGenomics_RTCG/relapse_germline/combinedFiles/" + sampleName + "_chrM.maf.fisher"
		outputFileName = sampleName + "_chrM.maf.fisher"
		#outputFileName = "/home/cwelsh1/projects/PCGP/Mito/COLO829/relapse_tumor/combinedFiles/" + sampleName + "_chrM.maf.fisher"
		##outputFileName = "/home/cwelsh1/projects/PCGP/Mito/COLO829/combinedFiles/" + sampleName + "_chrM.maf.fisher"
		#outputFileName = "/home/cwelsh1/projects/PCGP/Mito/single-cell/SJETV027/" + fileName + "_chrM.maf.fisher"
		#outputFileName = "/home/cwelsh1/projects/PCGP/Mito/single-cell/SJINF010/" + fileName + "_chrM.maf.fisher"
		#outputFileName = "/home/cwelsh1/projects/PCGP/Mito/rerunAll/combinedFiles/" + sampleName + "_chrM.maf.fisher"
		print fileName, outputFileName
		outputFile = open(outputFileName, 'w')
		for line in file:
			line = line.rstrip()
			cols = line.split('\t')
			newline = sampleName + '\t'
			if int(cols[3]) < 575 or int(cols[3]) >= 16022:
				newline += 'D-loop'
			elif cols[0].find('exonic') != -1:
				if cols[1].find("tRNA") != -1:
					newline += 'tRNA'
				elif cols[1].find("RNR") != -1:
					newline += 'rRNA'
			newline+= '\t'
			tumorMAF = 0
			normalMAF = 0
			for ind in xrange(0, len(cols)):
				col = cols[ind]
				if ind == 3 or ind == 4:
					#col = int(col) - 1	#no longer needed
					newline += str(col) + '\t'
				elif ind == 8 or ind == 9:
					mut_reads, all_reads = col.split('/')
					newline += mut_reads + '\t' + all_reads + '\t'
					if ind == 8:
						if int(all_reads) != 0:
							tumorMAF = float(mut_reads) / int(all_reads)
						else:
							tumorMAF = 0
						a = float(mut_reads)
						b = float(all_reads) - a
					else:
						if int(all_reads) != 0:
							normalMAF = float(mut_reads) / int(all_reads)
						else:
							normalMAF = 0
						c = float(mut_reads)
						d = float(all_reads) - c
				else:
					newline += col + '\t'
			oddsratio, pvalue = stats.fisher_exact([[a,b],[c,d]])
			newline += str(tumorMAF) + '\t' + str(normalMAF) + '\t' + str(pvalue) + '\n'
			outputFile.write(newline)
		file.close()
		outputFile.close()

if __name__ == "__main__":

	if len(sys.argv) > 1:
		fileName = sys.argv[1]
		calcMAF_fisher(fileName)
	else:
		print "Usage: calculateMAFfisherPvalue.py <variant_function_filename>"
