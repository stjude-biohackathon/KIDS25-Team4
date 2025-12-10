#!/usr/bin/env nextflow
/*******************************************************************************
 This Nextflow pipeline, is designed for mitochondria variant calling. User can customize to run portion or all of these analysis modules  

 Input: Common_Tumor pair sample file list. Each row has three columns that are separated by TAB and correspond to Sample_Name, Normal/Turmor, Normal_BAM/Tumor_BAM. Two rows are corresponds to one normal_tumor pair.  So, the samplesheet should be as " Sample1, normal1, norm1.bam; Sample1, tumor1, tumor1.bam; Sample2, normal2, norm2.bam; Sample2, tumor2, tumor2.bam;
 Output: Mitochondria variant calling result  
 Written by: Wenchao Zhang
 The Center for Applied Bioinformatics,St. Jude Children Research Hosptial
 Date: 10/24/2023
*********************************************************************************/

def DispConfig() {
// ContigName_Mito                         : ${params.ContigName_Mito}
 log.info """
Welocme to run Nextflow Pipeline Mitochondra_Variant_Call.nf  
Your configuration are the following:
  project                                 : ${params.project}
  Mitochondra_samplelist                  : ${params.Mitochondra_samplelist}
  Ref_MT                                  : ${params.Ref_MT}  
  Ref_Shifted_MT                          : ${params.Ref_Shifted_MT}
  Select_Variant_Call                     : ${params.Select_Variant_Call}
  Select_Variant_Combine_Filter           : ${params.Select_Variant_Combine_Filter}
  Select_Annotate_Mito_Variant            : ${params.Select_Annotate_Mito_Variant}
  Select_Annotate_MTDB_Parse_FinalizeTable: ${params.Select_Annotate_MTDB_Parse_FinalizeTable}
  Annovar_human_db_path                   : ${params.Annovar_human_db_path}
  ANNOVAR_SCript_Path                     : ${params.ANNOVAR_SCript_Path}
  Python_SCript_Path                      : ${params.Python_SCript_Path}
  Java_LIB_Path                           : ${params.Java_LIB_Path}
  outdir                                  : ${params.outdir}
 
 """  
//  exit 0    
}

def helpMessage() {
  //--ContigName_Mito                          Contig Name for Mitochondrial--ContigName_Mito                          Contig Name for Mitochondrial
  log.info """
        Welocme to run Nextflow Pipeline Mitochondra_Variant_Call.nf 
        Usage:
        A typical command for running the pipeline is as follows:
        nextflow run Mitochondra_Variant_Call.nf -profile standard 
        A command with more configurable arguments can be           
        nextflow run Mitochondra_Variant_Call.nf -profile standard -w ~/Nextflow_work --Mitochondra_samplelist .

        Configurable arguments:
        --project                                  Project name/atlas, a folder that used to distinguish other analysis 
        --Mitochondra_samplelist                   Long read WGS files list, unaligned BAM or Fastq file   
        --Ref_MT                                   Mitochondra Original reference Fasta file
        --Ref_Shifted_MT                           Mitochondra shifted reference Fasta file
        --Select_Variant_Call                      Whether select to call Mitochondra somatic variant calling 
        --Select_Variant_Combine_Filter            Whether select to combine and filter the called variants for downstream annotation  
        --Select_Annotate_Mito_Variant             Whether select to annotate Mito_Variant
        --Select_Annotate_MTDB_Parse_FinalizeTable Whether select to further annotate the MT variant against MT DataBase, Parse File by threshold, and finalize the final tab
        --Annovar_human_db_path                    The Path for ANNOVAR Database
        --ANNOVAR_SCript_Path                      The Path for ANNOVAR perl script 
        --Python_SCript_Path                       The Path for Python scripts 
        --Java_LIB_Path                            The Path for Java lib files     
        --outdir                                   The directory for the pipeline output              
        --help | h                                 This usage statement.
       Note:
         All the above arguments can be configured by the command line interface or in the nextflow.config (default)  
        """
}


// Show help message
if ((params.help) || (params.h)) {
    helpMessage()
    exit 0
}

//Parse the input Parameters and configs. 
//parse_config_parameters() 
// Display the configuration
DispConfig() 


// Set up the two external input channels for samples (fastq files) and the reference genome 
Normal_tumor_bam_ch    = Channel
            .fromPath(params.Mitochondra_samplelist)
            .splitText()
            .splitCsv(sep: '\t')

//Set up the two related geneome value channels 

/***************************************************************************************************************************************
* This Process is be responsible for extracting the MT.bam (Mitochondra mapped reads) and UM.bam(Unaligned bam) from the original normal or
* tumor aligned bam .       
***************************************************************************************************************************************/
process Extract_MT_UM_BAM {
   publishDir "${params.outdir}/${params.project}/Extract_MT_UM_BAM/${SampleName}/${Normal_Tumor}", mode: 'copy', overwrite: true
   
   input: 
   set SampleName, Normal_Tumor, bam  from Normal_tumor_bam_ch

   output:
   set SampleName, Normal_Tumor, file("${SampleName}_${Normal_Tumor}_MT_sort.bam"), file("${SampleName}_${Normal_Tumor}_UM.bam") into Extract_MT_BAM_ch
      
   script:
   """
     export TMPDIR=${params.TMP_DIR}
     ChrMCnt=`samtools view -H $bam|grep 'SN:chrM'|wc -l`
     MTCnt=`samtools view -H $bam|grep 'SN:MT'|wc -l`
     ContigName_Mito='chrM'

     if [ \${ChrMCnt} -ne 0 ]
     then
         ContigName_Mito='chrM'   
     elif [ \${MTCnt} -ne 0 ]
     then  
         ContigName_Mito='MT' 
     else
         echo -e 'Mitochondra contig was excluded in the input bam' 
         exit 1
     fi
      
     samtools view -b $bam \${ContigName_Mito} | samtools sort -@ 16 -T ${SampleName}_${Normal_Tumor}.tmp -o ${SampleName}_${Normal_Tumor}_MT_sort.bam 
     samtools index ${SampleName}_${Normal_Tumor}_MT_sort.bam   
     samtools view -b $bam -h -f 4 -o ${SampleName}_${Normal_Tumor}_UM.bam     
   """ 
}

/************************************************************************************************
* This Process is be responsible for Merginging the extracted MT.bam and UM.bam into one merged.bam. 
*************************************************************************************************/
process Merge_bam {
  publishDir "${params.outdir}/${params.project}/Merge_bam/${SampleName}/${Normal_Tumor}", mode: 'copy', overwrite: true
    input:
    set SampleName, Normal_Tumor, file(MT_bam), file(UM_bam) from Extract_MT_BAM_ch
    
    output:
    set SampleName, Normal_Tumor, file("${SampleName}_${Normal_Tumor}.merged.bam") into merge_bam_ch

    script:
    """
      export TMPDIR=${params.TMP_DIR}
      samtools merge --threads 16 -f ${SampleName}_${Normal_Tumor}.merged.bam $MT_bam $UM_bam  
    """
    
}

merge_bam_ch. into {merge_bam_ch1; merge_bam_ch2}
/***************************************************************************************************************************************
* This Process is be responsible for 1.) align the merged.bam against shiftedMT.fa, 2.)samtools sort and 3.) samtools index
***************************************************************************************************************************************/
process bwa_aln_shift {
    publishDir "${params.outdir}/${params.project}/bwa_aln_shift/$SampleName", mode: 'copy', overwrite: true

    memory { (16.GB + (8.GB * task.attempt)) } // First attempt 16GB, second 24GB, etc
    errorStrategy 'retry'
    maxRetries 3

    input:
    set SampleName, Normal_Tumor, file(merged_bam) from merge_bam_ch1

    output:
    set SampleName, Normal_Tumor, file("${SampleName}_${Normal_Tumor}_shift_merged_aln_sort.bam"), file("${SampleName}_${Normal_Tumor}_shift_merged_aln_sort.bam.bai") into bwa_aln_shift_ch
       
    script:
    """
      export TMPDIR=${params.TMP_DIR}
      bwa aln  ${params.Ref_Shifted_MT} -t 16 -f ${Normal_Tumor}_shift_merged.sai -b ${merged_bam}
      bwa samse ${params.Ref_Shifted_MT} ${Normal_Tumor}_shift_merged.sai ${merged_bam} | samtools sort -@ 16 -T {SampleName}_${Normal_Tumor}_shift_merged_aln.sort.tmp -o ${SampleName}_${Normal_Tumor}_shift_merged_aln_sort.bam 
      samtools index ${SampleName}_${Normal_Tumor}_shift_merged_aln_sort.bam         
    """
}

/***************************************************************************************************************************************
* This Process is be responsible pairing the shift aligned bams into ordered pairs as normal tumor , which can be sent to call variants
***************************************************************************************************************************************/
process Pair_Normal_Tumor_bwa_aln_shift {
    input:
    set SampleName, Normal_Tumors, file(aln_sort_bams), file(aln_sort_bam_bais) from bwa_aln_shift_ch.groupTuple()
    
    output:
    set SampleName, env(shift_org), file("${SampleName}_shift_normal_aln_sort.bam"), file("${SampleName}_shift_normal_aln_sort.bam.bai"),  file("${SampleName}_shift_tumor_aln_sort.bam"), file("${SampleName}_shift_tumor_aln_sort.bam.bai") into Pair_Normal_Tumor_bwa_aln_shift_ch

    script:
    if ( Normal_Tumors[0]== 'Tumor')
    """
     shift_org='shift'
     ln -s ${aln_sort_bams[1]} ${SampleName}_shift_normal_aln_sort.bam
     ln -s ${aln_sort_bam_bais[1]} ${SampleName}_shift_normal_aln_sort.bam.bai
     ln -s ${aln_sort_bams[0]} ${SampleName}_shift_tumor_aln_sort.bam
     ln -s ${aln_sort_bam_bais[0]} ${SampleName}_shift_tumor_aln_sort.bam.bai
    """
    else
    """
     shift_org='shift'
     ln -s ${aln_sort_bams[0]} ${SampleName}_shift_normal_aln_sort.bam
     ln -s ${aln_sort_bam_bais[0]} ${SampleName}_shift_normal_aln_sort.bam.bai
     ln -s ${aln_sort_bams[1]} ${SampleName}_shift_tumor_aln_sort.bam
     ln -s ${aln_sort_bam_bais[1]} ${SampleName}_shift_tumor_aln_sort.bam.bai
    """
}


/***************************************************************************************************************************************
* This Process is be responsible for 1.) align the merged.bam against MT.fa and 2.)samtools sort and 3.) samtools index
***************************************************************************************************************************************/
process bwa_aln_org {
	publishDir "${params.outdir}/${params.project}/bwa_aln_org/$SampleName", mode: 'copy', overwrite: true
    
    memory { (16.GB + (8.GB * task.attempt)) } // First attempt 16GB, second 24GB, etc
    errorStrategy 'retry'
    maxRetries 3
    
    input:
    set SampleName, Normal_Tumor, file(merged_bam) from merge_bam_ch2

    output:
    set SampleName, Normal_Tumor, file("${SampleName}_${Normal_Tumor}_org_merged_aln_sort.bam"), file("${SampleName}_${Normal_Tumor}_org_merged_aln_sort.bam.bai") into bwa_aln_org_ch
       
    script:
    """
      export TMPDIR=${params.TMP_DIR}
	  bwa aln ${params.Ref_MT} -t 16 -f ${Normal_Tumor}_merged_org.sai -b ${merged_bam}
	  bwa samse ${params.Ref_MT} ${Normal_Tumor}_merged_org.sai ${merged_bam} | samtools sort -@ 16 -T ${SampleName}_${Normal_Tumor}_org_merged_aln_sort.tmp -o ${SampleName}_${Normal_Tumor}_org_merged_aln_sort.bam
      samtools index ${SampleName}_${Normal_Tumor}_org_merged_aln_sort.bam          
    """
}


/***************************************************************************************************************************************
* This Process is be responsible pairing the original aligned bams into ordered pairs as normal tumor , which can be sent to call variants
***************************************************************************************************************************************/
process Pair_Normal_Tumor_bwa_aln_org {
    input:
    set SampleName, Normal_Tumors, file(aln_sort_bams), file(aln_sort_bam_bais) from bwa_aln_org_ch.groupTuple()
    
    output:
    set SampleName, env(Shift_Org),  file("${SampleName}_org_normal_aln_sort.bam"), file("${SampleName}_org_normal_aln_sort.bam.bai"),  file("${SampleName}_org_tumor_aln_sort.bam"), file("${SampleName}_org_tumor_aln_sort.bam.bai") into Pair_Normal_Tumor_bwa_aln_org_ch

    script:
    if ( Normal_Tumors[0]== 'Tumor')
    """
     Shift_Org='org'
     ln -s ${aln_sort_bams[1]} ${SampleName}_org_normal_aln_sort.bam
     ln -s ${aln_sort_bam_bais[1]} ${SampleName}_org_normal_aln_sort.bam.bai
     ln -s ${aln_sort_bams[0]} ${SampleName}_org_tumor_aln_sort.bam
     ln -s ${aln_sort_bam_bais[0]} ${SampleName}_org_tumor_aln_sort.bam.bai
    """
    else
    """
     Shift_Org='org'
     ln -s ${aln_sort_bams[0]} ${SampleName}_org_normal_aln_sort.bam
     ln -s ${aln_sort_bam_bais[0]} ${SampleName}_org_normal_aln_sort.bam.bai
     ln -s ${aln_sort_bams[1]} ${SampleName}_org_tumor_aln_sort.bam
     ln -s ${aln_sort_bam_bais[1]} ${SampleName}_org_tumor_aln_sort.bam.bai
    """
}

/*This process is responsible for calling the Mitochondra variants for one tumor and normal bam pair*/
Ace2_Variant_Call_In_Ch= (params.Select_Variant_Call == "N") ? Channel.empty() : Pair_Normal_Tumor_bwa_aln_org_ch.mix(Pair_Normal_Tumor_bwa_aln_shift_ch)
process Ace2_Variant_Call {
  publishDir "${params.outdir}/${params.project}/Ace2_Variant_Call/$SampleName", mode: 'copy', overwrite: true
  
  memory { (16.GB + (8.GB * task.attempt)) } // First attempt 16GB, second 24GB, etc
  errorStrategy 'retry'
  maxRetries 3
  
  input:
  set SampleName, Shift_Org, file(Normal_bam), file(Normal_bam_bai), file(Tumor_bam), file(Tumor_bam_bai) from Ace2_Variant_Call_In_Ch
  //set SampleName, file(Normal_bam), file(Normal_bam_bai), file(Tumor_bam), file(Tumor_bam_bai) from Ace2_Variant_Call_In_Ch
  
  output:
  set SampleName, Shift_Org, file("${SampleName}_${Shift_Org}_output_high.tab") into Ace2_Variant_Call_Ch

  // if ( Shift_Org=='org')  
  // -fasta ${params.Ref_MT} \
script:
   def Ref_fa=(Shift_Org=='org')? params.Ref_MT : params.Ref_Shifted_MT
   """
   export CLASSPATH="${params.Java_LIB_Path}/sam-1.65.jar:${params.Java_LIB_Path}/bambino-1.0.jar:${params.Java_LIB_Path}/third_party.jar:${params.Java_LIB_Path}/picard.jar:${params.Java_LIB_Path}/:mysql-connector-java-5.1.10.jar "

   java -Xmx8000m \
   Ace2.SAMStreamingSNPFinder \
   -bam  ${Normal_bam} \
   -tn N \
   -bam ${Tumor_bam} \
   -tn T \
   -of ${SampleName}_${Shift_Org}_output_high.tab \
   -no-dbsnp \
   -min-quality 20 \
   -min-flanking-quality 20 \
   -min-alt-allele-count 1 \
   -min-minor-frequency 0 \
   -broad-min-quality 10 \
   -mmf-max-hq-mismatches 4 \
   -mmf-min-quality 15 \
   -mmf-max-any-mismatches 6 \
   -unique-filter-coverage 2 \
   -min-mapq 1 \
   -optional-tags XT!=R \
   -fasta ${Ref_fa} \
   -fisher-strand-enable \
   -ignore-first-base-mismatches \
   -ignore-last-base-mismatches
   """
   
}

/*This process is responsible for Combine the Ace2 called Varaints against original reference(MT.fa) and shifted reference (MT_shifted.fa), and then filter by the Strand Bias  criteria */
Combine_Filter_In_Ch= (params.Select_Variant_Combine_Filter == "N") ? Channel.empty() : Ace2_Variant_Call_Ch.groupTuple()
process Combine_Filter_Mito_Variant {
  publishDir "${params.outdir}/${params.project}/Combine_Filter_Mito_Variant/$SampleName", mode: 'copy', overwrite: true  
  input:
  set SampleName, Shift_Orgs, file(output_high_tabs) from Combine_Filter_In_Ch
  
  output:
  set SampleName, file("${SampleName}_Combine_filter.tab") into Combine_Filter_Ch
  
  script:
  if ( Shift_Orgs[0]== 'org')  // [0]: org; [1] : shift
  """
    python ${params.Python_SCript_Path}/Combine_FilterbyNumberReadStrandBias.py \
    -org ${output_high_tabs[0]} \
    -shift ${output_high_tabs[1]} \
    -o ${SampleName}_Combine_filter.tab 
  """
  else  //[0]: shift;  [1] : org
  """
    python ${params.Python_SCript_Path}/Combine_FilterbyNumberReadStrandBias.py \
    -org ${output_high_tabs[1]} \
    -shift ${output_high_tabs[0]} \
    -o ${SampleName}_Combine_filter.tab  
  """
  
}

/*This process is responsible for Annotate the detected and filtered Mito Variant  */
Annotate_Mito_Variant_In_Ch = (params.Select_Annotate_Mito_Variant == "N" )? Channel.empty() : Combine_Filter_Ch
process Annotate_Mito_Variant {
  publishDir "${params.outdir}/${params.project}/Annotate_Mito_Variant/$SampleName", mode: 'copy', overwrite: true
  input:
  set SampleName, file(Combine_filter_tab) from Annotate_Mito_Variant_In_Ch

  output:
  set SampleName, file("${SampleName}_chrM.maf.fisher") into Annotate_Mito_Variant_Ch
  
  //conda '${baseDir}/annot_python2.yml' // Path to your YAML file
  // conda 'python=2.7'
   
  script:
  """
   python ${params.Python_SCript_Path}/bambino2annovarInput_MT.py ${Combine_filter_tab}
   ${params.ANNOVAR_SCript_Path}/annotate_variation.pl -geneanno ${Combine_filter_tab}.annovar_input ${params.Annovar_human_db_path} -buildver GRCh37_MT -dbtype ensGene -neargene 1000
   python ${params.Python_SCript_Path}/calculateMAFfisherPvalue.CAB.py ${Combine_filter_tab}.annovar_input.variant_function
  """ 
}

/*This process is responsible for Annotate MT_Database, ParseByThreshold, and Finalize Table of the upstream annotated Mito Variants*/
Annotate_MTDB_Parse_FinalizeTable_In_Ch  = (params.Select_Annotate_MTDB_Parse_FinalizeTable == "N" )? Channel.empty() :  Annotate_Mito_Variant_Ch
process Annotate_MTDB_Parse_FinalizeTable {  
  publishDir "${params.outdir}/${params.project}/Annotate_MTDB_Parse_FinalizeTable", mode: 'copy', overwrite: true
  input:
  set SampleNames, file(chrM_maf_fisher_Files) from Annotate_MTDB_Parse_FinalizeTable_In_Ch.groupTuple(by:3)

  output:
  file("all.maf.*") into AnnotateMTDB_Parse_FinalizeTable_Ch
   
  //conda '${baseDir}/annot_python2.yml' // Path to your YAML file
  //conda 'python=2.7'

  script:
  """
  cat ${chrM_maf_fisher_Files} >> all.maf.fisher
  python ${params.Python_SCript_Path}/annotateFiles.CAB.py 
  python ${params.Python_SCript_Path}/parseFilesByThreshold.CAB.py 
  python ${params.Python_SCript_Path}/finalizeTable.CAB.py
  """  
}

