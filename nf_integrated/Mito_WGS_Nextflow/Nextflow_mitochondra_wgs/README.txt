1) Install Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
sha256sum ./Miniconda3-latest-Linux-x86_64.sh 
dda3629462ba1cfa72eb74535214c2e315c77f1cfb0f02046537e99f1bf64abc
./Miniconda3-latest-Linux-x86_64.sh
source ~/.bashrc

2.)use the ../conda_envionment_yml/*.yml build two conda_env 
Notes.    Some python code are using python2,while some use python3 for long history reason.

conda env list
# conda environments:
#
base                   /home/wzhang/miniconda3
annot_python2          /home/wzhang/miniconda3/envs/annot_python2
nf_mitvar            * /home/wzhang/miniconda3/envs/nf_mitvar

3.) From Wenchao_Run,  copy the following files to your directory:
     nextflow pipeline main script: Mitochondra_Variant_Call.nf, 
     nextflow configure file:       nextflow.config   
     Samplesheet:                   Mit_WGS_Samplesheet.txt

4.)  In nextflow.config Configure the python2 module with your built  conda_python2
     e.g.
     process {
            withName: Annotate_Mito_Variant {
            conda = "/shared/Mito_WGS_Nextflow/conda_environment/envs/annot_python2"
            #conda = "/home/wzhang/miniconda3/envs/annot_python2"
            memory = '16GB'
            }
        }

        process {
            withName: Annotate_MTDB_Parse_FinalizeTable {
            conda = "/shared/Mito_WGS_Nextflow/conda_environment/envs/annot_python2"
            #conda = "/home/wzhang/miniconda3/envs/annot_python2"
            memory = '16GB'
            }
        }


5.)  In nextflow.config , specific your samplesheet, and projectname and output directory 
     ...
     Mitochondra_samplelist  = "./Mit_WGS_Samplesheet.txt" 
     //a txt file that recording the Mitochondra Normal and Tumor aligned BAMs. Each row has three columns that are separated by TAB and correspond to Sample_Name, Normal/Turmor, Normal_BAM/>

     project = "MT_VariantCall_HG008_Downsample_VM_Wencao" 
     //A project atlas for the current Mitochondra somatic variant calling

     outdir= "."

6.) Run nextflow
    conda activate /shared/Mito_WGS_Nextflow/conda_environment/envs/nf_mitvar
    nextflow run Mitochondra_Variant_Call.nf -profile VM

 
