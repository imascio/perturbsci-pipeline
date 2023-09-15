#!/bin/bash

#run as bash main_pipeline.sh > output.log 2> messages.log


## Location of the bashrc file to activate the conda environments.
#bashrc_location="/home/dalgarnoc/.bashrc"
#activate_location="/home/dalgarnoc/miniconda3/bin/activate"
#source $activate_location
conda activate EasySci

# Define script directory

script_path="/brahms/shared/EasySci/shell_pipeline/script_folder"

###################### DEFINE INPUT VARIABLES ######################

# define the project folder, this is where all output files will be written
project_folder="/brahms/shared/EasySci/test_rna1"

# define the fastq folder including all fastq files
fastq_folder="/brahms/shared/EasySci/fastq"
# change to the following line after troubleshooting phase is done
#fastq_folder=$project_folder/fastq

# define the PCR group sample id for each fastq file
sample_ID="/brahms/shared/EasySci/fastq/gex_sampleID.txt"

# define the output folder where all processing files will be written, inside the project folder
gex_processing_folder=$project_folder/gex_processing

# define the core number for parallele processing
core=16

# define the core number for reads filtering and sorting
samtools_core=16

# define the number of jobs for parallelizing
N_JOBS=16

# define the location of index files for reads alignment with STAR
# Human genome version 42
index="/brahms/shared/EasySci/human/STAR_original_pipeline"

# define the gtf file for gene counting
# Human genome version 42
gtf_file="/brahms/shared/EasySci/human/gencode.v42.primary_assembly.annotation.gtf"

# define the gtf file for exon counting
# Human genome version 42
gtf_file_exon="/brahms/shared/EasySci/human/gencode.v42.primary_assembly.annotation.exons.gtf"

# define the location of the ligation barcodes
ligation_barcode=$script_path/Ligation_barcodes_NextSeq.pickle2
# define the location of the RT barcodes
RT_barcode=$script_path/RT_barcodes.pickle2
# define the location of the randomN RT barcodes for the gene counting step
randomN_barcode_file=$script_path/RandomN_RT_barcodes.txt

# define the name of the final seuart object - should be .rds extension
seurat_object_name=ran1_seurat_obj.rds


############ BARCODE EXTRACTION ############

# the script take an input folder, a sample ID list, an output folder, the RT barcode list, the ligation barcode list and core number and a rendom hexamer barcode list. Then it extract the RT barcode from read1, the ligation barocde from read2, correct them to the nearest RT and ligation barcode (with edit distance <= 1), and attach the RT and ligation barcode and UMI sequence to the read name of read1 and read3. Reads with unmatched RT or ligation barcodes are discarded.

input_folder=$fastq_folder
output_folder=$gex_processing_folder/BC_attach
script=$script_path/barcode_extraction.py
echo "Changing the name of the fastq files..." >&2

echo $sample_ID >&2
for sample in $(cat $sample_ID); do echo changing name $sample; mv $input_folder/*$sample*R1*.fastq.gz $input_folder/$sample.R1.fastq.gz; mv $input_folder/*$sample*R2*.fastq.gz $input_folder/$sample.R2.fastq.gz; mv $input_folder/*$sample*R3*.fastq.gz $input_folder/$sample.R3.fastq.gz; done

echo "Attaching barcode and UMI...." >&2
mkdir -p $output_folder
parallel -j ${N_JOBS} --verbose python $script $input_folder {} $output_folder $ligation_barcode $RT_barcode $core $randomN_barcode_file :::: ${sample_ID}
echo "Barcode transformed and UMI attached." >&2

################# TRIM READ 2 #################

# This steps trims any poly A stretches present in read 2

echo "Trimming read 2..." >&2

input_folder=$gex_processing_folder/BC_attach
sampleID=$(head -n 1 $sample_ID)
output_folder=$gex_processing_folder/trimmed_fastq
mkdir -p $output_folder

parallel -j ${N_JOBS} --verbose bash $script_path/trimming.sh $input_folder {} $output_folder :::: ${sample_ID}

################# REMOVE EMPTY LINES #################

# This step removes any empty fastq lines resulting from trimming reads that are entirely poly A

echo "\nRemoving empty lines from fastqs" >&2

sampleID=$(head -n 1 $sample_ID)
input_folder=$gex_processing_folder/trimmed_fastq
output_folder=$gex_processing_folder/cleaned_fastq
mkdir -p $output_folder

parallel -j ${N_JOBS} --verbose bash ${script_path}/remove_empty_lines.sh {} $input_folder ${output_folder} :::: ${sample_ID}

echo "Empty lines removed" >&2

############ ALIGN ############


echo "Aligning cleaned fastq files..." >&2

input_folder=$gex_processing_folder/cleaned_fastq
output_folder=$gex_processing_folder/STAR_alignment
mkdir -p $output_folder

parallel -j ${N_JOBS} --verbose bash $script_path/align.sh $index {} $input_folder $output_folder :::: ${sample_ID}

echo "Alignment complete." >&2

############ FEATURE ############

# This step maps aligned reads to genes using the umi tools feature function

echo "Mapping reads to genes..." >&2

input_folder=$gex_processing_folder/STAR_alignment
output_folder=$gex_processing_folder/feature
mkdir -p $output_folder

parallel -j ${N_JOBS} --verbose bash $script_path/feature.sh $gtf_file {} $input_folder $output_folder :::: ${sample_ID}

echo "Mapping to genes complete." >&2

############ COUNT ############

# Remove UMI duplicates and generate a cell x gene counts matrix

echo "Generating a counts matrix! Almost done..." >&2

input_folder=$gex_processing_folder/feature
output_folder=$gex_processing_folder/count
mkdir -p $output_folder

parallel -j ${N_JOBS} --verbose bash $script_path/count.sh $input_folder $output_folder {} :::: ${sample_ID}

############ SEURAT ############

# Make a seurat object from the counts matrices

echo "Generating a seurat object! Last step..." >&2

input_folder=$gex_processing_folder/count
output_folder=$gex_processing_folder/seurat
mkdir -p $output_folder

Rscript $script_path/make_seurat_object.R $sampleID $input_folder $output_folder $seurat_object_name


echo "You're done ;)" >&2

