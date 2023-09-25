# perturbsci-pipeline
A processing pipeline that will take fastq reads produced from a PerturbSci library preparation with a NextSeq 550 75 cycle kit and output a counts matrix and a seurat object.

## Create the conda environment
Create a conda environment with the required packages from the provided YAML file.
```{bash}
conda env create --name EasySci --file=EasySci_conda_env.yaml
```

## Input parameters
### Sample ID files
Create two files, one for RNA samples and one for guide samples. Thesse are .txt files where each line is the sample name for the corresponding PCR barcode. Subsequently these are the sample names that are added to the fastq files after demultiplexing. An example for 4 PCR barcodes would be:
```
rna1
rna2
rna3
rna4
```

### Genome Reference
This pipeline uses STAR to align your reads. Build or download a reference such as the one from [10X Genomics](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest?). If you don't already have STAR indices built from the refernce, you will have to do that with:
```{bash}
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir path/to/star/output/reference --genomeFastaFiles path/to/fasta/refernce.tar.gz
```

### sgRNA Correction File and Info Table
You will need to provide a .pickle2 file of the guide sequences used in your experiment. This can be created with the `gRNA_whitelist_generation.py` script provided. To run this script you will need to provide a file with sgRNA information. This is a tab delimited .txt file with three columns. The file has a header with the columns as `gRNA_seq`, `names`, and `gene_symbol`. The first column is the `gRNA_seq` or sequence for your sgRNAs. Our sequencing structure with a 75 cycle kit only capture the first 19 out of 20 bp which is reflected in this column. The second column is the `names` or the name of your sgRNA (i.e. CD55g1). The third column is the `gene_symbol` or the gene symbol of the sgRNA target. Here is an example:
```
gRNA_seq        names   gene_symbol
GCGCGGCCCAGAGGCGGCA     CPSF6g1 CPSF6
GACCTGGGCTCTCACGTAC     CPSF6g2 CPSF6
GGCGGGCAGTAGCCGCTGA     NUDT21g1        NUDT21
GGCTACTGCCCGCCATTAA     NUDT21g2        NUDT21
```
Then change what your sgRNA info table is called in the `gRNA_pilot_screen_file` variable and what you want the .pickle2 file to be called with the `pickle_out` variable in the `gRNA_whitelist_generation.py` script before running it.

Both the info table .txt file and .pickle2 file are inputs for the main script file.


## Run the Pipeline
Run the pipeline as:
```{bash}
bash main_pipeline.sh > output.log 2> messages.log
```
Errors and other messages will be saved to messages.log and normal output will be written to outpt.log.

You will get fastq files and sam files ouput at each processing step in the `gex_processing` directory. The guide counts will be output to the `gdo_processing` directory. You will get a final seurat object with the guide and gene expression counts in the `gex_processing` directory. Proceed by checking the QC metrics (nCount_RNA, nFeature_RNA, nCount_GDO) and follow Seurat clustering workkflows to analyze your data, including the `MULTIseqDemux` function to assign guides. 
