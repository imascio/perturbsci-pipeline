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
