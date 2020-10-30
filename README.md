# FQ Utils
```
Mode 1 (Service)
Inputs:
#1 Job data
json file for job {"reference_genome_id": "1310806.3", \
                    "output_file": "rnaseq_baumanii_1505311", \
                    "recipe": ["FASTQC","TRIM","ALIGN"], "output_path": "/anwarren@patricbrc.org/home/test",\
                    "paired_end_libs": [{"read1": "/anwarren@patricbrc.org/home/rnaseq_test/MHB_R1.fq.gz",\
                    "read2": "/anwarren@patricbrc.org/home/rnaseq_test/MHB_R2.fq.gz"},\
                    {"read1": "/anwarren@patricbrc.org/home/rnaseq_test/MERO_75_R1.fq.gz",\
                    "read2": "/anwarren@patricbrc.org/home/rnaseq_test/MERO_75_R2.fq.gz"}]}
#2 Server string to get genome from using PATRIC ID system
#3 Override parameter string to govern number of processes used
#4 Output directory

Outputs:
#1 Trimmed Fastq files (2 files for paired samples)
#2 BAM/BAI file of aligned reads
#2.5 Samstat report for BAM/BAI file
#3 FastQC report
#4 A recipe specifying which of these steps to execute

This program will:
#1 Use the given SRA accession and/or fastq files to perform the designated operations:
#1.1 fastqc report
#1.2 fastq trimming
#1.3 read mapping

Dependencies:
Bowtie2
https://github.com/BenLangmead/bowtie2

Samstat
http://samstat.sourceforge.net

Fastqc
http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Cutadapt
https://github.com/marcelm/cutadapt

Trimgalore
https://github.com/FelixKrueger/TrimGalore

Samtools
http://www.htslib.org/

HISAT2
http://daehwankimlab.github.io/hisat2/
```
