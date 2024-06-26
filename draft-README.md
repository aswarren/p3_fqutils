# Fastq Utilities Service

## Overview

The Fastq Utilities Service makes available common operations for FASTQ files from high throughput sequencing, including: generating FastQC reports of base call quality; aligning reads to genomes using Bowtie2 to generate BAM files, saving unmapped reads and generating SamStat reports of the amount and quality of alignments; and trimming of adapters and low quality sequences using TrimGalore and CutAdapt. The Fastq Utiliites app allows the user to define a pipeline of activities to be performed to designated FASTQ files. The three components (trim, fastqc and align) can be used independently, or in any combination.These actions happen in the order in which they are specified. In the case of trimming, the action will replace untrimmed read files with trimmed ones as the target for all subsequent actions. FASTQ reads (paired-or single-end, long or short, zipped or not), as well as Sequence Read Archive accession numbers are supported.  



## About this module

This module is a component of the BV-BRC build system. It is designed to fit into the
`dev_container` infrastructure which manages development and production deployment of
the components of the BV-BRC. More documentation is available [here](https://github.com/BV-BRC/dev_container/tree/master/README.md).

This module provides the following application specfication(s):
* [FastqUtils](app_specs/FastqUtils.md)


## See also

* [Fastq Utilities Service Quick Reference](https://www.bv-brc.org/docs/quick_references/services/fastq_utilities_service.html)
* [Fastq Utilities Service](https://www.bv-brc.org/docs/https://bv-brc.org/app/FastqUtil.html)
* [Fastq Utilities Service Tutorial](https://www.bv-brc.org/docs//tutorial/fastq_utilities/fastq_utilities.html)



## References


1. Krueger, F., Trim Galore: a wrapper tool around Cutadapt and FastQC to consistently apply quality and adapter trimming to FastQ files, with some extra functionality for MspI-digested RRBS-type (Reduced Representation Bisufite-Seq) libraries. URL http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/. (Date of access: 28/04/2016), 2012.
2. Martin, M., Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet. journal, 2011. 17(1): p. 10-12.
3. Andrews, S., FastQC: a quality control tool for high throughput sequence data. 2010.
4. Langmead, B. and S.L. Salzberg, Fast gapped-read alignment with Bowtie 2. Nature methods, 2012. 9(4): p. 357.
5. Langmead, B., et al., Scaling read aligners to hundreds of threads on general-purpose processors. Bioinformatics, 2018. 35(3): p. 421-432.
6. Lassmann, T., Y. Hayashizaki, and C.O. Daub, SAMStat: monitoring biases in next generation sequencing data. Bioinformatics, 2010. 27(1): p. 130-131.
7. Bankevich, A., et al., SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. Journal of computational biology, 2012. 19(5): p. 455-477.
