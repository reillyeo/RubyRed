# RubyRed

Custom bioinformatics pipeline designed for high-throughput, read-by-read taxonomic classification of rRNA gene amplicons generated via Oxford Nanopore sequencing. RubyRed uses open-source tools and custom scripts to process raw fastq files and create taxonomically annotated feature tables.

Input data can be either raw or demultiplexed FASTQ files. Demultiplexing (if necessary) is done using Guppy Barcoder, trimming of primer binding regions with [Cutadapt](https://github.com/marcelm/cutadapt), quality and length filtering with [Chopper](https://github.com/wdecoster/chopper), and subsequent conversion to FASTA format with [VSEARCH](https://github.com/torognes/vsearch). To ensure data consistency, sequences with fewer than a minimum read count are discarded, and samples exceeding a certain user-definable read number threshold (default: 30,000 reads) are subsampled using [SeqKit](https://github.com/shenwei356/seqkit).

Next, all filtered reads are concatenated and imported into [QIIME2](https://docs.qiime2.org) as a single sequence artifact. A custom Python script is then employed to generate a feature table. Chimera removal is performed using VSEARCH’s uchime-ref algorithm against a reference database, and surviving sequences are reoriented with [RESCRIPt](https://github.com/bokulich-lab/RESCRIPt) to match reference strand orientation.

For taxonomic assignment, the pipeline supports three classification methods: sklearn, consensus-vsearch, and consensus-blast. The sklearn method requires a pre-trained classifier, but is considerably faster than either of the other methods. The final output includes a taxonomically annotated feature table (.biom and .tsv format) and an individual fasta file for each unique taxonomic classification, containing all sequences assigned to that taxon.


## Installation

RubyRed is developed for Linux/Unix. A conda installation is required prior to installation.
For out-of-the-box functionality, the following lines of code should be run:

    mkdir ~/my_scripts
    cd ~/my_scripts
    git clone https://github.com/reillyeo/RubyRed
    conda env create -n qiime2 --file  https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2024.5/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml
    conda activate qiime2
    conda install -c bioconda -c conda-forge chopper seqkit parallel
    echo 'export PATH="~/my_scripts/RubyRed/:$PATH"' >> ~/.bashrc
    echo 'export PATH="~/my_scripts/RubyRed/resources/ont-guppy/bin/:$PATH"' >> ~/.bashrc

## Usage

```text
Usage - RubyRed [OPTIONS]
                 
Options:                 
         -i      directory containing input FASTQ files (default: current directory)                
         -d      use this flag if data has already been demultiplexed                
         -q      minimum average read quality required to pass chopper quality filtering (default: 15)                
         -l      minimum read length allowed to pass chopper length filtering (default: 900)                
         -x      maximum read length allowed to pass chopper length filtering (default: 1200)                
         -p      number of threads to use for parallel processing (default: 20)                
         -m      minimum number of reads (post-filtering) to keep a file (default: 1)                
         -t      subsample fasta files with more than this number of reads (default: 30000)                
         -s      directory where python scripts are located (default: $HOME/my_scripts/RubyRed/scripts)                
         -w      path to the directory where resources (primer seqs, reference seqs, reference taxonomy, classifier) are located (default: $HOME/my_scripts/RubyRed/resources)                
         -u      path to fasta file containing primer sequences (default: $HOME/my_scripts/RubyRed/resources/UMI16s_primers.fasta)                
         -r      path to reference sequences for chimera filtering and reorientation (default: $HOME/my_scripts/RubyRed/resources/classifiers/MIMt_16s/MIMt_16s_refseqs.qza)                
         -y      path to reference taxonomy for vsearch/blast classifier (default: $HOME/my_scripts/RubyRed/resources/classifiers/MIMt_16s/MIMt_16s_taxonomy.qza)                
         -z      classification method to use (sklearn, vsearch, or blast) (default: sklearn)                
         -c      path to the classifier for taxonomy assignment (default: $HOME/my_scripts/RubyRed/resources/classifiers/MIMt_16s/MIMt_16s_classifier.qza)                
         -f      minimum frequency filter for taxonomic classifications (default: 2)                
         -o      directory name to save results (will be created if it doesn't exist). (default: outputs_{name of input directory} )                
         -h      display this help message and exit





