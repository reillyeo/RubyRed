# RubyRed

Custom bioinformatics pipeline designed for high-throughput, read-by-read taxonomic classification of rRNA gene amplicons generated via Oxford Nanopore sequencing. RubyRed uses open-source tools and custom scripts to process raw sequencing data into taxonomically annotated feature tables suitable for downstream data analysis.

Input data can be either raw or demultiplexed FASTQ files. Demultiplexing (if necessary) is done using Guppy Barcoder, trimming of primer binding regions with [Cutadapt](https://github.com/marcelm/cutadapt), quality and length filtering with [Chopper](https://github.com/wdecoster/chopper), and subsequently converted to FASTA format using [VSEARCH](https://github.com/torognes/vsearch). To ensure data consistency, sequences with fewer than a minimum read count are discarded, and samples exceeding a predefined threshold are subsampled using [SeqKit](https://github.com/shenwei356/seqkit).

Next, all filtered reads are concatenated and imported into [QIIME2](https://docs.qiime2.org) as a single sequence artifact. A custom Python script is then employed to generate a feature table. Chimera removal is performed using VSEARCH’s uchime-ref algorithm against a curated rRNA gene reference database, and surviving sequences are reoriented with [RESCRIPt](https://github.com/bokulich-lab/RESCRIPt) to match reference strand orientation.

For taxonomic assignment, the pipeline supports three classification methods: scikit-learn, VSEARCH, and BLAST. The scikit-learn method requires a pre-trained classifier, but is considerably faster than either of the other methods. Low abundance classifications (below a user-defined frequency threshold) are filtered to reduce noise in the final dataset. The end products include a taxonomically annotated feature table and representative sequences.


## Installation

Install RubyRed by cloning or forking this repository. (Note: the paths used by the default parameters assume that the RubyRed directory exists at the location $HOME/my_scripts/RubyRed). 

The fasta file containing primer sequences should be edited/replaced with whatever primer sequences you used.  

The Guppy barcoder binary is required to be downloaded and added to your $PATH if you want RubyRed to demultiplex your data. Otherwise, demultiplex prior to starting, and use -d flag.

QIIME2 amplicon distribution must be downloaded in a conda environment called qiime2:

    conda env create -n qiime2 --file https://data.qiime2.org/distro/amplicon/qiime2-amplicon-2024.10-py310-osx-conda.yml

Activate your qiime2 environment and install Chopper and SeqKit:

     conda activate qiime2
     conda install -c bioconda chopper seqkit

All other required packages are already included in QIIME2.

## Usage

```text
Usage - RubyRed [OPTIONS]
                 
Options:                 
         -i      directory containing input FASTQ files (default: current directory)                
         -d      use this flag if data has already been demultiplexed                
         -q      minimum average read quality required to pass chopper quality filtering (default: 15)                
         -l      minimum read length allowed to pass chopper length filtering (default: 800)                
         -x      maximum read length allowed to pass chopper length filtering (default: 1600)                
         -p      number of threads to use for parallel processing (default: 20)                
         -m      minimum number of reads (post-filtering) to keep a file (default: 100)                
         -t      subsample fasta files with more than this number of reads (default: 30000)                
         -s      directory where python scripts are located (default: $HOME/my_scripts/RubyRed/scripts)                
         -w      path to the directory where resources (primer seqs, reference seqs, reference taxonomy, classifier) are located (default: $HOME/my_scripts/RubyRed/resources)                
         -u      path to fasta file containing primer sequences (default: $HOME/my_scripts/RubyRed/resources/UMI16s_primers.fasta)                
         -r      path to reference sequences for chimera filtering and reorientation (default: $HOME/my_scripts/RubyRed/resources/classifiers/MIMt/MIMt_refseqs.qza)                
         -y      path to reference taxonomy for vsearch/blast classifier (default: $HOME/my_scripts/RubyRed/resources/classifiers/MIMt/MIMt_taxonomy.qza)                
         -z      classification method to use (sklearn, vsearch, or blast) (default: sklearn)                
         -c      path to the classifier for taxonomy assignment (default: $HOME/my_scripts/RubyRed/resources/classifiers/MIMt/MIMt_nb_classifier.qza)                
         -f      minimum frequency filter for taxonomic classifications (default: 2)                
         -o      directory name to save results (will be created if it doesn't exist). (default: outputs_{name of input directory} )                
         -h      display this help message and exit







