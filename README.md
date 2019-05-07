# miniBarcoder 

### To make it easier for installing dependencies, we would recommend setting up using conda (

get Miniconda here: https://docs.conda.io/en/latest/miniconda.html, this pipeline is written in Python 2.7

```
conda config --add channels bioconda
conda create -n mbconda python=2.7 mafft racon graphmap blast seqtk git fasta3
conda activate mbconda
conda install -c anaconda biopython 

git clone https://github.com/asrivathsan/miniBarcoder/
cd miniBarcoder
python setup.py install 

Please run the test files for the pipeline. In a few computers, we are experience issues with racon obtained from bioconda, if this happens: please compile directly from github
git clone --recursive https://github.com/isovic/racon.git racon
cd racon
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cp build/bin/racon ~/miniconda2/envs/mbconda/bin/racon

```
Alternatively, if dependencies are being installed separately, you can just use the scripts directly, or proceed from "git clone " onwards. 

### To test the pipeline, with the recent updates use the following. This works if setup.py install has been done. Else call each as a python script, i.e. python mb_parallel_demultiplex.py ... 

```
cd testfiles
mb_parallel_demultiplex.py -d demfile_2.txt -l 600 -o testout -f test.fasta
mb_parallel_consensus.py -i testout
mv testout/all_barcodes.fa test_mafft_barcode.fa
filter_by_Ns.py -n 6 -i test_mafft_barcode.fa
racon_consensus.py -i testout -d racon_out -o test_racon.fa -b test_mafft_barcode_Nfilter.fa
aacorrection.py -b test_mafft_barcode_Nfilter.fa -bo test_barcodes_megablast.txt -bf test_barcodes_megablast.fasta -o test_mafft_barcode_aacorr.fa
aacorrection.py -b test_racon.fa -bo test_barcodes_megablast.txt -bf test_barcodes_megablast.fasta -o test_racon_barcode_aacorr.fa
consolidate.py -m test_mafft_barcode_aacorr.fa -r test_racon_barcode_aacorr.fa -t con_temp -o test_consolidate.fa



If fasta file is not available:

seqtk seq -A fastqfile > fastafile

If BLAST output file not available for aacorrection.py
aacorrection.py -b test_mafft_barcode_Nfilter.fa -d /path/to/nt -o test_mafft_barcode_aacorr.fa

In a few computers, we are experience issues with racon obtained from bioconda, if this happens: please compile directly from github
git clone --recursive https://github.com/isovic/racon.git racon
cd racon
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cp build/bin/racon ~/miniconda2/envs/mbconda/bin/racon


```

##### Other details on scripts for quality assessments such as assess_corrbarcodes_wref.py are given below if you would like to compare minion barcodes with reference barcodes. 



##### LEGACY minibarcoder.py associated with doi: 10.1111/1755-0998.12890


##### Given the multiple steps in the pipeline, we have written out the steps for obtaining barcodes for datasets A,B, and C for our manuscript. Easiest way to get a handle on the pipeline is to test it out on one of these datasets as described in Steps to operate the pipeline document.

In order to run Racon correction, both FASTQ and FASTA files are required as Racon requires quality score information: the order of the sequences must be same in the 2 files. This can be simply obtained if fasta is generated from FASTQ using tools such as fq2fa.



##### REQUIREMENTS: This has been tested on Ubuntu and MacOS 10.12.6.

##### miniBarcoder.py (MAFFT barcodes)
1.	Python2 (tested with 2.7)
2.	MAFFT v7
3.	glsearch36 (https://github.com/wrpearson/fasta36) This is part of fasta36 suite and one of the installed programs is "glsearch36"
4.	numpy  
MAFFT and glsearch36 must be in path and callable by "mafft" and "glsearch36" on the terminal.

##### aacorrection.py (+AA barcodes)
1.	Python2 (tested with 2.7)
2.	Biopython (tested with v 1.69)
3.	MAFFT v7
4.	BLAST+ (tested with 2.2+, 2.6+, 2.7+)
5.	(If BLAST output and corresponding accession fasta does not exist) Preferably: a local copy of the NT database. Alternatively, you can use blastn -remote option to run commandline blast on NCBI server, but this may not be as fast. If relying on -remote option, then ensure your BLAST is updated. Also note that for frequent use Entrez email should be provided.  
“mafft”, “blastn” should be in path and callable from terminal. If local nt is used “blastdbcmd” should be callable from terminal, available as part of BLAST+

##### RACON correction: 
This worked in Ubuntu but we couldn’t get it to work in Mac OS 10.12.6  (we couldn’t install racon).  
 
1)	graphmap (https://github.com/isovic/graphmap)  (v0.5.2) (released currently with racon). v0.3 gave errors.
2)	racon (https://github.com/isovic/racon) : v1.3.1
Both can be installed by: 
```
git clone https://github.com/isovic/racon.git && cd racon && make modules && make tools && make -j  
export PATH=$PATH:$PWD/bin
export PATH=$PATH:$PWD/tools/graphmap/bin/Linux-x64
```
These must be installed and in path  and calling “graphmap” and “racon” in terminal should not throw an error.   


### DETAILS  

##### miniBarcoder.py : Script for obtaining MAFFT barcodes

This script performs the various steps from primer finding, demultiplexing, read alignment and consensus calling.  

For the purpose of the is pipeline, please avoid putting “;” in specimen name.  
Output is stored in outpudirectory/  

Basic usage for 658 bp barcode (100 reads used):  
```
python miniBarcoder.py –f inputfasta –d demultiplexingfile –o outputdirectory –l 600
```
Basic usage for 658 bp barcode (all reads used):
```
python miniBarcoder.py –f inputfasta –d demultiplexingfile –o outputdirectory –l 600 –D 0
```
Output is stored in outpudirectory/all_barcodes.fa
```
Arguments:
  -h, --help            show this help message and exit

  -f INFASTA, --infasta INFASTA
                        Path to input fasta file
  -d DEMFILE, --demfile DEMFILE
                        Path to demultiplexing file
  -o OUTDIR, --outdir OUTDIR
                        set output directory path
  -l MINLEN, --minlen MINLEN
                        exclude barcode sequences identified that are shorter
                        than specified length
  -m MODE, --mode MODE  run with unique tag mode (1) or dual tag mode (2),
                        mode 1 ignores mismatch setting and allows only 1 bp
                        mismatch, default 2
  -mm int, --mismatch int
                        number of mismatches allowed in tags, must be <=5,
                        default 2
  -e EVALUE, --evalue EVALUE
                        evalue for primer search using glsearch36,default 1e+6
  -g GAPS, --gaps GAPS  number of gaps allowed for primer identification,
                        default 5
  -D MAXDEPTH, --maxdepth MAXDEPTH
                        set max depth per coverage to improve speed, default
                        100, must be >2
  -t THREADS, --threads THREADS
                        number of threads for glsearch and mafft
  -bl BLEN, --barcodelength BLEN
                        estimated barcode length, used for unique tag mode
                        only, please keep it slightly shorter than actual barcode length, as                                insertion errors can make this run into primers, default 300

```
##### aacorrection.py : Script for correcting barcodes using conserved amino acids 
Basic usage for 658 bp COI  

MODE 1: Existing BLAST output and accession fasta file  
```
python aacorrection.py –b input_uncorrected_barcode_fasta –bf blast_accession_fasta_file –bo  blastout_output file –o outputfilename
``` 
MODE 2: Conducting BLAST to local NT 
```
python aacorrection.py –b input_uncorrected_barcode_fasta –d /path/to/nt/withprefix –o outputfilename
```
MODE 3: Conducting BLAST to remote NT (SLOW) 
```
python aacorrection.py –b input_uncorrected_barcode_fasta –d “nt –remote” –o outputfilename
```
Commands
```
arguments:
  -h, --help            show this help message and exit
  -b INFASTA, --barcodes INFASTA
                        Path to input barcode fasta file
  -o OUTFILE, --outfile OUTFILE
                        outfile file name
  -p THREADS, --threads THREADS
                        number of threads for BLAST,default=4
  -bo BLASTOUTFILE, --blastout BLASTOUTFILE
                        Path to blast output file, outputformat 6
  -bf BLASTACCFILE, --blastfasta BLASTACCFILE
                        Path to fasta file containing sequences of BLAST hits,
                        required if -bo or --blastout is given
  -d PATH_TO_DB, --db PATH_TO_DB
                        Path to nucleotide database with database prefix, if
                        local copy is unavailable you can try typing 'nt
                        -remote'. note that remote has not been extensively
                        tested and is slower
  -a NAMBS, --amb NAMBS
                        proportion of ambiguities allowed per barcode,
                        default=0.01
  -l MINLEN, --minlen MINLEN
                        exclude barcodes shorter than this length, default=640
  -L MAXLEN, --maxlen MAXLEN
                        exclude barcodes longer than this length, default=670
  -c CONGAPS, --congaps CONGAPS
                        exclude sequences containing any gap of length >=
                        value, default=5
  -n NAMINO, --namino NAMINO
                        number of flanking amino acids around the gap used for
                        correction, default=3
  -g GENCODE, --gencode GENCODE
                        genetic code https://www.ncbi.nlm.nih.gov/Taxonomy/Uti
                        ls/wprintgc.cgi, default=5, invertebrate mitochondrial
  -e EVALUE, --evalue EVALUE
                        e-value for BLAST search, default=1e-5
  -H HPLEN, --hplen HPLEN
                        minimum homopolymer length, default=2
  -s SUPPORT, --support SUPPORT
                        minimum support for indel in references. Reducing this
                        increases chances of errors and improper detection of
                        reading frames,default=2

```
Other scripts:  
##### run_racon_consensus.sh : Batch script performing fastq retrieval, graphmap and racon
Usage:   
```
sh racon_consensus.sh input_fastq_for_all_data input_fasta_for_all_data outputfolder_of_minibarcoder mafft_barcodes_obtained_by_minibarcoder outputdirectory
```
##### consolidate.py: Get consensus barcode of MAFFT+AA and RACON+AA barcode.  

usage: 
```
python consolidate.py –m mafft_aa_barcodefile –r racon_aa_barcodefile –o output_barcode_file
```
```
  -h, --help            show this help message and exit
  -m MAFFTB, --mafft MAFFTB
                        Path to input mafft corrected barcode fasta file
  -r RACONB, --racon RACONB
                        Path to input racon corrected barcode fasta file
  -o OUTFILE, --outfile OUTFILE
                        Path to output corrected barcode fasta file

  -o OUTDIR, --outdir OUTDIR
                        output directory
```

##### filter_by_Ns.py: filter barcode fasta file to remove sequences with lots of ambiguities
```
python filter_by_Ns.py -i input_barcode_fasta -n number of ambiguities
```
Output is stored in *Nfilter.fasta  
```
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        Path to input fasta file
  -n NAMBS, --nambs NAMBS
                        Number of ambiguities allowed

```
##### assess_corrbarcodes_wref.py: Compare corrected barcodes with reference barcodes to measure accuracy.
```
  -h, --help            show this help message and exit
  -m MINIONFASTA, --minionfasta MINIONFASTA
                        Path to input mafft corrected barcode fasta file
  -r REFFASTA, --reffasta REFFASTA
                        Path to input reference fasta file
  -t OUTDIR, --tempdir OUTDIR
                        Path to temporary directory
  -o OUTFILE, --outfile OUTFILE
                        Path to output file containing statistics
```

##### assess_corrbarcodes_wref.py: Compare uncorrected barcodes with reference barcodes to measure accuracy.
```
  -h, --help            show this help message and exit
  -m MINIONFASTA, --minionfasta MINIONFASTA
                        Path to input mafft corrected barcode fasta file
  -r REFFASTA, --reffasta REFFASTA
                        Path to input reference fasta file
  -t OUTDIR, --tempdir OUTDIR
                        Path to temporary directory
  -o OUTFILE, --outfile OUTFILE
                        Path to output file containing statistics
```
##### measure_ambs.py: Measures number of “N” in input non-interleaved barcode fastafile
```
python measure_ambs.py input_barcode_fasta
```

## Contact Information
For bugs, queries and other issues, please contact Amrita Srivathsan, asrivathsan@gmail.com

## Funding
This project was funded by Southeast Asian Biodiversity Genomics Centre (SEABIG), NUS and ASTAR

## Publication
The publication associated is currently in preprint server: https://doi.org/10.1101/253625
