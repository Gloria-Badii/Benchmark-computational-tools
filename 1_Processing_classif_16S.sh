#!/usr/bin/env

# Operative System: Debian GNU/Linux 11 (bullseye)
# programs and envs required: conda, qiime2-amplicon-2024.2, picrust2 2.4.2 (the latter is optional)




#################### SETTING VARIABLES ############################


### TO USE THIS SCRIPT AS A "PROGRAM"
shopt -s expand_aliases
source ~/miniconda3/etc/profile.d/conda.sh
# these rows load the environment variable "into the script"

### SILVA DATABASE
SILVA_trained_db=/media/matteo/SSD1/reference/SILVA_138_V3V4_BayesClassifier.qza

### THREADS TO USE
threads=60



##########################   PREPARING DIRECTORIES  ########################


mkdir QIIME
chmod +rwx QIIME




########################### PROCESSING RAW FILES #############################

##### FASTQ CHECKSUM

md5sum 16S_raw_fastq/* | sed 's+16S_raw_fastq/++g' | sed 's/  / \t/' > QIIME/checksum_FASTQ.tsv

dir 16S_raw_fastq/ | sed '/R2_001/d' >R1.txt
cat R1.txt | sed 's/R1_001/R2_001/' >R2.txt
paste R1.txt R2.txt --delimiters='\t' > QIIME/FASTQ_pairs.tsv
rm R1.txt R2.txt



##### OPENING QIIME ENVIRONMENT

conda activate qiime2-amplicon-2024.2

cd QIIME

qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path ../16S_raw_fastq --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path demux-paired-end.qza



##### PRIMER TRIMMING

qiime cutadapt trim-paired --i-demultiplexed-sequences demux-paired-end.qza --p-front-f CCTACGGGNBGCWSCAG --p-front-r GACTACNVGGGTWTCTAATCC --p-cores $threads --o-trimmed-sequences trimmed-seqs.qza --p-discard-untrimmed --p-minimum-length 100 --verbose >Log_di_cutadapt.txt

grep -i "Pairs written (passing filters)" Log_di_cutadapt.txt >Log_reads_with_primers.txt

cat Log_reads_with_primers.txt



##### DENOISING, QUALITY LENGTH TRIMMING and MERGING

qiime demux summarize --i-data trimmed-seqs.qza --o-visualization demux.qzv --p-n 500000 

# qiime tools view demux.qzv
#echo -e "\n put that demux.qzv in this site https://view.qiime2.org/ \n"


# the reads are cut to allowing an overlap up to 20nt (see https://bioconductor.org/packages/devel/bioc/vignettes/dada2/inst/doc/dada2-intro.html)
qiime dada2 denoise-paired --i-demultiplexed-seqs trimmed-seqs.qza --p-trim-left-f 0 --p-trim-left-r 0 --p-trunc-len-f 272 --p-trunc-len-r 170 --o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza --p-n-threads $threads

qiime tools export --input-path denoising-stats.qza --output-path ./      ### export and create the table of absolute and relative read abundances
cat stats.tsv | grep -v "q2:types" > DADA2_Initial_and_processed_reads.tsv
cat DADA2_Initial_and_processed_reads.tsv | sed 's/age of input /_/g'| sed 's/passed filter/filtered/' | cut -f 1,2,4,7,9


chmod +rwx stats.tsv




######################### TAXONOMIC CLASSIFICATION ##################################

 
### BAYESIAN CLASSIFICATION (using Scikit-learn retrained on V3-V4 SILVA 16S 138)
qiime feature-classifier classify-sklearn --i-classifier $SILVA_trained_db --i-reads rep-seqs.qza --o-classification taxonomy_SILVA.qza
qiime metadata tabulate --m-input-file taxonomy_SILVA.qza --o-visualization confidence_taxonomy_SILVA.qzv




############# EXPORTING THE ORIGINAL NUMBER OF READS ############################

qiime demux summarize --i-data demux-paired-end.qza --o-visualization original_reads.qzv 
qiime tools export --input-path original_reads.qzv  --output-path Raw_Reads_Info
mv Raw_Reads_Info/per-sample-fastq-counts.tsv Original_number_of_reads_for_sample.tsv
chmod +rw Original_number_of_reads_for_sample.tsv
chmod +rw Raw_Reads_Info -R
rm Raw_Reads_Info -R


############# CLEANING THE FOLDER FROM QIIME TEMP FILES #####################

rm denoising-stats.qza trimmed-seqs.qza demux-paired-end.qza stats.tsv original_reads.qzv 




