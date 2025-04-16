#!/user/bin/env bash

# Operative System: Ubuntu (Linux)
# required: kraken2, braken



###############      SETTING VARIABLES        ###########


### TO USE THIS SCRIPT AS A "PROGRAM"
shopt -s expand_aliases
source ~/miniconda3/etc/profile.d/conda.sh
# these rows load the environment variable "into the script"

### SETTING HOW THE FASTQ ARE NAMED
alias extract_name="sed 's/_R[1-2].*//g' | uniq"
# this command is used to modify and then to grep the unique file names from paired end files during the loops

### KRAKEN2 DATABASE
kraken_db='/media/matteo/SSD4/reference/kraken2_core_nt_db_v20241228/'

### THREADS TO USE
threads=60




##############     PREPARING THE DIRECTORIES    ###########


mkdir kraken_output




#########     STARTING THE ANALYSIS WITH KRAKEN2    ##########

conda activate kraken2

FOLDER_INPUT="processed_FASTQ"
SAMPLE_LIST=$(ls $FOLDER_INPUT | grep "CLEANED" | extract_name )


# Using Kraken2 on each file
for i in $SAMPLE_LIST
do

echo -e "\n\n\n\n *** Processing the sample $i with Kraken... *** \n"

TARGETS=$(ls $FOLDER_INPUT | grep $i)
FOR=$(echo "$TARGETS" | grep "_R1" )
REV=$(echo "$TARGETS" | grep "_R2" )

kraken2 --threads $threads --db $kraken_db --output kraken_output/${i}_output.txt --report kraken_output/${i}_report.txt --confidence 0.1 --paired $FOLDER_INPUT/$FOR $FOLDER_INPUT/$REV

bracken -r 150 -d $kraken_db -i kraken_output/${i}_report.txt  -l G  -o kraken_output/${i}_bracken_abund_Genus.tsv  # genus level
bracken -r 150 -d $kraken_db -i kraken_output/${i}_report.txt  -l S  -o kraken_output/${i}_bracken_abund_Species.tsv  # species level
# -r is the reads length

# Few entries have no "genus" name according to NCBI, hence they aren't featured in the genus report, even if they are abundant according to the 16S processing --> adding them manually from the species level report when possible
grep "Moraniibacteriota bacterium" kraken_output/${i}_bracken_abund_Species.tsv  > missing_genus_to_add_from_species.txt
grep "Kaiserbacteria bacterium" kraken_output/${i}_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt
grep "Nomurabacteria bacterium" kraken_output/${i}_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt
grep "SH-PL14" kraken_output/${i}_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt
#grep "Magasanik" kraken_output/${i}_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt   # none
grep "Saccharimonadia bacterium" kraken_output/${i}_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt
grep "Solirubrobacterales bacterium" kraken_output/${i}_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt   # 67-14
#grep -e "OPS17|OPS 17" kraken_output/${i}_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt # none
#grep -e "37-13" kraken_output/${i}_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt # none
#grep -e "Actinobacteria bacterium IMCC26207" kraken_output/${i}_bracken_abund_Species.tsv  >> missing_genus_to_add_from_species.txt # none (NB: different codes)
cat kraken_output/${i}_bracken_abund_Genus.tsv  missing_genus_to_add_from_species.txt > temp.tsv
mv temp.tsv kraken_output/${i}_bracken_abund_Genus.tsv
rm missing_genus_to_add_from_species.txt


done


# kraken2-inspect --db $kraken_db --threads 30 | cut -f 4,6 > FASTQ_check/Unclassified_Kraken/Actual_db.tsv


conda deactivate


# optimizing the required space to store files
rm kraken_output/*_output.txt




##########     EXTRA: QUALITY CONTROL AFTER THE CLEANING     ###########

# mkdir FASTQ_check/fastqc_processed_reads

# for i in ./processed_FASTQ/*CLEANED_* ; do fastqc $i -o FASTQ_check/fastqc_processed_reads --threads $threads; done

# multiqc FASTQ_check/fastqc_processed_reads -o FASTQ_check/fastqc_processed_reads

# rm FASTQ_check/fastqc_processed_reads/*fastqc.zip


# open FASTQ_check/fastqc_processed_reads/multiqc_report.html
# cat FASTQ_check/fastqc_processed_reads/multiqc_data/multiqc_general_stats.txt | cut -f 1,3,5,7,9 | sed 's/FastQC_mqc-generalstats-fastqc-//g' > FASTQ_check/fastqc_processed_reads/General_info_Sequences_recap_AFTER_THE_CLEANING.txt
# cat FASTQ_check/fastqc_processed_reads/General_info_Sequences_recap_AFTER_THE_CLEANING.txt



# fastqc does not work by remote (?) but I still need this info
for i in ./processed_FASTQ/*gz
do
echo "$i : $(zcat $i | grep "^@" -c)" >> FASTQ_check/remaining_number_of_seqs_after_process.txt
done
sed -e "s#./raw_FASTQ_dir/##g" -e "s/.gz//g" -i FASTQ_check/remaining_number_of_seqs_after_process.txt

