#!/user/bin/env bash

# required: bash, Kaiju v1.10 and R
# kaiju pre-indexed databases are available at https://bioinformatics-centre.github.io/kaiju/downloads.html

# PATH="/home/matteo/Desktop/programs/kaiju_v1_10_1/bin:$PATH"




#########    SETTING VARIABLES    #######


### TO USE THIS SCRIPT AS A "PROGRAM"
shopt -s expand_aliases

### SETTING HOW THE FASTQ ARE NAMED
alias extract_name="sed 's/_R[1-2].*//g' | uniq"
# this command is used to modify and then to grep the unique file names from paired end files during the loops

### WHERE KAIJU IS LOCATED
alias kaiju="/home/matteo/Desktop/programs/kaiju_v1_10_1/bin/kaiju"
alias kaiju2table="/home/matteo/Desktop/programs/kaiju_v1_10_1/bin/kaiju2table"

### KAIJU DATABASE
Kaiju_DB=/media/matteo/SSD4/reference/kaiju_db_nr_euk_2023-05-10

### THREADS TO USE
threads=45




##########    PREPARING THE DIRECTORIES    ########


mkdir Kaiju_output




####### STARTING THE ANALYSIS WITH KAIJU #######


FOLDER_INPUT="processed_FASTQ"
SAMPLE_LIST=$(ls $FOLDER_INPUT | grep ".fastq" | extract_name )

# Using Kaiju on each file ...
for i in $SAMPLE_LIST     # main loop
do
TARGETS=$(ls $FOLDER_INPUT | grep $i)
FOR=$(echo "$TARGETS" | grep "_R1" )
REV=$(echo "$TARGETS" | grep "_R2" )


# classification with Kaiju
kaiju -t $Kaiju_DB/nodes.dmp -f $Kaiju_DB/kaiju_db_*.fmi -i $FOLDER_INPUT/$FOR -j $FOLDER_INPUT/$REV  -o Kaiju_output/${i}_Kaiju_output.txt  -z $threads  -x -v
# E-value default: -E= 0.01
# score default: -s= 65
# matching length default: -m= 11 aa 
# -x enables SEG algorithm for filtering low-complexity region (suggested by the official tutorial)
# -v enables the addiction of other columns to the output
# nr_euk: Archaea, bacteria and viruses (nr) +  fungi and microbial eukaryotes

done


#################      SUMMARY OF RESULTS      ###############  
  
  
output_names=$(ls Kaiju_output | grep "output.txt" | sed 's#^#Kaiju_output/#g')  # List of inputs for the row below, with path added at the beginning

# inputs written as last argument
# -r choose the depth of taxonomic classification
# -m filters out taxa (of "r" level") with relative abundance below m threshold
# -l chooses which taxon level will be printed in the output
# -p adds the full taxonomic name (from phylum to species), not just the identified one (can't be used with l)



kaiju2table -t $Kaiju_DB/nodes.dmp -n $Kaiju_DB/names.dmp -r species -l species -o Kaiju_output/kaiju_classif_summary_species_level.tsv $output_names
kaiju2table -t $Kaiju_DB/nodes.dmp -n $Kaiju_DB/names.dmp -r genus -p -o Kaiju_output/kaiju_classif_summary_full.tsv Kaiju_output/kaiju_classif_summary.tsv $output_names


# Few entries are missing because they lack of the genus level in NCBI (but not in SILVA)! Taking these rows from the species level.
grep "Moranbacteria bacterium" Kaiju_output/kaiju_classif_summary_species_level.tsv  > missing_genus_to_add_from_species.txt
# NB: Moranii is Moranbacteria in Kaiju db, differently from Kraken2 nt_core
grep "Kaiserbacteria bacterium" Kaiju_output/kaiju_classif_summary_species_level.tsv  >> missing_genus_to_add_from_species.txt
grep "Nomurabacteria bacterium" Kaiju_output/kaiju_classif_summary_species_level.tsv  >> missing_genus_to_add_from_species.txt
grep "SH-PL14" Kaiju_output/kaiju_classif_summary_species_level.tsv  >> missing_genus_to_add_from_species.txt
grep "Magasanik" Kaiju_output/kaiju_classif_summary_species_level.tsv  >> missing_genus_to_add_from_species.txt   # none
grep "Saccharimonadia bacterium" Kaiju_output/kaiju_classif_summary_species_level.tsv  >> missing_genus_to_add_from_species.txt
grep "Solirubrobacterales bacterium" Kaiju_output/kaiju_classif_summary_species_level.tsv  >> missing_genus_to_add_from_species.txt   # 67-14
grep -e "OPS17|OPS 17" Kaiju_output/kaiju_classif_summary_species_level.tsv  >> missing_genus_to_add_from_species.txt # none
grep -e "37-13" Kaiju_output/kaiju_classif_summary_species_level.tsv  >> missing_genus_to_add_from_species.txt # none
grep -e "Actinobacteria bacterium IMCC26207" Kaiju_output/kaiju_classif_summary_species_level.tsv  >> missing_genus_to_add_from_species.txt # none (NB: different codes)

# Transforming the species name (various species) ensuring the manual parsing between SILVA and NCBI databases (which is coded in another script)
sed "s/Actinobacteria bacterium IMCC26207/IMCC26207/g" -i missing_genus_to_add_from_species.txt
sed "s/ bacterium.*/ bacterium/g" -i missing_genus_to_add_from_species.txt
sed "s/uncultured //g" -i missing_genus_to_add_from_species.txt
cat Kaiju_output/kaiju_classif_summary.tsv  missing_genus_to_add_from_species.txt > temp.tsv   # both of the files (original and updated) in an unique file
mv temp.tsv Kaiju_output/kaiju_classif_summary.tsv   # overwriting the original one
rm missing_genus_to_add_from_species.txt


# Reformatting Kaiju's output
sed "s/Candidatus /Candidatus_/g"   -i Kaiju_output/kaiju_classif_summary.tsv
sed "s/;$//g"   -i Kaiju_output/kaiju_classif_summary.tsv   # removing only the second ';' at the end of the lines
sed -i '/unclassified/d' Kaiju_output/kaiju_classif_summary.tsv   # removing unclassified reads
sed -i '/cannot be assigned to a (non-viral) genus/d' Kaiju_output/kaiju_classif_summary.tsv



  