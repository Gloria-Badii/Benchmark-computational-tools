#!/user/bin/env bash

# Operative System: Ubuntu (Linux)
# required: RiboFrame (downloaded from https://github.com/matteoramazzotti/riboFrame/tree/master )



##########    SETTING VARIABLES    ##############


### TO USE THIS SCRIPT AS A "PROGRAM"
shopt -s expand_aliases


### WHERE IS RIBOFRAME
alias riboTrap="/home/matteo/Desktop/programs/riboFrame/riboTrap.pl"
alias riboMap="/home/matteo/Desktop/programs/riboFrame/riboMap.pl"


### RDP VARIABLES
RDP_Folder="/home/matteo/Desktop/programs/rdp_classifier_2.14"  # rdp classifier
Retrained_RDP_Folder="/home/matteo/Desktop/programs/RDP_on_SILVA_138_2"
hmms_folder="/home/matteo/Desktop/programs/riboFrame/hmms" # hmmer3 16S models folder


### SETTING HOW THE FASTQ ARE NAMED
alias extract_name="sed 's/_R[1-2].*//g' | uniq"
# this command is used to modify and then grep the unique file names from paired end files in the loops

### THREADS TO USE
threads=62





################     PREPARING THE DIRECTORIES          ####################


mkdir RiboFr_temp
mkdir RiboFr_output




################    STARTING THE ANALYSIS WITH RIBOFRAME    #####################


# NB: RiboFrame searches the files according to their name --> do not change the file suffixes in the following lines


FOLDER_INPUT="processed_FASTQ"
SAMPLE_LIST=$(ls $FOLDER_INPUT | grep ".fastq" | extract_name )


for i in $SAMPLE_LIST    
do

echo -e "\n\n\n\n *** Preprocessing the sample $i with RiboFrame... *** \n"

# declaring the input for this cicle
TARGETS=$(ls $FOLDER_INPUT | grep $i)
FOR=$(echo "$TARGETS" | grep "_R1" )
REV=$(echo "$TARGETS" | grep "_R2" )


# transforming the fastq in fasta
zcat $FOLDER_INPUT/$FOR | grep "^@" -A 1 --no-group-separator | sed "s/@/>/g" | sed "s# 1.*#/1#g" > RiboFr_temp/$i.1.fasta
zcat $FOLDER_INPUT/$REV | grep "^@" -A 1 --no-group-separator | sed "s/@/>/g" | sed "s# 2.*#/2#g" > RiboFr_temp/$i.2.fasta
# NB: do not remove this fasta from the temp directory until riboTrap.pl is executed

# flagging which sequences are 16S through HMM
# NB: each forward and reverse has to be analysed with both the for and rev hmm model

# FORWARD
hmmsearch -E 0.00001 --domtblout RiboFr_temp/${i}.1.fwd.bact.ribosomal.table --noali --cpu $threads -o /dev/null $hmms_folder/16S_bact_for3.hmm   RiboFr_temp/$i.1.fasta
hmmsearch -E 0.00001 --domtblout RiboFr_temp/${i}.1.rev.bact.ribosomal.table --noali --cpu $threads -o /dev/null $hmms_folder/16S_bact_rev3.hmm   RiboFr_temp/$i.1.fasta
hmmsearch -E 0.00001 --domtblout RiboFr_temp/${i}.1.fwd.arch.ribosomal.table --noali --cpu $threads -o /dev/null $hmms_folder/16S_arch_for3.hmm   RiboFr_temp/$i.1.fasta
hmmsearch -E 0.00001 --domtblout RiboFr_temp/${i}.1.rev.arch.ribosomal.table --noali --cpu $threads -o /dev/null $hmms_folder/16S_arch_rev3.hmm   RiboFr_temp/$i.1.fasta
# REVERSE
hmmsearch -E 0.00001 --domtblout RiboFr_temp/${i}.2.fwd.bact.ribosomal.table --noali --cpu $threads -o /dev/null $hmms_folder/16S_bact_for3.hmm   RiboFr_temp/$i.2.fasta
hmmsearch -E 0.00001 --domtblout RiboFr_temp/${i}.2.rev.bact.ribosomal.table --noali --cpu $threads -o /dev/null $hmms_folder/16S_bact_rev3.hmm   RiboFr_temp/$i.2.fasta
hmmsearch -E 0.00001 --domtblout RiboFr_temp/${i}.2.fwd.arch.ribosomal.table --noali --cpu $threads -o /dev/null $hmms_folder/16S_arch_for3.hmm   RiboFr_temp/$i.2.fasta
hmmsearch -E 0.00001 --domtblout RiboFr_temp/${i}.2.rev.arch.ribosomal.table --noali --cpu $threads -o /dev/null $hmms_folder/16S_arch_rev3.hmm   RiboFr_temp/$i.2.fasta


# tagging the reads identified as 16S according to HMM
riboTrap RiboFr_temp/$i


mkdir RiboFr_output/${i/_L001/}   # output folder for each sample  (NB: the suffix "L001" can be modified according to the samples original names...)


# using RDP to classify the sequences
java -Xmx10g -jar $RDP_Folder/dist/classifier.jar  classify -c 0.8  -t $Retrained_RDP_Folder/rRNAClassifier.properties  -o RiboFr_output/${i/_L001/}/${i}.16S.rdp  -q RiboFr_temp/${i}.16S.fasta
# re-formatting for riboMap (required only for the re-trained classifier)
sed -e "s/Domain/domain/g" -e "s/Phylum/phylum/g" -e "s/Class/class/g" -e "s/Order/order/g" -e "s/Family/family/g" -e "s/Genus/genus/g" -i RiboFr_output/${i/_L001/}/${i}.16S.rdp


# extracting the results using all the hypervariable regions of 16S gene (var=full) and V3-V4 (var=V3,V4)
riboMap file=RiboFr_output/${i/_L001/}/${i}.16S.rdp var=full thr=0.8 cross=over out=RiboFr_output/${i/_L001/}/${i}_full_count
riboMap file=RiboFr_output/${i/_L001/}/${i}.16S.rdp var=V3,V4 thr=0.8 cross=over out=RiboFr_output/${i/_L001/}/${i}_V3_V4
# var = region, thr = confidence, percmin = abund_threshold
# NB: activating the plotting options of riboMap causes a bug in the program... do not turn these options on!


# rm  RiboFr_temp/*fasta   # big files, no more needed


done


rm RiboFr_temp -R
