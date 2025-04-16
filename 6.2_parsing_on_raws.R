library(reshape2)

#using the function "Import_from_Tot_Seq". If not in the environmrnt retrieve it from the preavious scritpt "5.1_parsing_on_percentages.txt"



### import from RiboFrame (full counts)
Import_from_Tot_Seq(input_folder="Obtained_Counts/", table_suffix="_full_count.genus.cnt", column_tax="Name", column_abund="Count" )
Counts_RiboF_full<-table_obtained
row.names(Counts_RiboF_full)<-Counts_RiboF_full$Genus

### import from RiboFrame (V3-V4)
Import_from_Tot_Seq(input_folder="Obtained_Counts/", table_suffix="_V3_V4.genus.cnt", column_tax="Name", column_abund="Count" )
Counts_RiboF_V3V4<-table_obtained
row.names(Counts_RiboF_V3V4)<-Counts_RiboF_V3V4$Genus

### import from Kraken2
Import_from_Tot_Seq(input_folder="Obtained_Counts/", table_suffix="_bracken_abund_Genus", column_tax="name", column_abund="new_est_reads" )
Counts_Kraken<-table_obtained
row.names(Counts_Kraken)<-Counts_Kraken$Genus







#################   TAXONOMY ADJUSTMENTS




### working on 16S feature table
true_16_count<- read.table(file="Obtained_Counts/true_16S_counts.tsv", sep="\t", row.names=1, header=T)
true_16_tax<- read.table(file="Obtained_Counts/true_16S_taxonomy.tsv", sep="\t", header=T)
true_16_tax<-true_16_tax[ , c(1,2)]   # removing the confidence column
true_16_tax$Taxon<-gsub("; s__.*","",true_16_tax$Taxon)   # removing info about species
no_genus_classif<- ! grepl("g__", true_16_tax$Taxon )    # annotating where the classifier did not reach the genus level
no_fam_classif<- ! grepl("f__", true_16_tax$Taxon )    # and where did not reach the family either
phylum<-gsub(".*p__","", true_16_tax$Taxon)
phylum<-gsub(";.*","", phylum)
family<-gsub(".*f__","", true_16_tax$Taxon)
family<-gsub(";.*","", family)
genus <-  true_16_tax$Taxon
genus[ no_genus_classif ] <- "NA"
genus<-gsub(".*g__","", genus)
genus<-gsub(";.*","", genus)
# adding more info where required
genus_new <- genus
genus_new [ genus %in% c("NA","uncultured") ]<- paste0(genus[genus %in% c("NA","uncultured")] , "_f_", family[genus %in% c("NA","uncultured")] )
# cleaning where there is not info at family level either (e.g. "NA_f_NA")
genus_new [ grepl("_f_uncultured", genus_new) | no_fam_classif ] <- gsub("_f_.*","", genus_new [ grepl("_f_uncultured", genus_new) | no_fam_classif ] )
# if the genus name is just "NA" or "uncultered" again then the family name was also missing... then adding the phylum name instead
genus_new [ genus %in% c("NA","uncultured") ]<- paste0(genus[genus %in% c("NA","uncultured")] , "_p_", phylum[genus %in% c("NA","uncultured")] )
true_16_tax$Taxon<-genus_new   # updating
# creating an unique 16S table
true_16_count <- true_16_count[ true_16_tax$Feature.ID , ]   # NB: same observation order
Counts_true_16<-cbind.data.frame(true_16_count , "Genus"=true_16_tax$Taxon )
Counts_true_16 <- aggregate( .~Genus, Counts_true_16, FUN= sum ) # aggregating at Genus level
row.names(Counts_true_16)<- Counts_true_16$Genus



# working on Kaiju table
kaiju_count <- read.table(file="Kaiju_output/kaiju_classif_summary.tsv", sep="\t", header=T)
kaiju_count <- kaiju_count[ , c(1,3,5)]
kaiju_wide <- dcast(kaiju_count, taxon_name ~ file, value.var = "reads", fun.aggregate = function(x) x[1])
kaiju_wide[is.na(kaiju_wide)] <- 0 
row.names(kaiju_wide) <- kaiju_wide$taxon_name





# translating names between SILVA and NCBI (at least the most abundants)
total_16S <- rowSums( Counts_true_16[ , colnames(Counts_true_16)!="Genus" ])
total_16S <- sort(total_16S, decreasing=TRUE)
top_16S <- names( total_16S [1:50] )
Counts_Kraken$Genus <- gsub("Candidatus ","Candidatus_", Counts_Kraken$Genus)
top_Kraken <- names( sort (rowSums( Counts_Kraken[, colnames(Counts_Kraken)!="Genus"] ), decreasing=TRUE ) ) [1:50]
top_Kraken <-gsub(" ", "_", top_Kraken)
unique_names_Kraken <- top_Kraken [!top_Kraken %in% top_16S ]
unique_names_SILVA <- top_16S [!top_16S %in% top_Kraken ] 
# different "tops" between the two database ... but, maybe, the few of these names are still featured in the other object, only with smaller counts...
unique_names_Kraken <- unique_names_Kraken[ ! unique_names_Kraken%in% Counts_true_16$Genus ]
unique_names_SILVA<- unique_names_SILVA[ ! unique_names_SILVA%in% Counts_Kraken$Genus ]
# unique_names_Kraken
# unique_names_SILVA

# solving the different names where possible (through gsub), or adding from the species level output in the previous script...
Counts_Kraken$Genus <- gsub("Nostocoides","Tetrasphaera", Counts_Kraken$Genus ) 
Counts_Kraken$Genus <- gsub("Candidatus_Moraniibacteriota bacterium","Candidatus_Moranbacteria", Counts_Kraken$Genus )
Counts_Kraken$Genus <- gsub("Candidatus_Kaiserbacteria bacterium","Candidatus_Kaiserbacteria", Counts_Kraken$Genus )
Counts_Kraken$Genus <- gsub("Candidatus_Nomurabacteria bacterium","Candidatus_Nomurabacteria", Counts_Kraken$Genus )
# the SILVA genus SH-PL14 has an NCBI taxID 1632864 according to MIDAS, which corresponds to the "missing" genus of the species "Planctomyces sp. SH-PL14"
Counts_Kraken$Genus <- gsub("Planctomyces sp. SH-PL14","SH-PL14", Counts_Kraken$Genus )
# the SILVA genus IMCC26207 has the NCBI code "1641811" according to MIDAS --> corresponds to the "missing" genus of the species Actinobacteria bacterium IMCC26207
# the SILVA genus 67-14 corresponds to the "missing" genus of the species Solirubrobacterales bacterium 67-14 on NBCI
Counts_Kraken$Genus <- gsub("Solirubrobacterales bacterium","67-14", Counts_Kraken$Genus )
# the SILVA genus 37-13 corresponds to the "missing" genus of the species Bacteroidetes bacterium 37-13
# Bdellovibrio, Lentimicrobium and Prosthecobacter are actually also in NCBI, they are just not featured in the Kraken dataset 
# "Pir4 lineage" (of Pirellulaceae) is feautured in both SILVA and MIDAS but it is not in NCBI at all!
# "SM1A02" does not exist on NCBI either...
# codes such as "1-20"  are not found in NCBI...
Counts_Kraken$Genus <- gsub("Gracilibacteria","Altimarinota", Counts_Kraken$Genus ) 




# solving the different names where possible (through gsub), in Kaiju

kaiju_wide$taxon_name <- gsub("Candidatus_Moranbacteria bacterium","Candidatus_Moranbacteria", kaiju_wide$taxon_name )
kaiju_wide$taxon_name <- gsub("Candidatus_Kaiserbacteria bacterium","Candidatus_Kaiserbacteria", kaiju_wide$taxon_name )
kaiju_wide$taxon_name <- gsub("Candidatus_Nomurabacteria bacterium","Candidatus_Nomurabacteria", kaiju_wide$taxon_name )
kaiju_wide$taxon_name <- gsub("Planctomyces sp. SH-PL14","SH-PL14", kaiju_wide$taxon_name )
kaiju_wide$taxon_name <- gsub("Solirubrobacterales bacterium","67-14", kaiju_wide$taxon_name )














################ DELETION OF EUCARIOTES, VIRUSES AND ARCHAEA + FILTERING IS PERFORMED AFTERWARDS




### New colnames
codes<- c("X1510023F1498733"="16Strue_Lab", "X1621223F1607793"="16Strue_T1",
          "X1621235F1607805"="16Strue_T2", "X1623017F1607817"="16Strue_T3",
          "X1623038F1607838"="16Strue_T4")
colnames(Counts_true_16)<- c("Genus", codes [ colnames(Counts_true_16)[colnames(Counts_true_16)!="Genus"] ] )

codes<- c("1623909_S50.tsv"="Kraken_Lab", "1623905_S46.tsv"="Kraken_T1",
          "1623906_S47.tsv"="Kraken_T2", "1623907_S48.tsv"="Kraken_T3",
          "1623908_S49.tsv"="Kraken_T4",
          "mock_AGS_DNAseq.tsv"="Kraken_Mock")
colnames(Counts_Kraken)<- c("Genus", codes [ colnames(Counts_Kraken)[colnames(Counts_Kraken)!="Genus"] ] )

codes<- c("1623909_S50"="RiboFull_Lab",
          "mock_AGS_DNAseq"="RiboFull_Mock",  
          "1623905_S46"="RiboFull_T1",
          "1623906_S47"="RiboFull_T2", "1623907_S48"="RiboFull_T3",
          "1623908_S49"="RiboFull_T4"
)
colnames(Counts_RiboF_full) <- c( "Genus", codes [ colnames(Counts_RiboF_full)[colnames(Counts_RiboF_full)!="Genus"] ] )

codes<- c("1623909_S50"="RiboV3V4_Lab",
          "mock_AGS_DNAseq"="RiboV3V4_Mock",
          "1623905_S46"="RiboV3V4_T1",
          "1623906_S47"="RiboV3V4_T2", "1623907_S48"="RiboV3V4_T3",
          "1623908_S49"="RiboV3V4_T4"
)
colnames(Counts_RiboF_V3V4) <- c( "Genus", codes [ colnames(Counts_RiboF_V3V4)[colnames(Counts_RiboF_V3V4)!="Genus"] ] )

codes<- c("1623909_S50_Kaiju"="Kaiju_Lab",
          "mock_AGS_Kaiju"="Kaiju_Mock",
          "1623905_S46_Kaiju"="Kaiju_T1",
          "1623906_S47_Kaiju"="Kaiju_T2", "1623907_S48_Kaiju"="Kaiju_T3",
          "1623908_S49_Kaiju"="Kaiju_T4"
)
colnames(kaiju_wide) <- c( "taxon_name", codes [ colnames(kaiju_wide)[colnames(kaiju_wide)!="taxon_name"] ] )








### parsing between programs (while adding the known "true" abundances of the mock)

# importing the mock abundances
table_mock <- read.delim(file="true_abundances_of_the_DNAmock.tsv", sep=" ", col.names= c("Genus", "Known_mock") , header= F )

#fixing colname for Kaiju table before parsing
colnames(kaiju_wide)[colnames(kaiju_wide) == "taxon_name"] <- "Genus"


files_list<- list( Counts_Kraken, Counts_RiboF_full, Counts_RiboF_V3V4, Counts_true_16 , kaiju_wide, table_mock )
names(files_list) <-  c("Counts_Kraken", "Counts_RiboF_full", "Counts_RiboF_V3V4", "Counts_true_16" , "kaiju_wide", "table_mock" )
every_entry<-NULL
for (x in 1:length(files_list)) {
  every_entry <- c(every_entry, files_list[[x]][["Genus"]])
}
# from list (each set_pipeline is repeated along a column) to feature table (each set_pipeline is a column)
feature_table<-cbind("Temp_column"=unique(every_entry))
# each bacterium in a row --> every bacterium in this table --> searching across datasets
for( s in 1:length(files_list) ) {
  set_pipeline<- names(files_list)[[s]]  # 1 set_pipeline at time
  table_4_that_set_pipeline<- files_list[[s]]
  abundances<-NULL  # resets the object, ready for the next bacterium
  # an empty column to be filled with each sample...
  temp_set_pipeline<- matrix( "temp", length(unique(every_entry)) , length(colnames(table_4_that_set_pipeline)) )
  temp_set_pipeline <- as.data.frame( temp_set_pipeline )
  temp_set_pipeline[ ,1] <- unique(every_entry)
  for( b in unique(every_entry)) {
    # scan each bacterium in order: if found in that set_pipeline then appends the correspective abundances, if absent paste a 0
    if( b %in% table_4_that_set_pipeline$Genus ){
      abundances<- as.numeric(table_4_that_set_pipeline[ table_4_that_set_pipeline$Genus==b, !colnames(table_4_that_set_pipeline)%in%"Genus"  ])   # NB: the entire ROW will be captured this time
    } else {
      abundances<- rep(0, length( colnames(table_4_that_set_pipeline)[!colnames(table_4_that_set_pipeline)%in%"Genus" ]) )  # that bacterium is absent --> a series of zero
    }
    temp_set_pipeline[ temp_set_pipeline$V1 == b , ] <- c(b , abundances)  # overwrite the "temp"s row with actual values
  }
  colnames(temp_set_pipeline)<- colnames(table_4_that_set_pipeline)
  feature_table <-cbind.data.frame(feature_table, temp_set_pipeline)  # this adds the sample to the complete table
}


feature_table <- feature_table[ , ! colnames(feature_table) %in% "Genus" ]
beasts <- feature_table$Temp_column  # this column correspond to another "Genus" column
# NB: "apply" will delete the colnames, then the genera names will be assigned as row.names afterwards
feature_table$Temp_column <- NULL
feature_table<-as.data.frame(apply(feature_table, MARGIN=2, FUN= as.numeric))
row.names(feature_table) <- beasts


#filtering the genera with low abundances
filtered_taxa <- read.table(file="Counts_from_every_pipeline.tsv", sep="\t", row.names = 1) 
raw_counts <- feature_table[ row.names(feature_table) %in% row.names(filtered_taxa) , ]

#rounding the raw counts from riboFrame
raw_counts[ , grepl("Ribo", names(raw_counts)) ] <- round(raw_counts[ , grepl("Ribo", names(raw_counts)) ])



write.table(raw_counts, file="Raw_counts_from_every_pipeline.tsv", sep = "\t", quote=F, row.names = T, col.names = T)


