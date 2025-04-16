

############################# STARTING WITH THE ANALYSIS ##########################

# launching the necessary packages
library(ggplot2)
library(vegan) 
library(ecodist) 
library(ggvenn)
library(RColorBrewer)
library(reshape2) 
library(ggrepel)
library(dplyr)
library(gridExtra)

# Colors for the stacked bar plots 
fill_color_15 <- c("brown3","springgreen2","wheat","darkmagenta","coral",
                   "yellow3","magenta","pink3", "blue2","darkslategray3","darkgreen","gold",
                   "gray","chartreuse3","violet", "deepskyblue2")

fill_color_10 <- c("springgreen2","orange", "darkmagenta","coral","yellow2","blue2","firebrick3","chartreuse4","gray","deepskyblue2","violet")







##################################           MOCK   COMMUNITY     ##############################

######################                identified mock GENERA         ######################

# loading the dataset
feature_table<- read.table(file="Counts_from_every_pipeline.tsv", sep="\t", row.names = 1)
mock_true <- feature_table[feature_table$Known_mock != 0 , "Known_mock", drop= FALSE]

# preparing the dataset to plot
mock_table <- feature_table[row.names(feature_table) %in% row.names(mock_true) , grepl("ock", names(feature_table))] # data frame with ALL the mocks
mock_table$genera <- row.names(mock_table)
others_mock_table <- feature_table[!(row.names(feature_table) %in% row.names(mock_true)) , grepl("ock", names(feature_table)) , ]   
others_mock_table$genera <- "others"   # dataframe of "others"

# saving the mock abundances 
write.table(mock_table, file = "mock_abundances_top", sep="/t", quote = FALSE)

# input dataframe to plot
mock_to_plot <- melt(mock_table, id.vars = "genera", variable.name = "Type", value.name = "Relative_abund")
others_mock_to_plot <-melt(others_mock_table, id.vars= "genera", variable.name = "Type", value.name = "Relative_abund")
table_mock_to_plot <- rbind(mock_to_plot, others_mock_to_plot)

# transforming into factor to make sure that "others" is at the bottom of the bar plot
table_mock_to_plot$genera <- factor(table_mock_to_plot$genera, levels = c(unique(table_mock_to_plot$genera)[!unique(table_mock_to_plot$genera) %in% "Others"], "Others"))

# transforming into factor to use facet in ggplot2
table_mock_to_plot$Type <- factor(table_mock_to_plot$Type, levels = c("Known_mock", "Kraken_Mock", "RiboFull_Mock", "RiboV3V4_Mock", "Kaiju_Mock"))

# editing the labels
table_mock_to_plot$Type <- gsub("_Mock", "", table_mock_to_plot$Type)
table_mock_to_plot$Type <- gsub("Full", "Frame Full", table_mock_to_plot$Type)
table_mock_to_plot$Type <- gsub("V3V4", "Frame V3-V4", table_mock_to_plot$Type)
table_mock_to_plot$Type <- gsub("Known_mock", "Mock", table_mock_to_plot$Type)

# plotting
ggplot(data = table_mock_to_plot, aes(x = Type, y = Relative_abund, fill = genera)) +
  geom_bar(stat = "identity", position = "stack") + 
  theme_classic(base_size = 8.5) +
  theme_set(theme_minimal())  +
  scale_fill_manual(values = fill_color_10) +
  scale_y_continuous(expand = c(0, 1)) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.x = element_text(angle = 42, vjust = 1, hjust = 1, size = 10),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14, vjust = 1, margin = margin(t = 10)),  
    strip.text = element_text(size = 14),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width = unit(0.38, "cm"),
    legend.text = element_text(size = 10),
    legend.position = "bottom",                                
    legend.margin = margin(9, 28, 0, 1),
    plot.margin = margin(4, 1, 4, 1)
  ) +
  guides(fill = guide_legend(nrow = 3)) +
  labs(
    x = "Processing methods", 
    y = "Percentual abundance",
    title = "Microbial composition of the mock \nat genera level",
    fill = "",
    caption = "'Others' includes every genus not present in the mock sequencing"
  )

dev.off()


####################                    REAL TOP 10               ######################

# top 10 Kraken2
mock_table_Kraken <- feature_table[, grepl("Kraken_Mock", names(feature_table)), drop = FALSE]
top_Kraken <- mock_table_Kraken[order(mock_table_Kraken[[1]], decreasing = TRUE), , drop = FALSE][1:10, , drop = FALSE]

# top 10 Kaiju
mock_table_Kaiju <- feature_table[, grepl("Kaiju_Mock", names(feature_table)), drop = FALSE]
top_Kaiju <- mock_table_Kaiju[order(mock_table_Kaiju[[1]], decreasing = TRUE), , drop = FALSE][1:10, , drop = FALSE]

# top 10 riboFull
mock_table_riboFull <- feature_table[, grepl("RiboFull_Mock", names(feature_table)), drop = FALSE]
top_riboFull <- mock_table_riboFull[order(mock_table_riboFull[[1]], decreasing = TRUE), , drop = FALSE][1:10, , drop = FALSE]

# top 10 riboV3V4
mock_table_riboV3V4 <- feature_table[, grepl("RiboV3V4_Mock", names(feature_table)), drop = FALSE]
top_riboV3V4 <- mock_table_riboV3V4[order(mock_table_riboV3V4[[1]], decreasing = TRUE), , drop = FALSE][1:10, , drop = FALSE]


top_10 <- data.frame(Kaiju = rownames(top_Kaiju), Kraken = rownames(top_Kraken), riboFull = rownames(top_riboFull),riboV3V4 = rownames(top_riboV3V4))

write.table(top_10, file = "real_top_10", sep = "/t", quote = FALSE)


######################                   VENN DIAGRAM mock        ######################

# THREE-SET Venn diagram 
venn_data <- list("Mock" = row.names(mock_true),      # genera chr vectors 
                  "Kraken" =  row.names(feature_table[feature_table$Kraken_Mock != 0 , "Kraken_Mock", drop= FALSE]) , 
                  "Kaiju" = row.names(feature_table[ feature_table$Kaiju_Mock != 0 , "Kaiju_Mock" ,drop = FALSE ]))
# plotting
ggvenn(venn_data,
       fill_color = brewer.pal(3, "Set2"),
       stroke_size = 0.5,
       set_name_size = 6,  
       label_sep = "\n",
       show_percentage = FALSE,
       text_size = 6) +  
  ggtitle("Venn diagram of the mock community") +
  theme(plot.title = element_text(vjust=4, hjust = 0.5, size = 18, face = "bold"),
        text = element_text(size = 15),
        plot.margin = margin(2, 2, 2, 2, "cm")  
  )




# THREE-SET Venn diagram 
venn_data <- list("Mock" = row.names(mock_true),      # genera chr vectors 
                  "Riboframe Full " = row.names(feature_table[feature_table$RiboFull_Mock != 0 , "RiboFull_Mock", drop= FALSE]) , 
                  "Riboframe V3-V4" = row.names(feature_table[ feature_table$RiboV3V4_Mock != 0 , "RiboV3V4_Mock" ,drop = FALSE ]))


# plotting
ggvenn(venn_data,
       fill_color = brewer.pal(3, "Set2"),
       stroke_size = 0.5,
       set_name_size = 6,  
       label_sep = "\n",
       show_percentage = FALSE,
       text_size = 6) +  
  ggtitle("Venn diagram of the mock community") +
  theme(plot.title = element_text(vjust=4, hjust = 0.5, size = 18, face = "bold"),
        text = element_text(size = 15),
        plot.margin = margin(2, 2, 2, 2, "cm")  
  )





###################     ALPHA DIVERSITY MOCK      ##############


#################        SHANNON & EVENNESS      ######

# loading dataset
feature_table<- read.table(file="Raw_counts_from_every_pipeline.tsv", sep="\t", row.names = 1)
table_mock_to_plot <- feature_table[, (grepl("ock", names(feature_table) ))]

# defining the function to calculate the shannon index
shannon_index <- function(x) {
  p <- x / sum(x)  
  p <- p[p > 0]   
  -sum(p * log(p)) 
}

# shannon table to plot
shan_mock_to_plot <- data.frame("Shannon_Index" = apply(table_mock_to_plot, 2, shannon_index), "Sample" = names(table_mock_to_plot))

# editing the labels
shan_mock_to_plot$Sample <- gsub("_Mock", "", shan_mock_to_plot$Sample)
shan_mock_to_plot$Sample <- gsub("Full", "Frame \nFull", shan_mock_to_plot$Sample)
shan_mock_to_plot$Sample <- gsub("V3V4", "Frame \nV3-V4", shan_mock_to_plot$Sample)
shan_mock_to_plot$Sample <- gsub("Known_mock", "Mock", shan_mock_to_plot$Sample)

# saving the results
write.table(shan_mock_to_plot, file = "shannon_mock", sep = "/t", quote = FALSE)

# defining the function to calculate the evenness index
evenness_index <- function(x) {
  H <- shannon_index(x)  
  S <- length(x[x > 0]) 
  H_max <- log(S)  
  E <- H / H_max 
  return(E)
}

# evenness table to plot
even_table_to_plot <- data.frame("Evenness" = apply(table_mock_to_plot, 2, evenness_index), "Sample" = names(table_mock_to_plot))

# editing the labels
even_table_to_plot$Sample <- gsub( "_Mock", ""  , even_table_to_plot$Sample)
even_table_to_plot$Sample <- gsub("Full", "Frame \nFull", even_table_to_plot$Sample)
even_table_to_plot$Sample <- gsub("V3V4", "Frame \nV3-V4", even_table_to_plot$Sample)
even_table_to_plot$Sample <- gsub("Known_mock", "Mock", even_table_to_plot$Sample)

# saving the results
write.table(even_table_to_plot, file = "eveness_mock", sep = "/t", quote = FALSE)


# plotting the Shannon Index
shannon_plot <-ggplot(shan_mock_to_plot, aes(x = Sample, y = Shannon_Index)) +
  geom_boxplot(outlier.shape = NA) +   
  geom_jitter(width = 0.2, size = 2, aes(color = Sample)) +  
  theme_minimal() +
  labs(title = "Shannon Index at genera level",
       x = "",  
       y = "Shannon Index") +  
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),  
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# plotting the Evenness Index
evenness_plot <-ggplot(even_table_to_plot, aes(x = Sample, y = Evenness)) +
  geom_boxplot(outlier.shape = NA) +   
  geom_jitter(width = 0.2, size = 2, aes(color = Sample)) +  
  theme_minimal() +
  labs(title = "Evenness Index at genera level",
       x = "",  
       y = "Evenness Index") +  
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),  
    legend.title = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# to couple the images
grid.arrange(shannon_plot, evenness_plot, ncol = 2, widths = c(4, 4), heights = c(2, 2))  

dev.off()



######################              BETA DIVERSITY mock        ######################

# dataset preparation
feature_table<- read.table(file="Counts_from_every_pipeline.tsv", sep="\t", row.names = 1)
feature_table <- feature_table[, (grepl("ock", names(feature_table) ))]

names(feature_table) <- gsub( "_Mock", ""  , names(feature_table))
names(feature_table) <- gsub("Full", "Frame Full", names(feature_table))
names(feature_table) <- gsub("V3V4", "Frame V3-V4", names(feature_table))
names(feature_table) <- gsub("Known_mock", "Mock", names(feature_table))


# PCoA analysis
this <- t(feature_table)
this <- sqrt(this) 
dist <- vegdist(this, method = "euclidean")
pcoa <- cmdscale(dist, eig = TRUE)

# plotting
# defining groups and levels
group <- factor(c("Kraken", "RiboFrame Full", "RiboFrame V3-V4", "Kaiju" , "Mock"))


# defining colors
group_colors <- c("Kaiju" = "purple",
                  "Kraken" = "red", 
                  "Mock" = "darkgreen",
                  "RiboFrame Full" = "blue", 
                  "RiboFrame V3-V4" = "orange" 
)


# PCoA site scores extraction 
site_scores <- scores(pcoa, display = "sites")

# adjusting labels
par(mgp = c(2.7, 1.5, 0))  
# basic plot
ordiplot(pcoa, display = 'sites', type = 'n',
         xlab = paste("PC1 (", round(pcoa$eig[1]/sum(pcoa$eig)*100, 2), "%)", sep=""),
         ylab = paste("PC2 (", round(pcoa$eig[2]/sum(pcoa$eig)*100, 2), "%)", sep=""), cex.lab = 1.5)

# adding title and legend
title(main = "Beta diversity using Hellinger distance \nat genera level", font.main = 2, cex.main = 1.5)
legend("bottomleft", legend = levels(group), col = group_colors, 
       pch = c(16, 16, 18, 16),  
       title = "Processing method", pt.cex = 1.5)

# adding points
points(site_scores[group != "Mock", 1], site_scores[group != "Mock", 2], 
       col = group_colors[group[group != "Mock"]], pch = 19, cex = 1.5)  

points(site_scores[group == "Mock", 1], site_scores[group == "Mock", 2], 
       col = group_colors["Mock"], pch = 18, cex = 2.5) 









########################          WASTEWATER SAMPLES        #############################




#########################          TOP 15 GENERA            ################

# loading the dataset 
feature_table<- read.table(file="Counts_from_every_pipeline.tsv", sep="\t", row.names = 1)

# preparing the dataset to plot
feature_table <- feature_table[, !grepl("ock", names(feature_table))]   # remove the mock
feature_table$sum <- apply(feature_table, 1, sum)  # total abundance across samples for each taxon
top_feature_table <- feature_table[order(feature_table$sum, decreasing = TRUE), ][1:15, ]   # data frame with 15 top
top_feature_table$sum <- NULL
top_feature_table$genera <- row.names(top_feature_table)                                      
others_feature_table <- feature_table[!(row.names(feature_table) %in% row.names(top_feature_table)) , ]   # data frame with "others"
others_feature_table$sum <- NULL
others_feature_table$genera <- "others"

# input dataframe to plot
top_to_plot <- melt(top_feature_table, id.vars = "genera", variable.name = "Type", value.name = "Relative_abund")
others_to_plot <-melt(others_feature_table, id.vars= "genera", variable.name = "Type", value.name = "Relative_abund")
table_to_plot <- rbind(top_to_plot, others_to_plot)

# transforming into factor to make sure that "others" is at the bottom of the bar plot
table_to_plot$genera <- factor(table_to_plot$genera, levels = c(unique(table_to_plot$genera)[!unique(table_to_plot$genera) %in% "Others"], "Others"))

# transforming into factor to use facet in ggplot2
table_to_plot$Type <- factor(table_to_plot$Type, levels = unique(table_to_plot$Type))


#################################   STACKED BAR PLOT TOP GENERA   #####################

# editing the labels
table_to_plot$Sample <- table_to_plot$Type
table_to_plot$Sample <- gsub(".*_", "", table_to_plot$Sample)
table_to_plot$Type <- gsub("_.*", "", table_to_plot$Type)
table_to_plot$Type <- gsub("Full", "Frame \nFull", table_to_plot$Type)
table_to_plot$Type <- gsub("V3V4", "Frame \nV3-V4", table_to_plot$Type)
table_to_plot$Type <- gsub("X16Strue", "16S amplicon", table_to_plot$Type)

# plotting
ggplot(data = table_to_plot, aes(x = Sample, y = Relative_abund, fill = genera)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(~ Type, scales = "free_x", space = "free_x") +
  theme_classic(base_size = 9) +
  # theme_set(theme_minimal())  +
  scale_fill_manual(values = fill_color_15) +
  scale_y_continuous(expand = c(0, 1)) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.x = element_text(angle = 42, vjust = 1, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14, margin = margin(t = 10)),  
    strip.text = element_text(size = 14),
    legend.key.height = unit(0.35, "cm"),
    legend.key.width = unit(0.28, "cm"),
    legend.text = element_text(size = 8.3),
    legend.position = "bottom",                                 
    legend.margin = margin(9, 38, 0, 1),
    plot.margin = margin(4, 1, 4, 1)
  ) +
  guides(fill = guide_legend(nrow = 4)) +
  labs(
    x = "Samples", 
    y = "Percentual abundance",
    title = "Microbial composition at genera level",
    fill = "",
    caption = "'Others' includes every genus below rank 15"
  )

dev.off()



############################       ALPHA DIVERSITY        ##############################



############################        RICHNESS

# preparing the dataframe to plot
feature_table <- read.table(file="Raw_counts_from_every_pipeline.tsv", sep="\t", row.names = 1)
feature_table <- feature_table[ , !grepl( "ock", names(feature_table)) ]   # remove mocks
feature_table_richness <- colSums(feature_table > 0)                       # numeric vector with richness
richness_to_plot <- data.frame("Type" = names(feature_table_richness), "Richness" = feature_table_richness, "Sample"= names(feature_table_richness))

# editing the labels
richness_to_plot$Sample <- gsub( ".*_", ""  , richness_to_plot$Sample)
richness_to_plot$Type <- gsub( "_.*", ""  , richness_to_plot$Type)
richness_to_plot$Type <- gsub("Full", "Frame Full", richness_to_plot$Type)
richness_to_plot$Type <- gsub("V3V4", "Frame V3-V4", richness_to_plot$Type)
richness_to_plot$Type <- gsub("X16Strue", "16S amplicon", richness_to_plot$Type)

# defining values for error bars
rich_summary <- richness_to_plot %>%
  group_by(Type) %>%
  summarise(
    mean_Richness = mean(Richness),  # mean
    sd_Richness = sd(Richness),  # standard deviation
    n = n(),  # number of samples
    se = sd_Richness / sqrt(n)  # Standard error
  )

# saving the results
write.table(rich_summary, file= "rich_summary_wws", sep = "/t", quote = FALSE)

# plotting
ggplot(richness_to_plot, aes(x = Type, y = Richness)) +
  geom_boxplot(outlier.shape = NA) +   
  geom_jitter(width = 0.2, size = 3, aes(color = Type)) +  
  geom_text_repel(aes(label = Sample), size = 4, color = "black", max.overlaps = Inf) +          
  geom_errorbar(data = rich_summary, 
                aes(x = Type, y = mean_Richness, ymin = mean_Richness - se, ymax = mean_Richness + se), 
                width = 0.2, color = "black") +
  geom_point(data = rich_summary, aes(x = Type, y = mean_Richness), size = 3, color = "black") + 
  scale_y_continuous(breaks = seq(0, max(richness_to_plot$Richness), by = 100)) +
  theme_minimal() +
  labs(x = "Processing method",
       y = "Richness",
       title = "Observed richness at genera level",
       color = "") +  
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),  
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  scale_color_manual(values = c(
    "Kaiju" = "purple",
    "Kraken" = "firebrick3",           
    "RiboFrame Full" = "blue2",   
    "RiboFrame V3-V4" = "gold",  
    "16S amplicon" = "springgreen2"          
  ))

#ggsave(file = "Rplots/richness_boxplot.png", width = 7, height = 5.4, dpi = 300)
dev.off()

############################        SHANNON  

# loading dataset
raw_counts <- read.table(file="Raw_counts_from_every_pipeline.tsv", sep="\t", row.names = 1)
raw_counts <- raw_counts[ , !(grepl("ock", names(raw_counts)))]

# using the function previously defined
# table to plot
shan_table_to_plot <- data.frame("Type" = names(raw_counts), "Shannon_Index" = apply(raw_counts, 2, shannon_index), "Sample" = names(raw_counts))

# editing the labels
shan_table_to_plot$Sample <- gsub( ".*_", ""  , shan_table_to_plot$Sample)
shan_table_to_plot$Type <- gsub( "_.*", ""  , shan_table_to_plot$Type)
shan_table_to_plot$Type <- gsub("Full", "Frame Full", shan_table_to_plot$Type)
shan_table_to_plot$Type <- gsub("V3V4", "Frame V3-V4", shan_table_to_plot$Type)
shan_table_to_plot$Type <- gsub("X16Strue", "16S amplicon", shan_table_to_plot$Type)

# defining values for error bars
shan_summary <- shan_table_to_plot %>%
  group_by(Type) %>%
  summarise(
    mean_Shannon = mean(Shannon_Index),
    sd_Shannon = sd(Shannon_Index),
    n = n(),  
    se = sd_Shannon / sqrt(n)  
  )

# saving the results
write.table(shan_summary, file= "shan_summary_wws", sep = "/t", quote = FALSE)

# plotting
ggplot(shan_table_to_plot, aes(x = Type, y = Shannon_Index)) +
  geom_boxplot(outlier.shape = NA) +   
  geom_jitter(width = 0.2, size = 3, aes(color = Type)) +  
  geom_text_repel(aes(label = Sample), size = 4, color = "black", max.overlaps = Inf) +          
  geom_errorbar(data = shan_summary, 
                aes(x = Type, y = mean_Shannon, ymin = mean_Shannon - se, ymax = mean_Shannon + se), 
                width = 0.2, color = "black") +
  geom_point(data = shan_summary, aes(x = Type, y = mean_Shannon), size = 3, color = "black") + 
  scale_y_continuous(breaks = seq(0, max(shan_table_to_plot$Shannon_Index), by = 0.5)) +
  theme_minimal() +
  labs(x = "Processing method",
       y = "Shannon Index",
       title = "Shannon Index at genera level",
       color = "") +  
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),  
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  ) +
  scale_color_manual(values = c(
    "Kaiju" = "purple",
    "Kraken" = "#FF00FF",           
    "RiboFrame Full" = "#1E90FF",   
    "RiboFrame V3-V4" = "#00FA9A",  
    "16S amplicon" = "#DAA520"     
  ))

dev.off()


#########################        EVENNESS

# using the function previously defined
# table to plot
even_table_to_plot <- data.frame("Type" = names(raw_counts), "Eveness" = apply(raw_counts, 2, evenness_index), "Sample" = names(raw_counts))

# editing the labels
even_table_to_plot$Sample <- gsub( ".*_", ""  , even_table_to_plot$Sample)
even_table_to_plot$Type <- gsub( "_.*", ""  , even_table_to_plot$Type)
even_table_to_plot$Type <- gsub("Full", "Frame Full", even_table_to_plot$Type)
even_table_to_plot$Type <- gsub("V3V4", "Frame V3-V4", even_table_to_plot$Type)
even_table_to_plot$Type <- gsub("X16Strue", "16S amplicon", even_table_to_plot$Type)

# defining values for error bars
even_summary <- even_table_to_plot %>%
  group_by(Type) %>%
  summarise(
    mean_Eveness = mean(Eveness),
    sd_Eveness = sd(Eveness),
    n = n(),  
    se = sd_Eveness / sqrt(n)  
  )

# saving the results
write.table(even_summary, file= "even_summary_wws", sep = "/t", quote = FALSE)

# plotting
ggplot(even_table_to_plot, aes(x = Type, y = Eveness)) +
  geom_boxplot(outlier.shape = NA) +   
  geom_jitter(width = 0.2, size = 3, aes(color = Type)) +  
  geom_text_repel(aes(label = Sample), size = 4, color = "black", max.overlaps = Inf) +          
  geom_errorbar(data = even_summary, 
                aes(x = Type, y = mean_Eveness, ymin = mean_Eveness - se, ymax = mean_Eveness + se), 
                width = 0.2, color = "black") +
  geom_point(data = even_summary, aes(x = Type, y = mean_Eveness), size = 3, color = "black") + 
  scale_y_continuous(breaks = seq(0, max(shan_table_to_plot$Shannon_Index), by = 0.1)) +
  theme_minimal() +
  labs(x = "Processing method",
       y = "Evenness Index",
       title = "Evenness Index at genera level",
       color = "") +  
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18),  
    legend.title = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)  
  ) +
  scale_color_manual(values = c(
    "Kaiju" = "orange",
    "Kraken" = "green",           
    "RiboFrame Full" = "darkred",   
    "RiboFrame V3-V4" = "darkblue",  
    "16S amplicon" = "purple"      
  )) 

dev.off()









########################         BETA DIVERSITY        ##########################




####################      ONLY WASTEWATER SAMPLES  


# dataset preparation
feature_table<- read.table(file="Counts_from_every_pipeline.tsv", sep="\t", row.names = 1)
feature_table <- feature_table[, !(grepl("ock", names(feature_table) ))]

# labels editing
names(feature_table) <- sapply(names(feature_table), function(x) {
  x <- gsub("RiboFull_", "", x)  
  x <- gsub("RiboV3V4_", "", x)
  x <- gsub("X16Strue_", "", x)   
  x <- gsub("Kraken_", "", x) 
  x <- gsub("Kaiju_", "", x)
  return(x)  
})

# PCoA analysis
this <- t(feature_table)
this <- sqrt(this) 
dist <- vegdist(this, method = "euclidean")
pcoa <- cmdscale(dist, eig = TRUE)

# plotting
# defining groups and levels
group <- factor(c(rep("Kraken", 5),
                  rep("RiboFrame Full", 5),
                  rep("RiboFrame V3-V4", 5),
                  rep("Kaiju", 5),
                  rep("16S amplicon", 5)),  
                levels = c("Kraken", "RiboFrame Full", "RiboFrame V3-V4", "Kaiju", "16S amplicon"))

# defining colors
group_colors <- c("Kraken" = "red", 
                  "RiboFrame Full" = "blue", 
                  "RiboFrame V3-V4" = "orange",
                  "Kaiju" = "purple",
                  "16S amplicon" = "green")

# PCoA site scores extraction 
site_scores <- scores(pcoa, display = "sites")

# adjusting margins and labels
par(mar = c(5, 5, 3, 1), mgp = c(3, 1.5, 0))

# basic plot
ordiplot(pcoa, display = 'sites', type = 'n',
         xlab = paste("PC1 (", round(pcoa$eig[1]/sum(pcoa$eig)*100, 2), "%)", sep=""),
         ylab = paste("PC2 (", round(pcoa$eig[2]/sum(pcoa$eig)*100, 2), "%)", sep=""), cex.lab = 1.5
)

# adding points
points(site_scores[,1], site_scores[,2], col = group_colors[group], pch = 19, cex = 1.5)

# adding title and legend
title(main = "Beta diversity using Hellinger distance \nat genera level", font.main = 2, cex.main = 1.5)
legend("bottomleft", legend = levels(group), col = group_colors, pch = 16, title = "Processing method", pt.cex = 1.5)

# adding text
text(site_scores[,1], site_scores[,2] - 0.05,  
     labels = names(feature_table),               
     pos = 2, , cex = 1)

# adding confidence ellipses
for (grp in levels(group)) {
  grp_idx <- group == grp
  grp_scores <- site_scores[grp_idx, ]
  if (is.null(dim(grp_scores))) {
    grp_scores <- matrix(grp_scores, nrow = 1)
  }
  ordiellipse(grp_scores, groups = rep(grp, sum(grp_idx)), 
              kind = "sd", draw = "polygon", 
              col = adjustcolor(group_colors[grp], alpha.f = 0.2), 
              border = group_colors[grp])
}



dev.off()




################          ONLY AGS COMMON

feature_table<- read.table(file="Counts_from_every_pipeline.tsv", sep="\t", row.names = 1)   
feature_table <- feature_table[, !(grepl("ock", names(feature_table) ))]

# extract target columns
feature_table_16S <- feature_table[, grepl("16S", names(feature_table))]                     
feature_table_Kraken <- feature_table[, grepl("Kraken", names(feature_table))]
feature_table_Ribo_Full <- feature_table[ , grepl("Full", names(feature_table))]
feature_table_Ribo_V3V4 <- feature_table[, grepl("V3V4", names(feature_table))]
feature_table_Kaiju <- feature_table[, grepl("aiju", names(feature_table))]


# data frames with the only present bacteria for each dataset
names_16S <- row.names(feature_table_16S[ names(which(rowSums(feature_table_16S) > 0)) , ])    
names_Kraken <- row.names(feature_table_Kraken[ names(which(rowSums(feature_table_Kraken) > 0)) , ])
names_Ribo_Full <- row.names(feature_table_Ribo_Full[ names(which(rowSums(feature_table_Ribo_Full) > 0)) , ])
names_Ribo_V3V4 <- row.names(feature_table_Ribo_V3V4[ names(which(rowSums(feature_table_Ribo_V3V4) > 0)) , ])
names_Kaiju <- row.names(feature_table_Kaiju[ names(which(rowSums(feature_table_Kaiju) > 0)) , ])


# common bacteria
common_names <- Reduce(intersect, list(names_16S, names_Kraken, names_Ribo_Full, names_Ribo_V3V4, names_Kaiju))
common_feature_table <- feature_table[ row.names(feature_table) %in% common_names, !(grepl("ock", names(feature_table) ))]

# editing names
names(common_feature_table) <- sapply(names(common_feature_table), function(x) {
  if (grepl("ock", x)) {
    x <- gsub("_", " ", x)  
  } else {
    x <- gsub("RiboFull_", "", x)  
    x <- gsub("RiboV3V4_", "", x)
    x <- gsub("X16Strue_", "", x)   
    x <- gsub("Kraken_", "", x)  
    x <- gsub("Kaiju_", "", x)
    return(x)  
  }
})

# compute pcoa
this <- t(common_feature_table)
this <- sqrt(this) 
dist <- vegdist(this, method = "euclidean")
pcoa <- cmdscale(dist, eig = TRUE)

# plotting
# defining groups and levels
group <- factor(c(rep("Kraken", 5),
                  rep("RiboFrame Full", 5),
                  rep("RiboFrame V3-V4", 5),
                  rep("Kaiju", 5),
                  rep("16S amplicon", 5)),  
                levels = c("Kraken", "RiboFrame Full", "RiboFrame V3-V4", "Kaiju", "16S amplicon"))

# defining colors
group_colors <- c("Kraken" = "red", 
                  "RiboFrame Full" = "blue", 
                  "RiboFrame V3-V4" = "orange",
                  "Kaiju" = "purple",
                  "16S amplicon" = "green")

# PCoA site scores extraction 
site_scores <- scores(pcoa, display = "sites")

# adjusting margins and labels
par(mar = c(5, 5, 3, 1), mgp = c(3, 1.5, 0))

# basic plot
ordiplot(pcoa, display = 'sites', type = 'n',
         xlab = paste("PC1 (", round(pcoa$eig[1]/sum(pcoa$eig)*100, 2), "%)", sep=""),
         ylab = paste("PC2 (", round(pcoa$eig[2]/sum(pcoa$eig)*100, 2), "%)", sep=""), cex.lab = 1.5) 

# adding title and legend
title(main = "Beta diversity using Hellinger distance \non common genera", font.main = 2, cex.main = 1.5)
legend("bottomleft", legend = levels(group), col = group_colors, pch = 16, title = "Processing method", pt.cex = 1.5)

# adding points
points(site_scores[,1], site_scores[,2], col = group_colors[group], pch = 19, cex = 1.5)

# adding labels  
text(site_scores[, 1], site_scores[, 2],
     label = names(common_feature_table),
     pos = 4, cex = 1)

# confidence ellipses
for (grp in levels(group)) {
  # selection of the points corresponding to each group 
  grp_idx <- group == grp
  grp_scores <- site_scores[grp_idx, ]
  
  # if there is only one cpd, make sure that 'grp_scores' is a 2d matrix
  if (is.null(dim(grp_scores))) {
    grp_scores <- matrix(grp_scores, nrow = 1)
  }
  
  # draw the confidence ellipse for each group
  ordiellipse(grp_scores, groups = rep(grp, sum(grp_idx)), 
              kind = "sd", draw = "polygon", 
              col = adjustcolor(group_colors[grp], alpha.f = 0.2), 
              border = group_colors[grp])
}


devoff()






