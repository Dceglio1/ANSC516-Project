#Install the packages, IF YOU NEED TO :)
#install.packages("tidyverse")
#install.packages("vegan")
#install.packages("devtools")
library(devtools)
#devtools::install_github("jbisanz/qiime2R")

#Load the packages. Everyone needs to do this.
library(tidyverse)
library(vegan)
library(qiime2R)
library(ggpubr)

##################################################
#Set working directory to C:/Users/Jvlan/Downloads/AliyaRawData
##################################################

#How to load a file into R
metadata2 <- read.delim("AliyaMetaDataMonthFlush.tsv", sep = "\t", header = T, quote = "", stringsAsFactors = F)

#Now the qiime2R method
metadata<-read_q2metadata("AliyaMetaDataMonthFlush.tsv")
#levels(metadata$`Material`)
#colnames(metadata)[7] <- "Material"
#colnames(metadata)[8] <- "Wall"
#colnames(metadata)[15] <- "Stagnation"
#colnames(metadata)[16] <- "WallLocation"
#colnames(metadata)[17] <- "WallType"
#colnames(metadata)[18] <- "FlushType"
#colnames(metadata)[19] <- "StagType"

row.names(metadata) <- metadata[,1]
row.names(metadata) <- metadata$SampleID
#metadata <- metadata[,-1]

Wall_colors <- c("aquamarine", "coral", "black", "darkorchid", "gold3", "gold3", "royalblue", "royalblue" )

# Read QZA
bc_PCoA <-read_qza("Qiime Outputs/Monthly Flush/bray_curtis_pcoa_results.qza")
wUF_PCoA <- read_qza("Qiime Outputs/Monthly Flush/weighted_unifrac_pcoa_results.qza")
uwUF_PCoA <- read_qza("Qiime Outputs/Monthly Flush/unweighted_unifrac_pcoa_results.qza")
jac_PCoA <- read_qza("Qiime Outputs/Monthly Flush/jaccard_pcoa_results.qza")

### Bray Curtis
bc_meta <- bc_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

bc_centroids <- aggregate(cbind(PC1,PC2)~get("Stagnation"),bc_meta,mean)
colnames(bc_centroids)[1] <- "Stagnation"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get("Stagnation"))) +
  geom_point(aes(shape=Material), size = 3) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  ggtitle("Bray Curtis") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=Wall_colors, name = "Stagnation")
ggsave(paste0("output/Monthly Flush/bc-ellipse", "Stagnation","-subject.pdf"), height=4, width=4.5, device="pdf")

### Weighted UniFrac
wUF_meta <- wUF_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

wUF_centroids <- aggregate(cbind(PC1,PC2)~get("Stagnation"),wUF_meta,mean)
colnames(wUF_centroids)[1] <- "Stagnation"

ggplot(wUF_meta, aes(x=PC1, y=PC2, color=get("Stagnation"))) +
  geom_point(aes(shape= Material), size = 3) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  ggtitle("Weighted UniFrac") +
  xlab(paste0("PC1 (", round(100*wUF_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*wUF_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=Wall_colors, name = "Stagnation")
ggsave(paste0("output/Monthly Flush/wUF-ellipse", "Stagnation","-subject.pdf"), height=4, width=4.5, device="pdf")

### UnWeighted UniFrac
uwUF_meta <- uwUF_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

uwUF_centroids <- aggregate(cbind(PC1,PC2)~get("Stagnation"),bc_meta,mean)
colnames(uwUF_centroids)[1] <- "Stagnation"

ggplot(uwUF_meta, aes(x=PC1, y=PC2, color=get("Stagnation"))) +
  geom_point(aes(shape= Material), size = 3) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  ggtitle("UnWeighted UniFrac") +
  xlab(paste0("PC1 (", round(100*uwUF_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*uwUF_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=Wall_colors, name = "Stagnation")
ggsave(paste0("output/Monthly Flush/uwUF-ellipse", "Stagnation","-subject.pdf"), height=4, width=4.5, device="pdf")

### Jaccard
jac_meta <- jac_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "SampleID"))

jac_centroids <- aggregate(cbind(PC1,PC2)~get("Stagnation"),bc_meta,mean)
colnames(jac_centroids)[1] <- "Stagnation"

ggplot(jac_meta, aes(x=PC1, y=PC2, color=get("Stagnation"))) +
  geom_point(aes(shape= Material), size = 3) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  ggtitle("Jaccard") +
  xlab(paste0("PC1 (", round(100*jac_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*jac_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=Wall_colors, name = "Stagnation")
ggsave(paste0("output/Monthly Flush/jac-ellipse", "Stagnation","-subject.pdf"), height=4, width=4.5, device="pdf")

####################################################
#Alpha Diversity
####################################################

evenness = read_qza("Qiime Outputs/Monthly Flush/evenness_vector.qza")
evenness<-evenness$data %>% rownames_to_column("SampleID") 

observed_features = read_qza("Qiime Outputs/Monthly Flush/observed_features_vector.qza")
observed_features<-observed_features$data %>% rownames_to_column("SampleID") 

shannon = read_qza("Qiime Outputs/Monthly Flush/shannon_vector.qza")
shannon<-shannon$data %>% rownames_to_column("SampleID") 

faith_pd = read_qza("Qiime Outputs/Monthly Flush/faith_pd_vector.qza")
faith_pd<-faith_pd$data %>% rownames_to_column("SampleID")
faith_pd <- faith_pd[,-1]
colnames(faith_pd) <- c('SampleID', 'faith_pd')

alpha_diversity = merge(x=faith_pd, y=evenness, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, observed_features, by.x = "SampleID", by.y = "SampleID")
alpha_diversity = merge(alpha_diversity, shannon, by.x = "SampleID", by.y = "SampleID")
metadata = merge(metadata, alpha_diversity, by.x = "SampleID", by.y = "SampleID")
row.names(metadata) <- metadata$SampleID

###Graphs
hist(metadata$shannon_entropy, main="Shannon diversity", xlab="", breaks=10)
hist(metadata$faith_pd, main="Faith phylogenetic diversity", xlab="", breaks=10)
hist(metadata$pielou_e, main="Evenness", xlab="", breaks=10)
hist(as.numeric(metadata$observed_features), main="Observed Features", xlab="", breaks=10)

ggqqplot(metadata$shannon_entropy, title = "Shannon")
ggqqplot(metadata$faith_pd, title = "Faith PD")
ggqqplot(metadata$pielou_e, title = "Evenness")
ggqqplot(metadata$observed_features, title = "Observed Features")

shapiro.test(metadata$shannon)
shapiro.test(metadata$faith_pd)
shapiro.test(metadata$pielou_e)
shapiro.test(metadata$observed_features)

levels(metadata$StagType)
######################################################
#Evenness
######################################################

metadata <- filter(metadata, !StagType %in% c(''))

evenness_summary <- metadata %>% # the names of the new data frame and the data frame to be summarised
  group_by(StagType) %>%   # the grouping variable
  summarise(mean_evenness = mean(pielou_evenness),  # calculates the mean of each group
            sd_evenness = sd(pielou_evenness), # calculates the standard deviation of each group
            n_evenness = n(),  # calculates the sample size per group
            se_evenness = sd(pielou_evenness)/sqrt(n()))

boxplot(pielou_evenness ~ StagType, data=metadata, ylab="Peilou evenness")

evenness_boxplot <- ggplot(metadata, aes(x = StagType, pielou_evenness)) + 
  geom_boxplot() + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("output/Monthly Flush/evenness_boxplot.pdf", evenness_boxplot, height = 3, width = 3)

evenness_se <- ggplot(evenness_summary, aes(x = StagType, mean_evenness, fill = StagType)) +
  geom_col() + 
  geom_errorbar(aes(ymin = mean_evenness - se_evenness, ymax = mean_evenness + se_evenness), width=0.2) + 
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Pielou's evenness  ± s.e.", x = "") 

ggsave("output/Monthly Flush/evenness_se.pdf", evenness_se, height = 2.5, width = 3)

######################################################
#Faith
######################################################

faith_summary <- metadata %>% # the names of the new data frame and the data frame to be summarised
  group_by(StagType) %>%   # the grouping variable
  summarise(mean_faith = mean(faith_pd),  # calculates the mean of each group
            sd_faith = sd(faith_pd), # calculates the standard deviation of each group
            n_faith = n(),  # calculates the sample size per group
            se_faith = sd(faith_pd)/sqrt(n()))

boxplot(faith_pd ~ StagType, data=metadata, ylab="Faith PD")

faith_boxplot <- ggplot(metadata, aes(StagType, faith_pd)) + 
  geom_boxplot() + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("output/Monthly Flush/faith_boxplot.pdf", faith_boxplot, height = 3, width = 3)

faith_se <- ggplot(faith_summary, aes(StagType, mean_faith, fill = StagType)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = mean_faith - se_faith, ymax = mean_faith + se_faith), width=0.2) + 
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Faith PD  ± s.e.", x = "") 

ggsave("output/Monthly Flush/faith_se.pdf", faith_se, height = 2.5, width = 3)

######################################################
#Observed
######################################################

obs_summary <- metadata %>% # the names of the new data frame and the data frame to be summarised
  group_by(StagType) %>%   # the grouping variable
  summarise(mean_obs = mean(observed_features),  # calculates the mean of each group
            sd_obs = sd(observed_features), # calculates the standard deviation of each group
            n_obs = n(),  # calculates the sample size per group
            se_obs = sd(observed_features)/sqrt(n()))

boxplot(observed_features ~ StagType, data=metadata, ylab="Observed features")

obs_boxplot <- ggplot(metadata, aes(StagType, observed_features)) + 
  geom_boxplot() + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("output/Monthly Flush/obs_boxplot.pdf", obs_boxplot, height = 3, width = 3)

obs_se <- ggplot(obs_summary, aes(StagType, mean_obs, fill = StagType)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = mean_obs - se_obs, ymax = mean_obs + se_obs), width=0.2) + 
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Observed features  ± s.e.", x = "") 

ggsave("output/Monthly Flush/obs_se.pdf", obs_se, height = 2.5, width = 3)

######################################################
#Shannon
######################################################

shannon_summary <- metadata %>% # the names of the new data frame and the data frame to be summarised
  group_by(StagType) %>%   # the grouping variable
  summarise(mean_shannon = mean(shannon_entropy),  # calculates the mean of each group
            sd_shannon = sd(shannon_entropy), # calculates the standard deviation of each group
            n_shannon = n(),  # calculates the sample size per group
            se_shannon = sd(shannon_entropy)/sqrt(n()))

boxplot(shannon_entropy ~ StagType, data=metadata, ylab="Shannon")

shannon_boxplot <- ggplot(metadata, aes(StagType, shannon_entropy)) + 
  geom_boxplot() + 
  #ylim(c(0.5,1)) +
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("output/Monthly Flush/shannon_boxplot.pdf", shannon_boxplot, height = 3, width = 3)

shannon_se <- ggplot(shannon_summary, aes(StagType, mean_shannon, fill = StagType)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = mean_shannon - se_shannon, ymax = mean_shannon + se_shannon), width=0.2) + 
  theme_q2r() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  theme(legend.title = element_blank()) +
  labs(y="Shannon  ± s.e.", x = "") 

ggsave("output/Monthly Flush/shannon_se.pdf", shannon_se, height = 2.5, width = 3)