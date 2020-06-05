#### These libraries are needed to run this script if any are missing you will need to install them
# install.packages("tidyverse", "ggforce", "reshape2", "readxl")
library(tidyverse)
library(ggforce)
library(reshape2)
library(readxl)
library(ggthemes)
##### End libraries needed for this script

##### Variables section - Update these to reflect your files.

pheno_file <- "CAAPA2_Ethiopian_Saliva_Samples_Received_10282019_JWM.xlsx"
rdata_file <- "CAAPA2_Ethiopia_MEGA_012220_updated_callrate_passing_QC_ld_prune_50kb_5_2_pcair_pcrel.Rdata"
# rdata_file <- "CAAPA2_Ethiopia_MEGA_012220_pcair_pcrel.Rdata"
# This is what you want your output files named.
prefix <- "Ethiopia_old_GENESIS_ldpruned"
#prefix <- "Ethiopia"

##### End Variables section 

pheno <- read_excel(path = pheno_file, sheet = "Sheet2")  %>%
  rename(`Gender (M/F/U)` = "Gender (", Ethnicity = "Ethnic group", `Sample ID` = "sample_id") %>%
  drop_na("Sample ID")

load(rdata_file)

### You will need to adjust this section to get the data in the correct format
pcs <- mypcair$vectors %>% as_tibble(rownames = "IID") %>%
  separate(IID, into = c(NA, NA, NA, "Sample ID")) %>% select(1:11)
plot_pca <- inner_join(pheno, pcs, by = "Sample ID") %>%
  select("Sample ID", "Ethnicity", "Language group", num_range("V", 1:10))

# Create and save the scree plot
plot_scree <- tibble(Eigenvalues = mypcair$values, Percent = mypcair$values / mypcair$sum.values * 100)
ggplot(plot_scree, aes(x = as.integer(row.names(plot_scree)))) + 
  geom_bar(aes(y = Percent), stat = "identity", color = plot_scree$Percent > 1.25) +
  xlab("PC Number") + ylab("Percent Explained") + scale_color_ptol()
ggsave(paste(prefix, "screeplot.png", sep = "_"))

# Create and save PC1 vs PC2 plot
ggplot(plot_pca, aes(x = V1, y = V2, shape = `Language group`, color = Ethnicity)) + 
  geom_point(size = 2.5) + scale_color_ptol()
ggsave(paste(prefix, "PC1vPC2.png", sep = "_"))

# Create and save PC plot matrix, PC's 1-5
ggplot(plot_pca, aes(x = .panel_x, y = -.panel_y, shape = `Language group`, color = Ethnicity)) + 
  geom_point(size = 3) + geom_autodensity(aes(fill = Ethnicity)) + geom_density2d() +
  facet_matrix(rows = vars(num_range("V", 1:5)), switch = "y", layer.diag = 2, layer.upper = 3) +
  scale_color_ptol() + scale_fill_ptol()
ggsave(paste(prefix, "PCA_plots_PC1-5.png", sep = "_"), height = 14, width = 14)


# Gather the data for the kinship plot - You may need to adjust this.
plot_kin <- mypcrel$kinship %>% melt() %>% 
  separate("Var1", sep = "_", into = c(NA, NA, "ID1")) %>%
  separate("Var2", sep = "_", into = c(NA, NA, "ID2")) %>%
  rename(Kinship = "value") %>%
  filter(Kinship > 2^(-11/2), ID1 %in% plot_pca$`Sample ID`, ID1 != ID2) %>%
  distinct(ID1, ID2, .keep_all = TRUE)

# Find first degree related samples and write to a file 
rel_1 <- plot_kin %>% filter(Kinship > 0.2) %>% select(-Kinship) %>% as.list() %>% flatten_chr() %>% unique()
non_1_degree_related <- plot_pca %>% filter(! `Sample ID` %in% rel_1) %>% select(`Sample ID`)
write_csv(non_1_degree_related, paste(prefix, "non_1st_Degree_related.csv", sep = "_"))

# Find second degree related samples and write to a file 
rel_2 <- plot_kin %>% filter(Kinship > 0.125 & Kinship < 0.2) %>% select(-Kinship) %>% flatten_chr() %>% unique()
non_1_or_2_degree_related <- plot_pca %>% filter(! `Sample ID` %in% c(rel_1, rel_2)) %>% select(`Sample ID`)
write_csv(non_1_or_2_degree_related, paste(prefix, "non_1st_or_2nd_Degree_related.csv", sep = "_"))

# Create and save Kinship plot
ggplot(plot_kin, aes(x = "kinship", y = Kinship)) +
  geom_boxplot() +
  coord_flip() +
  geom_violin(fill = "grey80",alpha=0.4, position = position_dodge(width = .75),size=1,color="#CCCCCC") +
  ggtitle("Kinship for Ethiopia") +
  geom_hline(yintercept = 0.25, linetype="dashed", color = "blue") +
  annotate(x = 1.5, y = 0.27, geom = "text", label = "1st Degree", color = "blue",) +
  annotate(x = 1.43, y = 0.26, geom = "text", label = length(rel_1), color = "blue") +
  geom_hline(yintercept = 0.125, linetype="dashed", color = "darkgreen") +
  annotate(x = 1.5, y = 0.145, geom = "text", label = "2nd Degree", color = "darkgreen") +
  annotate(x = 1.43, y = 0.135, geom = "text", label = length(rel_2), color = "darkgreen") +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave(paste(prefix, "Ethiopia_kinship_plot.png", sep = "_"))