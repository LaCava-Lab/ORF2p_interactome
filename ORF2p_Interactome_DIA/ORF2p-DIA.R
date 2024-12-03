library(magrittr)
library(dplyr)
library(tidyr)
library(stringr)
library(protti)
library(forcats)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(ggtext)
library(eulerr)

## Spectronaut

sprect_data <- read_protti(filename = "../Data/20240311_011036_34-ORS-A1_Report.csv")
head(sprect_data)

df_test <- sprect_data %>% 
  filter(str_detect(pg_protein_accessions, "Q5T7N2")) %>% 
  group_by(r_condition, r_file_name, pep_stripped_sequence) %>% 
  summarise(Med_fg_quantity = median(fg_quantity)) %>% 
  mutate(r_condition = factor(r_condition, levels = c("WT", "p427", "p428", "p430")),
         r_file_name = factor(r_file_name, levels = as.character(unique(df_test[order(df_test$r_condition),]$r_file_name)))) 

data_names <- data.frame(
  r_file_name = as.character(unique(df_test[order(df_test$r_condition),]$r_file_name)),
  r_fixed_name = c("WT1", "WT2", "WT3", "scramble1", "scramble2", "scramble3", "shRNA1-1", "shRNA1-2", "shRNA1-3", "shRNA2-1", "shRNA2-2", "shRNA2-3")
)

df_test %>% 
  left_join(data_names, by = "r_file_name") %>% 
  mutate(r_fixed_name = factor(r_fixed_name, levels = c("WT1", "WT2", "WT3", 
                                                        "scramble1", "scramble2", "scramble3", 
                                                        "shRNA1-1", "shRNA1-2", "shRNA1-3", 
                                                        "shRNA2-1", "shRNA2-2", "shRNA2-3"))) %>% 
  ggplot(aes(x = r_fixed_name, y = log2(Med_fg_quantity),
             color = r_condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2),
             size = 0.4, alpha = 0.2) +
  labs(x = "", y = "log2(PEP Quantity)") +
  scale_color_manual(name = "condition", 
                     values = c("WT"="#C9441E","p427"="#6ECBB6","p428"= "#619CFF","p430"= "#E69F00"),
                     labels = c("WT"="WT","p427"="Scramble","p428"= "shRNA1","p430"= "shRNA2"))+
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )


df_test <- sprect_data %>% 
  filter(str_detect(pg_protein_accessions, "Q9UN81")) %>% 
  group_by(r_condition, r_file_name, pep_stripped_sequence) %>% 
  summarise(Med_fg_quantity = median(fg_quantity)) %>% 
  mutate(r_condition = factor(r_condition, levels = c("WT", "p427", "p428", "p430")),
         r_file_name = factor(r_file_name, levels = as.character(unique(df_test[order(df_test$r_condition),]$r_file_name)))) 

df_test %>% 
  left_join(data_names, by = "r_file_name") %>% 
  mutate(r_fixed_name = factor(r_fixed_name, levels = c("WT1", "WT2", "WT3", 
                                                        "scramble1", "scramble2", "scramble3", 
                                                        "shRNA1-1", "shRNA1-2", "shRNA1-3", 
                                                        "shRNA2-1", "shRNA2-2", "shRNA2-3"))) %>% 
  ggplot(aes(x = r_fixed_name, y = log2(Med_fg_quantity),
             color = r_condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2),
             size = 0.4, alpha = 0.2) +
  labs(x = "", y = "log2(PEP Quantity)") +
  scale_color_manual(name = "condition", 
                     values = c("WT"="#C9441E","p427"="#6ECBB6","p428"= "#619CFF","p430"= "#E69F00"),
                     labels = c("WT"="WT","p427"="Scramble","p428"= "shRNA1","p430"= "shRNA2"))+
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

## Proteins

data_proteins <- sprect_data %>% 
  filter(pep_used_for_protein_group_quantity == "TRUE") %>% 
  select(r_condition, r_file_name, pg_protein_accessions, pg_quantity) %>% 
  filter(!is.nan(pg_quantity)) %>% 
  group_by(r_condition, r_file_name, pg_protein_accessions) %>% 
  summarise(med_quantity = median(pg_quantity))  %>% 
  mutate(r_condition = factor(r_condition, levels = c("WT", "p427", "p428", "p430")),
         r_file_name = factor(r_file_name, levels = as.character(unique(df_test[order(df_test$r_condition),]$r_file_name)))
  )

data_proteins %>% 
  group_by(r_condition, r_file_name) %>% 
  summarize(Proteins = n()) %>%
  left_join(data_names, by = "r_file_name") %>% 
  mutate(r_fixed_name = factor(r_fixed_name, levels = c("WT1", "WT2", "WT3", 
                                                        "scramble1", "scramble2", "scramble3", 
                                                        "shRNA1-1", "shRNA1-2", "shRNA1-3", 
                                                        "shRNA2-1", "shRNA2-2", "shRNA2-3"))) %>% 
  ggplot() +
  geom_bar(aes(y = r_fixed_name, x = Proteins, fill = r_condition),
           stat = "identity", show.legend = FALSE) +
  geom_text(aes(y = r_fixed_name, x = Proteins, label = Proteins),
            nudge_x = 50, size = 3) +
  scale_fill_manual(name = "condition", 
                    values = c("WT"="#C9441E","p427"="#6ECBB6","p428"= "#619CFF","p430"= "#E69F00"),
                    labels = c("WT"="WT","p427"="Scramble","p428"= "shRNA1","p430"= "shRNA2"))+
  labs(x = "Proteins", y = "") +
  theme_classic()


x <- data_proteins %>% 
  left_join(data_names, by = "r_file_name") %>% 
  mutate(r_fixed_name = factor(r_fixed_name, levels = c("WT1", "WT2", "WT3", 
                                                        "scramble1", "scramble2", "scramble3", 
                                                        "shRNA1-1", "shRNA1-2", "shRNA1-3", 
                                                        "shRNA2-1", "shRNA2-2", "shRNA2-3")))
x %>% 
  ggplot(aes(x = r_fixed_name, y = log2(med_quantity),
             color = r_condition)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.2),
             size = 0.4, alpha = 0.1) +
  # L1TD1
  geom_point(data = x %>% filter(str_detect(pg_protein_accessions, "Q5T7N2")),
             color = "black") +
  geom_text_repel(data = x %>% filter(str_detect(pg_protein_accessions, "Q5T7N2")),
                  label = "L1TD1",
                  color = "black", size = 1.5)+
  # ORF1 protein
  geom_point(data = x %>% filter(str_detect(pg_protein_accessions, "Q9UN81")),
             color = "darkred") +
  geom_text_repel(data = x %>% filter(str_detect(pg_protein_accessions, "Q9UN81")),
                  label = "L1RE1",
                  color = "black", size = 1.5)+
  labs(x = "", y = "log2(PG Quantity)") +
  scale_color_manual(name = "condition", 
                     values = c("WT"="#C9441E","p427"="#6ECBB6","p428"= "#619CFF","p430"= "#E69F00"),
                     labels = c("WT"="WT","p427"="Scramble","p428"= "shRNA1","p430"= "shRNA2"))+
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )

## PCA

pca_data <- read.delim("../Data/Chart Data.tsv",
                       stringsAsFactors = FALSE)

pca_data %>% 
  ggplot() +
  geom_point(aes(x = PC1, 
                 y = PC2,
                 color = Data),
             size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
  scale_color_manual(name = "Condition", 
                     values = c("WT"="#C9441E","scramble"="#6ECBB6","shRNA1"= "#619CFF","shRNA2"= "#E69F00"))+
  labs(x = "PC1 (30.8%)", y = "PC2 (20.1%)") +
  theme_classic()

## Candidates
candidates_csv <- read.delim("../Data/34-ORS-A1_candidates.tsv", 
                             stringsAsFactors = FALSE)
head(candidates_csv)
names(candidates_csv)
dim(candidates_csv)

## Load DDA results
dda_set <- readRDS("../../ORF2_IP/Output/ORF2_IP_General_MatchRuns_Global_Significant.rds")

## STRAP data
STRAP_Sig <- dda_set$Significant_Tables$ORF2_vs_IgG_Hua_N2102EP

head(STRAP_Sig)

dda_significant <- STRAP_Sig %>% 
  filter(log_FC_Post >= 2, p.adj_Post <= 0.05)
dda_significant

## IDIRT
IDIRT <- read.csv("../../ORF2_IP/ExternalData/CustomAnnotations/DIRT_LINE1.csv", 
                  stringsAsFactors = FALSE)
head(IDIRT)

#Pizarro
other_data <- read.csv("../../ORF2_IP/ExternalData/CustomAnnotations/Cristofari_data.csv",
                       stringsAsFactors = FALSE, na.strings = "") %>% 
  fill(From, .direction = "down") %>% 
  mutate(Datasets = "LINE1-Interactors") %>% 
  select(-Protein) %>% 
  pull(UniprotID) %>% 
  unique()
head(other_data)

df <- candidates_csv %>% 
  filter(Comparison..group1.group2. == "p427 / WT") %>% 
  mutate(
    ## Add DDA significant match
    DDA_sig = if_else(UniProtIds %in% dda_significant$Protein.ID,
                      "Yes", "No"),
    ## Add IDIRT data
    IDIRT = if_else(UniProtIds %in% IDIRT$Uniprot.symbol, "i-DIRT", "No"),
    
    ## Add Datasets
    LINE1_interactors = if_else(UniProtIds %in% other_data, 
                                "LINE1-interactors", "No"),
  )
head(df)

ggplot() +
  geom_point(data = df %>% filter(abs(AVG.Log2.Ratio) > 1,
                                  Qvalue < 0.05, DDA_sig != "Yes"),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)),
             alpha = 0.7, size = 0.7)  +
  ## add DDA sig data
  geom_point(data = df %>% filter(abs(AVG.Log2.Ratio) > 1,
                                  Qvalue < 0.05, DDA_sig == "Yes"),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue),
                 color = DDA_sig), size = 0.7, show.legend = FALSE) +
  ## Not significant
  geom_point(data = df %>% filter(!(abs(AVG.Log2.Ratio) > 1 & Qvalue < 0.05)),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)), alpha = 0.2, color = "grey80", size = 0.7)+ 
  
  ## Add ORF1
  geom_point(data = df %>% filter(UniProtIds %in% c("Q9UN81", "Q5T7N2")),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)),
             color = "darkred", size = 2) +
  ## Add labels
  geom_text_repel(data = df %>% filter(abs(AVG.Log2.Ratio) > 1, Qvalue < 0.05),
                  aes(x = AVG.Log2.Ratio,
                      y = -log10(Qvalue),
                      label = Genes),
                  size = 1.5, box.padding = 0.25, min.segment.length = 0.1,
                  #max.overlaps = Inf
  ) +
  
  #   # Legend
  scale_color_manual(values = c("#C9441E"), labels = c("DDA")) +
  geom_vline(xintercept = c(-1, 1), alpha = 0.2, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), alpha = 0.2, linetype = "dashed") +
  labs(x = "log<sub>2</sub> Fold-Change", y = "-log<sub>10</sub> p-value adj") +
  theme_classic() +
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(),
    legend.title = element_markdown(),
    legend.position = "bottom"
  )

### Shared


sig_proteins <- list(
  "DDA" = dda_significant$Protein.ID,
  "DIA-scramble" = df %>% filter(AVG.Log2.Ratio > 1, Qvalue < 0.05) %>% pull(UniProtIds),
  "DIA-WT" = df %>% filter(AVG.Log2.Ratio < -1, Qvalue < 0.05) %>% pull(UniProtIds)
)

shared_pro  <- euler(sig_proteins, shape = "ellipse")
plot(shared_pro,
     quantities = TRUE,
     labels = list(font = 1, cex = 0.95))


### L1 Interactors highlighted



ggplot() +
  geom_point(data = df %>% filter(abs(AVG.Log2.Ratio) > 1,
                                  Qvalue < 0.05,  LINE1_interactors== "No"),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)),
             alpha = 0.2, size = 0.5)  +
  ## add DDA sig data
  geom_point(data = df %>% filter(abs(AVG.Log2.Ratio) > 1,
                                  Qvalue < 0.05,  LINE1_interactors!= "No"),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue),
                 color = LINE1_interactors), size = 2,
             show.legend = FALSE) +
  ## Not significant
  geom_point(data = df %>% filter(!(abs(AVG.Log2.Ratio) > 1 & Qvalue < 0.05)),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)), alpha = 0.2, color = "grey80", size = 0.7)+ 
  
  ## Add ORF1
  geom_point(data = df %>% filter(UniProtIds %in% c("Q9UN81", "Q5T7N2")),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)),
             color = "darkred", size = 2) +
  ## Add labels
  geom_text_repel(data = df %>% filter(abs(AVG.Log2.Ratio) > 1, 
                                       Qvalue < 0.05,
                                       LINE1_interactors!= "No"),
                  aes(x = AVG.Log2.Ratio,
                      y = -log10(Qvalue),
                      label = Genes),
                  size = 3, box.padding = 0.25, min.segment.length = 0.1, 
                  fontface = "bold",
                  #max.overlaps = Inf
  ) +
  
  #   # Legend
  #scale_color_manual(values = c("#C9441E"), labels = c("LINE-1 Interactors")) +
  geom_vline(xintercept = c(-1, 1), alpha = 0.2, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), alpha = 0.2, linetype = "dashed") +
  labs(x = "log<sub>2</sub> Fold-Change", y = "-log<sub>10</sub> p-value adj") +
  theme_classic() +
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(),
    legend.title = element_markdown(),
    legend.position = "bottom"
  )

ggsave("ORF2p_IP_Volcano_comp1_L1Interactorsmark.pdf")


### IDIRT


ggplot() +
  geom_point(data = df %>% filter(abs(AVG.Log2.Ratio) > 1,
                                  Qvalue < 0.05,  IDIRT!= "i-DIRT"),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)),
             alpha = 0.4, size = 0.7)  +
  ## add DDA sig data
  geom_point(data = df %>% filter(abs(AVG.Log2.Ratio) > 1,
                                  Qvalue < 0.05, IDIRT == "i-DIRT"),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue),
                 color = IDIRT), size = 1,
             show.legend = FALSE) +
  ## Not significant
  geom_point(data = df %>% filter(!(abs(AVG.Log2.Ratio) > 1 & Qvalue < 0.05)),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)), alpha = 0.2, color = "grey80", size = 0.7)+ 
  
  ## Add ORF1
  geom_point(data = df %>% filter(UniProtIds %in% c("Q9UN81", "Q5T7N2")),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)),
             color = "darkred", size = 2) +
  ## Add labels
  geom_text_repel(data = df %>% filter(abs(AVG.Log2.Ratio) > 1, 
                                       Qvalue < 0.05,
                                       IDIRT == "i-DIRT"),
                  aes(x = AVG.Log2.Ratio,
                      y = -log10(Qvalue),
                      label = Genes),
                  size = 3, box.padding = 0.25, min.segment.length = 0.1, 
                  fontface = "bold",
                  #max.overlaps = Inf
  ) +
  
  #   # Legend
  scale_color_manual(values = c("#6ECBB6"), labels = c("i-DIRT")) +
  geom_vline(xintercept = c(-1, 1), alpha = 0.2, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), alpha = 0.2, linetype = "dashed") +
  labs(x = "log<sub>2</sub> Fold-Change", y = "-log<sub>10</sub> p-value adj") +
  theme_classic() +
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(),
    legend.title = element_markdown(),
    legend.position = "bottom"
  )
ggsave("ORF2p_IP_Volcano_comp1_IDIRTHigh.pdf")




## "p428 / p427"


df <- candidates_csv %>% 
  filter(Comparison..group1.group2. == "p428 / p427") %>% 
  mutate(
    ## Add DDA significant match
    DDA_sig = if_else(UniProtIds %in% dda_significant$Protein.ID,
                      "Yes", "No"),
    ## Add IDIRT data
    IDIRT = if_else(UniProtIds %in% IDIRT$Uniprot.symbol, "i-DIRT", "No"),
    
    ## Add Datasets
    LINE1_interactors = if_else(UniProtIds %in% other_data, 
                                "LINE1-interactors", "No"),
  )




df %>% 
  filter(AVG.Log2.Ratio < -10)




ggplot() +
  geom_point(data = df %>% filter(abs(AVG.Log2.Ratio) > 1,
                                  Qvalue < 0.05, DDA_sig != "Yes"),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)),
             alpha = 0.7, size = 0.7)  +
  ## add DDA sig data
  geom_point(data = df %>% filter(abs(AVG.Log2.Ratio) > 1,
                                  Qvalue < 0.05, DDA_sig == "Yes"),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue),
                 color = DDA_sig), size = 0.7, show.legend = FALSE) +
  ## Not significant
  geom_point(data = df %>% filter(!(abs(AVG.Log2.Ratio) > 1 & Qvalue < 0.05)),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)), alpha = 0.2, color = "grey80", size = 0.7)+ 
  
  ## Add ORF1
  geom_point(data = df %>% filter(UniProtIds %in% c("Q9UN81", "Q5T7N2")),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)),
             color = "darkred", size = 2) +
  ## Add labels
  geom_text_repel(data = df %>% filter(abs(AVG.Log2.Ratio) > 1, Qvalue < 0.05),
                  aes(x = AVG.Log2.Ratio,
                      y = -log10(Qvalue),
                      label = Genes),
                  size = 1.5, box.padding = 0.25, min.segment.length = 0.1,
                  #max.overlaps = Inf
  ) +
  
  #   # Legend
  scale_color_manual(values = c("#C9441E"), labels = c("DDA")) +
  scale_x_continuous(limits = c(-8, 7)) +
  geom_vline(xintercept = c(-1, 1), alpha = 0.2, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), alpha = 0.2, linetype = "dashed") +
  labs(x = "log<sub>2</sub> Fold-Change", y = "-log<sub>10</sub> p-value adj") +
  theme_classic() +
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(),
    legend.title = element_markdown(),
    legend.position = "bottom"
  )

ggsave("ORF2p_IP_Volcano_comp2_DDAmark.pdf")



### Shared

How many proteins are shared between DDA and DIA.


sig_proteins <- list(
  "DDA" = dda_significant$Protein.ID,
  "DIA-shRNA1" = df %>% filter(AVG.Log2.Ratio > 1, Qvalue < 0.05) %>% pull(UniProtIds),
  "DIA-scramble" = df %>% filter(AVG.Log2.Ratio < -1, Qvalue < 0.05) %>% pull(UniProtIds)
)

shared_pro  <- euler(sig_proteins, shape = "ellipse")
plot(shared_pro,
     quantities = TRUE,
     labels = list(font = 1, cex = 0.95))


### L1 Interactors highlighted



ggplot() +
  geom_point(data = df %>% filter(abs(AVG.Log2.Ratio) > 1,
                                  Qvalue < 0.05,  LINE1_interactors== "No"),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)),
             alpha = 0.2, size = 0.5)  +
  ## add DDA sig data
  geom_point(data = df %>% filter(abs(AVG.Log2.Ratio) > 1,
                                  Qvalue < 0.05,  LINE1_interactors!= "No"),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue),
                 color = LINE1_interactors), size = 2,
             show.legend = FALSE) +
  ## Not significant
  geom_point(data = df %>% filter(!(abs(AVG.Log2.Ratio) > 1 & Qvalue < 0.05)),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)), alpha = 0.2, color = "grey80", size = 0.7)+ 
  
  ## Add ORF1
  geom_point(data = df %>% filter(UniProtIds %in% c("Q9UN81", "Q5T7N2")),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)),
             color = "darkred", size = 2) +
  ## Add labels
  geom_text_repel(data = df %>% filter(abs(AVG.Log2.Ratio) > 1, 
                                       Qvalue < 0.05,
                                       LINE1_interactors!= "No"),
                  aes(x = AVG.Log2.Ratio,
                      y = -log10(Qvalue),
                      label = Genes),
                  size = 3, box.padding = 0.25, min.segment.length = 0.1, 
                  fontface = "bold",
                  #max.overlaps = Inf
  ) +
  
  #   # Legend
  #scale_color_manual(values = c("#C9441E"), labels = c("LINE-1 Interactors")) +
  scale_x_continuous(limits = c(-8, 7)) +
  geom_vline(xintercept = c(-1, 1), alpha = 0.2, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), alpha = 0.2, linetype = "dashed") +
  labs(x = "log<sub>2</sub> Fold-Change", y = "-log<sub>10</sub> p-value adj") +
  theme_classic() +
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(),
    legend.title = element_markdown(),
    legend.position = "bottom"
  )

ggsave("ORF2p_IP_Volcano_comp2_L1Interactorsmark.pdf")


### IDIRT


ggplot() +
  geom_point(data = df %>% filter(abs(AVG.Log2.Ratio) > 1,
                                  Qvalue < 0.05,  IDIRT!= "i-DIRT"),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)),
             alpha = 0.4, size = 0.7)  +
  ## add DDA sig data
  geom_point(data = df %>% filter(abs(AVG.Log2.Ratio) > 1,
                                  Qvalue < 0.05, IDIRT == "i-DIRT"),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue),
                 color = IDIRT), size = 1,
             show.legend = FALSE) +
  ## Not significant
  geom_point(data = df %>% filter(!(abs(AVG.Log2.Ratio) > 1 & Qvalue < 0.05)),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)), alpha = 0.2, color = "grey80", size = 0.7)+ 
  
  ## Add ORF1
  geom_point(data = df %>% filter(UniProtIds %in% c("Q9UN81", "Q5T7N2")),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)),
             color = "darkred", size = 2) +
  ## Add labels
  geom_text_repel(data = df %>% filter(abs(AVG.Log2.Ratio) > 1, 
                                       Qvalue < 0.05,
                                       IDIRT == "i-DIRT"),
                  aes(x = AVG.Log2.Ratio,
                      y = -log10(Qvalue),
                      label = Genes),
                  size = 3, box.padding = 0.25, min.segment.length = 0.1, 
                  fontface = "bold",
                  #max.overlaps = Inf
  ) +
  
  #   # Legend
  scale_color_manual(values = c("#6ECBB6"), labels = c("i-DIRT")) +
  scale_x_continuous(limits = c(-8, 7)) +
  geom_vline(xintercept = c(-1, 1), alpha = 0.2, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), alpha = 0.2, linetype = "dashed") +
  labs(x = "log<sub>2</sub> Fold-Change", y = "-log<sub>10</sub> p-value adj") +
  theme_classic() +
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(),
    legend.title = element_markdown(),
    legend.position = "bottom"
  )
ggsave("ORF2p_IP_Volcano_comp2_IDIRTHigh.pdf")



## "p430 / p427"


df <- candidates_csv %>% 
  filter(Comparison..group1.group2. == "p428 / p427") %>% 
  mutate(
    ## Add DDA significant match
    DDA_sig = if_else(UniProtIds %in% dda_significant$Protein.ID,
                      "Yes", "No"),
    ## Add IDIRT data
    IDIRT = if_else(UniProtIds %in% IDIRT$Uniprot.symbol, "i-DIRT", "No"),
    
    ## Add Datasets
    LINE1_interactors = if_else(UniProtIds %in% other_data, 
                                "LINE1-interactors", "No"),
  )



ggplot() +
  geom_point(data = df %>% filter(abs(AVG.Log2.Ratio) > 1,
                                  Qvalue < 0.05, DDA_sig != "Yes"),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)),
             alpha = 0.7, size = 0.7)  +
  ## add DDA sig data
  geom_point(data = df %>% filter(abs(AVG.Log2.Ratio) > 1,
                                  Qvalue < 0.05, DDA_sig == "Yes"),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue),
                 color = DDA_sig), size = 0.7, show.legend = FALSE) +
  ## Not significant
  geom_point(data = df %>% filter(!(abs(AVG.Log2.Ratio) > 1 & Qvalue < 0.05)),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)), alpha = 0.2, color = "grey80", size = 0.7)+ 
  
  ## Add ORF1
  geom_point(data = df %>% filter(UniProtIds %in% c("Q9UN81", "Q5T7N2")),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)),
             color = "darkred", size = 2) +
  ## Add labels
  geom_text_repel(data = df %>% filter(abs(AVG.Log2.Ratio) > 1, Qvalue < 0.05),
                  aes(x = AVG.Log2.Ratio,
                      y = -log10(Qvalue),
                      label = Genes),
                  size = 1.5, box.padding = 0.25, min.segment.length = 0.1,
                  #max.overlaps = Inf
  ) +
  
  #   # Legend
  scale_color_manual(values = c("#C9441E"), labels = c("DDA")) +
  geom_vline(xintercept = c(-1, 1), alpha = 0.2, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), alpha = 0.2, linetype = "dashed") +
  labs(x = "log<sub>2</sub> Fold-Change", y = "-log<sub>10</sub> p-value adj") +
  theme_classic() +
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(),
    legend.title = element_markdown(),
    legend.position = "bottom"
  )

ggsave("ORF2p_IP_Volcano_comp3_DDAmark.pdf")


### Shared

How many proteins are shared between DDA and DIA.


sig_proteins <- list(
  "DDA" = dda_significant$Protein.ID,
  "DIA-shRNA2" = df %>% filter(AVG.Log2.Ratio > 1, Qvalue < 0.05) %>% pull(UniProtIds),
  "DIA-scramble" = df %>% filter(AVG.Log2.Ratio < -1, Qvalue < 0.05) %>% pull(UniProtIds)
)

shared_pro  <- euler(sig_proteins, shape = "ellipse")
plot(shared_pro,
     quantities = TRUE,
     labels = list(font = 1, cex = 0.95))


### L1 Interactors highlighted



ggplot() +
  geom_point(data = df %>% filter(abs(AVG.Log2.Ratio) > 1,
                                  Qvalue < 0.05,  LINE1_interactors== "No"),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)),
             alpha = 0.2, size = 0.5)  +
  ## add DDA sig data
  geom_point(data = df %>% filter(abs(AVG.Log2.Ratio) > 1,
                                  Qvalue < 0.05,  LINE1_interactors!= "No"),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue),
                 color = LINE1_interactors), size = 2,
             show.legend = FALSE) +
  ## Not significant
  geom_point(data = df %>% filter(!(abs(AVG.Log2.Ratio) > 1 & Qvalue < 0.05)),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)), alpha = 0.2, color = "grey80", size = 0.7)+ 
  
  ## Add ORF1
  geom_point(data = df %>% filter(UniProtIds %in% c("Q9UN81", "Q5T7N2")),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)),
             color = "darkred", size = 2) +
  ## Add labels
  geom_text_repel(data = df %>% filter(abs(AVG.Log2.Ratio) > 1, 
                                       Qvalue < 0.05,
                                       LINE1_interactors!= "No"),
                  aes(x = AVG.Log2.Ratio,
                      y = -log10(Qvalue),
                      label = Genes),
                  size = 3, box.padding = 0.25, min.segment.length = 0.1, 
                  fontface = "bold",
                  #max.overlaps = Inf
  ) +
  
  #   # Legend
  #scale_color_manual(values = c("#C9441E"), labels = c("LINE-1 Interactors")) +
  #scale_x_continuous(limits = c(-8, 7)) +
  geom_vline(xintercept = c(-1, 1), alpha = 0.2, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), alpha = 0.2, linetype = "dashed") +
  labs(x = "log<sub>2</sub> Fold-Change", y = "-log<sub>10</sub> p-value adj") +
  theme_classic() +
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(),
    legend.title = element_markdown(),
    legend.position = "bottom"
  )

ggsave("ORF2p_IP_Volcano_comp3_L1Interactorsmark.pdf")


### IDIRT


ggplot() +
  geom_point(data = df %>% filter(abs(AVG.Log2.Ratio) > 1,
                                  Qvalue < 0.05,  IDIRT!= "i-DIRT"),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)),
             alpha = 0.4, size = 0.7)  +
  ## add DDA sig data
  geom_point(data = df %>% filter(abs(AVG.Log2.Ratio) > 1,
                                  Qvalue < 0.05, IDIRT == "i-DIRT"),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue),
                 color = IDIRT), size = 1,
             show.legend = FALSE) +
  ## Not significant
  geom_point(data = df %>% filter(!(abs(AVG.Log2.Ratio) > 1 & Qvalue < 0.05)),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)), alpha = 0.2, color = "grey80", size = 0.7)+ 
  
  ## Add ORF1
  geom_point(data = df %>% filter(UniProtIds %in% c("Q9UN81", "Q5T7N2")),
             aes(x = AVG.Log2.Ratio,
                 y = -log10(Qvalue)),
             color = "darkred", size = 2) +
  ## Add labels
  geom_text_repel(data = df %>% filter(abs(AVG.Log2.Ratio) > 1, 
                                       Qvalue < 0.05,
                                       IDIRT == "i-DIRT"),
                  aes(x = AVG.Log2.Ratio,
                      y = -log10(Qvalue),
                      label = Genes),
                  size = 3, box.padding = 0.25, min.segment.length = 0.1, 
                  fontface = "bold",
                  #max.overlaps = Inf
  ) +
  
  #   # Legend
  scale_color_manual(values = c("#6ECBB6"), labels = c("i-DIRT")) +
  #scale_x_continuous(limits = c(-8, 7)) +
  geom_vline(xintercept = c(-1, 1), alpha = 0.2, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), alpha = 0.2, linetype = "dashed") +
  labs(x = "log<sub>2</sub> Fold-Change", y = "-log<sub>10</sub> p-value adj") +
  theme_classic() +
  theme(
    axis.title.x = element_markdown(),
    axis.title.y = element_markdown(),
    legend.title = element_markdown(),
    legend.position = "bottom"
  )
ggsave("ORF2p_IP_Volcano_comp3_IDIRTHigh.pdf")


