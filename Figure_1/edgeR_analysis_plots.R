library(edgeR)
library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(ggrepel)

BioCC_input = read.table(file = "/data2/lackey_lab/randazza/microexons/short_read/merge_7to12/counts.txt", sep = "\t", header = T)
#change row names to gene id
rownames(BioCC_input) <- BioCC_input$Geneid
#delete unnecessary columns
extra <- BioCC_input[,1:6]
BioCC_input[,1:6] <- NULL
colnames(BioCC_input) <- c("x5B1", "x5B2", "x5B3", "x5B4", "x5B5", "x5B6", "x5H1", "x5H2", "x5H3", "x5H4",
                           "x5H5", "x5H6", "x7B1", "x7B2", "x7B3", "x7B4", "x7B5", "x7B6", "x7H1", "x7H3", "x7H4",
                           "x7H7", "x7H9", "x7H10", "x9B10",
                           "x9B11", "x9B13", "x9B5", "x9B6", "x9H4", "x9H5", "x9H13", "x9H14", "x9H15", "x9H9",
                           "S11B12", "S11B13", "S11B14", "S11B15",
                           "S11B1", "S11B4", "S11B19", "S11H12", "S11H13", "S11H14", "S11H1", "S11H4", "S11H5",
                           "S11H19", "S15B1", "S15B4", "S15B5", "S15B6", "S15B7", "S15B8", "S15B9", "S15H1", "S15H4",
                           "S15H5", "S15H6", "S15H7", "S15H8", "S15H9", "9B15", "9H22", "11H20")

BioCC_input <- BioCC_input[, c("x5B1", "x5B2", "x5B3", "x5B4", "x5B5", "x5B6", "x5H1", "x5H2", "x5H3", "x5H4",
                               "x5H5", "x5H6", "x7B1", "x7B2", "x7B3", "x7B4", "x7B5", "x7B6", "x7H1", "x7H3", "x7H4",
                               "x7H7", "x7H9", "x7H10", "x9B10", 
                               "x9B11", "x9B13", "x9B5", "x9B6", "9B15", "x9H4", "x9H5", "x9H13", "x9H14", "x9H15", "x9H9", "9H22",
                               "S11B12", "S11B13", "S11B14", "S11B15",
                               "S11B1", "S11B4", "S11B19", "S11H12", "S11H13", "S11H14", "S11H1", "S11H4", "S11H5",
                               "S11H19", "11H20", "S15B1", "S15B4", "S15B5", "S15B6", "S15B7", "S15B8", "S15B9", "S15H1", "S15H4",
                               "S15H5", "S15H6", "S15H7", "S15H8", "S15H9")]
#assign samples to groups
Groups <- as.factor(c("A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B", "C", "C", "C", "C", "C",
                      "C", "D", "D", "D", "D", "D", "D", "E", "E", "E", "E", "E", "E", "F", "F", "F", "F", "F", "F", "F", "G", "G", "G", "G", "G", "G",
                      "G", "H", "H", "H", "H", "H", "H", "H", "H", "I", "I", "I", "I", "I", "I", "I", "J", "J", "J", "J", "J",
                      "J", "J"))

#create DGEList, counts = numeric matrix of read counts, group = factor giving the conditions of each sample
rpk.norm.g <- DGEList(counts = BioCC_input, group = Groups)
keep <- filterByExpr(rpk.norm.g, group = Groups) 
sum(keep) #20649 
cpm_data <- cpm(rpk.norm.g)

rpk.norm.g <- rpk.norm.g[keep, , keep.lib.sizes= FALSE]

###TPM calculations
counts <- rpk.norm.g$counts
gene_lengths <- extra[, c("Geneid", "Length")]
gene_lengths$KB <- gene_lengths$Length / 1000

gene_len_vec <- setNames(
  gene_lengths$KB,
  gene_lengths$Geneid
)
common_genes <- intersect(
  rownames(counts),
  names(gene_len_vec)
)

length(common_genes) #20649
counts <- counts[common_genes, , drop = FALSE]
gene_len_aligned <- gene_len_vec[rownames(counts)]
# Calculate RPK
rpk <- sweep(
  counts,
  1,
  gene_len_aligned,
  FUN = "/"
)

# Calculate TPM
tpm <- sweep(
  rpk,
  2,
  colSums(rpk, na.rm = TRUE) / 1e6,
  FUN = "/"
)

colSums(tpm) 
tpm_data <- as.data.frame(tpm)
colnames(tpm_data) <- c("Male_5B1", "Female_5B1", "Female_5B2", "Male_5B2", "Female_5B3", "Male_5B3",
                        "Male_5H1", "Female_5H1", "Female_5H2", "Female_5H3", "Male_5H2", "Male_5H3",
                        "Female_7B1", "Male_7B1", "Female_7B2", "Female_7B3", "Male_7B2", "Male_7B3",
                        "Female_7H1", "Female_7H2", "Male_7H1", "Female_7H3", "Male_7H2", "Male_7H3",
                        "Female_9B1", "Male_9B1", "Male_9B2", "Female_9B2", "Male_9B3", "Female_9B3", 
                        "Female_9H1", "Female_9H2", "Female_9H3", "Male_9H1", "Male_9H2", "Female_9H4", "Male_9H3", 
                        "Male_11B1", "Male_11B2", "Female_11B1", "Female_11B2", "Male_11B3", "Male_11B4", "Female_11B3",
                        "Male_11H1", "Male_11H2", "Female_11H1", "Male_11H3", "Male_11H4", "Male_11H5", "Female_11H2", "Female_11H3",
                        "Female_15B1", "Female_15B2", "Male_15B1", "Male_15B2", "Male_15B3", "Female_15B3", "Male_15B4",
                        "Female_15H1", "Female_15H2", "Female_15H3", "Male_15H1", "Male_15H2", "Female_15H4", "Male_15H3")
tpm_data$IDs <- rownames(tpm_data)
gene2id <- read.delim("/data2/lackey_lab/DownloadedSequenceData/randazza/new_chick/gene_id_name.txt", sep = "\t", header = FALSE)
colnames(gene2id) <- c("IDs", "gene")
tpm_data <- merge(tpm_data, gene2id, by = "IDs")
write.table(tpm_data, "/data2/lackey_lab/randazza/microexons/short_read/merge_7to12/all_tpms.txt", sep = "\t")

gois_micro <- c("ENSGALG00000004045", "ENSGALG00000047813", "ENSGALG00000016419", "ENSGALG00000007167",
          "ENSGALG00000033512", "ENSGALG00000035156", "ENSGALG00000010967", "ENSGALG00000003969",
          "ENSGALG00000017075", "ENSGALG00000036677",
          "ENSGALG00000032366", "ENSGALG00000037423", "ENSGALG00000015511")

limited <- tpm_data %>%
  filter(IDs %in% gois_micro)
limited_long <- limited %>%
  pivot_longer(cols = c("Male_5B1", "Female_5B1", "Female_5B2", "Male_5B2", "Female_5B3", "Male_5B3",
                        "Male_5H1", "Female_5H1", "Female_5H2", "Female_5H3", "Male_5H2", "Male_5H3",
                        "Female_7B1", "Male_7B1", "Female_7B2", "Female_7B3", "Male_7B2", "Male_7B3",
                        "Female_7H1", "Female_7H2", "Male_7H1", "Female_7H3", "Male_7H2", "Male_7H3",
                        "Female_9B1", "Male_9B1", "Male_9B2", "Female_9B2", "Male_9B3",
                        "Female_9H1", "Female_9H2", "Female_9H3", "Male_9H1", "Male_9H2", "Female_9H4",
                        "Male_11B1", "Male_11B2", "Female_11B1", "Female_11B2", "Male_11B3", "Male_11B4", "Female_11B3",
                        "Male_11H1", "Male_11H2", "Female_11H1", "Male_11H3", "Male_11H4", "Male_11H5", "Female_11H2",
                        "Female_15B1", "Female_15B2", "Male_15B1", "Male_15B2", "Male_15B3", "Female_15B3", "Male_15B4",
                        "Female_15H1", "Female_15H2", "Female_15H3", "Male_15H1", "Male_15H2", "Female_15H4", "Male_15H3",
                        "Female_9B3", "Male_9H3", "Female_11H3"),
               names_to = "Sample", values_to = "TPMs")
tpm_long <- limited_long %>%
  mutate(
    Sex = str_extract(Sample, "Male|Female"),
    Stage = str_extract(Sample, "(?<=_)[0-9]+"),
    Tissue = sub(".*_[0-9]+([A-Z])[0-9]+$", "\\1", Sample),
    Replicate = str_extract(Sample, "[0-9]+$")
  )
names_micro <- data.frame(
  IDs = c("ENSGALG00000004045", "ENSGALG00000047813", "ENSGALG00000016419", "ENSGALG00000007167",
          "ENSGALG00000033512", "ENSGALG00000035156", "ENSGALG00000010967", "ENSGALG00000003969",
          "ENSGALG00000017075", "ENSGALG00000036677",
          "ENSGALG00000032366", "ENSGALG00000037423", "ENSGALG00000015511"),
  Gene = c("AGAP1", "APBB2", "ASAP2", "CLEC16A", "CPEB4", "DCTN1", "DOCK7", "EMC1", "FRY", "ITSN1",
           "MEF2A", "PTPRK", "ROBO1")
)

tpm_long <- merge(tpm_long, names_micro, 
                  by.x = "IDs", by.y = "IDs", 
                  all.x = TRUE)
tpm_summary <- tpm_long %>%
  group_by(IDs, Sex, Tissue, Stage) %>%
  summarise(mean_TPM = mean(TPMs, na.rm = TRUE), .groups = "drop")

tpm_long %>%
  filter(Tissue == "B") %>%
  ggplot(aes(x = Stage, y = TPMs, fill = Sex)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8)) +
  geom_point(aes(group = Sex), color = "black",
             position = position_jitterdodge(
               jitter.width = 0.15, dodge.width = 0.8), size = 1.5, alpha = 0.7) +
  facet_wrap(~Gene, scales = "free_y") +
  theme_bw() +
  labs(title = "Expression in tissue", y = "TPM", x = "Stage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

tpm_long %>%
  filter(Tissue == "H") %>%
  ggplot(aes(x = Stage, y = TPMs, fill = Sex)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8)) +
  geom_point(aes(group = Sex), color = "black",
             position = position_jitterdodge(
               jitter.width = 0.15, dodge.width = 0.8), size = 1.5, alpha = 0.7) +
  facet_wrap(~Gene, scales = "free_y") +
  theme_bw() +
  labs(title = "Expression in tissue", y = "TPM", x = "Stage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

wilcox_by_group <- function(
    df,
    stage_col = "Stage",
    tissue_col = "Tissue",
    gene_col = "Gene",
    sex_col = "Sex",
    value_col = "TPMs"
) {
  
  required_cols <- c(
    stage_col, tissue_col, gene_col,
    sex_col, value_col
  )
  
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  dat <- df[
    !is.na(df[[value_col]]) &
      !is.na(df[[sex_col]]) &
      df[[sex_col]] %in% c("Male", "Female"),
    ,
    drop = FALSE
  ]
  
  groups <- unique(dat[c(stage_col, tissue_col, gene_col)])
  
  results <- lapply(seq_len(nrow(groups)), function(i) {
    
    subgroup <- dat[
      dat[[stage_col]] == groups[[stage_col]][i] &
        dat[[tissue_col]] == groups[[tissue_col]][i] &
        dat[[gene_col]] == groups[[gene_col]][i],
      ,
      drop = FALSE
    ]
    
    male <- subgroup[subgroup[[sex_col]] == "Male", value_col]
    female <- subgroup[subgroup[[sex_col]] == "Female", value_col]
    
    if (length(male) < 1 || length(female) < 1) {
      return(data.frame(
        Stage = groups[[stage_col]][i],
        Tissue = groups[[tissue_col]][i],
        Gene = groups[[gene_col]][i],
        n_male = length(male),
        n_female = length(female),
        statistic = NA,
        p_value = NA
      ))
    }
    
    test <- wilcox.test(
      male,
      female,
      exact = FALSE
    )
    
    data.frame(
      Stage = groups[[stage_col]][i],
      Tissue = groups[[tissue_col]][i],
      Gene = groups[[gene_col]][i],
      n_male = length(male),
      n_female = length(female),
      statistic = unname(test$statistic),
      p_value = test$p.value
    )
  })
  
  results <- do.call(rbind, results)
  results$FDR <- p.adjust(results$p_value, method = "BH")
  
  results
}
results <- wilcox_by_group(tpm_long)

## RBP expression levels
gois_rbp <- c("ENSGALG00000064136", "ENSGALG00000066499", "ENSGALG00000012540", "ENSGALG00000032889", "ENSGALG00000034325", "ENSGALG00000011305", "ENSGALG00000009288", "ENSGALG00000001962")

limited <- tpm_data %>%
  filter(IDs %in% gois_rbp)
limited_long <- limited %>%
  pivot_longer(cols = c("Male_5B1", "Female_5B1", "Female_5B2", "Male_5B2", "Female_5B3", "Male_5B3",
                        "Male_5H1", "Female_5H1", "Female_5H2", "Female_5H3", "Male_5H2", "Male_5H3",
                        "Female_7B1", "Male_7B1", "Female_7B2", "Female_7B3", "Male_7B2", "Male_7B3",
                        "Female_7H1", "Female_7H2", "Male_7H1", "Female_7H3", "Male_7H2", "Male_7H3",
                        "Female_9B1", "Male_9B1", "Male_9B2", "Female_9B2", "Male_9B3",
                        "Female_9H1", "Female_9H2", "Female_9H3", "Male_9H1", "Male_9H2", "Female_9H4",
                        "Male_11B1", "Male_11B2", "Female_11B1", "Female_11B2", "Male_11B3", "Male_11B4", "Female_11B3",
                        "Male_11H1", "Male_11H2", "Female_11H1", "Male_11H3", "Male_11H4", "Male_11H5", "Female_11H2",
                        "Female_15B1", "Female_15B2", "Male_15B1", "Male_15B2", "Male_15B3", "Female_15B3", "Male_15B4",
                        "Female_15H1", "Female_15H2", "Female_15H3", "Male_15H1", "Male_15H2", "Female_15H4", "Male_15H3",
                        "Female_9B3", "Male_9H3", "Female_11H3"),
               names_to = "Sample", values_to = "TPMs")
tpm_long <- limited_long %>%
  mutate(
    Sex = str_extract(Sample, "Male|Female"),
    Stage = str_extract(Sample, "(?<=_)[0-9]+"),
    Tissue = sub(".*_[0-9]+([A-Z])[0-9]+$", "\\1", Sample),
    Replicate = str_extract(Sample, "[0-9]+$")
  )

names_rbp <- data.frame(
  IDs = c("ENSGALG00000064136", "ENSGALG00000066499", "ENSGALG00000012540", "ENSGALG00000032889", "ENSGALG00000034325", "ENSGALG00000011305", "ENSGALG00000009288", "ENSGALG00000001962"),
  Gene = c("NOVA1", "SRRM4", "RBFOX2", "RBFOX3", "RBFOX1", "SRSF11", "RNPS1", "PTBP1")
)

tpm_long <- merge(tpm_long, names_rbp, 
                  by.x = "IDs", by.y = "IDs", 
                  all.x = TRUE)
tpm_summary <- tpm_long %>%
  group_by(IDs, Sex, Tissue, Stage) %>%
  summarise(mean_TPM = mean(TPMs, na.rm = TRUE), .groups = "drop")

tpm_long %>%
  filter(Tissue == "B") %>%
  ggplot(aes(x = Stage, y = TPMs, fill = Sex)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8)) +
  geom_point(aes(group = Sex), color = "black",
             position = position_jitterdodge(
               jitter.width = 0.15, dodge.width = 0.8), size = 1.5, alpha = 0.7) +
  facet_wrap(~Gene, scales = "free_y") +
  theme_bw() +
  labs(title = "Expression in tissue", y = "TPM", x = "Stage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

tpm_long %>%
  filter(Tissue == "H") %>%
  ggplot(aes(x = Stage, y = TPMs, fill = Sex)) +
  geom_bar(stat = "summary", fun = "mean", position = position_dodge(width = 0.8)) +
  geom_point(aes(group = Sex), color = "black",
             position = position_jitterdodge(
               jitter.width = 0.15, dodge.width = 0.8), size = 1.5, alpha = 0.7) +
  facet_wrap(~Gene, scales = "free_y") +
  theme_bw() +
  labs(title = "Expression in tissue", y = "TPM", x = "Stage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

wilcox_by_group <- function(
    df,
    stage_col = "Stage",
    tissue_col = "Tissue",
    gene_col = "Gene",
    sex_col = "Sex",
    value_col = "TPMs"
) {
  
  required_cols <- c(
    stage_col, tissue_col, gene_col,
    sex_col, value_col
  )
  
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  dat <- df[
    !is.na(df[[value_col]]) &
      !is.na(df[[sex_col]]) &
      df[[sex_col]] %in% c("Male", "Female"),
    ,
    drop = FALSE
  ]
  
  groups <- unique(dat[c(stage_col, tissue_col, gene_col)])
  
  results <- lapply(seq_len(nrow(groups)), function(i) {
    
    subgroup <- dat[
      dat[[stage_col]] == groups[[stage_col]][i] &
        dat[[tissue_col]] == groups[[tissue_col]][i] &
        dat[[gene_col]] == groups[[gene_col]][i],
      ,
      drop = FALSE
    ]
    
    male <- subgroup[subgroup[[sex_col]] == "Male", value_col]
    female <- subgroup[subgroup[[sex_col]] == "Female", value_col]
    
    if (length(male) < 1 || length(female) < 1) {
      return(data.frame(
        Stage = groups[[stage_col]][i],
        Tissue = groups[[tissue_col]][i],
        Gene = groups[[gene_col]][i],
        n_male = length(male),
        n_female = length(female),
        statistic = NA,
        p_value = NA
      ))
    }
    
    test <- wilcox.test(
      male,
      female,
      exact = FALSE
    )
    
    data.frame(
      Stage = groups[[stage_col]][i],
      Tissue = groups[[tissue_col]][i],
      Gene = groups[[gene_col]][i],
      n_male = length(male),
      n_female = length(female),
      statistic = unname(test$statistic),
      p_value = test$p.value
    )
  })
  
  results <- do.call(rbind, results)
  results$FDR <- p.adjust(results$p_value, method = "BH")
  
  results
}
results <- wilcox_by_group(tpm_long)
