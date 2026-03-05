library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggeasy)
library(readxl)
library(writexl)
library(cluster)
library(pheatmap)

ev_info <- read.table("/Users/allierandazza/R_research/VastDB/EVENT_INFO-hg38.tab", sep="\t", header=TRUE)
psi <- read.table("/Users/allierandazza/R_research/VastDB/PSI_TABLE-hg38.tab", sep="\t", header=TRUE)
combo <- merge(psi, event, by = "EVENT")

#filter by length; microexons <= 27, midexons > 100 & <= 200
micro <- combo[combo$LENGTH <= 27, ]
#filter for exon skipping events
filter_1 <- micro %>%
  filter(grepl("HsaEX", EVENT))
odd_col <- filter_1[, seq(7, ncol(filter_1), by = 2)]
filter <- odd_col[apply(odd_col >= 5 & odd_col <=95, 1, function(row) sum(row, na.rm = TRUE) >= 4),]
filter_a <- filter_1[rownames(filter),]

#filter out genes that have N for score in more than 1/2 of the tissues (filter_b)
even_cols <- seq(8, ncol(filter_a), by = 2)
N_counts <- rowSums(sapply(filter_a[, even_cols], function(x) grepl("N,", x)))
filter_b <- filter_a[N_counts <= length(even_cols)/2,]
#filter for genes that have OK or SOK for scores in 1/5 of the tissues (filter_c)
ok_counts <- rowSums(sapply(filter_b[, even_cols], function(x) grepl("SOK,|OK,", x)))
filter_c <- filter_b[ok_counts >= length(even_cols)/5,]

#sort into neural and non-neural
micro <- filter_c
micro_psi <- cbind(micro$EVENT, micro[, seq(7, ncol(micro), by = 2)])
micro_psi[ ,298:304] <- NULL
neural <- c( "Cerebellum_a",
             "Cerebellum_c", "Cortex", "Embr_Cortex_13_17wpc",
             "Embr_Forebrain_9_12wpc", "Embr_Forebrain_St13_14",
             "Embr_Forebrain_St17_20", "Embr_Forebrain_St22_23",
             "Frontal_Gyrus_old", "Frontal_Gyrus_young", 
             "NCC_Cranial_c", "NCC_default_b", "NCC_Enteric_b",
             "NCC_Melano_biased_b", "Neuroblastoma", "Neurons",
             "Neurons_Cortex_KCl_0h", "Neurons_Cortex_KCl_1h",
             "Neurons_Cortex_KCl_6h", 
             "Retina_a", "Retina_macular", "Retina_peripheral",
             "Sup_Temporal_Gyrus", 
             "Whole_Brain_b")
not_neural <- c("Brain_Endoth", "Microglia", "Astrocytes", "Oligodendrocytes", "Adipose_b", "Adipose_c", "Adipose_d", "Adrenal_b", "Adrenal_c",
                "Amnion", "Bladder_a", "Bone_marrow_a", "Bone_marrow_b", "Bone_marrow_c",
                "Breast_a", "Breast_Epith_a", "Chorion", "CL_293T", "CL_EndoBH1C_a",
                "CL_Gm12878", "CL_HeLa", "CL_K562", "CL_LP1", "CL_MB231", "CL_MCF7",
                "CL_PNT2", "CL_SHSY5Y_noRA", "CL_SHSY5Y_RA", "Colon_b", "Colon_sigmoid",
                "Colon_transverse", "Decidua", "Embr_2C_a_A", "Embr_2C_a_B", "Embr_4C_a_A",
                "Embr_4C_a_B", "Embr_4C_a_C", "Embr_8C_a_A", "Embr_8C_a_B", "Embr_8C_a_C",
                "Embr_8C_a_D", "Embr_ICM_a_A", "Embr_Morula_a_A", "Embr_Morula_a_B", "Embr_TE_a_A",
                "Embr_TE_a_B", "EndomStromCells", "EndothCells", "EpithelialCells", "ESC_H1_a",
                "ESC_H1_b", "ESC_H1_c", "ESC_H1_d", "ESC_H9_a", "ESC_H9_b", "Fibroblasts",
                "GLS_cells", "Heart_a", "Heart_b", "Heart_c", "HFDPC", "HMEpC_a", "iPS_a", "iPS_b",
                "Kidney_d", "Kidney_e", "Kidney_f", "Liver_a", "Liver_b", "Liver_c", "Lung_b",
                "Lung_e", "Lung_f", "Lymph_node_b", "Lymph_node_c", "Melanocytes", "MonoNucCells",
                "MSC", "Muscle_b", "Muscle_d", "Muscle_e", "NCC_Cranial_c", "NCC_default_b", "NCC_Enteric_b",
                "NCC_Melano_biased_b", "NPC_a", "NPC_b", "Oocyte_a_A", "Ovary_a", "Ovary_b", "Pancreas_Acinar_Old",
                "Pancreas_Acinar_Young", "Pancreas_Alpha_Old", "Pancreas_Alpha_Young", "Pancreas_Beta_Old",
                "Pancreas_Beta_Young", "Pancreas_Duct_Old", "Pancreas_Duct_Young", "Placenta_a", "Placenta_b",
                "Placenta_c", "Placenta_Epith", "Prostate_b", "Prostate_c", "Prostate_d", "Skin_a", "Skin_b",
                "Skin_c", "Skin_d", "Small_intestine", "Spleen_a", "Spleen_b", "Stomach_a", "Stomach_b",
                "Testis_a", "Testis_b", "Testis_c", "Thymus_a", "Thymus_b", "Thyroid_b", "Thyroid_c",
                "Thyroid_d", "WBC_d", "WBC_e", "Zygote_a_A")

psi_neural <- micro_psi[, neural, drop = FALSE]
psi_non_neural <- micro_psi[, not_neural, drop = FALSE]

p_values <- numeric(nrow(micro_psi))
# Loop over each splicing event (row) to perform a t-test
for (i in seq_len(nrow(micro_psi))) {
  # Extract PSI values for the current splicing event
  neural_values <- as.numeric(psi_neural[i, ])
  non_neural_values <- as.numeric(psi_non_neural[i, ])
  
  # Perform a t-test
  t_test_result <- t.test(neural_values, non_neural_values, 
                          alternative = "greater", var.equal = FALSE)  # Welch's t-test
  p_values[i] <- t_test_result$p.value
}
micro_psi$p_value <- p_values

# Identify significant splicing events
significant_events <- micro_psi[micro_psi$p_value < 0.001, ]

micro_psi$adjusted_p_value <- p.adjust(micro_psi$p_value, method = "fdr")
significant_events_fdr_g <- micro_psi[micro_psi$adjusted_p_value < 0.001, ]

for (i in seq_len(nrow(micro_psi))) {
  # Extract PSI values for the current splicing event
  neural_values <- as.numeric(psi_neural[i, ])
  non_neural_values <- as.numeric(psi_non_neural[i, ])
  
  # Perform a t-test
  t_test_result <- t.test(neural_values, non_neural_values, 
                          alternative = "less", var.equal = FALSE)  # Welch's t-test
  p_values[i] <- t_test_result$p.value
}
micro_psi$p_value_less <- p_values

# Identify significant splicing events
significant_events_less <- micro_psi[micro_psi$p_value_less < 0.001, ]

micro_psi$adjusted_p_value_less <- p.adjust(micro_psi$p_value_less, method = "fdr")
significant_events_fdr_l <- micro_psi[micro_psi$adjusted_p_value_less < 0.001, ]

eliminated_g <- significant_events %>%
  filter(!significant_events$`micro$EVENT` %in% significant_events_fdr_g$`micro$EVENT`)
eliminated_g[ ,2:148] <- NULL
eliminated_l <- significant_events_less %>%
  filter(!significant_events_less$`micro$EVENT` %in% significant_events_fdr_l$`micro$EVENT`)
eliminated_l[ ,2:150] <- NULL
eliminated <- rbind(eliminated_g, eliminated_l)

fdr_g <- micro %>%
  filter(EVENT %in% significant_events_fdr_g$`micro$EVENT`) %>%
  mutate(tissue = "fdr_g")
fdr_l <- micro %>%
  filter(EVENT %in% significant_events_fdr_l$`micro$EVENT`) %>%
  mutate(tissue = "fdr_l")

fail <- micro %>%
  filter(!EVENT %in% significant_events_fdr_g$`micro$EVENT`) %>%
  filter(!EVENT %in% significant_events_fdr_l$`micro$EVENT`)

failed_psi <- micro_psi %>%
  filter(micro$EVENT %in% fail$EVENT)
failed_psi[ ,147:151] <- NULL

results <- failed_psi %>%
  rowwise() %>%
  mutate(
    mean_psi = mean(c_across(!starts_with("micro")), na.rm = TRUE),
    variance = var(c_across(!starts_with("micro")), na.rm = TRUE)
  )

rare <- results[results$mean_psi <= 25, ]
common <- results[results$mean_psi >= 75, ]

rare_set <- micro %>%
  filter(EVENT %in% rare$`micro$EVENT`) %>%
  mutate(tissue = "rare")

common_set <- micro %>%
  filter(EVENT %in% common$`micro$EVENT`) %>%
  mutate(tissue = "common")

other <- micro %>%
  filter(!EVENT %in% fdr_g$EVENT) %>%
  filter(!EVENT %in% fdr_l$EVENT) %>%
  filter(!EVENT %in% rare_set$EVENT) %>%
  filter(!EVENT %in% common_set$EVENT) %>%
  mutate(tissue = "other")

tissue_sets <- rbind(rare_set, common_set)
fdr_set <- rbind(fdr_g, fdr_l, other)
full <- rbind(tissue_sets, fdr_set)
full[ ,7:296] <- NULL

write.table(full, file = "/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/psi_grouped.txt", sep = "\t")

## add column for sizes and in-frame/out-of-frame
full <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/psi_grouped.txt", sep = "\t")
size_label <- full %>%
  mutate(size = case_when(
    LENGTH <10 ~ "1to9",
    LENGTH >= 10 & LENGTH <= 18 ~ "10to18",
    LENGTH >= 19 ~ "19to27"
  ))
final <- size_label %>%
  mutate(frame = ifelse(LENGTH %% 3 == 0, "in-frame", "out-of-frame"))
final$size <- factor(final$size, levels = c("1to9", "10to18", "19to27"))

ggplot(final, aes(x = size, fill = frame)) +
  geom_bar() +
  labs(
    title = "In-Frame Microexons per Size Class",
    x = "Microexon Length",
    y = "Count"
  ) +
  ggeasy::easy_center_title()

final$tissue <- factor(final$tissue, levels = c("fdr_g", "fdr_l", "rare", "common", "other"))

ggplot(final, aes(x = tissue, fill = size)) +
  geom_bar() +
  labs(
    title = "Tissue Expression of Microexon Size Classes",
    x = "Tissue Expression",
    y = "Count"
  ) +
  ggeasy::easy_center_title()

ggplot(final, aes(x = tissue, fill = frame)) +
  geom_bar() +
  labs(
    title = "In-Frame Microexons per Tissue",
    x = "Tissue",
    y = "Count"
  ) +
  ggeasy::easy_center_title()

write.table(final, file = "/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/detailed_list.txt", sep = "\t")
