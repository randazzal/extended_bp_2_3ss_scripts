library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(pheatmap)
library(purrr)

setwd("/data2/lackey_lab/randazza/microexons/short_read/merge_7to12/")

b11vs15 <- read.delim("HH11vsHH15_brain/SE.MATS.JC.txt")

b15vsday5 <- read.delim("HH15vsHH27_brain/SE.MATS.JC.txt")

bday5vsday7 <- read.delim("HH27vsHH31_brain/SE.MATS.JC.txt")

bday7vsday9 <- read.delim("HH31vsHH35_brain/SE.MATS.JC.txt")

h11vs15 <- read.delim("HH11vsHH15_heart/SE.MATS.JC.txt")

h15vsday5 <- read.delim("HH15vsHH27_heart/SE.MATS.JC.txt")

hday5vsday7 <- read.delim("HH27vsHH31_heart/SE.MATS.JC.txt")

hday7vsday9 <- read.delim("HH31vsHH35_heart/SE.MATS.JC.txt")


HH11 <- b11vs15[ ,c("GeneID", "geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncLevel1")]
colnames(HH11) <- c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
HH15_2 <- b11vs15[ ,c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel2")]
colnames(HH15_2) <- c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
HH15 <- b15vsday5[ ,c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncLevel1")]
colnames(HH15) <- c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
day5_2 <- b15vsday5[ ,c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel2")]
colnames(day5_2) <- c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
day5 <- bday5vsday7[ ,c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncLevel1")]
colnames(day5) <- c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
day7_2 <- bday5vsday7[ ,c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel2")]
colnames(day7_2) <- c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
day7 <- bday7vsday9[ ,c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncLevel1")]
colnames(day7) <- c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
day9 <- bday7vsday9[ ,c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel2")]
colnames(day9) <- c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
HH11h <- h11vs15[ ,c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncLevel1")]
colnames(HH11h) <- c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
HH15_2h <- h11vs15[ ,c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel2")]
colnames(HH15_2h) <- c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
HH15h <- h15vsday5[ ,c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncLevel1")]
colnames(HH15h) <- c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
day5_2h <- h15vsday5[ ,c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel2")]
colnames(day5_2h) <- c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
day5h <- hday5vsday7[ ,c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncLevel1")]
colnames(day5h) <- c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
day7_2h <- hday5vsday7[ ,c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel2")]
colnames(day7_2h) <- c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
day7h <- hday7vsday9[ ,c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncLevel1")]
colnames(day7h) <- c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
day9h <- hday7vsday9[ ,c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel2")]
colnames(day9h) <- c("GeneID","geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")

merge_HH15 <- full_join(HH15_2, HH15, by = c("GeneID","chr", "strand", "exonStart_0base", "exonEnd", "geneSymbol", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel"))
merge_day5 <- full_join(day5_2, day5, by = c("GeneID","chr", "strand", "exonStart_0base", "exonEnd", "geneSymbol", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel"))
merge_day7 <- full_join(day7_2, day7, by = c("GeneID","chr", "strand", "exonStart_0base", "exonEnd", "geneSymbol", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel"))
merge_HH15h <- full_join(HH15_2h, HH15h, by = c("GeneID","chr", "strand", "exonStart_0base", "exonEnd", "geneSymbol", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel"))
merge_day5h <- full_join(day5_2h, day5h, by = c("GeneID","chr", "strand", "exonStart_0base", "exonEnd", "geneSymbol", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel"))
merge_day7h <- full_join(day7_2h, day7h, by = c("GeneID","chr", "strand", "exonStart_0base", "exonEnd", "geneSymbol", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel"))

HH11 <- HH11 %>%
  rowwise() %>%
  mutate(
    IncLevel_mean_HH11 = mean(as.numeric(strsplit(IncLevel, ",")[[1]]), na.rm = TRUE),
    IJC_sum_npc = sum(as.numeric(strsplit(IJC_SAMPLE, ",")[[1]]), na.rm = TRUE),
    SJC_sum_npc = sum(as.numeric(strsplit(SJC_SAMPLE, ",")[[1]]), na.rm = TRUE)
  ) %>%
  ungroup()
merge_HH15 <- merge_HH15 %>%
  rowwise() %>%
  mutate(
    IncLevel_mean_HH15 = mean(as.numeric(strsplit(IncLevel, ",")[[1]]), na.rm = TRUE),
    IJC_sum_npc = sum(as.numeric(strsplit(IJC_SAMPLE, ",")[[1]]), na.rm = TRUE),
    SJC_sum_npc = sum(as.numeric(strsplit(SJC_SAMPLE, ",")[[1]]), na.rm = TRUE)
  ) %>%
  ungroup()
merge_day5 <- merge_day5 %>%
  rowwise() %>%
  mutate(
    IncLevel_mean_day5 = mean(as.numeric(strsplit(IncLevel, ",")[[1]]), na.rm = TRUE),
    IJC_sum_npc = sum(as.numeric(strsplit(IJC_SAMPLE, ",")[[1]]), na.rm = TRUE),
    SJC_sum_npc = sum(as.numeric(strsplit(SJC_SAMPLE, ",")[[1]]), na.rm = TRUE)
  ) %>%
  ungroup()
merge_day7 <- merge_day7 %>%
  rowwise() %>%
  mutate(
    IncLevel_mean_day7 = mean(as.numeric(strsplit(IncLevel, ",")[[1]]), na.rm = TRUE),
    IJC_sum_npc = sum(as.numeric(strsplit(IJC_SAMPLE, ",")[[1]]), na.rm = TRUE),
    SJC_sum_npc = sum(as.numeric(strsplit(SJC_SAMPLE, ",")[[1]]), na.rm = TRUE)
  ) %>%
  ungroup()
day9 <- day9 %>%
  rowwise() %>%
  mutate(
    IncLevel_mean_day9 = mean(as.numeric(strsplit(IncLevel, ",")[[1]]), na.rm = TRUE),
    IJC_sum_npc = sum(as.numeric(strsplit(IJC_SAMPLE, ",")[[1]]), na.rm = TRUE),
    SJC_sum_npc = sum(as.numeric(strsplit(SJC_SAMPLE, ",")[[1]]), na.rm = TRUE)
  ) %>%
  ungroup()
HH11h <- HH11h %>%
  rowwise() %>%
  mutate(
    IncLevel_mean_HH11h = mean(as.numeric(strsplit(IncLevel, ",")[[1]]), na.rm = TRUE),
    IJC_sum_npc = sum(as.numeric(strsplit(IJC_SAMPLE, ",")[[1]]), na.rm = TRUE),
    SJC_sum_npc = sum(as.numeric(strsplit(SJC_SAMPLE, ",")[[1]]), na.rm = TRUE)
  ) %>%
  ungroup()
merge_HH15h <- merge_HH15h %>%
  rowwise() %>%
  mutate(
    IncLevel_mean_HH15h = mean(as.numeric(strsplit(IncLevel, ",")[[1]]), na.rm = TRUE),
    IJC_sum_npc = sum(as.numeric(strsplit(IJC_SAMPLE, ",")[[1]]), na.rm = TRUE),
    SJC_sum_npc = sum(as.numeric(strsplit(SJC_SAMPLE, ",")[[1]]), na.rm = TRUE)
  ) %>%
  ungroup()
merge_day5h <- merge_day5h %>%
  rowwise() %>%
  mutate(
    IncLevel_mean_day5h = mean(as.numeric(strsplit(IncLevel, ",")[[1]]), na.rm = TRUE),
    IJC_sum_npc = sum(as.numeric(strsplit(IJC_SAMPLE, ",")[[1]]), na.rm = TRUE),
    SJC_sum_npc = sum(as.numeric(strsplit(SJC_SAMPLE, ",")[[1]]), na.rm = TRUE)
  ) %>%
  ungroup()
merge_day7h <- merge_day7h %>%
  rowwise() %>%
  mutate(
    IncLevel_mean_day7h = mean(as.numeric(strsplit(IncLevel, ",")[[1]]), na.rm = TRUE),
    IJC_sum_npc = sum(as.numeric(strsplit(IJC_SAMPLE, ",")[[1]]), na.rm = TRUE),
    SJC_sum_npc = sum(as.numeric(strsplit(SJC_SAMPLE, ",")[[1]]), na.rm = TRUE)
  ) %>%
  ungroup()
day9h <- day9h %>%
  rowwise() %>%
  mutate(
    IncLevel_mean_day9h = mean(as.numeric(strsplit(IncLevel, ",")[[1]]), na.rm = TRUE),
    IJC_sum_npc = sum(as.numeric(strsplit(IJC_SAMPLE, ",")[[1]]), na.rm = TRUE),
    SJC_sum_npc = sum(as.numeric(strsplit(SJC_SAMPLE, ",")[[1]]), na.rm = TRUE)
  ) %>%
  ungroup()
HH11$total <- HH11$IJC_sum_npc + HH11$SJC_sum_npc
HH11 <- HH11 %>%
  filter(total >= 15)
HH11h$total <- HH11h$IJC_sum_npc + HH11h$SJC_sum_npc
HH11h <- HH11h %>%
  filter(total >= 15)
merge_HH15$total <- merge_HH15$IJC_sum_npc + merge_HH15$SJC_sum_npc
merge_HH15 <- merge_HH15 %>%
  filter(total >= 15)
merge_HH15h$total <- merge_HH15h$IJC_sum_npc + merge_HH15h$SJC_sum_npc
merge_HH15h <- merge_HH15h %>%
  filter(total >= 15)
merge_day5$total <- merge_day5$IJC_sum_npc + merge_day5$SJC_sum_npc
merge_day5 <- merge_day5 %>%
  filter(total >= 15)
merge_day5h$total <- merge_day5h$IJC_sum_npc + merge_day5h$SJC_sum_npc
merge_day5h <- merge_day5h %>%
  filter(total >= 15)
merge_day7$total <- merge_day7$IJC_sum_npc + merge_day7$SJC_sum_npc
merge_day7 <- merge_day7 %>%
  filter(total >= 15)
merge_day7h$total <- merge_day7h$IJC_sum_npc + merge_day7h$SJC_sum_npc
merge_day7h <- merge_day7h %>%
  filter(total >= 15)
day9$total <- day9$IJC_sum_npc + day9$SJC_sum_npc
day9 <- day9 %>%
  filter(total >= 15)

HH11 <- HH11[ ,c("GeneID","chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd", "IncLevel_mean_HH11")]
HH11_avg <- HH11 %>%
  group_by(GeneID, geneSymbol, exonStart_0base, exonEnd, chr, strand) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
merge_HH15 <- merge_HH15[ ,c("GeneID","chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd", "IncLevel_mean_HH15")]
HH15_avg <- merge_HH15 %>%
  group_by(GeneID,geneSymbol, exonStart_0base, exonEnd, chr, strand) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
merge_day5 <- merge_day5[ ,c("GeneID","chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd", "IncLevel_mean_day5")]
day5_avg <- merge_day5 %>%
  group_by(GeneID,geneSymbol, exonStart_0base, exonEnd, chr, strand) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
merge_day7 <- merge_day7[ ,c("GeneID","chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd", "IncLevel_mean_day7")]
day7_avg <- merge_day7 %>%
  group_by(GeneID,geneSymbol, exonStart_0base, exonEnd, chr, strand) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
day9 <- day9[ ,c("GeneID","chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd", "IncLevel_mean_day9")]
day9_avg <- day9 %>%
  group_by(GeneID,geneSymbol, exonStart_0base, exonEnd, chr, strand) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
HH11h <- HH11h[ ,c("GeneID","chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd", "IncLevel_mean_HH11h")]
HH11h_avg <- HH11h %>%
  group_by(GeneID,geneSymbol, exonStart_0base, exonEnd, chr, strand) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
merge_HH15h <- merge_HH15h[ ,c("GeneID","chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd", "IncLevel_mean_HH15h")]
HH15h_avg <- merge_HH15h %>%
  group_by(GeneID,geneSymbol, exonStart_0base, exonEnd, chr, strand) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
merge_day5h <- merge_day5h[ ,c("GeneID","chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd", "IncLevel_mean_day5h")]
day5h_avg <- merge_day5h %>%
  group_by(GeneID,geneSymbol, exonStart_0base, exonEnd, chr, strand) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
merge_day7h <- merge_day7h[ ,c("GeneID","chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd", "IncLevel_mean_day7h")]
day7h_avg <- merge_day7h %>%
  group_by(GeneID,geneSymbol, exonStart_0base, exonEnd, chr, strand) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
day9h <- day9h[ ,c("GeneID","chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd", "IncLevel_mean_day9h")]
day9h_avg <- day9h %>%
  group_by(GeneID,geneSymbol, exonStart_0base, exonEnd, chr, strand) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

combo <- full_join(HH11_avg, HH15_avg, by = c("GeneID","chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd"))
combo <- full_join(day5_avg, combo, by = c("GeneID","chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd"))
combo <- full_join(day7_avg, combo, by = c("GeneID","chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd"))
combo <- full_join(day9_avg, combo, by = c("GeneID","chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd"))
combo <- full_join(HH11h_avg, combo, by = c("GeneID","chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd"))
combo <- full_join(HH15h_avg, combo, by = c("GeneID","chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd"))
combo <- full_join(day5h_avg, combo, by = c("GeneID","chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd"))
combo <- full_join(day7h_avg, combo, by = c("GeneID","chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd"))
combo <- full_join(day9h_avg, combo, by = c("GeneID","chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd"))
write.table(combo, "all_psi_mandf.txt", sep = "\t")

megacombo <- combo
megacombo$exonEnd <- as.character(megacombo$exonEnd)
megacombo$exonStart_0base <- as.character(megacombo$exonStart_0base)
megacombo$name <- paste(megacombo$geneSymbol, megacombo$exonStart_0base, megacombo$exonEnd, sep = "_")
megacombo$exonEnd <- as.numeric(megacombo$exonEnd)
megacombo$exonStart_0base <- as.numeric(megacombo$exonStart_0base)

megacombo$length <- megacombo$exonEnd - megacombo$exonStart_0base
megacombo <- megacombo[megacombo$length <= 27, ]
#write.table(megacombo, "megacombo_sizelimited.txt", sep = "\t")

heatmap <- megacombo[ ,c("name", "IncLevel_mean_HH11", "IncLevel_mean_HH15", "IncLevel_mean_day5",
                         "IncLevel_mean_day7", "IncLevel_mean_day9", "IncLevel_mean_HH11h",
                         "IncLevel_mean_HH15h", "IncLevel_mean_day5h", "IncLevel_mean_day7h",
                         "IncLevel_mean_day9h")]
heatmap <- as.data.frame(heatmap)
rownames(heatmap) <- heatmap$name
heatmap$name <- NULL

heatmap <- heatmap[rowSums(!is.na(heatmap)) >= 5, ]
heatmap[] <- lapply(heatmap, function(col) {
  col <- as.numeric(col)
  col[is.na(col) | is.nan(col) | is.infinite(col)] <- 0
  return(col)
})

heatmap <- as.matrix(heatmap)
pheatmap(heatmap,
         cluster_cols = FALSE,
         clustering_method = "complete",
         clustering_distance_rows = "manhattan",
         na_col = "black",
         na_row = "black",
         fontsize_col = 5,
         display_numbers = FALSE,
         show_rownames = F)
p <- pheatmap(heatmap,
              cluster_cols = FALSE,
              clustering_method = "complete",
              clustering_distance_rows = "manhattan",
              cutree_rows = 10)
row_clusters <- cutree(p$tree_row, k = 10)
cluster_df <- data.frame(
  Event = names(row_clusters),
  Cluster = row_clusters
)
write.table(cluster_df, "/data2/lackey_lab/randazza/microexons/short_read/merge_7to12/combo_mf_psi_clusters_k13.txt", sep = "\t")

matches <- c("AGAP1", "ROBO1", "CPEB4", "EMC1", "DCTN1", "FRY", "DOCK7", "MEF2A", "ASAP2", "PTPRK", "CLEC16A", "APBB2", "ITSN1")
search <- c("ENSGALG00000016419", "ENSGALG00000007167", "ENSGALG00000016006")
less <- megacombo %>%
  filter(geneSymbol %in% matches)
lesser <- megacombo %>%
  filter(GeneID %in% search)
less <- rbind(less, lesser)
less$exonEnd <- as.numeric(less$exonEnd)
less$exonStart_0base <- as.numeric(less$exonStart_0base)
less <- less %>%
  filter((exonEnd - exonStart_0base) <= 27)
less_1 <- less[1, ]
less_2 <- less[3:4, ]
less_3 <- less[6, ]
less_4 <- less[9:11, ]
less_5 <- less[15:16, ]
less_6 <- less[19, ]
less_7 <- less[25:26, ]
less <- rbind(less_1, less_2, less_3, less_4, less_5, less_6, less_7)
less$name <- c("EMC1","AGAP1","DOCK7","ROBO1","FRY","MEF2A","CPEB4","DCTN1","PTPRK","APBB2","CLEC16A", "ASAP2")
heatmap <- less[ ,c("geneSymbol", "IncLevel_mean_HH11", "IncLevel_mean_HH15", "IncLevel_mean_day5",
                    "IncLevel_mean_day7", "IncLevel_mean_day9", "IncLevel_mean_HH11h",
                    "IncLevel_mean_HH15h", "IncLevel_mean_day5h", "IncLevel_mean_day7h",
                    "IncLevel_mean_day9h")]
heatmap$geneSymbol <- NULL
heatmap[] <- lapply(heatmap, function(col) {
  col <- as.numeric(col)
  col[is.na(col) | is.nan(col) | is.infinite(col)] <- -0.25
  return(col)
})
heatmap <- as.matrix(heatmap)
pheatmap(heatmap,
         cluster_cols = FALSE,
         clustering_method = "complete",
         clustering_distance_rows = "manhattan",
         na_col = "black",
         na_row = "black",
         fontsize_col = 5,
         display_numbers = TRUE)
