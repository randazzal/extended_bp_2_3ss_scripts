library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(pheatmap)
library(purrr)

setwd("/data2/lackey_lab/randazza/microexons/short_read/merge_7to12/RMATS/")
data <- read.delim("/data2/lackey_lab/randazza/microexons/devo/chicken/gga_classes.txt")
data2 <- read.delim("/data2/lackey_lab/randazza/microexons/devo/chicken/chicken_data.txt")
data2$tissue <- "none"
data <- full_join(data, data2, by = "EVENT")
bed_file <- read.delim("/data2/lackey_lab/randazza/microexons/devo/chicken/bed_files/fixed_gga_class_6.bed", header = FALSE)
bed_file2 <- read.delim("/data2/lackey_lab/randazza/microexons/devo/chicken/bed_files/fixed_micro_gal6.bed", header = FALSE)
bed_file <- full_join(bed_file, bed_file2, by = "V4")
colnames(bed_file) <- c("chr", "start", "end", "EVENT")
dctn1 <- c("4", "91039774", "91039788", "GgaEX1018908")
agap1 <- c("7", "5429561", "5429572", "GgaEX0005391")
asap2 <- c("3", "96151176", "96151184", "GgaEXoo29181")
fry <- c("1", "175973227", "175973235", "GgaEX0029190")
bed_file <- rbind(bed_file, dctn1, agap1, asap2, fry)
data <- merge(data, bed_file, by = "EVENT")
rm(data2)
rm(bed_file2)

b11vs15 <- read.delim("rmats_output/post/combo_HH11vsHH15_brain/SE.MATS.JC.txt")
b11vs15 <- b11vs15 %>%
  filter((exonEnd - exonStart_0base) <= 27)

b15vsday5 <- read.delim("rmats_output/post/combo_HH15_vs_day5_brain/SE.MATS.JC.txt")
b15vsday5 <- b15vsday5 %>%
  filter((exonEnd - exonStart_0base) <= 27)

bday5vsday7 <- read.delim("rmats_output/post/combo_day5_vs_day7_brain/SE.MATS.JC.txt")
bday5vsday7 <- bday5vsday7 %>%
  filter((exonEnd - exonStart_0base) <= 27)

bday7vsday9 <- read.delim("rmats_output/post/combo_day7_vs_day9_brain/SE.MATS.JC.txt")
bday7vsday9 <- bday7vsday9 %>%
  filter((exonEnd - exonStart_0base) <= 27)

h11vs15 <- read.delim("rmats_output/post/combo_HH11vsHH15_heart/SE.MATS.JC.txt")
h11vs15 <- h11vs15 %>%
  filter((exonEnd - exonStart_0base) <= 27)

h15vsday5 <- read.delim("rmats_output/post/combo_HH15_vs_day5_heart/SE.MATS.JC.txt")
h15vsday5 <- h15vsday5 %>%
  filter((exonEnd - exonStart_0base) <= 27)

hday5vsday7 <- read.delim("rmats_output/post/combo_day5_vs_day7_heart/SE.MATS.JC.txt")
hday5vsday7 <- hday5vsday7 %>%
  filter((exonEnd - exonStart_0base) <= 27)

hday7vsday9 <- read.delim("rmats_output/post/combo_day7_vs_day9_heart/SE.MATS.JC.txt")
hday7vsday9 <- hday7vsday9 %>%
  filter((exonEnd - exonStart_0base) <= 27)

HH11 <- b11vs15[ ,c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncLevel1")]
colnames(HH11) <- c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
HH15_2 <- b11vs15[ ,c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel2")]
colnames(HH15_2) <- c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
HH15 <- b15vsday5[ ,c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncLevel1")]
colnames(HH15) <- c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
day5_2 <- b15vsday5[ ,c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel2")]
colnames(day5_2) <- c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
day5 <- bday5vsday7[ ,c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncLevel1")]
colnames(day5) <- c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
day7_2 <- bday5vsday7[ ,c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel2")]
colnames(day7_2) <- c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
day7 <- bday7vsday9[ ,c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncLevel1")]
colnames(day7) <- c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
day9 <- bday7vsday9[ ,c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel2")]
colnames(day9) <- c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
HH11h <- h11vs15[ ,c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncLevel1")]
colnames(HH11h) <- c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
HH15_2h <- h11vs15[ ,c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel2")]
colnames(HH15_2h) <- c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
HH15h <- h15vsday5[ ,c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncLevel1")]
colnames(HH15h) <- c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
day5_2h <- h15vsday5[ ,c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel2")]
colnames(day5_2h) <- c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
day5h <- hday5vsday7[ ,c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncLevel1")]
colnames(day5h) <- c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
day7_2h <- hday5vsday7[ ,c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel2")]
colnames(day7_2h) <- c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
day7h <- hday7vsday9[ ,c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IncLevel1")]
colnames(day7h) <- c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")
day9h <- hday7vsday9[ ,c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncLevel2")]
colnames(day9h) <- c("geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel")

merge_HH15 <- full_join(HH15_2, HH15, by = c("chr", "strand", "exonStart_0base", "exonEnd", "geneSymbol", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel"))
merge_day5 <- full_join(day5_2, day5, by = c("chr", "strand", "exonStart_0base", "exonEnd", "geneSymbol", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel"))
merge_day7 <- full_join(day7_2, day7, by = c("chr", "strand", "exonStart_0base", "exonEnd", "geneSymbol", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel"))
merge_HH15h <- full_join(HH15_2h, HH15h, by = c("chr", "strand", "exonStart_0base", "exonEnd", "geneSymbol", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel"))
merge_day5h <- full_join(day5_2h, day5h, by = c("chr", "strand", "exonStart_0base", "exonEnd", "geneSymbol", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel"))
merge_day7h <- full_join(day7_2h, day7h, by = c("chr", "strand", "exonStart_0base", "exonEnd", "geneSymbol", "IJC_SAMPLE", "SJC_SAMPLE", "IncLevel"))

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

HH11 <- HH11[ ,c("chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd", "IncLevel_mean_HH11")]
HH11_avg <- HH11 %>%
  group_by(geneSymbol, exonStart_0base, exonEnd, chr, strand) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
merge_HH15 <- merge_HH15[ ,c("chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd", "IncLevel_mean_HH15")]
HH15_avg <- merge_HH15 %>%
  group_by(geneSymbol, exonStart_0base, exonEnd, chr, strand) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
merge_day5 <- merge_day5[ ,c("chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd", "IncLevel_mean_day5")]
day5_avg <- merge_day5 %>%
  group_by(geneSymbol, exonStart_0base, exonEnd, chr, strand) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
merge_day7 <- merge_day7[ ,c("chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd", "IncLevel_mean_day7")]
day7_avg <- merge_day7 %>%
  group_by(geneSymbol, exonStart_0base, exonEnd, chr, strand) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
day9 <- day9[ ,c("chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd", "IncLevel_mean_day9")]
day9_avg <- day9 %>%
  group_by(geneSymbol, exonStart_0base, exonEnd, chr, strand) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
HH11h <- HH11h[ ,c("chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd", "IncLevel_mean_HH11h")]
HH11h_avg <- HH11h %>%
  group_by(geneSymbol, exonStart_0base, exonEnd, chr, strand) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
merge_HH15h <- merge_HH15h[ ,c("chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd", "IncLevel_mean_HH15h")]
HH15h_avg <- merge_HH15h %>%
  group_by(geneSymbol, exonStart_0base, exonEnd, chr, strand) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
merge_day5h <- merge_day5h[ ,c("chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd", "IncLevel_mean_day5h")]
day5h_avg <- merge_day5h %>%
  group_by(geneSymbol, exonStart_0base, exonEnd, chr, strand) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
merge_day7h <- merge_day7h[ ,c("chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd", "IncLevel_mean_day7h")]
day7h_avg <- merge_day7h %>%
  group_by(geneSymbol, exonStart_0base, exonEnd, chr, strand) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")
day9h <- day9h[ ,c("chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd", "IncLevel_mean_day9h")]
day9h_avg <- day9h %>%
  group_by(geneSymbol, exonStart_0base, exonEnd, chr, strand) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop")

combo <- full_join(HH11_avg, HH15_avg, by = c("chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd"))
combo <- full_join(day5_avg, combo, by = c("chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd"))
combo <- full_join(day7_avg, combo, by = c("chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd"))
combo <- full_join(day9_avg, combo, by = c("chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd"))
combo <- full_join(HH11h_avg, combo, by = c("chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd"))
combo <- full_join(HH15h_avg, combo, by = c("chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd"))
combo <- full_join(day5h_avg, combo, by = c("chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd"))
combo <- full_join(day7h_avg, combo, by = c("chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd"))
combo <- full_join(day9h_avg, combo, by = c("chr", "strand", "geneSymbol", "exonStart_0base", "exonEnd"))

megacombo <- combo
megacombo$exonEnd <- as.character(megacombo$exonEnd)
megacombo$exonStart_0base <- as.character(megacombo$exonStart_0base)
megacombo$name <- paste(megacombo$geneSymbol, megacombo$exonStart_0base, megacombo$exonEnd, sep = "_")
write.table(megacombo, "../megacombo_sizelimited.txt", sep = "\t")
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
         display_numbers = FALSE)
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

matches <- c("AGAP1", "ROBO1", "CPEB4", "EMC1", "DCTN1", "FRY", "DOCK7", "MEF2A", "ASAP2", "PTPRK", "CLEC16A", "APBB2", "ITSN1")
less <- megacombo %>%
  filter(geneSymbol %in% matches)
less_1 <- less[2, ]
less_2 <- less[4:5, ]
less_3 <- less[9:10, ]
less_4 <- less[12, ]
less_5 <- less[14:16, ]
less_6 <- less[18, ]
less <- rbind(less_1, less_2, less_3, less_4, less_5, less_6)
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
