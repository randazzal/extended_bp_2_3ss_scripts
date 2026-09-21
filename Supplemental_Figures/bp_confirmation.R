library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(ggeasy)
library(rstatix)

## look at results of lariat sequencing and compare to bpfinder predictions
########## midexons
#use bedtools intersect for overlap with exons of interest and those found in CoLa-seq
lariats <- read.delim("/data2/lackey_lab/DownloadedSequenceData/randazza/lariat_seq/Zeng_colaseq/ab_avg_intersect.bed", sep = "\t", header = FALSE)
lariats <- lariats %>%
  separate(V10, into = c("GENE.x", "EVENT", "exon"), sep = "_", remove = TRUE)
lariats <- lariats[c("V1", "V2", "V3", "V5", "V6", "V8", "V9","GENE.x", "EVENT")]
avg <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/avg_regroup.txt", sep = "\t")
avg <- avg[,c("GENE.x", "EVENT", "COORD", "tissue")]
less <- merge(lariats, avg, by = c("GENE.x", "EVENT"))
less <- less[,c("GENE.x", "EVENT", "V1", "V2", "V3","V5","V6", "COORD", "tissue")]
colnames(less) <- c("GENE.x", "EVENT", "bp_chr", "bp_first", "bp_last","score","strand", "COORD", "tissue")
less <- less %>%
  separate(COORD, into = c("mic_chr", "range"), sep = ":", remove = TRUE) %>%
  separate(range, into = c("mic_first", "mic_last"), sep = "-", remove = TRUE)
filtered <- less %>%
  filter(
    (strand == "+" & bp_first < mic_first & bp_last < mic_last) |
      (strand == "-" & bp_first > mic_first & bp_last > mic_last)
  )

boundaries <- function(df) {
  df$distance <- ifelse(df$strand == "+", df$mic_first - df$bp_last, df$bp_first - df$mic_last)
  return(df)
}
filtered$mic_first <- as.integer(filtered$mic_first)
filtered$mic_last <- as.integer(filtered$mic_last)
filtered <- boundaries(filtered)
filtered <- filtered %>%
  filter(abs(distance) < 100)
filtered$name <- paste(filtered$EVENT,filtered$GENE.x, sep = "_")

rm(avg)
rm(lariats)

#load SVM-bfpinder predicted branchpoints
neu_mid <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/SVM-BPfinder/neu_avg_intron_results.txt", sep = "\t")
neu_mid <- neu_mid %>%
  separate(seq_id, into = c("name","coord"), sep = "::", remove = TRUE)
neu_mid <- neu_mid[,c("name","coord", "ss_dist", "bp_seq", "svm_scr")]
nn_mid <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/SVM-BPfinder/nn_avg_intron_results.txt", sep = "\t")
nn_mid <- nn_mid %>%
  separate(seq_id, into = c("name","coord"), sep = "::", remove = TRUE)
nn_mid <- nn_mid[,c("name","coord", "ss_dist", "bp_seq", "svm_scr")]

neu_filter <- merge(filtered, neu_mid, by = "name")
nn_filter <- merge(filtered, nn_mid, by = "name")
#neu = 92 nn = 2133
neu_filter$new_name <- paste(neu_filter$name, neu_filter$bp_first, neu_filter$bp_last, sep = "_")
nn_filter$new_name <- paste(nn_filter$name, nn_filter$bp_first, nn_filter$bp_last, sep = "_")

ranked_neu <- neu_filter %>%
  group_by(new_name) %>%
  mutate(rank = min_rank(desc(svm_scr))) %>%
  ungroup()
ranked_nn <- nn_filter %>%
  group_by(new_name) %>%
  mutate(rank = min_rank(desc(svm_scr))) %>%
  ungroup()

ranked_neu <- ranked_neu %>%
  filter(abs(distance - ss_dist) <= 5)
ranked_nn <- ranked_nn %>%
  filter(abs(distance - ss_dist) <= 5)

clean_ranked_neu <- ranked_neu %>%
  mutate(distance_diff = abs(distance - ss_dist)) %>%
  group_by(new_name) %>%
  slice_min(distance_diff, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-distance_diff)

clean_ranked_nn <- ranked_nn %>%
  mutate(distance_diff = abs(distance - ss_dist)) %>%
  group_by(new_name) %>%
  slice_min(distance_diff, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-distance_diff)

clean_ranked_neu <- clean_ranked_neu %>%
  group_by(name) %>%
  distinct(ss_dist, .keep_all = TRUE) %>%
  ungroup()
clean_ranked_nn <- clean_ranked_nn %>%
  group_by(name) %>%
  distinct(ss_dist, .keep_all = TRUE) %>%
  ungroup()

clean_ranked_neu <- clean_ranked_neu %>%
  group_by(name) %>%
  filter(
    n_distinct(new_name) <= 1 |
      rank == min(rank, na.rm = TRUE)
  ) %>%
  ungroup()
clean_ranked_nn <- clean_ranked_nn %>%
  group_by(name) %>%
  filter(
    n_distinct(new_name) <= 1 |
      rank == min(rank, na.rm = TRUE)
  ) %>%
  ungroup()

combo <- rbind(clean_ranked_neu, clean_ranked_nn)

combo <- combo %>%
  mutate(
    new_rank = factor(
      if_else(rank >= 3, "other", as.character(rank)),
      levels = c("1", "2", "other")
    )
  )

combo %>%
  ggplot(aes(x = new_rank)) +
  geom_bar() +
  theme_minimal()
sum(combo$new_rank == "other", na.rm = TRUE)
break_down <- data.frame(rank = c("1", "2", "Other", "No Prediction"),
                         count = as.integer(c("1166", "397", "259", "403")))
break_down$fraction <- break_down$count/2225
break_down$rank <- factor(break_down$rank, levels = c("1", "2", "Other", "No Prediction"))
break_down %>%
  ggplot(aes(x = rank, y = fraction)) +
  geom_col() +
  theme_minimal()

########## microexon's turn
lariats <- read.delim("/data2/lackey_lab/DownloadedSequenceData/randazza/lariat_seq/Zeng_colaseq/ab_micro_intersect.bed", sep = "\t", header = FALSE)
lariats <- lariats %>%
  separate(V4, into = c("GENE.x", "EVENT", "exon"), sep = "_", remove = TRUE)
micro <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/detailed_list.txt", sep = "\t")
less <- merge(lariats, micro, by = c("GENE.x", "EVENT"))
less <- less[,c("GENE.x", "EVENT", "V1", "V2", "V3","V6", "COORD")]
colnames(less) <- c("GENE.x", "EVENT", "bp_chr", "bp_first", "bp_last","strand", "COORD")
less <- less %>%
  separate(COORD, into = c("mic_chr", "range"), sep = ":", remove = TRUE) %>%
  separate(range, into = c("mic_first", "mic_last"), sep = "-", remove = TRUE)
filtered <- less %>%
  filter(
    (strand == "+" & bp_first < mic_first & bp_last < mic_last) |
      (strand == "-" & bp_first > mic_first & bp_last > mic_last)
  )

boundaries <- function(df) {
  df$distance <- ifelse(df$strand == "+", df$mic_first - df$bp_last, df$bp_first - df$mic_last)
  return(df)
}
filtered$mic_first <- as.integer(filtered$mic_first)
filtered$mic_last <- as.integer(filtered$mic_last)
filtered <- boundaries(filtered)
filtered <- filtered %>%
  filter(abs(distance) < 100)
filtered$name <- paste(filtered$EVENT,filtered$GENE.x, sep = "_")

rm(micro)
rm(lariats)

#load SVM-bfpinder predicted branchpoints
neu_mid <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/SVM-BPfinder/neu_micro_intron_results.txt", sep = "\t")
neu_mid <- neu_mid %>%
  separate(seq_id, into = c("name","coord"), sep = "::", remove = TRUE)
neu_mid <- neu_mid[,c("name","coord", "ss_dist", "bp_seq", "svm_scr")]
nn_mid <- read.delim("/data2/lackey_lab/randazza/microexons/gblocks/group_troubleshoot/polyU_tract/SVM-BPfinder/nn_micro_intron_results.txt", sep = "\t")
nn_mid <- nn_mid %>%
  separate(seq_id, into = c("name","coord"), sep = "::", remove = TRUE)
nn_mid <- nn_mid[,c("name","coord", "ss_dist", "bp_seq", "svm_scr")]

neu_filter <- merge(filtered, neu_mid, by = "name")
nn_filter <- merge(filtered, nn_mid, by = "name")
#neu = 9 nn = 39
neu_filter$new_name <- paste(neu_filter$name, neu_filter$bp_first, neu_filter$bp_last, sep = "_")
nn_filter$new_name <- paste(nn_filter$name, nn_filter$bp_first, nn_filter$bp_last, sep = "_")

ranked_neu <- neu_filter %>%
  group_by(new_name) %>%
  mutate(rank = min_rank(desc(svm_scr))) %>%
  ungroup()
ranked_nn <- nn_filter %>%
  group_by(new_name) %>%
  mutate(rank = min_rank(desc(svm_scr))) %>%
  ungroup()

ranked_neu <- ranked_neu %>%
  filter(abs(distance - ss_dist) <= 5)
ranked_nn <- ranked_nn %>%
  filter(abs(distance - ss_dist) <= 5)

clean_ranked_neu <- ranked_neu %>%
  mutate(distance_diff = abs(distance - ss_dist)) %>%
  group_by(new_name) %>%
  slice_min(distance_diff, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-distance_diff)

clean_ranked_nn <- ranked_nn %>%
  mutate(distance_diff = abs(distance - ss_dist)) %>%
  group_by(new_name) %>%
  slice_min(distance_diff, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-distance_diff)

clean_ranked_neu <- clean_ranked_neu %>%
  group_by(name) %>%
  distinct(ss_dist, .keep_all = TRUE) %>%
  ungroup()
clean_ranked_nn <- clean_ranked_nn %>%
  group_by(name) %>%
  distinct(ss_dist, .keep_all = TRUE) %>%
  ungroup()

clean_ranked_neu <- clean_ranked_neu %>%
  group_by(name) %>%
  filter(
    n_distinct(new_name) <= 1 |
      rank == min(rank, na.rm = TRUE)
  ) %>%
  ungroup()
clean_ranked_nn <- clean_ranked_nn %>%
  group_by(name) %>%
  filter(
    n_distinct(new_name) <= 1 |
      rank == min(rank, na.rm = TRUE)
  ) %>%
  ungroup()

combo <- rbind(clean_ranked_neu, clean_ranked_nn)

combo <- combo %>%
  mutate(
    new_rank = factor(
      if_else(rank >= 3, "other", as.character(rank)),
      levels = c("1", "2", "other")
    )
  )

combo %>%
  ggplot(aes(x = new_rank)) +
  geom_bar() +
  theme_minimal()
sum(combo$new_rank == "2", na.rm = TRUE)
break_down <- data.frame(rank = c("1", "2", "Other", "No Prediction"),
                         count = as.integer(c("34", "8", "2", "4")))
break_down$fraction <- break_down$count/48
break_down$rank <- factor(break_down$rank, levels = c("1", "2", "Other", "No Prediction"))
break_down %>%
  ggplot(aes(x = rank, y = fraction)) +
  geom_col() +
  theme_minimal()

