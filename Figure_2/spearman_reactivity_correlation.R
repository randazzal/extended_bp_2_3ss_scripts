library(dplyr)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
library(zoo)
library(beeswarm)

setwd("/data2/lackey_lab/randazza/microexons/SHAPE/")
dock <- read.delim("combined_5NIA/human_v_chicken_images/DOCK7_aligned_reactivities.txt")
dock$nucleotide <- 1:502
dock_hsa <- dock[ ,1:3]
dock_hsa$Nucleotide <- 1:502
dock_gga <- dock[ ,4:6]
dock_gga$Nucleotide <- 1:502
mef2a <- read.delim("combined_5NIA/human_v_chicken_images/MEF2A_aligned_reactivities.txt")
mef2a$nucleotide <- 1:530
mef2a_hsa <- mef2a[ ,1:3]
mef2a_hsa$Nucleotide <- 1:530
mef2a_gga <- mef2a[ ,4:6]
mef2a_gga$Nucleotide <- 1:530
dctn1 <- read.delim("combined_5NIA/human_v_chicken_images/DCTN1_aligned_reactivities.txt")
dctn1$nucleotide <- 1:545
dctn1_hsa <- dctn1[ ,1:3]
dctn1_hsa$Nucleotide <- 1:545
dctn1_gga <- dctn1[ ,4:6]
dctn1_gga$Nucleotide <- 1:545
asap2 <- read.delim("combined_5NIA/human_v_chicken_images/ASAP2_aligned_reactivities.txt")
asap2$nucleotide <- 1:517
asap2_hsa <- asap2[ ,1:3]
asap2_hsa$Nucleotide <- 1:517
asap2_gga <- asap2[ ,4:6]
asap2_gga$Nucleotide <- 1:517
cpeb4 <- read.delim("combined_5NIA/human_v_chicken_images/CPEB4_aligned_reactivities.txt")
cpeb4$nucleotide <- 1:536
cpeb4_hsa <- cpeb4[ ,1:3]
cpeb4_hsa$Nucleotide <- 1:536
cpeb4_gga <- cpeb4[ ,4:6]
cpeb4_gga$Nucleotide <- 1:536
agap1 <- read.delim("combined_5NIA/human_v_chicken_images/trimmed_AGAP1_aligned_reactivities.txt")
agap1$nucleotide <- 1:498
agap1_hsa <- agap1[ ,1:3]
agap1_hsa$Nucleotide <- 1:498
agap1_gga <- agap1[ ,4:6]
agap1_gga$Nucleotide <- 1:498
emc1 <- read.delim("combined_5NIA/human_v_chicken_images/trimmed_EMC1_aligned_reactivities.txt")
emc1$nucleotide <- 1:531
emc1_hsa <- emc1[ ,1:3]
emc1_hsa$Nucleotide <- 1:531
emc1_gga <- emc1[ ,4:6]
emc1_gga$Nucleotide <- 1:531
robo1 <- read.delim("combined_5NIA/human_v_chicken_images/trimmed_ROBO1_aligned_reactivities.txt")
robo1$nucleotide <- 1:499
robo1_hsa <- robo1[ ,1:3]
robo1_hsa$Nucleotide <- 1:499
robo1_gga <- robo1[ ,4:6]
robo1_gga$Nucleotide <- 1:499
ptprk <- read.delim("combined_5NIA/human_v_chicken_images/trimmed_PTPRK_aligned_reactivities.txt")
ptprk$nucleotide <- 1:529
ptprk_hsa <- ptprk[ ,1:3]
ptprk_hsa$Nucleotide <- 1:529
ptprk_gga <- ptprk[ ,4:6]
ptprk_gga$Nucleotide <- 1:529
fry <- read.delim("combined_5NIA/human_v_chicken_images/trimmed_FRY_aligned_reactivities.txt")
fry$nucleotide <- 1:524
fry_hsa <- fry[ ,1:3]
fry_hsa$Nucleotide <- 1:524
fry_gga <- fry[ ,4:6]
fry_gga$Nucleotide <- 1:524
clec <- read.delim("combined_5NIA/human_v_chicken_images/trimmed_CLEC16A_aligned_reactivities.txt")
clec$nucleotide <- 1:516
clec_hsa <- clec[ ,1:3]
clec_hsa$Nucleotide <- 1:516
clec_gga <- clec[ ,4:6]
clec_gga$Nucleotide <- 1:516
apbb <- read.delim("combined_5NIA/human_v_chicken_images/trimmed_APBB2_aligned_reactivities.txt")
apbb$nucleotide <- 1:532
apbb_hsa <- apbb[ ,1:3]
apbb_hsa$Nucleotide <- 1:532
apbb_gga <- apbb[ ,4:6]
apbb_gga$Nucleotide <- 1:532
itsn <- read.delim("combined_5NIA/human_v_chicken_images/trimmed_ITSN1_aligned_reactivities.txt")
itsn$nucleotide <- 1:530
itsn_hsa <- itsn[ ,1:3]
itsn_hsa$Nucleotide <- 1:530
itsn_gga <- itsn[ ,4:6]
itsn_gga$Nucleotide <- 1:530

run_spearman <- function(x, y, title) {
  df <- data.frame(x = x, y = y) %>% na.omit()
  
  # Cap values: below 0 → 0, above 4 → 4
  df$x <- pmin(pmax(df$x, 0), 4)
  df$y <- pmin(pmax(df$y, 0), 4)
  
  # Calculate Spearman correlation (suppress warning for ties)
  cor_result <- suppressWarnings(cor.test(df$x, df$y, method = "spearman"))
  rho <- cor_result$estimate
  pval <- cor_result$p.value
  
  # Plot
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_jitter(width = 0.05, height = 0.05, alpha = 0.5) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +  # optional trend line
    xlim(0, 4) + ylim(0, 4) +
    labs(x = "Human", y = "Chicken") +
    ggtitle(paste0(
      title, "\nSpearman rho = ", round(rho, 3),
      " | p = ", formatC(pval, format = "e", digits = 2)
    ))
  
  print(p)
  return(list(rho = rho, p.value = pval))
}

#exclude control hairpin
temp <- run_spearman(agap1_hsa[25:498, c("React1")], agap1_gga[25:498, c("React2")], "AGAP1")
temp <- run_spearman(dock_hsa[25:502, c("React1")], dock_gga[25:502, c("React2")], "DOCK7")
temp <- run_spearman(apbb_hsa[25:532, c("React1")], apbb_gga[25:532, c("React2")], "APBB2")
temp <- run_spearman(asap2_hsa[25:517, c("React1")], asap2_gga[25:517, c("React2")], "ASAP2")
temp <- run_spearman(clec_hsa[25:516, c("React1")], clec_gga[25:516, c("React2")], "CLEC16A")
temp <- run_spearman(cpeb4_hsa[25:536, c("React1")], cpeb4_gga[25:536, c("React2")], "CPEB4")
temp <- run_spearman(dctn1_hsa[25:545, c("React1")], dctn1_gga[25:545, c("React2")], "DCTN1")
temp <- run_spearman(emc1_hsa[25:531, c("React1")], emc1_gga[25:531, c("React2")], "EMC1")
temp <- run_spearman(fry_hsa[25:524, c("React1")], fry_gga[25:524, c("React2")], "FRY")
temp <- run_spearman(itsn_hsa[25:530, c("React1")], itsn_gga[25:530, c("React2")], "ITSN1")
temp <- run_spearman(mef2a_hsa[25:530, c("React1")], mef2a_gga[25:530, c("React2")], "MEF2A")
temp <- run_spearman(ptprk_hsa[25:529, c("React1")], ptprk_gga[25:529, c("React2")], "PTPRK")
temp <- run_spearman(robo1_hsa[25:499, c("React1")], robo1_gga[25:499, c("React2")], "ROBO1")